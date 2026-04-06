#!/usr/bin/env python3
# ============================================================================
#  Step 4: 差分発現解析（Differential Expression Analysis）— Figure 2 の再現
# ============================================================================
#
#  【概要】
#  このスクリプトは、前処理済みのプロテオミクスデータに対して
#  差分発現解析（Differential Expression Analysis）を行い、
#  腫瘍（Tumor）と正常（Normal）の間で発現量が大きく異なるタンパク質を
#  統計的に特定します。
#
#  【処理の流れ】
#    1. Welch の t 検定（Welch's t-test）で各タンパク質の群間差を検定
#    2. p値 < 0.05 かつ fold change > 2 または < 0.5 でフィルタリング
#       （log2 fold change で ±1.0 に相当）
#    3. Volcano プロットで結果を可視化
#    4. 有意差のあるトップ N タンパク質で階層的クラスタリングと PCA を実行
#
#  【Welch の t 検定とは？】
#    通常の t 検定（Student's t-test）は「2群の分散が等しい」という仮定が必要
#    ですが、Welch の t 検定はこの仮定を緩和したバージョンです。
#    生物学的データでは群間で分散が異なることが多いため、
#    Welch の t 検定のほうが適切です。
#    scipy では equal_var=False を指定することで Welch 版になります。
#
#  【Volcano プロットとは？】
#    Volcano プロット（火山プロット）は、差分発現解析の結果を
#    一目で把握するための散布図です。
#    - X軸: log2 fold change（発現量の変化の大きさ）
#    - Y軸: -log10(p値)（統計的有意性の高さ）
#    火山のような形になることが多いためこの名前がついています。
#    右上のドットは「腫瘍で有意に増加」、左上は「腫瘍で有意に減少」を示します。
#
#  【入力ファイル】
#    - results/preprocessed_data.csv   : 前処理済みのタンパク質発現量行列
#    - results/sample_info.csv         : 各サンプルの条件（Normal / Tumor）情報
#
#  【出力ファイル】
#    - results/tables/differential_proteins.csv  : 全タンパク質の検定結果テーブル
#    - results/figures/fig_bonus_volcano.png      : Volcano プロット（ボーナス図）
#    - results/figures/fig2a_heatmap_all.png      : 全有意差タンパク質のヒートマップ
#    - results/figures/fig2b_pca_all.png          : 全有意差タンパク質の PCA
#    - results/figures/fig2c_clustering_top50.png : トップ 50 の階層的クラスタリング
#    - results/figures/fig2d_clustering_top100.png: トップ 100 の階層的クラスタリング
#    - results/figures/fig2e_clustering_top200.png: トップ 200 の階層的クラスタリング
#    - results/figures/fig2f_pca_top50.png        : トップ 50 の PCA プロット
#    - results/figures/fig2g_pca_top100.png       : トップ 100 の PCA プロット
#    - results/figures/fig2h_pca_top200.png       : トップ 200 の PCA プロット
# ============================================================================

# ============================================================================
#  ライブラリのインポート
# ============================================================================

# os: ファイルパスの結合やディレクトリ作成など、OS レベルの操作に使用
import os

# numpy: 数値計算ライブラリ。log10 や配列操作など数学的な処理に使用
import numpy as np

# pandas: データフレーム（表形式データ）の読み込み・操作・保存に使用
#   CSV の読み書き、行列のスライス、統計量の計算などを担う
import pandas as pd

# matplotlib.pyplot: グラフ描画ライブラリ。散布図（Volcano プロット）や
#   PCA プロットなどの図を作成・保存するために使用
import matplotlib.pyplot as plt

# matplotlib.patches.Ellipse: 信頼楕円（confidence ellipse）を描画するために使用
from matplotlib.patches import Ellipse

# matplotlib.transforms: 楕円のアフィン変換（回転・拡大・平行移動）に使用
import matplotlib.transforms as transforms

# seaborn: matplotlib をベースにした統計可視化ライブラリ。
#   clustermap（階層的クラスタリング付きヒートマップ）の描画に使用
import seaborn as sns

# scipy.stats: 科学計算用ライブラリ scipy の統計モジュール。
#   ttest_ind（独立2標本の t 検定）を使って Welch's t-test を実行するために使用
from scipy import stats

# sklearn.decomposition.PCA: scikit-learn の主成分分析（PCA）クラス。
#   高次元のタンパク質発現データを 2次元に圧縮して可視化するために使用
from sklearn.decomposition import PCA

# ============================================================================
#  パス設定（プロジェクト構造に基づくディレクトリ定義）
# ============================================================================

# このスクリプト自身が置かれているディレクトリの絶対パスを取得
# __file__ はこのスクリプトのパス、abspath で絶対パスに変換し、dirname で親ディレクトリを取得
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# プロジェクトルートディレクトリ（scripts/ の1つ上の階層）
PROJECT_DIR = os.path.join(SCRIPT_DIR, "..")

# 解析結果全般を格納するディレクトリ
RESULTS_DIR = os.path.join(PROJECT_DIR, "results")

# 図を保存するディレクトリ（results/figures/）
FIG_DIR = os.path.join(RESULTS_DIR, "figures")

# 表（CSV）を保存するディレクトリ（results/tables/）
TABLE_DIR = os.path.join(RESULTS_DIR, "tables")

# 出力先ディレクトリが存在しなければ作成する（exist_ok=True で既存でもエラーにならない）
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

# ============================================================================
#  定数・閾値の定義
# ============================================================================

# 【p値の閾値】
#   統計的に有意とみなす p 値の上限。0.05 は「5%の確率で偶然こうなる」という意味。
#   p < 0.05 なら「偶然ではなく本当に差がある可能性が高い」と判断する。
P_VALUE_THRESHOLD = 0.05

# 【Log2 Fold Change の閾値】
#   log2(2) = 1.0 なので、LOG2FC_THRESHOLD = 1.0 は fold change > 2 に対応する。
#   つまり、腫瘍群の平均発現量が正常群の2倍以上（または1/2以下）のタンパク質を
#   「生物学的に意味のある変化」として抽出する。
#   - log2FC > 1.0  → 腫瘍で2倍以上に増加（Up-regulated）
#   - log2FC < -1.0 → 腫瘍で1/2以下に減少（Down-regulated）
LOG2FC_THRESHOLD = 1.0  # log2(2) = 1.0 → fold change > 2

# ============================================================================
#  カラーパレットの定義（図の配色）
# ============================================================================

# 正常サンプルの色（青）— PCA・クラスタリングのサイドバーに使用
COLOR_NORMAL = "#3498DB"

# 腫瘍サンプルの色（赤）— PCA・クラスタリングのサイドバーに使用
COLOR_TUMOR = "#E74C3C"

# Volcano プロットで「腫瘍で発現上昇（Up-regulated）」を示す色（赤）
COLOR_UP = "#E74C3C"

# Volcano プロットで「腫瘍で発現低下（Down-regulated）」を示す色（青）
COLOR_DOWN = "#3498DB"

# Volcano プロットで「有意差なし（Not Significant）」を示す色（灰色）
COLOR_NS = "#CCCCCC"


# ============================================================================
#  関数定義: 95%信頼楕円の描画
# ============================================================================

def confidence_ellipse(x, y, ax, n_std=2.0, **kwargs):
    """
    散布図上に信頼楕円を描画する関数。

    【何をするか】
      2次元データ (x, y) の共分散行列を計算し、データの散らばりを表す
      楕円を描画する。n_std=2.0 は約95%の信頼区間に対応する。

    【引数】
      x     : array-like — X座標のデータ
      y     : array-like — Y座標のデータ
      ax    : matplotlib.axes.Axes — 楕円を描画する軸オブジェクト
      n_std : float — 標準偏差の倍率（2.0 で約95%信頼区間）
      **kwargs : Ellipse に渡す追加パラメータ（色、透明度など）

    【戻り値】
      matplotlib.patches.Ellipse — 描画された楕円パッチ
    """
    if len(x) < 2:
        return
    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, **kwargs)
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)
    transf = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)
    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


# ============================================================================
#  関数定義: データ読み込み
# ============================================================================

def load_data():
    """
    前処理済みデータとサンプル情報を CSV から読み込む関数。

    【何をするか】
      Step 3（前処理）で出力された2つの CSV ファイルを pandas の DataFrame として
      メモリに読み込みます。

    【なぜ必要か】
      差分発現解析には、タンパク質ごとの発現量データと、
      各サンプルがどの群（Normal / Tumor）に属するかの情報が必要なため。

    【引数】
      なし（ファイルパスはグローバル定数 RESULTS_DIR から自動構築）

    【戻り値】
      df          : pandas.DataFrame — 行がタンパク質、列がサンプルの発現量行列。
                    index_col=0 により最初の列（タンパク質名）が行インデックスになる。
      sample_info : pandas.DataFrame — 各サンプルの条件（"Normal" or "Tumor"）を含む表。
                    "Sample" 列と "Condition" 列を持つ。
    """
    # preprocessed_data.csv: 行=タンパク質, 列=サンプル の発現量行列
    # index_col=0 で最初の列（タンパク質名）を行インデックスとして読み込む
    df = pd.read_csv(os.path.join(RESULTS_DIR, "preprocessed_data.csv"), index_col=0)

    # sample_info.csv: 各サンプルが Normal / Tumor のどちらかを記録した表
    sample_info = pd.read_csv(os.path.join(RESULTS_DIR, "sample_info.csv"))

    return df, sample_info


# ============================================================================
#  関数定義: Welch's t-test（差分発現解析の中核）
# ============================================================================

def welch_ttest(df, normal_samples, tumor_samples):
    """
    全タンパク質に対して Welch の t 検定を実行し、発現変動の有意性を判定する関数。

    【何をするか】
      1. 各タンパク質について Normal 群と Tumor 群の発現量を取り出す
      2. Welch の t 検定で p 値を計算する
      3. Log2 Fold Change（発現量の変化量）を計算する
      4. p値と fold change の閾値に基づいて Up / Down / NS に分類する

    【Welch の t 検定の詳細説明】
      - 帰無仮説（H0）: 「2群の母平均は等しい」（つまり差がない）
      - 対立仮説（H1）: 「2群の母平均は異なる」（つまり差がある）
      - p 値が小さいほど「差がない」という仮説が棄却され、「差がある」と判断できる
      - equal_var=False を指定する理由:
          通常の t 検定（Student's t-test）は2群の分散が等しいと仮定するが、
          生物学的データでは腫瘍群と正常群で分散が異なることが一般的。
          equal_var=False にすると Welch の補正が適用され、
          不等分散でも正確な検定ができる。自由度が Welch-Satterthwaite の式で
          近似的に計算される。

    【Log2 Fold Change の意味】
      - データが既に log2 スケールの場合、差（引き算）がそのまま log2 FC になる
      - log2FC = mean(Tumor) - mean(Normal)
      - log2FC > 0 → 腫瘍で発現上昇、log2FC < 0 → 腫瘍で発現低下

    【引数】
      df             : pandas.DataFrame — 行=タンパク質, 列=サンプル の発現量行列
      normal_samples : list of str — Normal 群に属するサンプル名のリスト
      tumor_samples  : list of str — Tumor 群に属するサンプル名のリスト

    【戻り値】
      result_df : pandas.DataFrame — 以下の列を持つ検定結果テーブル:
        - Protein      : タンパク質名
        - Mean_Normal  : Normal 群の平均発現量
        - Mean_Tumor   : Tumor 群の平均発現量
        - Log2FC       : Log2 Fold Change (Tumor / Normal)
        - T_statistic  : t 統計量（検定統計量）
        - P_value      : p 値（有意確率）
        - Neg_log10_P  : -log10(p値)（Volcano プロットの Y 軸に使用）
        - Significant  : "Up" / "Down" / "NS" の分類ラベル

    論文: p < 0.05, fold change > 2 or < 0.5
    """
    # 各タンパク質の検定結果を格納するリスト
    results = []

    # --- 全タンパク質に対してループ処理で検定を実行 ---
    for protein in df.index:
        # このタンパク質の Normal 群のデータを取得し、欠損値（NaN）を除外
        normal_vals = df.loc[protein, normal_samples].dropna()
        # このタンパク質の Tumor 群のデータを取得し、欠損値（NaN）を除外
        tumor_vals = df.loc[protein, tumor_samples].dropna()

        # 各群に最低2サンプル必要（t 検定には最低2つのデータ点が必要なため）
        # 1つ以下だと分散が計算できず t 統計量が定義できない
        if len(normal_vals) < 2 or len(tumor_vals) < 2:
            continue

        # 【Welch の t 検定を実行】
        #   stats.ttest_ind(): 独立2標本の t 検定
        #   - 第1引数: tumor_vals（腫瘍群の値）
        #   - 第2引数: normal_vals（正常群の値）
        #   - equal_var=False: Welch 補正を使用（分散が等しいと仮定しない）
        #   戻り値:
        #   - t_stat: t 統計量（正なら腫瘍群のほうが高い、負なら低い）
        #   - p_val : p 値（0に近いほど有意差が大きい）
        t_stat, p_val = stats.ttest_ind(tumor_vals, normal_vals, equal_var=False)

        # 【Log2 Fold Change の計算】
        #   データは既に log2 スケールに変換されていると想定。
        #   log2 スケールでは: log2(A/B) = log2(A) - log2(B)
        #   よって平均値の差がそのまま log2 fold change になる。
        #   正の値 → 腫瘍で発現上昇、負の値 → 腫瘍で発現低下
        log2fc = tumor_vals.mean() - normal_vals.mean()

        # 検定結果を辞書形式でリストに追加
        results.append({
            "Protein": protein,                                # タンパク質名
            "Mean_Normal": normal_vals.mean(),                 # Normal群の平均値
            "Mean_Tumor": tumor_vals.mean(),                   # Tumor群の平均値
            "Log2FC": log2fc,                                  # Log2 Fold Change
            "T_statistic": t_stat,                             # t統計量
            "P_value": p_val,                                  # p値（有意確率）
            # -log10(p値): Volcano プロットの Y 軸に使用
            # p値が小さいほど -log10(p) は大きくなる（上に位置する）
            # max(p_val, 1e-300) は p_val=0 のときに log10(0) = -inf を避けるための安全策
            "Neg_log10_P": -np.log10(max(p_val, 1e-300)),
        })

    # リストを DataFrame に変換（行が各タンパク質、列が統計量の表）
    result_df = pd.DataFrame(results)

    # ========================================================================
    #  有意性の分類（Significant カラムの付与）
    # ========================================================================
    #  各タンパク質を以下の3カテゴリに分類する:
    #    "Up"   : 腫瘍で有意に発現上昇（p < 0.05 かつ log2FC > 1.0）
    #    "Down" : 腫瘍で有意に発現低下（p < 0.05 かつ log2FC < -1.0）
    #    "NS"   : 有意差なし（Not Significant）— 上記に該当しないもの

    # まず全タンパク質を "NS"（有意差なし）で初期化
    result_df["Significant"] = "NS"

    # 【Up-regulated の条件マスク】
    #   p値が閾値未満 AND log2FC が正の閾値より大きい
    #   → 腫瘍で統計的に有意に発現量が増加しているタンパク質
    up_mask = (result_df["P_value"] < P_VALUE_THRESHOLD) & (result_df["Log2FC"] > LOG2FC_THRESHOLD)

    # 【Down-regulated の条件マスク】
    #   p値が閾値未満 AND log2FC が負の閾値より小さい
    #   → 腫瘍で統計的に有意に発現量が減少しているタンパク質
    down_mask = (result_df["P_value"] < P_VALUE_THRESHOLD) & (result_df["Log2FC"] < -LOG2FC_THRESHOLD)

    # 条件に合致する行の "Significant" 列を "Up" または "Down" に更新
    result_df.loc[up_mask, "Significant"] = "Up"
    result_df.loc[down_mask, "Significant"] = "Down"

    return result_df


# ============================================================================
#  関数定義: Volcano プロット描画（ボーナス図）
# ============================================================================

def plot_volcano(result_df):
    """
    差分発現解析の結果を Volcano プロット（火山プロット）として可視化する関数。

    【Volcano プロットとは】
      差分発現解析の結果を俯瞰するための散布図。
      - X軸: Log2 Fold Change（発現量の変化の大きさと方向）
        - 右に行くほど腫瘍で発現上昇
        - 左に行くほど腫瘍で発現低下
      - Y軸: -Log10(P-value)（統計的有意性の高さ）
        - 上に行くほど統計的に有意（p値が小さい）
      閾値線で区切られた右上が「有意に上昇」、左上が「有意に低下」を示す。
      名前の由来は、ドットの分布が火山の噴火に見えることから。

    【引数】
      result_df : pandas.DataFrame — welch_ttest() が返した検定結果テーブル。
                  "Log2FC", "Neg_log10_P", "Significant" の列が必要。

    【戻り値】
      なし（PNG ファイルとして保存）
    """
    print("  Volcanoプロット（ボーナス図）")

    # --- 図と軸（Axes）オブジェクトの作成 ---
    # figsize=(8, 6): 図のサイズを 幅8インチ x 高さ6インチ に設定
    fig, ax = plt.subplots(figsize=(8, 6))

    # --- 3つのカテゴリごとにドットを描画 ---
    # NS（灰色・半透明）→ Up（赤）→ Down（青）の順で描画
    # 後から描画したほうが前面に来るため、有意なドットが目立つようにする
    for sig, color, alpha in [
        ("NS", COLOR_NS, 0.3),     # Not Significant: 灰色、透明度70%（alpha=0.3）
        ("Up", COLOR_UP, 0.6),     # Up-regulated: 赤色、透明度40%（alpha=0.6）
        ("Down", COLOR_DOWN, 0.6), # Down-regulated: 水色、透明度40%（alpha=0.6）
    ]:
        # 現在のカテゴリに属するタンパク質のみを選択するブールマスク
        mask = result_df["Significant"] == sig
        ax.scatter(
            result_df.loc[mask, "Log2FC"],       # X軸: Log2 Fold Change
            result_df.loc[mask, "Neg_log10_P"],  # Y軸: -Log10(P-value)
            c=color,                              # ドットの色
            s=10,                                 # ドットのサイズ（ポイント単位）
            alpha=alpha,                          # 透明度（0=完全透明, 1=不透明）
            label=f"{sig} ({mask.sum()})",        # 凡例ラベル（カテゴリ名と個数）
        )

    # --- 閾値線の描画（点線） ---
    # 水平線: p 値の閾値（-log10(0.05) ≈ 1.3）
    # この線より上のドットは統計的に有意
    ax.axhline(
        -np.log10(P_VALUE_THRESHOLD),  # Y座標: -log10(0.05) ≈ 1.301
        color="gray",                  # 線の色（灰色）
        linestyle="--",                # 線のスタイル（破線）
        linewidth=0.5,                 # 線の太さ（0.5ポイント）
    )
    # 垂直線（右側）: log2FC の正の閾値（+1.0 = fold change 2倍）
    ax.axvline(
        LOG2FC_THRESHOLD,              # X座標: +1.0
        color="gray",
        linestyle="--",
        linewidth=0.5,
    )
    # 垂直線（左側）: log2FC の負の閾値（-1.0 = fold change 0.5倍）
    ax.axvline(
        -LOG2FC_THRESHOLD,             # X座標: -1.0
        color="gray",
        linestyle="--",
        linewidth=0.5,
    )

    # --- 軸ラベルとタイトルの設定 ---
    ax.set_xlabel("Log2 Fold Change (Tumor / Normal)")              # X軸のラベル
    ax.set_ylabel("-Log10(P-value)")                                 # Y軸のラベル
    ax.set_title("Differential Protein Abundance: Tumor vs Non-tumor")  # 図のタイトル

    # --- 凡例の設定 ---
    # frameon=False: 凡例の枠線を非表示にしてスッキリさせる
    # loc="upper right": 凡例を右上に配置
    ax.legend(frameon=False, loc="upper right")

    # --- 枠線（spines）の非表示 ---
    # 上と右の枠線を消すことで、すっきりした学術論文風のデザインにする
    ax.spines["top"].set_visible(False)     # 上辺の枠線を非表示
    ax.spines["right"].set_visible(False)   # 右辺の枠線を非表示

    # --- 図の保存 ---
    # ボーナス図として保存（論文の Figure 2 には含まれない補足図）
    filepath = os.path.join(FIG_DIR, "fig_bonus_volcano.png")
    fig.savefig(
        filepath,
        dpi=150,                # 解像度: 150 DPI（Dots Per Inch）。高いほど高精細
        bbox_inches="tight",    # 余白を最小限にトリミングして保存
    )
    # メモリ解放のためにプロットを閉じる（大量の図を生成する場合に重要）
    plt.close()
    print(f"    保存: {filepath}")


# ============================================================================
#  関数定義: 全有意差タンパク質のヒートマップ（Figure 2a）
# ============================================================================

def plot_heatmap_all(df, result_df, sample_info):
    """
    全ての有意差タンパク質（2642個）を使ったクラスタリング付きヒートマップを描画する関数。
    論文の Figure 2a に対応する。

    【何をするか】
      有意差のある全タンパク質（Up + Down）の発現量データを用いて、
      行=タンパク質、列=サンプル の階層的クラスタリング付きヒートマップを作成する。
      行と列の両方をクラスタリングし、サンプルの条件（Normal/Tumor）と
      タンパク質の変化方向（Up/Down）をサイドカラーで表示する。

    【引数】
      df          : pandas.DataFrame — 行=タンパク質, 列=サンプル の発現量行列
      result_df   : pandas.DataFrame — welch_ttest() が返した検定結果テーブル
      sample_info : pandas.DataFrame — サンプルの条件情報（Sample, Condition列）

    【戻り値】
      なし（PNG ファイルとして保存）
    """
    print("  Figure 2a: 全有意差タンパク質のヒートマップ")

    # sample_info から Condition（Normal/Tumor）を引くための Series を作成
    conditions = sample_info.set_index("Sample")["Condition"]

    # --- 有意差のあるタンパク質のみを抽出 ---
    sig_df = result_df[result_df["Significant"] != "NS"].copy()

    if len(sig_df) == 0:
        print("    スキップ: 有意なタンパク質が0個")
        return

    # 有意差タンパク質の発現量行列を取得（行=タンパク質、列=サンプル）
    sig_proteins = sig_df["Protein"].values
    sig_data = df.loc[df.index.isin(sig_proteins)]

    # --- 列（サンプル）のカラーバーを作成 ---
    # Normal → 青、Tumor → 赤
    col_colors = conditions.reindex(df.columns).map({
        "Normal": COLOR_NORMAL,
        "Tumor": COLOR_TUMOR,
    })

    # --- 列（タンパク質）のカラーバーを作成（上部に表示） ---
    # 【論文準拠】タンパク質が列なので、上部のカラーバーで Up/Down を表示
    protein_direction = sig_df.set_index("Protein")["Log2FC"]
    col_color_map = sig_data.index.map(
        lambda p: COLOR_UP if protein_direction.get(p, 0) > 0 else COLOR_DOWN
    )
    protein_colors = pd.Series(col_color_map, index=sig_data.index, name="Expression")

    # --- 行（サンプル）のカラーバーを作成（左側に表示） ---
    sample_row_colors = conditions.reindex(sig_data.columns).map({
        "Normal": COLOR_NORMAL,
        "Tumor": COLOR_TUMOR,
    })

    # 【論文準拠の向き】行=サンプル、列=タンパク質
    g = sns.clustermap(
        sig_data.T,                # 転置: 行=サンプル, 列=タンパク質
        method="ward",             # クラスタリング手法: Ward法
        cmap="RdBu_r",            # カラーマップ: 赤=高発現, 青=低発現
        z_score=1,                 # 列方向（タンパク質ごと）で Z スコア正規化
        row_cluster=True,          # 行（サンプル）をクラスタリング
        col_cluster=True,          # 列（タンパク質）をクラスタリング
        row_colors=sample_row_colors,   # 左のカラーバー: Condition
        col_colors=protein_colors,      # 上のカラーバー: Expression
        xticklabels=False,         # X軸（タンパク質名）は非表示
        yticklabels=True,          # Y軸（サンプル名）を表示
        figsize=(10, 8),           # 論文と同じプロポーション
        vmin=-3, vmax=3,           # Z-scoreの色範囲を固定
    )

    # Y軸ラベル（サンプル名）のフォントサイズを調整
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=7)

    # ヒートマップを PNG として保存
    filepath = os.path.join(FIG_DIR, "fig2a_heatmap_all.png")
    g.savefig(filepath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"    保存: {filepath}")


# ============================================================================
#  関数定義: 全有意差タンパク質の PCA（Figure 2b）
# ============================================================================

def plot_pca_all(df, result_df, sample_info):
    """
    全ての有意差タンパク質を使った PCA プロットを描画する関数。
    論文の Figure 2b に対応する。

    【引数】
      df          : pandas.DataFrame — 行=タンパク質, 列=サンプル の発現量行列
      result_df   : pandas.DataFrame — welch_ttest() が返した検定結果テーブル
      sample_info : pandas.DataFrame — サンプルの条件情報（Sample, Condition列）

    【戻り値】
      なし（PNG ファイルとして保存）
    """
    print("  Figure 2b: 全有意差タンパク質の PCA")

    # sample_info から Condition（Normal/Tumor）を引くための Series を作成
    conditions = sample_info.set_index("Sample")["Condition"]

    # --- 有意差のあるタンパク質のみを抽出 ---
    sig_df = result_df[result_df["Significant"] != "NS"].copy()

    if len(sig_df) == 0:
        print("    スキップ: 有意なタンパク質が0個")
        return

    # 有意差タンパク質の発現量行列を取得
    sig_proteins = sig_df["Protein"].values
    sig_data = df.loc[df.index.isin(sig_proteins)]

    # PCA オブジェクトを作成（2成分に圧縮）
    pca = PCA(n_components=2)
    # fit_transform(): データに PCA を適用し、変換後の座標を返す
    # sig_data.T: 転置して行=サンプル, 列=タンパク質にする
    scores = pca.fit_transform(sig_data.T)

    # --- PCA 散布図の描画 ---
    fig, ax = plt.subplots(figsize=(8, 6))

    # Normal と Tumor をそれぞれ異なる色で描画（両方とも丸マーカー "o"）
    for condition, color in [
        ("Normal", COLOR_NORMAL),   # Normal: 青の丸
        ("Tumor", COLOR_TUMOR),     # Tumor: 赤の丸
    ]:
        # 現在の条件に属するサンプルのブールマスクを作成
        mask = conditions.reindex(df.columns) == condition
        ax.scatter(
            scores[mask, 0],          # X座標: PC1 のスコア
            scores[mask, 1],          # Y座標: PC2 のスコア
            c=color,                  # ドットの色
            marker="o",               # マーカーの形（丸）
            s=80,                     # マーカーのサイズ（80ポイント）
            alpha=0.8,                # 透明度
            label=condition,          # 凡例ラベル
            edgecolors="white",       # マーカーの縁の色
        )
        # 95%信頼楕円を描画
        confidence_ellipse(
            scores[mask, 0], scores[mask, 1], ax,
            n_std=2.0, facecolor=color, alpha=0.15, edgecolor=color, linewidth=1.5,
        )

    # --- 軸ラベルに分散寄与率を表示 ---
    ax.set_xlabel(f"Component 1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"Component 2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title("PCA using All Differentially Abundant Proteins")

    # 凡例と枠線の設定
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # 薄いグリッド線を追加
    ax.grid(True, alpha=0.3, linestyle="--")

    # PCA プロットを PNG として保存
    filepath = os.path.join(FIG_DIR, "fig2b_pca_all.png")
    fig.savefig(filepath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"    保存: {filepath}")


# ============================================================================
#  関数定義: Top N タンパク質での解析（クラスタリング + PCA）
# ============================================================================

def plot_top_n_analysis(df, result_df, sample_info, n_top=50):
    """
    有意差のあるタンパク質のうち、変化量の大きい上位 N 個を用いて
    階層的クラスタリングと PCA を実行し、可視化する関数（Figure 2c-h の再現）。

    【何をするか】
      1. 有意差のあるタンパク質を |log2FC| の大きい順にソート
      2. 上位 50, 100, 200 個について以下を実行:
         a. 階層的クラスタリング（Ward法）付きヒートマップ
         b. PCA（主成分分析）による2次元プロット

    【階層的クラスタリングとは】
      似たパターンを持つサンプルやタンパク質をグループ化する手法。
      Ward 法は、クラスタ内の分散を最小化する方法で、
      バランスの取れたクラスタが得られやすい。

    【PCA（主成分分析）とは】
      高次元データ（数千のタンパク質）を、情報をなるべく保持したまま
      2次元や3次元に圧縮する手法。第1主成分（PC1）がデータの分散を
      最も説明し、第2主成分（PC2）がその次に説明する。
      NormalとTumorがPCA上で分離すれば、選んだタンパク質セットが
      2群を区別する能力を持つことを意味する。

    【引数】
      df          : pandas.DataFrame — 行=タンパク質, 列=サンプル の発現量行列
      result_df   : pandas.DataFrame — welch_ttest() が返した検定結果テーブル
      sample_info : pandas.DataFrame — サンプルの条件情報（Sample, Condition列）
      n_top       : int — デフォルト値50。実際にはループで50, 100, 200を試行する

    【戻り値】
      なし（PNG ファイルとして保存）
    """
    # sample_info から Condition（Normal/Tumor）を引くための Series を作成
    # set_index("Sample") で Sample 列をインデックスにし、Condition 列だけ取り出す
    # → conditions["SampleA"] = "Normal" のようにアクセスできるようになる
    conditions = sample_info.set_index("Sample")["Condition"]

    # --- 有意差のあるタンパク質のみを抽出 ---
    # "NS"（有意差なし）以外 = "Up" または "Down" のタンパク質
    sig_df = result_df[result_df["Significant"] != "NS"].copy()

    # --- Up / Down を分離 ---
    # 論文では「top 50 proteins with increased AND decreased abundance」と記載。
    # つまり上昇Top N + 低下Top N を別々に選択して結合する。
    # 例: "top 50" = 上昇Top50 + 低下Top50 = 合計100タンパク質
    sig_up = sig_df[sig_df["Significant"] == "Up"].copy()
    sig_down = sig_df[sig_df["Significant"] == "Down"].copy()

    # 各方向で |fold change| の降順にソート
    sig_up["AbsLog2FC"] = sig_up["Log2FC"].abs()
    sig_down["AbsLog2FC"] = sig_down["Log2FC"].abs()
    sig_up = sig_up.sort_values("AbsLog2FC", ascending=False)
    sig_down = sig_down.sort_values("AbsLog2FC", ascending=False)

    # --- Figure 2c-e のパネル番号マッピング ---
    # Top 50 → fig2c, Top 100 → fig2d, Top 200 → fig2e
    heatmap_labels = {50: "c", 100: "d", 200: "e"}
    # Top 50 → fig2f, Top 100 → fig2g, Top 200 → fig2h
    pca_labels = {50: "f", 100: "g", 200: "h"}

    # --- Top 50, 100, 200 それぞれで解析を実行 ---
    for n in [50, 100, 200]:
        # 有意差のあるタンパク質が n 個未満の場合は、あるだけ使う
        if len(sig_df) < n:
            n = len(sig_df)
            # 有意差のあるタンパク質が0個の場合は解析をスキップ
            if n == 0:
                print(f"  スキップ: 有意なタンパク質が0個")
                return

        # 【論文準拠】上昇Top N + 低下Top N を別々に選択して結合
        # 例: n=50 → 上昇50 + 低下50 = 合計100タンパク質
        n_up = min(n, len(sig_up))
        n_down = min(n, len(sig_down))
        top_up = sig_up.head(n_up)["Protein"].values
        top_down = sig_down.head(n_down)["Protein"].values
        top_proteins = np.concatenate([top_up, top_down])

        # 発現量行列から選択したタンパク質の行だけを抽出
        top_data = df.loc[df.index.isin(top_proteins)]
        print(f"    （上昇{n_up} + 低下{n_down} = 合計{len(top_data)} タンパク質）")

        # ================================================================
        #  階層的クラスタリング付きヒートマップの描画（Figure 2c-e）
        # ================================================================
        print(f"  Figure 2{heatmap_labels.get(n, '')}: Top {n} タンパク質: 階層的クラスタリング")

        # --- 列（サンプル）のカラーバーを作成 ---
        # 各サンプルの条件を色にマッピング（Normal→青, Tumor→赤）
        # これはヒートマップのサイドバーに使う色情報
        col_colors = conditions.reindex(df.columns).map({
            "Normal": COLOR_NORMAL,
            "Tumor": COLOR_TUMOR,
        })

        # --- 列（タンパク質）のカラーバーを作成（上部に表示） ---
        # 論文Figure 2c-eでは、タンパク質が列、サンプルが行。
        # 上部のカラーバーで Up（赤）/ Down（青）を色分けする。
        protein_direction = pd.concat([sig_up, sig_down]).set_index("Protein")["Log2FC"]
        col_color_map = top_data.index.map(
            lambda p: COLOR_UP if protein_direction.get(p, 0) > 0 else COLOR_DOWN
        )
        protein_colors = pd.Series(col_color_map, index=top_data.index, name="Expression")

        # --- 行（サンプル）のカラーバーを作成（左側に表示） ---
        # Normal=青, Tumor=赤 でサンプルを色分け
        sample_row_colors = conditions.reindex(df.columns).map({
            "Normal": COLOR_NORMAL,
            "Tumor": COLOR_TUMOR,
        })

        # 【論文準拠の向き】行=サンプル、列=タンパク質
        # top_data.T で転置: タンパク質(行)×サンプル(列) → サンプル(行)×タンパク質(列)
        g = sns.clustermap(
            top_data.T,                # 転置: 行=サンプル, 列=タンパク質（論文と同じ向き）
            method="ward",             # クラスタリング手法: Ward法
            cmap="RdBu_r",            # カラーマップ: 赤=高発現, 青=低発現
            z_score=1,                 # 列方向（タンパク質ごと）で Z スコア正規化
            row_cluster=True,          # 行（サンプル）をクラスタリング
            col_cluster=True,          # 列（タンパク質）をクラスタリング
            row_colors=sample_row_colors,   # 左のカラーバー: Condition（Normal=青, Tumor=赤）
            col_colors=protein_colors,      # 上のカラーバー: Expression（Up=赤, Down=青）
            figsize=(8, 10),           # 縦長の図（論文と同じプロポーション）
            xticklabels=False,         # X軸（タンパク質名）は非表示（多すぎるため）
            yticklabels=True,          # Y軸（サンプル名）を表示
            vmin=-3, vmax=3,           # Z-score の色範囲を固定
        )
        # ヒートマップの Y 軸ラベル（サンプル名）のフォントサイズを7ptに縮小
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=7)

        # クラスタリングヒートマップを PNG として保存
        panel = heatmap_labels.get(n, "")
        filepath = os.path.join(FIG_DIR, f"fig2{panel}_clustering_top{n}.png")
        g.savefig(filepath, dpi=150, bbox_inches="tight")
        plt.close()  # メモリ解放
        print(f"    保存: {filepath}")

        # ================================================================
        #  PCA（主成分分析）プロットの描画（Figure 2f-h）
        # ================================================================
        panel_pca = pca_labels.get(n, "")
        print(f"  Figure 2{panel_pca}: Top {n} タンパク質: PCA")

        # PCA オブジェクトを作成（2成分に圧縮）
        # n_components=2: 第1主成分（PC1）と第2主成分（PC2）の2次元に削減
        pca = PCA(n_components=2)

        # fit_transform(): データに PCA を適用し、変換後の座標を返す
        # top_data.T: 転置して行=サンプル, 列=タンパク質にする
        #   （PCA は行を観測点、列を特徴量として扱うため）
        # scores: 各サンプルの PC1, PC2 座標を格納した 2D配列（shape: サンプル数 x 2）
        scores = pca.fit_transform(top_data.T)

        # --- PCA 散布図の描画 ---
        # figsize=(8, 6): 幅8 x 高さ6 インチ
        fig, ax = plt.subplots(figsize=(8, 6))

        # Normal と Tumor をそれぞれ異なる色で描画（両方とも丸マーカー "o"）
        for condition, color in [
            ("Normal", COLOR_NORMAL),   # Normal: 青の丸（circle）
            ("Tumor", COLOR_TUMOR),     # Tumor: 赤の丸（circle）
        ]:
            # 現在の条件に属するサンプルのブールマスクを作成
            # reindex で発現量行列の列（サンプル）の順序に合わせる
            mask = conditions.reindex(df.columns) == condition
            ax.scatter(
                scores[mask, 0],          # X座標: PC1 のスコア
                scores[mask, 1],          # Y座標: PC2 のスコア
                c=color,                  # ドットの色
                marker="o",               # マーカーの形（丸 "o" のみ使用）
                s=80,                     # マーカーのサイズ（80ポイント）
                alpha=0.8,                # 透明度（0.8 = やや透明）
                label=condition,          # 凡例ラベル
                edgecolors="white",       # マーカーの縁の色（白で囲むと見やすくなる）
            )
            # 95%信頼楕円を描画
            confidence_ellipse(
                scores[mask, 0], scores[mask, 1], ax,
                n_std=2.0, facecolor=color, alpha=0.15, edgecolor=color, linewidth=1.5,
            )

        # --- 軸ラベルに分散寄与率を表示 ---
        # explained_variance_ratio_ には各主成分が全分散の何%を説明するかが格納されている
        # 例: PC1 が 45.3%, PC2 が 12.1% なら → "Component 1 (45.3%)", "Component 2 (12.1%)"
        ax.set_xlabel(f"Component 1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        ax.set_ylabel(f"Component 2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        ax.set_title(f"PCA using Top {n} Differentially Abundant Proteins")

        # frameon=False: 凡例の外枠を非表示にする
        ax.legend(frameon=False)

        # 上辺と右辺の枠線を非表示にして学術論文風のすっきりした見た目にする
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # 薄いグリッド線を追加
        ax.grid(True, alpha=0.3, linestyle="--")

        # PCA プロットを PNG として保存
        filepath = os.path.join(FIG_DIR, f"fig2{panel_pca}_pca_top{n}.png")
        fig.savefig(filepath, dpi=150, bbox_inches="tight")
        plt.close()  # メモリ解放
        print(f"    保存: {filepath}")


# ============================================================================
#  メイン関数: 全処理の実行フロー
# ============================================================================

def main():
    """
    差分発現解析のメインエントリーポイント。

    【何をするか】
      1. 前処理済みデータとサンプル情報を読み込む
      2. Normal 群と Tumor 群のサンプルリストを作成
      3. Welch の t 検定を全タンパク質に対して実行
      4. 結果サマリーを表示し、CSV に保存
      5. Volcano プロットを描画・保存（ボーナス図）
      6. 全有意差タンパク質のヒートマップと PCA を描画（Figure 2a, 2b）
      7. Top N タンパク質で階層的クラスタリングと PCA を実行・保存（Figure 2c-h）

    【引数】
      なし

    【戻り値】
      なし（結果はファイルとして保存され、進捗はコンソールに表示される）
    """
    # 処理開始のヘッダーを表示（"=" x 60 の区切り線で囲む）
    print("=" * 60)
    print("Step 4: 差分発現解析（Figure 2 の再現）")
    print("=" * 60)

    # --- データの読み込み ---
    df, sample_info = load_data()

    # sample_info から条件ごとのサンプル名リストを作成
    # Condition列が "Normal" の行の Sample列の値をリストとして取得
    normal_samples = sample_info[sample_info["Condition"] == "Normal"]["Sample"].tolist()
    # Condition列が "Tumor" の行の Sample列の値をリストとして取得
    tumor_samples = sample_info[sample_info["Condition"] == "Tumor"]["Sample"].tolist()

    # ====================================================================
    #  [1/5] Welch's t-test の実行
    # ====================================================================
    print("\n[1/5] Welch's t-test")
    result_df = welch_ttest(df, normal_samples, tumor_samples)

    # 結果のサマリーを表示（有意差ありの個数を集計）
    n_up = (result_df["Significant"] == "Up").sum()       # 発現上昇タンパク質の個数
    n_down = (result_df["Significant"] == "Down").sum()   # 発現低下タンパク質の個数
    n_total = n_up + n_down                                # 有意差ありの合計
    print(f"  検定タンパク質数: {len(result_df)}")
    print(f"  有意差あり: {n_total} （↑{n_up}, ↓{n_down}）")
    # 論文で報告されている値と比較するための参考情報
    print(f"  （論文値: 2642, ↑1475, ↓1167）")

    # 検定結果テーブルを CSV として保存
    result_df.to_csv(os.path.join(TABLE_DIR, "differential_proteins.csv"), index=False)

    # ====================================================================
    #  [2/5] Volcano プロットの描画（ボーナス図）
    # ====================================================================
    print("\n[2/5] Volcanoプロット（ボーナス図）")
    plot_volcano(result_df)

    # ====================================================================
    #  [3/5] 全有意差タンパク質のヒートマップ（Figure 2a）
    # ====================================================================
    print("\n[3/5] Figure 2a: 全有意差タンパク質のヒートマップ")
    plot_heatmap_all(df, result_df, sample_info)

    # ====================================================================
    #  [4/5] 全有意差タンパク質の PCA（Figure 2b）
    # ====================================================================
    print("\n[4/5] Figure 2b: 全有意差タンパク質の PCA")
    plot_pca_all(df, result_df, sample_info)

    # ====================================================================
    #  [5/5] Top N タンパク質での解析（クラスタリング + PCA）
    # ====================================================================
    print("\n[5/5] Top N タンパク質の解析（Figure 2c-h）")
    plot_top_n_analysis(df, result_df, sample_info)

    # 処理完了メッセージの表示
    print("\n" + "=" * 60)
    print("差分発現解析完了！")
    print("=" * 60)


# ============================================================================
#  スクリプト実行時のエントリーポイント
# ============================================================================
# このファイルが直接実行された場合のみ main() を呼び出す。
# 他のスクリプトから import された場合は main() は実行されない。
if __name__ == "__main__":
    main()
