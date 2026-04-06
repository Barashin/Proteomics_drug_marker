#!/usr/bin/env python3
# ============================================================================
#  Step 6: ステージ別解析（Figure 3 の再現）
# ============================================================================
#
#  このスクリプトは、がんの【臨床ステージ（I〜IV期）】ごとに
#  タンパク質の発現パターンがどう変化するかを統計的に解析し、
#  可視化するためのものです。
#
#  【処理の流れ】
#    1. 前処理済みデータと臨床情報を読み込む
#    2. 各タンパク質について【One-way ANOVA】（一元配置分散分析）を実行
#       → ステージ間で発現量に有意な差があるタンパク質を特定
#    3. 【FDR補正】（偽発見率の制御）により多重比較の問題を解決
#    4. 有意なタンパク質を【階層的クラスタリング】で30群に分類
#    5. 【Z-score正規化】で発現量を標準化し、プロファイルプロットを作成
#
#  【One-way ANOVA（一元配置分散分析）とは？】
#    - 3群以上のグループ間で平均値に差があるかを検定する手法です。
#    - t検定は「2群」の比較しかできませんが、ANOVAは「3群以上」を
#      一度に比較できます。
#    - 例: Normal vs Stage I vs Stage II vs Stage III vs Stage IV の
#      5群を同時に比較する場合、t検定だと10回（5C2通り）の検定が
#      必要ですが、ANOVAなら1回の検定で済みます。
#    - F統計量が大きいほど、グループ間の差が大きいことを意味します。
#
#  【FDR補正（偽発見率の制御）とは？】
#    - 数千個のタンパク質に対して検定を行うと、偶然だけでも
#      「有意」と判定されるものが多数出てしまいます（多重比較問題）。
#    - 例: 10,000個のタンパク質を p < 0.05 で検定すると、
#      差がなくても約500個が偶然「有意」になります。
#    - FDR補正は「有意と判定したもののうち、本当は差がないものの
#      割合（偽発見率）」を制御する手法です。
#    - Benjamini-Hochberg法は最も広く使われるFDR補正法で、
#      p値を小さい順に並べ、順位に応じた閾値で判定します。
#
#  【Z-score正規化とは？】
#    - 各タンパク質の発現値を「平均0、標準偏差1」に変換する手法です。
#    - 計算式: Z = (値 - 平均) / 標準偏差
#    - これにより、発現量のスケールが異なるタンパク質同士を
#      同じ尺度で比較・可視化できるようになります。
#    - Z > 0 は平均より高発現、Z < 0 は平均より低発現を意味します。
#
#  入力ファイル:
#    - results/preprocessed_data.csv   （前処理済み発現データ）
#    - results/sample_info.csv         （サンプル情報）
#    - results/clinical_info.csv       （臨床情報、存在する場合）
#
#  出力ファイル:
#    - results/tables/anova_results.csv         （ANOVA検定結果）
#    - results/tables/cluster_assignments.csv   （クラスター割り当て）
#    - results/figures/fig3_stage_heatmap.png    （Figure 3: ヒートマップ+ラインプロット複合パネル）
# ============================================================================


# ============================================================================
#  ライブラリのインポート
# ============================================================================

import os                       # ファイルパス操作・ディレクトリ作成に使用
import numpy as np              # 数値計算ライブラリ（配列操作、ソートなど）
import pandas as pd             # データフレーム操作（CSV読み書き、データ整形）
import matplotlib.pyplot as plt  # グラフ描画の基本ライブラリ（Figure, Axesの管理）
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # サブプロットのレイアウト制御
from matplotlib.colors import ListedColormap  # カスタムカラーマップ作成
import seaborn as sns           # 統計グラフの高機能描画（ヒートマップなど）

from scipy import stats         # 科学技術計算ライブラリの統計モジュール（ANOVA実行に使用）

from scipy.cluster.hierarchy import (
    linkage,      # 【階層的クラスタリング】の連結行列を計算する関数
    fcluster,     # 連結行列からクラスター割り当てを取得する関数
    dendrogram,   # 樹形図（デンドログラム）を描画する関数（本スクリプトでは未使用だが将来利用可能）
)

from scipy.spatial.distance import pdist  # 観測間のペアワイズ距離を計算（クラスタリングの前段階）

from statsmodels.stats.multitest import multipletests  # 【多重検定補正】（FDR補正など）を実行する関数


# ============================================================================
#  定数・設定
# ============================================================================

# --- ディレクトリパスの設定 ---
# このスクリプトが置かれているディレクトリの絶対パスを取得
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# プロジェクトのルートディレクトリ（scriptsの1つ上の階層）
PROJECT_DIR = os.path.join(SCRIPT_DIR, "..")

# 解析結果を保存するディレクトリ
RESULTS_DIR = os.path.join(PROJECT_DIR, "results")

# 図を保存するディレクトリ
FIG_DIR = os.path.join(RESULTS_DIR, "figures")

# 表（CSV）を保存するディレクトリ
TABLE_DIR = os.path.join(RESULTS_DIR, "tables")

# 出力ディレクトリが存在しない場合は自動的に作成する
# exist_ok=True により、既にディレクトリが存在してもエラーにならない
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

# --- 論文に準拠した解析パラメータ ---

# 【クラスター数】: 階層的クラスタリングで分割するクラスターの数
# 論文ではタンパク質を30個のグループに分類して発現パターンを分析している
N_CLUSTERS = 30

# 【FDR閾値】: この値より小さいFDR値を持つタンパク質を「有意」と判定する
# 0.01 = 1%、つまり有意と判定したタンパク質のうち、偽陽性の割合が1%以下
# 一般的な閾値は 0.05（5%）だが、本論文ではより厳しい 0.01 を採用
FDR_THRESHOLD = 0.01

# 【ステージ別カラーマップ】（論文のFigure 3に準拠した配色）
# がんの進行度に合わせて直感的な色を割り当てている:
#   - Normal（正常組織）: 青 → 健康・正常を連想させる寒色系
#   - Stage I（早期がん）: 緑 → まだ初期段階であることを示す
#   - Stage II（中期がん）: 黄 → 注意レベルの上昇を示す
#   - Stage III（進行がん）: オレンジ → 危険度の上昇を示す
#   - Stage IV（末期がん）: 赤 → 最も深刻な状態を示す暖色系
# この信号機的な配色により、ステージの進行を視覚的に把握しやすくなる
STAGE_COLORS = {
    "Normal": "#4EAED1",   # 青（正常組織）
    "I": "#66BB6A",        # 緑（Stage I: 早期）
    "II": "#FFD54F",       # 黄（Stage II: 中期）
    "III": "#FFA726",      # オレンジ（Stage III: 進行期）
    "IV": "#E8524A",       # 赤（Stage IV: 末期）
}


# ============================================================================
#  データ読み込み関数
# ============================================================================

def load_data():
    """
    前処理済みデータとサンプル情報をCSVファイルから読み込む。

    【この関数の役割】
    Step 5までで作成された前処理済みデータ（欠損値処理・正規化済み）を
    読み込み、以降の解析で使えるDataFrame形式にする。

    【戻り値】
    - df : pd.DataFrame
        前処理済みのタンパク質発現データ。
        行 = タンパク質名、列 = サンプル名、値 = 発現量。
    - sample_info : pd.DataFrame
        サンプルの属性情報（サンプル名、条件[Normal/Tumor]など）。
    """
    # index_col=0: 最初の列（タンパク質名）をインデックスとして読み込む
    df = pd.read_csv(os.path.join(RESULTS_DIR, "preprocessed_data.csv"), index_col=0)

    # サンプル情報（各サンプルがNormalかTumorか、どの患者由来かなど）
    sample_info = pd.read_csv(os.path.join(RESULTS_DIR, "sample_info.csv"))

    return df, sample_info


# ============================================================================
#  ステージ情報の付与
# ============================================================================

def assign_stages(sample_info):
    """
    各サンプルに【臨床ステージ】情報を付与する。

    【この関数の役割】
    がんの臨床ステージ（I〜IV期）はがんの進行度を示す指標で、
    ステージが上がるほど進行していることを意味します。
    この関数は、サンプル情報にステージ列を追加します。

    【ステージの意味】
    - Normal : 正常組織（がんではない）
    - Stage I : 早期がん（腫瘍が小さく、浸潤が限定的）
    - Stage II : 中期がん（腫瘍がやや大きい、またはリンパ節転移あり）
    - Stage III : 進行がん（周囲組織への浸潤が進んでいる）
    - Stage IV : 末期がん（遠隔転移あり）

    【処理の流れ】
    1. clinical_info.csv が存在する場合:
       → そこからステージ情報を取得し、各サンプルに割り当てる
    2. clinical_info.csv が存在しない場合:
       → Condition列（Normal/Tumor）をそのままステージとして使用

    【引数】
    - sample_info : pd.DataFrame
        サンプル情報のデータフレーム（Sample列、Condition列を含む）

    【戻り値】
    - sample_info : pd.DataFrame
        Stage列が追加されたサンプル情報のデータフレーム
    """
    # 臨床情報ファイルのパスを構築
    clinical_path = os.path.join(RESULTS_DIR, "clinical_info.csv")

    if os.path.exists(clinical_path):
        # 臨床情報ファイルが存在する場合: 患者ごとのステージ情報を利用
        clinical = pd.read_csv(clinical_path)

        # サンプル名 → ステージのマッピング辞書を作成
        stage_map = {}
        for _, row in clinical.iterrows():
            # 正常組織サンプル（Sample_N）には "Normal" を割り当て
            stage_map[row["Sample_N"]] = "Normal"
            # 腫瘍サンプル（Sample_T）には臨床ステージ（I, II, III, IV）を割り当て
            stage_map[row["Sample_T"]] = row["Stage"]

        # マッピング辞書を使って、sample_infoにStage列を追加
        # map()はSample列の各値をstage_mapで変換する
        # _dupN サフィックス（重複Gene Symbol処理で付与）を除去してからマッピング
        import re
        sample_info["Stage"] = sample_info["Sample"].apply(
            lambda x: stage_map.get(re.sub(r"_dup\d+$", "", x))
        )
    else:
        # 臨床情報ファイルが存在しない場合:
        # Condition列（"Normal" または "Tumor"）をステージの代わりに使用
        # この場合、ステージ別の詳細な解析はできないが、
        # Normal vs Tumor の2群比較は可能
        sample_info["Stage"] = sample_info["Condition"]

    return sample_info


# ============================================================================
#  One-way ANOVA（一元配置分散分析）
# ============================================================================

def run_anova(df, sample_info):
    """
    各タンパク質についてOne-way ANOVA（一元配置分散分析）を実行する。

    【One-way ANOVAの仕組み（初心者向け解説）】
    ANOVAは「グループ間のばらつき」と「グループ内のばらつき」を比較します。

    ・【グループ間のばらつき】: ステージごとの平均値のばらつき
      → 大きいほど「ステージ間で発現量が異なる」証拠になる
    ・【グループ内のばらつき】: 同じステージ内のサンプル間のばらつき
      → 個体差による自然なばらつき

    F統計量 = グループ間のばらつき / グループ内のばらつき
      → F値が大きい = グループ間の差が個体差より大きい = 有意な差がある可能性

    【ANOVAとt検定の違い】
    ・t検定: 2群の平均値の差を検定（例: Normal vs Tumor）
    ・ANOVA: 3群以上の平均値の差を一度に検定
      （例: Normal vs Stage I vs Stage II vs Stage III vs Stage IV）
    ・3群以上でt検定を繰り返すと多重比較の問題が生じるため、ANOVAを使う

    【この関数の処理の流れ】
    1. 各タンパク質について、ステージごとの発現値を抽出
    2. scipy.stats.f_oneway() でF検定を実行し、F統計量とp値を取得
    3. 全タンパク質のp値に対してBenjamini-Hochberg法でFDR補正を適用
    4. FDR < 0.01 のタンパク質を「有意」と判定

    【引数】
    - df : pd.DataFrame
        タンパク質発現データ（行=タンパク質、列=サンプル）
    - sample_info : pd.DataFrame
        サンプル情報（Stage列を含む）

    【戻り値】
    - result_df : pd.DataFrame
        各タンパク質のANOVA結果。以下の列を含む:
        - Protein: タンパク質名
        - F_statistic: F統計量（大きいほどグループ間差が大きい）
        - P_value: p値（小さいほど有意）
        - FDR: FDR補正後のp値（多重比較を考慮した値）
        - Significant: FDR < 閾値 かどうかのブール値
    """
    # ステージを生物学的順序で並べる
    stage_order = ["Normal", "I", "II", "III", "IV"]
    available = [s for s in sample_info["Stage"].unique() if pd.notna(s)]
    stages_with_normal = [s for s in stage_order if s in available]

    # 各タンパク質のANOVA結果を格納するリスト
    results = []

    # --- 全タンパク質に対してループ処理でANOVAを実行 ---
    for protein in df.index:
        # 各ステージの発現値をリストとして収集
        groups = []
        for stage in stages_with_normal:
            # 現在のステージに属するサンプル名を取得
            samples = sample_info[sample_info["Stage"] == stage]["Sample"].tolist()

            # 発現データに実際に存在するサンプルのみを抽出し、欠損値を除外
            vals = df.loc[protein, [s for s in samples if s in df.columns]].dropna()

            # サンプル数が2以上の場合のみグループに追加
            # （分散を計算するには最低2つのデータ点が必要）
            if len(vals) >= 2:
                groups.append(vals.values)

        # 比較可能なグループが2つ未満の場合はスキップ
        # （ANOVAには最低2グループが必要）
        if len(groups) < 2:
            continue

        # 【One-way ANOVAの実行】
        # *groups でリストを展開し、各ステージのデータを引数として渡す
        # 戻り値: F統計量（f_stat）とp値（p_val）
        # p値が小さいほど「ステージ間で発現量に差がある」可能性が高い
        f_stat, p_val = stats.f_oneway(*groups)

        # 結果を辞書としてリストに追加
        results.append({
            "Protein": protein,        # タンパク質名
            "F_statistic": f_stat,     # F統計量
            "P_value": p_val,          # p値（未補正）
        })

    # 結果リストをDataFrameに変換
    result_df = pd.DataFrame(results)

    # --- 【FDR補正】（Benjamini-Hochberg法）---
    # 多重比較問題への対処:
    # 数千個のタンパク質を同時に検定すると、偶然「有意」になるものが多発する。
    # Benjamini-Hochberg法は、p値を小さい順にランク付けし、
    # 各p値に (総検定数 / ランク) を掛けて調整済みp値（FDR）を計算する。
    # これにより、「有意と判定した中の偽陽性の割合」を制御できる。
    if len(result_df) > 0:
        # multipletests関数の戻り値:
        #   reject: 帰無仮説を棄却するかどうかのブール配列
        #   fdr: FDR補正後のp値（調整済みp値）
        #   _: Sidak補正の値（ここでは不使用）
        #   _: Bonferroni補正の値（ここでは不使用）
        reject, fdr, _, _ = multipletests(result_df["P_value"], method="fdr_bh")

        # FDR補正後のp値を列として追加
        result_df["FDR"] = fdr

        # FDR < 閾値（0.01）なら「有意」と判定
        result_df["Significant"] = fdr < FDR_THRESHOLD

    return result_df


# ============================================================================
#  ステージ別中央値の計算
# ============================================================================

def compute_stage_medians(df, sample_info):
    """
    各ステージごとにタンパク質発現量の【中央値】を計算する。

    【この関数の役割】
    各ステージに属する複数のサンプルの発現量を1つの代表値（中央値）に
    集約することで、ステージ間の発現パターンを比較しやすくします。

    【なぜ中央値を使うのか？】
    - 中央値は外れ値の影響を受けにくい（ロバスト）
    - 平均値は極端に大きい/小さい値に引きずられやすいが、
      中央値はデータの真ん中の値なので安定している
    - プロテオミクスデータには外れ値が含まれることが多いため、
      中央値が代表値として適している

    【引数】
    - df : pd.DataFrame
        タンパク質発現データ（行=タンパク質、列=サンプル）
    - sample_info : pd.DataFrame
        サンプル情報（Stage列を含む）

    【戻り値】
    - pd.DataFrame
        行 = タンパク質名、列 = ステージ名、値 = 中央値
        例: Normal, I, II, III, IV の5列を持つデータフレーム
    """
    # 【重要】ステージを生物学的順序で並べる（Non-tumor → Stage I → II → III → IV）
    # sorted() だとアルファベット順（I, II, III, IV, Normal）になり、
    # Normalが最後に来てしまう。これだとZ-score計算やクラスタリングの結果が変わる。
    stage_order = ["Normal", "I", "II", "III", "IV"]
    available = [s for s in sample_info["Stage"].unique() if pd.notna(s)]
    stages = [s for s in stage_order if s in available]

    # ステージ名 → 中央値Series のマッピング辞書
    medians = {}

    for stage in stages:
        # 現在のステージに属するサンプル名を取得
        samples = sample_info[sample_info["Stage"] == stage]["Sample"].tolist()

        # 実際に発現データに存在するサンプルのみを抽出
        valid_samples = [s for s in samples if s in df.columns]

        if valid_samples:
            # axis=1: 行方向（=各タンパク質について）中央値を計算
            # → 各タンパク質の、このステージにおける代表発現値
            medians[stage] = df[valid_samples].median(axis=1)

    # 辞書からDataFrameを作成（列名がステージ名になる）
    return pd.DataFrame(medians)


# ============================================================================
#  階層的クラスタリングとプロファイルプロット用データの作成
# ============================================================================

def cluster_and_plot(median_df, anova_df, n_clusters=N_CLUSTERS):
    """
    ANOVAで有意と判定されたタンパク質を【階層的クラスタリング】で分類する。

    【階層的クラスタリングとは？（初心者向け解説）】
    - 似た発現パターンを持つタンパク質同士をグループ（クラスター）にまとめる手法
    - 「ボトムアップ」方式: 最初は各タンパク質が1つのクラスター、
      最も似ているペアを順次結合していく
    - Ward法: クラスター内の分散が最小になるように結合する方法
      → コンパクトで均一なクラスターが作られやすい

    【Z-score正規化の必要性】
    - クラスタリングの前に、各タンパク質の発現値をZ-scoreに変換する
    - これにより「発現量の絶対値」ではなく「発現パターンの形」で
      タンパク質を分類できる
    - 例: 発現量が100→200に増加するタンパク質と、
      発現量が1→2に増加するタンパク質は、Z-scoreにすると同じパターン

    【引数】
    - median_df : pd.DataFrame
        ステージ別中央値データ（行=タンパク質、列=ステージ）
    - anova_df : pd.DataFrame
        ANOVA結果（Significant列を含む）
    - n_clusters : int
        分割するクラスター数（デフォルト: 30）

    【戻り値】
    - タプル (z_data, clusters, actual_clusters) または None
        - z_data : pd.DataFrame - Z-score正規化後のデータ
        - clusters : np.ndarray - 各タンパク質のクラスター番号
        - actual_clusters : int - 実際のクラスター数
    """
    # --- 有意なタンパク質の抽出 ---
    # FDR < 0.01 を満たすタンパク質のみを対象とする
    sig_proteins = anova_df[anova_df["Significant"]]["Protein"].tolist()

    if not sig_proteins:
        # 有意なタンパク質が0個の場合: FDR閾値を自動的に緩和
        # FDR値が小さい上位200個（またはデータ数が200未満ならその全数）を使用
        print("  注意: FDR < 0.01 のタンパク質がありません。閾値を緩めます。")
        sig_proteins = anova_df.nsmallest(min(200, len(anova_df)), "FDR")["Protein"].tolist()

    # 中央値データから有意なタンパク質のみを抽出
    sig_data = median_df.loc[median_df.index.isin(sig_proteins)]

    # --- 【Z-score正規化】（行方向 = 各タンパク質ごとに正規化）---
    # 計算式: Z = (x - mean) / std
    # axis=1: 行方向（各タンパク質の全ステージにわたる値で平均・標準偏差を計算）
    # lambda関数で各行に対して「(値 - 行平均) / 行標準偏差」を適用
    # → 各タンパク質の発現パターンが平均0、標準偏差1に標準化される
    z_data = sig_data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

    # 標準偏差が0（= 全ステージで同じ値）のタンパク質はNaNになるので除外
    z_data = z_data.dropna()

    if len(z_data) == 0:
        # Z-score変換後に有効なデータがない場合はスキップ
        print("  スキップ: 有効なデータがありません")
        return None

    # --- 【階層的クラスタリング】の実行 ---
    # クラスター数をデータ数で制限（データ数が30未満ならデータ数をクラスター数とする）
    actual_clusters = min(n_clusters, len(z_data))

    # linkage関数: 階層的クラスタリングの連結行列Zを計算
    # method="ward": Ward法（クラスター内分散の増加を最小化する方法）
    # z_data.values: Z-score正規化後の数値行列を入力
    # 出力のZ行列は (n-1) x 4 の行列で、各行が1つの結合ステップを表す:
    #   [結合するクラスター1, クラスター2, 距離, 新クラスターのサンプル数]
    Z = linkage(z_data.values, method="ward")

    # fcluster関数: 連結行列から指定数のクラスターに分割
    # t=actual_clusters: 目標クラスター数
    # criterion="maxclust": 最大クラスター数で切断する方法
    # 戻り値: 各タンパク質のクラスター番号（1〜actual_clusters）の配列
    clusters = fcluster(Z, t=actual_clusters, criterion="maxclust")

    # --- クラスター割り当て結果をCSVに保存 ---
    cluster_df = pd.DataFrame({
        "Protein": z_data.index,     # タンパク質名
        "Cluster": clusters,         # 割り当てられたクラスター番号
    })
    cluster_df.to_csv(os.path.join(TABLE_DIR, "cluster_assignments.csv"), index=False)

    return z_data, clusters, actual_clusters


# ============================================================================
#  Figure 3: ヒートマップ + ラインプロットの複合パネル描画
# ============================================================================

def load_paper_clusters(raw_dir):
    """
    論文の補足データ（Table S9-S12）からクラスター割り当てを読み込む関数。

    【論文のクラスター】
      - Cluster 3: 34 タンパク質（ステージ進行で一貫して増加）
      - Cluster 14: 1324 タンパク質（ステージ進行で増加）
      - Cluster 20: 16 タンパク質（ステージ進行で一貫して減少）
      - Cluster 25: 1062 タンパク質（ステージ進行で減少）

    戻り値: list of dict（各クラスターの情報）
    """
    cluster_files = [
        ("Supplementary Table S9 Fig3A_Clustre3_34 protein.xlsx", 3, 34, "Increased"),
        ("Supplementary Table S10 Fig3B_Clustre14_1324 protein.xlsx", 14, 1324, "Increased"),
        ("Supplementary Table S11 Fig3C_Clustre20_16 protein.xlsx", 20, 16, "Decreased"),
        ("Supplementary Table S12 Fig3D_Clustre25_1062 protein.xlsx", 25, 1062, "Decreased"),
    ]

    paper_clusters = []
    for fname, cl_num, expected_n, direction in cluster_files:
        filepath = os.path.join(raw_dir, fname)
        if not os.path.exists(filepath):
            print(f"  注意: {fname} が見つかりません。スキップします。")
            continue
        cluster_df = pd.read_excel(filepath, header=0)
        gene_col = [c for c in cluster_df.columns if "gene" in c.lower() or "symbol" in c.lower()]
        if gene_col:
            genes = cluster_df[gene_col[0]].dropna().tolist()
        else:
            genes = cluster_df.iloc[:, 5].dropna().tolist()
        paper_clusters.append({
            "cluster": cl_num,
            "n": len(genes),
            "genes": genes,
            "direction": direction,
        })
        print(f"  Cluster {cl_num}: {len(genes)} タンパク質 ({direction})")

    return paper_clusters


def plot_figure3(df, sample_info, z_data, clusters):
    """
    Toyota et al. 2025 の Figure 3 を再現する。
    論文の補足データ（Table S9-S12）のクラスターを使用して完全再現。

    【Figure 3 の構成】
    2x2 の4パネル構成で、各パネルは以下の2つの要素で構成される:
      - 上段: ヒートマップ（サンプル×タンパク質のZ-score発現量マップ）
      - 下段: ラインプロット（ステージ別の中央値プロファイル）

    パネル (a), (b): 疾患進行とともに発現が【増加】するクラスター
    パネル (c), (d): 疾患進行とともに発現が【減少】するクラスター

    【ヒートマップの詳細】
    - X軸 = サンプル（ステージ順に並べ替え: Normal → I → II → III → IV）
    - Y軸 = タンパク質（階層的クラスタリングで並べ替え）
    - 色: RdBu_r カラーマップ（赤 = 高発現、青 = 低発現）
    - 上部に2つのカラーバーを表示:
      - "Histology": 青 = Non-tumor、赤 = Tumor
      - "Clinical stage": 緑(I) → 黄(II) → オレンジ(III) → 赤(IV)

    【ラインプロットの詳細】
    - X軸: 臨床ステージ（Normal, I, II, III, IV）
    - Y軸: 平均タンパク質発現量（Z-score中央値の中央値）
    - 各ステージに対応する色のドットを表示
    - クラスター全体の中央値トレンドを1本の折れ線で表示

    【クラスター選択の基準】
    - 全クラスターのうち、タンパク質数が5以上のものを対象
    - Normal vs Stage IV のZ-score中央値の差（trend）を算出
    - trend が最も大きい2クラスター（増加）と最も小さい2クラスター（減少）を選択

    【引数】
    - df : pd.DataFrame
        前処理済みのタンパク質発現データ（行=タンパク質、列=サンプル）
    - sample_info : pd.DataFrame
        サンプル情報（Sample列、Condition列、Stage列を含む）
    - z_data : pd.DataFrame
        Z-score正規化後のステージ別中央値データ（行=タンパク質、列=ステージ）
    - clusters : np.ndarray
        各タンパク質のクラスター番号
    """
    print("  Figure 3: ヒートマップ + ラインプロットの複合パネル")

    # === 1. 論文のクラスターデータを読み込む ===
    # 自前でクラスタリングするのではなく、論文の補足データ（Table S9-S12）の
    # クラスター割り当てをそのまま使用して完全再現する。
    RAW_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "raw")
    paper_clusters = load_paper_clusters(RAW_DIR)

    if len(paper_clusters) == 4:
        # 論文のクラスターデータが4つ揃っている場合: そのまま使用
        # 順序: (a) Cluster 3, (b) Cluster 14, (c) Cluster 20, (d) Cluster 25
        selected = []
        for pc in paper_clusters:
            # 論文のGene Symbolリストから、前処理済みデータに存在するものを抽出
            matched_genes = [g for g in pc["genes"] if g in df.index]
            selected.append({
                "cluster": pc["cluster"],
                "n": len(matched_genes),
                "genes": matched_genes,
                "direction": pc["direction"],
            })
            print(f"    Cluster {pc['cluster']}: {len(matched_genes)}/{len(pc['genes'])} マッチ")
    else:
        # 補足データがない場合: 自前クラスタリング結果からフォールバック
        print("  注意: 論文のクラスターデータが不完全。自前クラスタリングを使用。")
        unique_clusters = np.unique(clusters)
        cluster_info = []
        for cl in unique_clusters:
            mask = clusters == cl
            n_proteins = mask.sum()
            if n_proteins < 5:
                continue
            profile = z_data.iloc[mask].mean(axis=0)
            normal_val = profile.get("Normal", profile.iloc[0])
            last_val = profile.get("IV", profile.iloc[-1])
            trend = last_val - normal_val
            genes = z_data.index[mask].tolist()
            direction = "Increased" if trend > 0 else "Decreased"
            cluster_info.append({"cluster": cl, "n": n_proteins, "genes": genes, "direction": direction, "trend": trend})

        inc = sorted([c for c in cluster_info if c["direction"] == "Increased"], key=lambda x: x["trend"], reverse=True)
        dec = sorted([c for c in cluster_info if c["direction"] == "Decreased"], key=lambda x: x["trend"])
        selected = (inc[:2] if len(inc) >= 2 else inc) + (dec[:2] if len(dec) >= 2 else dec)
    panel_labels = ["(a)", "(b)", "(c)", "(d)"]

    # === 2. サンプル順序の決定 ===
    # ステージ順にサンプルを並べ替える: Normal → I → II → III → IV
    stage_order = ["Normal", "I", "II", "III", "IV"]
    ordered_samples = []
    sample_stage_map = {}  # サンプル名 → ステージのマッピング

    for stage in stage_order:
        # 各ステージに属するサンプル名を取得
        stage_samples = sample_info[sample_info["Stage"] == stage]["Sample"].tolist()
        # 実際のデータに存在するサンプルのみ抽出
        valid = [s for s in stage_samples if s in df.columns]
        ordered_samples.extend(valid)
        for s in valid:
            sample_stage_map[s] = stage

    # === 3. サンプルごとのZ-score計算（全サンプルレベル）===
    # ヒートマップには個別サンプルのデータが必要なため、
    # dfから対象タンパク質のデータを抽出し、行方向にZ-score正規化する
    # （z_dataはステージ別中央値のZ-scoreなので、サンプルレベルのデータが別途必要）

    # === 4. カラーバー用のデータ作成 ===
    # Histology（組織型）: Non-tumor（青）/ Tumor（赤）
    histology_colors_map = {"Normal": "#4EAED1", "Tumor": "#E8524A"}

    # === 5. Figure作成 ===
    # 2x2 のパネル配置、各パネルが上下2段（ヒートマップ + ラインプロット）
    fig = plt.figure(figsize=(16, 20))
    outer_gs = GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.3)

    for idx, info in enumerate(selected):
        # --- 各パネルの内部レイアウト ---
        # 上段: ヒートマップ（height_ratio=3）、下段: ラインプロット（height_ratio=1）
        inner_gs = GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer_gs[idx],
            height_ratios=[3, 1], hspace=0.15
        )
        ax_heat = fig.add_subplot(inner_gs[0])
        ax_line = fig.add_subplot(inner_gs[1])

        # --- クラスターデータの抽出 ---
        cl = info["cluster"]
        if "genes" in info and info["genes"]:
            cluster_proteins = [g for g in info["genes"] if g in df.index]
        else:
            mask = clusters == cl
            cluster_proteins = z_data.index[mask].tolist()
            cluster_proteins = [p for p in cluster_proteins if p in df.index]

        # 【論文準拠】各ステージのグループ中央値を計算（5列のヒートマップ）
        # 論文Methods: "the median value for each group was calculated
        # and converted to a normalized z-score"
        stage_median_data = {}
        for stage in stage_order:
            stage_samples = sample_info[sample_info["Stage"] == stage]["Sample"].tolist()
            valid_s = [s for s in stage_samples if s in df.columns]
            if valid_s:
                stage_median_data[stage] = df.loc[cluster_proteins, valid_s].median(axis=1)

        # ステージ中央値のDataFrame（行=タンパク質、列=5ステージ）
        median_matrix = pd.DataFrame(stage_median_data)

        # 行方向（タンパク質ごと）にZ-score正規化
        cluster_z = median_matrix.apply(
            lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x * 0, axis=1
        )

        # タンパク質の階層的クラスタリングによる並べ替え
        if len(cluster_z) > 1:
            from scipy.cluster.hierarchy import leaves_list
            Z_linkage = linkage(cluster_z.values, method="ward")
            leaf_order = leaves_list(Z_linkage)
            cluster_z = cluster_z.iloc[leaf_order]

        # --- ヒートマップの描画（5列 = ステージ中央値） ---
        sns.heatmap(
            cluster_z,
            cmap="RdYlGn_r",             # 【論文準拠カラーマップ】赤=高, 黄=中, 緑=低
            center=0,
            vmin=-2, vmax=2,
            ax=ax_heat,
            xticklabels=True,             # 5列なのでステージ名を表示
            yticklabels=False,
            cbar_kws={"shrink": 0.5, "label": "Z-score"},
        )
        ax_heat.set_xticklabels(
            ["Non-tumor" if s == "Normal" else f"Stage {s}" for s in cluster_z.columns],
            fontsize=7, rotation=45, ha="right"
        )

        # --- 上部カラーバー（1段: Histology/Clinical stage 統合） ---
        heat_pos = ax_heat.get_position()
        bar_height = heat_pos.height * 0.04
        ax_bar = fig.add_axes([
            heat_pos.x0, heat_pos.y1 + 0.005,
            heat_pos.width * 0.85, bar_height,
        ])
        stage_colors_list = [STAGE_COLORS.get(s, "#888") for s in cluster_z.columns]
        for i, color in enumerate(stage_colors_list):
            ax_bar.axvspan(i, i + 1, facecolor=color, edgecolor="white", linewidth=0.5)
        ax_bar.set_xlim(0, len(stage_colors_list))
        ax_bar.set_xticks([])
        ax_bar.set_yticks([])

        # --- タイトル ---
        trend_label = info.get("direction", "Increased" if info.get("trend", 0) > 0 else "Decreased")
        ax_heat.set_title(
            f"{panel_labels[idx]} Cluster {cl}: {info['n']} proteins ({trend_label})",
            fontsize=11, fontweight="bold", pad=15,
        )

        # --- ラインプロット（スパゲッティプロット） ---
        # 論文と同じく、個々のタンパク質の線を薄く描画し、平均線を太く描画
        stage_labels_display = ["Non-tumor" if s == "Normal" else s for s in cluster_z.columns]
        x_pos = range(len(stage_labels_display))

        # 個々のタンパク質のプロファイル（薄い灰色の線）
        for _, row in cluster_z.iterrows():
            ax_line.plot(x_pos, row.values, color="gray", alpha=0.15, linewidth=0.3)

        # 全タンパク質の平均プロファイル（太い黒線）
        mean_profile = cluster_z.mean(axis=0)
        ax_line.plot(x_pos, mean_profile.values, "k-", linewidth=2.5)

        # ステージ色のドット
        for i, (xp, ym) in enumerate(zip(x_pos, mean_profile.values)):
            ax_line.plot(xp, ym, "o", color=stage_colors_list[i],
                         markersize=8, markeredgecolor="black", markeredgewidth=0.5, zorder=5)

        ax_line.set_xticks(list(x_pos))
        ax_line.set_xticklabels(stage_labels_display, fontsize=8)
        ax_line.set_ylabel("Z-score", fontsize=8)
        ax_line.axhline(0, color="gray", linestyle="--", linewidth=0.5, alpha=0.5)

        # Z-score = 0 の水平線を破線で描画（基準線）
        ax_line.axhline(0, color="gray", linestyle="--", linewidth=0.5)

        # 上辺と右辺の枠線を非表示にする（すっきりした見た目）
        ax_line.spines["top"].set_visible(False)
        ax_line.spines["right"].set_visible(False)

    # --- 凡例の追加 ---
    # Figure下部にHistologyとClinical stageの凡例を配置
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=histology_colors_map["Normal"], label="Non-tumor"),
        Patch(facecolor=histology_colors_map["Tumor"], label="Tumor"),
        Patch(facecolor="white", edgecolor="white", label="   "),  # スペーサー
        Patch(facecolor=STAGE_COLORS["I"], label="Stage I"),
        Patch(facecolor=STAGE_COLORS["II"], label="Stage II"),
        Patch(facecolor=STAGE_COLORS["III"], label="Stage III"),
        Patch(facecolor=STAGE_COLORS["IV"], label="Stage IV"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=7,
        fontsize=9,
        frameon=False,
        bbox_to_anchor=(0.5, -0.01),
    )

    # 図をPNGファイルとして保存
    filepath = os.path.join(FIG_DIR, "fig3_stage_heatmap.png")
    fig.savefig(
        filepath,
        dpi=150,                            # 解像度: 150 DPI
        bbox_inches="tight"                 # 余白を最小限にして保存
    )
    plt.close()                             # メモリ節約のため図を閉じる
    print(f"    保存: {filepath}")


# ============================================================================
#  メイン処理
# ============================================================================

def main():
    """
    ステージ別解析のメイン処理を実行する。

    【処理の全体フロー】
    1. データ読み込み → 前処理済みデータとサンプル情報を取得
    2. ステージ情報の付与 → 各サンプルにNormal/I/II/III/IVを割り当て
    3. ステージ別中央値の計算 → 各ステージの代表値を算出
    4. One-way ANOVA → ステージ間で有意に変動するタンパク質を特定
    5. 階層的クラスタリング → 有意なタンパク質を30群に分類
    6. 可視化 → ヒートマップとプロファイルプロットを作成

    【引数】 なし
    【戻り値】 なし（結果はファイルに保存される）
    """
    # === ヘッダー表示 ===
    print("=" * 60)
    print("Step 6: ステージ別解析（Figure 3 の再現）")
    print("=" * 60)

    # === データ読み込みとステージ情報の付与 ===
    df, sample_info = load_data()         # 前処理済みデータとサンプル情報を読み込み
    sample_info = assign_stages(sample_info)  # 各サンプルにステージ情報を付与

    # データの概要を表示（デバッグ・確認用）
    print(f"\nデータ: {df.shape[0]} タンパク質 × {df.shape[1]} サンプル")
    print(f"ステージ分布:")
    # value_counts(): 各ステージのサンプル数をカウントして表示
    print(sample_info["Stage"].value_counts().to_string())

    # === 1. ステージ別中央値の計算 ===
    print("\n[1/4] ステージ別中央値の計算")
    median_df = compute_stage_medians(df, sample_info)
    print(f"  ステージ: {list(median_df.columns)}")

    # === 2. One-way ANOVA（一元配置分散分析）===
    print("\n[2/4] One-way ANOVA")
    anova_df = run_anova(df, sample_info)

    # 有意なタンパク質の数を取得（Significant列が存在する場合のみ）
    n_sig = anova_df["Significant"].sum() if "Significant" in anova_df.columns else 0
    print(f"  検定タンパク質数: {len(anova_df)}")
    print(f"  有意（FDR < {FDR_THRESHOLD}）: {n_sig}")

    # ANOVA結果をCSVファイルに保存
    anova_df.to_csv(os.path.join(TABLE_DIR, "anova_results.csv"), index=False)

    # === 3. 階層的クラスタリング ===
    print("\n[3/4] 階層的クラスタリング")
    result = cluster_and_plot(median_df, anova_df)

    if result is not None:
        # クラスタリングが成功した場合
        z_data, clusters, n_cl = result  # Z-scoreデータ、クラスター番号、クラスター数
        print(f"  {n_cl} クラスターに分割")

        # === 4. 可視化（Figure 3: ヒートマップ + ラインプロットの複合パネル）===
        print("\n[4/4] 可視化")
        plot_figure3(df, sample_info, z_data, clusters)  # Figure 3 の再現

    # === フッター表示 ===
    print("\n" + "=" * 60)
    print("ステージ別解析完了！")
    print("=" * 60)


# ============================================================================
#  スクリプト実行のエントリーポイント
# ============================================================================
# このファイルが直接実行された場合（python step_06_stage_analysis.py）にのみ
# main()を実行する。他のスクリプトからimportされた場合は実行されない。
if __name__ == "__main__":
    main()
