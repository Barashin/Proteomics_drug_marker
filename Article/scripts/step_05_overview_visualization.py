#!/usr/bin/env python3
"""
Step 3: 全体像の可視化（Figure 1 の再現）
========================================
前処理済みのプロテオミクスデータを使い、データ全体の傾向を可視化するスクリプトです。
論文 Toyota et al. 2025 の Figure 1 に相当する3種類の図を作成します。

【作成する図】
  - (a) 相関行列ヒートマップ（サンプル間の類似度を色で表現）
  - (b) 教師なし階層的クラスタリング（似たサンプルを自動でグループ化）
  - (c) PCA（主成分分析）（高次元データを2次元に圧縮して全体像を把握）

【なぜ全体像の可視化が重要か？】
  - 差次発現解析（Step 4）の前に、データの品質やバッチ効果を確認できる
  - Normal群とTumor群が自然に分離するかどうかを目視で確認できる
  - 外れ値サンプル（異常なサンプル）を早期発見できる
  - 論文のレビュアーに対して「データの妥当性」を示す根拠になる

【3つの手法の違い】
  - 相関行列: サンプル間の「発現パターンの類似度」を見る（ピアソン相関）
  - クラスタリング: 似たサンプルを自動でツリー状に分類する（教師なし＝ラベル不使用）
  - PCA: 数千次元のデータを2次元に要約し、群の分離を確認する

入力: results/preprocessed_data.csv, results/sample_info.csv
出力: results/figures/fig1a_correlation.png
      results/figures/fig1b_clustering.png
      results/figures/fig1c_pca.png
"""

# ============================================================
# ライブラリの読み込み
# ============================================================
import os                                    # ファイルパス操作・ディレクトリ作成に使用
import numpy as np                           # 数値計算ライブラリ（累積和の計算等に使用）
import pandas as pd                          # データフレーム操作ライブラリ（CSV読み込み・表形式データの操作）
import matplotlib.pyplot as plt              # グラフ描画の基盤ライブラリ（図の作成・保存に使用）
import seaborn as sns                        # 統計データ可視化ライブラリ（ヒートマップ・クラスタマップの作成に使用）
from scipy.cluster.hierarchy import linkage, dendrogram  # 階層的クラスタリングの計算とデンドログラム描画（scipy: 科学計算ライブラリ）
from scipy.spatial.distance import pdist     # サンプル間の距離行列を計算する関数（pdist = pairwise distance）
from sklearn.decomposition import PCA        # 主成分分析（PCA）の実装（scikit-learn: 機械学習ライブラリ）
from matplotlib.patches import Ellipse       # PCAプロットに95%信頼楕円を描画するために使用
import matplotlib.transforms as transforms   # 楕円のアフィン変換（回転・拡大・平行移動）に使用


# ============================================================
# 設定（パスとカラーの定義）
# ============================================================

# --- ディレクトリパスの設定 ---
# __file__ はこのスクリプト自身のパスを指す
# os.path.abspath() で絶対パスに変換し、os.path.dirname() でディレクトリ部分を取得
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ".." は1つ上のディレクトリを意味する（scripts/ → Article/）
PROJECT_DIR = os.path.join(SCRIPT_DIR, "..")

# 前処理済みデータの保存先ディレクトリ
RESULTS_DIR = os.path.join(PROJECT_DIR, "results")

# 図の保存先ディレクトリ
FIG_DIR = os.path.join(RESULTS_DIR, "figures")

# ディレクトリが存在しなければ作成（exist_ok=True で既に存在してもエラーにならない）
os.makedirs(FIG_DIR, exist_ok=True)

# --- カラー設定（論文準拠） ---
# 論文の図と同じ色を使用して、再現性を高める
COLOR_NORMAL = "#4EAED1"    # 青色: Non-tumor（正常組織）を示す
COLOR_TUMOR = "#E8524A"     # 赤色: Tumor（腫瘍組織）を示す


# ============================================================
# データ読み込み関数
# ============================================================
def load_data():
    """
    前処理済みデータとサンプル情報を読み込む関数。

    【読み込むファイル】
      - preprocessed_data.csv: Step 2 で前処理した タンパク質（行）× サンプル（列）の定量マトリクス
        - 値は Log2 変換済み・欠損値補完済み
      - sample_info.csv: 各サンプルの群情報（Normal / Tumor）

    引数:
        なし

    戻り値:
        df (pd.DataFrame): タンパク質定量データ（行=タンパク質, 列=サンプル）
        sample_info (pd.DataFrame): サンプル情報（Sample列, Condition列）
    """
    # index_col=0: 最初の列（タンパク質名）をインデックス（行名）として読み込む
    df = pd.read_csv(os.path.join(RESULTS_DIR, "preprocessed_data.csv"), index_col=0)

    # サンプル情報を読み込む（Sample列とCondition列を含む）
    sample_info = pd.read_csv(os.path.join(RESULTS_DIR, "sample_info.csv"))

    return df, sample_info


# ============================================================
# (a) 相関行列ヒートマップの作成関数
# ============================================================
def plot_correlation_matrix(df, sample_info):
    """
    Figure 1a: サンプル間の相関行列をヒートマップとして可視化する関数。

    【相関行列とは？】
      各サンプル同士の「発現パターンの類似度」をピアソン相関係数で数値化したもの。
      - 相関係数 = 1.0: 完全に同じパターン（完全一致）
      - 相関係数 = 0.0: 全く関係がない
      - 相関係数 = -1.0: 完全に逆のパターン
      NormalとTumorが異なる発現パターンを持つなら、群間の相関が群内より低くなる。

    【ピアソン相関係数】
      2つのサンプルの発現値ベクトル間の線形相関を測る指標。
      値が1に近いほど発現パターンが似ていることを意味する。

    【クラスタマップ（clustermap）とは？】
      ヒートマップに加えて、行・列を自動的に並び替え（クラスタリング）する。
      似たサンプル同士が隣り合うように並ぶため、群の構造が見やすくなる。

    引数:
        df (pd.DataFrame): タンパク質定量データ（行=タンパク質, 列=サンプル）
        sample_info (pd.DataFrame): サンプル情報（Sample列, Condition列）

    戻り値:
        なし（図をファイルに保存）
    """
    print("  (a) 相関行列ヒートマップ")

    # ----------------------------------------------------------
    # ピアソン相関係数の計算
    # ----------------------------------------------------------
    # df.corr() はデフォルトでピアソン相関係数を計算する
    # method="pearson" を明示的に指定（他に "spearman", "kendall" もある）
    # 結果はサンプル数 × サンプル数の正方行列になる
    corr = df.corr(method="pearson")

    # ----------------------------------------------------------
    # サンプルの色分け（行・列に沿ったカラーバー）
    # ----------------------------------------------------------
    # sample_info の "Condition" 列を "Sample" 列をキーにして辞書化
    # → 各サンプル名に対応する条件（Normal/Tumor）を取得
    conditions = sample_info.set_index("Sample")["Condition"]

    # 条件名を色に変換（Normal → 青, Tumor → 赤）
    # .map() は値を辞書に基づいて置換するメソッド
    colors = conditions.map({"Normal": COLOR_NORMAL, "Tumor": COLOR_TUMOR})

    # 相関行列のインデックス順に色を並べ替え
    # reindex() で相関行列のサンプル順に合わせる
    row_colors = colors.reindex(corr.index)

    # ----------------------------------------------------------
    # クラスタマップの作成
    # ----------------------------------------------------------
    # sns.clustermap() は seaborn のヒートマップ＋クラスタリング関数
    g = sns.clustermap(
        corr,                   # 描画するデータ（相関行列）
        method="average",       # 【クラスタリング手法】average = UPGMA法（群平均法）
                                #   2つのクラスタ間の距離を「全ペアの平均距離」で定義
                                #   他に "ward"（ウォード法）, "complete"（最遠隣法）等がある
        metric="correlation",   # 【距離指標】correlation = 1 - ピアソン相関係数
                                #   相関が高い（似ている）ほど距離が短くなる
        cmap="Reds",            # 【カラーマップ】Reds: 白=低い相関、濃い赤=高い相関
                                #   論文原図と同じ赤単色グラデーション
                                #   対角線（相関1.0）が最も濃い赤、群間（低相関）が薄い色になる
        vmin=0.7,               # 【色の下限値】0.7以下の相関はすべて同じ色（薄黄）で表示
                                #   → 論文では0.744程度の値も見られるため、0.7を下限に設定
        vmax=1.0,               # 【色の上限値】1.0（完全一致）が最も赤くなる
        figsize=(10, 10),       # 図のサイズ（横10インチ × 縦10インチ）
        row_colors=row_colors,  # 行の左に表示するカラーバー（Normal=青, Tumor=赤）
        col_colors=row_colors,  # 列の上に表示するカラーバー（同じ色分け）
        linewidths=0,           # セル間の線の太さ（0 = 線なし）
        xticklabels=True,       # X軸のサンプル名を表示する
        yticklabels=True,       # Y軸のサンプル名を表示する
    )

    # X軸ラベルのフォントサイズと回転角度を調整
    # fontsize=6: サンプル名が多いため小さめのフォント
    # rotation=90: ラベルを90度回転（縦書き）して重ならないようにする
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=6, rotation=90)

    # Y軸ラベルのフォントサイズを調整
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=6)

    # ----------------------------------------------------------
    # 凡例（レジェンド）の追加
    # ----------------------------------------------------------
    # matplotlib.patches.Patch を使って色付き四角の凡例を手動作成
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=COLOR_NORMAL, label="Non-tumor"),   # 青い四角 + ラベル
        Patch(facecolor=COLOR_TUMOR, label="Tumor"),        # 赤い四角 + ラベル
    ]

    # ヒートマップの右上に凡例を配置
    # bbox_to_anchor=(1.05, 1.0): 図の右端外側（x=1.05）の上端（y=1.0）に配置
    # frameon=False: 凡例の枠線を表示しない
    g.ax_heatmap.legend(
        handles=legend_elements, loc="upper left",
        bbox_to_anchor=(1.05, 1.0), frameon=False
    )

    # ----------------------------------------------------------
    # 図の保存
    # ----------------------------------------------------------
    filepath = os.path.join(FIG_DIR, "fig1a_correlation.png")
    # dpi=150: 解像度150ドット/インチ（論文品質には300が推奨だが、ファイルサイズ節約）
    # bbox_inches="tight": 余白を自動調整して図全体がきれいに収まるようにする
    g.savefig(filepath, dpi=150, bbox_inches="tight")

    # メモリ節約のため図を閉じる（plt.close() しないとメモリリークの原因になる）
    plt.close()
    print(f"    保存: {filepath}")

    # ----------------------------------------------------------
    # 相関係数の数値データもCSVとして保存（後で数値確認に使える）
    # ----------------------------------------------------------
    corr.to_csv(os.path.join(RESULTS_DIR, "tables", "correlation_matrix.csv"))


# ============================================================
# (b) 教師なし階層的クラスタリングの作成関数
# ============================================================
def plot_hierarchical_clustering(df, sample_info):
    """
    Figure 1b: 教師なし階層的クラスタリングのヒートマップを作成する関数。

    【階層的クラスタリングとは？】
      データの「似ている度合い」に基づいて、サンプルやタンパク質を
      ツリー状（デンドログラム）に分類する手法。
      「教師なし」= Normal/Tumorのラベルを使わずに、データだけで分類する。
      もしNormalとTumorが正しくクラスタとして分離すれば、
      データに生物学的な差が確かに存在する証拠になる。

    【ウォード法（Ward's method）とは？】
      クラスタリングの手法の1つ。2つのクラスタを結合したときに
      「クラスタ内の分散増加量が最小になるペア」を優先的に結合する。
      比較的コンパクトで均等な大きさのクラスタを作る傾向がある。
      生物学データでよく使用される。

    【ユークリッド距離とは？】
      2点間の直線距離（ピタゴラスの定理）。
      発現プロファイルが近いサンプルほど距離が短くなる。

    【z_score（Zスコア標準化）とは？】
      各行（サンプル）の値を平均0、標準偏差1に変換する処理。
      サンプル間で発現レベルの絶対値が異なっても、
      「相対的な高低パターン」で比較できるようになる。
      計算式: z = (値 - 平均) / 標準偏差

    引数:
        df (pd.DataFrame): タンパク質定量データ（行=タンパク質, 列=サンプル）
        sample_info (pd.DataFrame): サンプル情報（Sample列, Condition列）

    戻り値:
        なし（図をファイルに保存）
    """
    print("  (b) 階層的クラスタリング")

    # ----------------------------------------------------------
    # サンプルの色分け
    # ----------------------------------------------------------
    # sample_info から条件（Normal/Tumor）を取得し、色に変換
    conditions = sample_info.set_index("Sample")["Condition"]
    sample_colors = conditions.map({"Normal": COLOR_NORMAL, "Tumor": COLOR_TUMOR})

    # ----------------------------------------------------------
    # タンパク質ごとのfold change方向を計算（列カラーバー用）
    # ----------------------------------------------------------
    # 論文のFigure 1bでは、タンパク質（列）の上にカラーバーがあり、
    # 赤 = 腫瘍で上昇（Up-regulated）、青 = 腫瘍で低下（Down-regulated）を示す
    normal_samples_list = sample_info[sample_info["Condition"] == "Normal"]["Sample"].tolist()
    tumor_samples_list = sample_info[sample_info["Condition"] == "Tumor"]["Sample"].tolist()

    # 各タンパク質のlog2 fold changeを計算（Tumor - Normal の平均差）
    normal_mean = df[[s for s in normal_samples_list if s in df.columns]].mean(axis=1)
    tumor_mean = df[[s for s in tumor_samples_list if s in df.columns]].mean(axis=1)
    fc = tumor_mean - normal_mean

    # Up/Down を色に変換
    protein_direction_colors = fc.apply(
        lambda x: "#E74C3C" if x > 0 else "#3498DB"  # 赤=Up, 青=Down
    )

    # ----------------------------------------------------------
    # クラスタマップの作成
    # ----------------------------------------------------------
    # df.T で転置: サンプル（行）× タンパク質（列）
    g = sns.clustermap(
        df.T,                   # 転置データ: サンプル（行）× タンパク質（列）
        method="ward",          # 【クラスタリング手法】ウォード法
        metric="euclidean",     # 【距離指標】ユークリッド距離
        cmap="RdBu_r",          # 【カラーマップ】Red-Blue (reversed): 赤=高発現、青=低発現
        center=0,               # 【色の中心値】0を白色にする
        figsize=(14, 8),        # 図のサイズ
        row_colors=sample_colors.reindex(df.columns),  # 行（サンプル）のカラーバー: Histology
        col_colors=protein_direction_colors,  # 列（タンパク質）のカラーバー: Expression pattern
                                # 赤=Up-regulated in tumor, 青=Down-regulated in tumor
        col_cluster=True,       # 列（タンパク質）もクラスタリングする
        row_cluster=True,       # 行（サンプル）もクラスタリングする
        xticklabels=False,      # X軸ラベル（タンパク質名）は非表示
        yticklabels=True,       # Y軸ラベル（サンプル名）は表示する
        z_score=1,              # 【Zスコア標準化】1 = 列方向（タンパク質方向）で標準化
                                #   各タンパク質の平均=0, SD=1に変換
                                #   タンパク質間の発現レベル差を除去し、
                                #   サンプル間の相対的な発現パターンを可視化する
        vmin=-3, vmax=3,        # 色の範囲を固定（Z-scoreの±3）
    )

    # Y軸（サンプル名）のフォントサイズを調整
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=7)

    # 軸ラベルの設定
    g.ax_heatmap.set_xlabel("Proteins")   # X軸: タンパク質
    g.ax_heatmap.set_ylabel("Samples")    # Y軸: サンプル

    # ----------------------------------------------------------
    # 凡例の追加（Protein expression pattern）
    # ----------------------------------------------------------
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#E74C3C", label="Up-regulated in tumor"),
        Patch(facecolor="#3498DB", label="Down-regulated in tumor"),
    ]
    g.ax_heatmap.legend(
        handles=legend_elements, loc="lower right",
        bbox_to_anchor=(1.3, -0.15), frameon=False, fontsize=8
    )

    # ----------------------------------------------------------
    # 図の保存
    # ----------------------------------------------------------
    filepath = os.path.join(FIG_DIR, "fig1b_clustering.png")
    g.savefig(filepath, dpi=150, bbox_inches="tight")
    plt.close()  # メモリ節約のため図を閉じる
    print(f"    保存: {filepath}")


# ============================================================
# 95%信頼楕円の描画ヘルパー関数
# ============================================================
def confidence_ellipse(x, y, ax, n_std=2.0, **kwargs):
    """
    2変量データの信頼楕円を描画するヘルパー関数。

    【信頼楕円とは？】
      2変量正規分布を仮定したとき、データの散らばりを楕円で表現したもの。
      n_std=2.0 のとき約95%の信頼領域に相当する。
      PCAプロットで群の広がりと重なりを視覚的に示すために使用する。

    引数:
        x (array-like): X座標の値
        y (array-like): Y座標の値
        ax (matplotlib.axes.Axes): 描画先のAxesオブジェクト
        n_std (float): 標準偏差の倍数（2.0で約95%信頼領域）
        **kwargs: Ellipse に渡す追加キーワード引数（facecolor, edgecolor, alpha 等）

    戻り値:
        matplotlib.patches.Ellipse または None（データ不足時）
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


# ============================================================
# (c) 主成分分析（PCA）の作成関数
# ============================================================
def plot_pca(df, sample_info):
    """
    Figure 1c: PCA（主成分分析）のスコアプロットを作成する関数。

    【PCA（主成分分析）とは？】
      高次元データ（ここでは数百〜数千タンパク質の発現値）を、
      情報をなるべく失わずに少ない次元（ここでは2次元）に圧縮する手法。
      - PC1（第1主成分）: データの分散が最も大きい方向
      - PC2（第2主成分）: PC1と直交する方向で、次に分散が大きい方向
      PC1とPC2の2軸で散布図を描くことで、数千次元のデータを2次元平面で可視化できる。

    【寄与率（Explained Variance Ratio）とは？】
      各主成分が元データの分散をどれだけ説明しているかの割合。
      例: PC1の寄与率=40% → PC1だけで元データの40%の情報を説明できる
      PC1+PC2の寄与率が高いほど、2Dプロットでの可視化が元データをよく反映している。

    【期待される結果】
      NormalサンプルとTumorサンプルがPCA空間上で明確に分離すれば、
      2群の発現パターンに系統的な違いがあることの強い証拠になる。

    引数:
        df (pd.DataFrame): タンパク質定量データ（行=タンパク質, 列=サンプル）
        sample_info (pd.DataFrame): サンプル情報（Sample列, Condition列）

    戻り値:
        なし（図をファイルに保存、寄与率テーブルもCSVに保存）
    """
    print("  (c) 主成分分析（PCA）")

    # ----------------------------------------------------------
    # 条件情報の取得
    # ----------------------------------------------------------
    # sample_info から各サンプルの条件（Normal/Tumor）を取得
    conditions = sample_info.set_index("Sample")["Condition"]

    # ----------------------------------------------------------
    # PCA の実行
    # ----------------------------------------------------------
    # PCA(n_components=2): 上位2主成分（PC1, PC2）を計算
    # n_components を増やすと PC3, PC4... も計算できるが、2Dプロットには2つで十分
    pca = PCA(n_components=2)

    # fit_transform() でPCAモデルの学習と変換を同時に行う
    # df.T: 転置して サンプル（行）× タンパク質（列）の形にする
    #   PCAの入力は「サンプルが行、変数が列」の形式が必要
    # scores: サンプル数 × 2 の配列（各サンプルのPC1, PC2スコア）
    scores = pca.fit_transform(df.T)

    # ----------------------------------------------------------
    # 散布図の作成
    # ----------------------------------------------------------
    # fig: 図全体のオブジェクト、ax: 1つのプロット領域（Axes）
    # figsize=(8, 6): 横8インチ × 縦6インチ
    fig, ax = plt.subplots(figsize=(8, 6))

    # --- 論文準拠のカラー設定 ---
    # 論文原図に合わせ、Normal=#3498DB（青）、Tumorはステージ別に色分け
    # 論文: Non-tumor=青, Stage I=緑, Stage IIA=黄, Stage IIIB=橙, Stage IVC=紫
    COLOR_NORMAL_PCA = "#3498DB"  # 青（Non-tumor）

    # 臨床ステージ情報を読み込み（存在する場合）
    import re
    clinical_path = os.path.join(RESULTS_DIR, "clinical_info.csv")
    stage_map = {}
    if os.path.exists(clinical_path):
        clinical = pd.read_csv(clinical_path)
        for _, row in clinical.iterrows():
            stage_map[row["Sample_N"]] = "Normal"
            stage_map[row["Sample_T"]] = row["Stage"]

    # ステージ別カラー（論文準拠）
    stage_colors = {
        "Normal": "#3498DB",   # 青（Non-tumor tissue）
        "I":      "#66BB6A",   # 緑（Stage I）
        "II":     "#FFD54F",   # 黄（Stage II / IIA）
        "III":    "#FFA726",   # オレンジ（Stage III / IIIB）
        "IV":     "#AB47BC",   # 紫（Stage IV / IVC）
    }

    # 各サンプルのステージを割り当て
    sample_stages = []
    for s in df.columns:
        clean = re.sub(r"_dup\d+$", "", s)
        if clean in stage_map:
            sample_stages.append(stage_map[clean])
        elif "-N" in s:
            sample_stages.append("Normal")
        else:
            sample_stages.append("Tumor")  # ステージ不明のTumor用フォールバック

    # サンプルごとにステージ色でプロット
    plotted_labels = set()
    for i, (sample, stage) in enumerate(zip(df.columns, sample_stages)):
        color = stage_colors.get(stage, "#E74C3C")
        # 凡例にはステージごとに1回だけ表示
        if stage == "Normal":
            label_text = "Non-tumor tissue"
        else:
            label_text = f"Stage {stage}"
        label = label_text if label_text not in plotted_labels else None
        plotted_labels.add(label_text)

        ax.scatter(
            scores[i, 0],              # X座標: PC1スコア
            scores[i, 1],              # Y座標: PC2スコア
            c=color,                   # 点の色（ステージ別）
            marker="o",                # 点の形状（丸）
            s=120,                     # 点のサイズ（論文原図に合わせて大きめ）
            alpha=0.85,                # 透過度
            label=label,               # 凡例ラベル
            edgecolors="white",        # 点の縁の色
            linewidth=0.8,             # 点の縁の太さ
            zorder=3,                  # 描画順（グリッドの上に表示）
        )

    # 信頼楕円は描画しない（論文原図に合わせて除去）

    # ----------------------------------------------------------
    # 軸ラベルの設定（寄与率を％で表示）
    # ----------------------------------------------------------
    # pca.explained_variance_ratio_ にPC1, PC2の寄与率が格納されている
    # ×100 で百分率に変換、:.1f で小数第1位まで表示
    # 論文原図に合わせ「Component 1」「Component 2」表記を使用
    ax.set_xlabel(f"Component 1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"Component 2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")

    # 図のタイトル
    ax.set_title("PCA of Non-tumor and Tumor Tissues")

    # 凡例を表示（frameon=False: 凡例の枠線なし）
    ax.legend(frameon=False)

    # グリッド線を追加（論文原図に合わせ、薄いグレーで背景に表示）
    ax.grid(True, color="gray", alpha=0.3, linestyle="-", linewidth=0.5)

    # 上と右の枠線（スパイン）を非表示にする
    # 学術論文では左と下の2本だけ残すスタイルが一般的
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ----------------------------------------------------------
    # 図の保存
    # ----------------------------------------------------------
    filepath = os.path.join(FIG_DIR, "fig1c_pca.png")
    fig.savefig(filepath, dpi=150, bbox_inches="tight")
    plt.close()  # メモリ節約のため図を閉じる
    print(f"    保存: {filepath}")

    # ----------------------------------------------------------
    # 寄与率テーブルの保存
    # ----------------------------------------------------------
    # 各主成分の寄与率と累積寄与率をDataFrameにまとめてCSVに保存
    # これにより、PC1+PC2でデータ全体のどれだけを説明できているか確認できる
    var_df = pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(len(pca.explained_variance_ratio_))],  # PC番号（PC1, PC2）
        "Variance_Ratio": pca.explained_variance_ratio_,      # 各主成分の寄与率
        "Cumulative": np.cumsum(pca.explained_variance_ratio_),  # 累積寄与率
        # np.cumsum() は累積和を計算: [0.4, 0.2] → [0.4, 0.6]
    })
    var_df.to_csv(os.path.join(RESULTS_DIR, "tables", "pca_variance.csv"), index=False)


# ============================================================
# メイン処理
# ============================================================
def main():
    """
    メイン関数: Step 3 の全処理を順番に実行する。

    【処理の流れ】
      1. 前処理済みデータの読み込み（Step 2 の出力）
      2. 相関行列ヒートマップの作成（Figure 1a）
      3. 階層的クラスタリングの作成（Figure 1b）
      4. PCA散布図の作成（Figure 1c）
    """
    print("=" * 60)
    print("Step 3: 全体像の可視化（Figure 1 の再現）")
    print("=" * 60)

    # テーブル保存用ディレクトリを作成
    os.makedirs(os.path.join(RESULTS_DIR, "tables"), exist_ok=True)

    # ----------------------------------------------------------
    # データ読み込み
    # ----------------------------------------------------------
    df, sample_info = load_data()

    # データの概要を表示（行=タンパク質数、列=サンプル数）
    print(f"\nデータ: {df.shape[0]} タンパク質 × {df.shape[1]} サンプル")

    # ----------------------------------------------------------
    # 3つの図を順番に作成
    # ----------------------------------------------------------
    print("\n--- Figure 1 の再現 ---")
    plot_correlation_matrix(df, sample_info)          # (a) 相関行列ヒートマップ
    plot_hierarchical_clustering(df, sample_info)     # (b) 階層的クラスタリング
    plot_pca(df, sample_info)                         # (c) PCA散布図

    # ----------------------------------------------------------
    # 完了メッセージ
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("可視化完了！ results/figures/ を確認してください。")
    print("=" * 60)


# ============================================================
# スクリプトの実行エントリーポイント
# ============================================================
# このファイルを直接実行した場合のみ main() を呼び出す。
# 他のスクリプトから import した場合は実行されない。
if __name__ == "__main__":
    main()
