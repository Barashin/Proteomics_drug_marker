---
title: "Welch's t-testとVolcanoプロットで差分発現タンパク質を見つける【論文再現シリーズ #8】"
emoji: "🌋"
type: "tech"
topics: ["proteomics", "python", "statistics", "labcode"]
published: false
---

# Welch's t-testとVolcanoプロットで差分発現タンパク質を見つける

## はじめに

この記事では、腫瘍 vs 正常組織でタンパク質の発現量が有意に異なるものを統計的に特定します。Toyota et al. 2025 では**2,642タンパク質**が差分発現として報告されています（↑1,475、↓1,167）。本書の解析（sage + 32ファイル、2,081 タンパク質）では **1,055 タンパク質（↑867、↓188）** が有意差として検出されました。

:::message
**この記事で行う処理**
各タンパク質についてNormal群とTumor群の平均発現量をWelchのt検定で比較し、統計的に有意（p<0.05）かつ発現変化が大きい（fold change > 2倍）タンパク質を特定します。結果をVolcanoプロットで可視化し、腫瘍で増加・減少したタンパク質を一目で把握します。さらにTop 50〜200のタンパク質でクラスタリングとPCAを行い、少数のマーカーだけでも群分離が可能かを検証します。
:::

## 前提

- [#7 全体像の可視化](article-07-visualization.md) が完了していること
- **対応Notebook**: [`notebooks/step_08.ipynb`](../notebooks/step_08.ipynb) — この記事のコードをセルごとに実行できます

## コード全文（対応Notebook: step_08.ipynb）

### ライブラリと設定

```python
import numpy as np           # 数値計算ライブラリ（配列操作・数学関数に使用）
import pandas as pd          # データフレーム操作ライブラリ（表形式データの読み込み・加工に使用）
import matplotlib.pyplot as plt  # グラフ描画ライブラリ（Volcanoプロット等の作成に使用）
from matplotlib.patches import Ellipse  # 楕円パッチ（PCAプロットに95%信頼楕円を描画するため）
import matplotlib.transforms as transforms  # 座標変換（楕円の回転・拡縮・移動に使用）
import seaborn as sns        # 統計的可視化ライブラリ（ヒートマップ等の高水準プロットに使用）
from scipy import stats      # 統計検定ライブラリ（Welch's t-testの実行に使用）
from sklearn.decomposition import PCA  # 主成分分析（PCA）クラス（次元削減による群分離の可視化に使用）

# --- パス ---
RESULTS  = "../results"              # 解析結果の保存先ディレクトリへのパス
FIG_DIR  = f"{RESULTS}/figures"      # 図の保存先ディレクトリへのパス
TABLE_DIR = f"{RESULTS}/tables"      # テーブル（CSV等）の保存先ディレクトリへのパス

# --- 閾値 ---
# p値の閾値: 0.05未満なら「統計的に有意」と判定する（慣例的に広く使われる基準）
P_THRESHOLD    = 0.05
# log2 fold changeの閾値: 1.0以上なら2倍以上の発現変化があることを意味する
LOG2FC_THRESHOLD = 1.0

# --- カラー ---
NORMAL, TUMOR = "#3498DB", "#E74C3C"  # Normal群を青、Tumor群を赤で表示する色コード
UP, DOWN, NS  = "#E74C3C", "#3498DB", "#CCCCCC"  # 上昇(赤)・低下(青)・非有意(灰)の色コード
# 条件名から色への対応辞書（グラフ描画時に群ごとの色を自動で割り当てるため）
COND_MAP = {"Normal": NORMAL, "Tumor": TUMOR}
```

### データ読み込み

```python
# データ読み込み & Normal/Tumor サンプル分離
# 前処理済みタンパク質発現データを読み込む（行=タンパク質、列=サンプル、値=log2発現量）
df = pd.read_csv(f"{RESULTS}/preprocessed_data.csv", index_col=0)
# サンプル情報（各サンプルがNormalかTumorか）を読み込む
sample_info = pd.read_csv(f"{RESULTS}/sample_info.csv")
# サンプル名をインデックスにして条件列だけを取り出す（後でサンプルから条件を引くための辞書的Series）
conditions = sample_info.set_index("Sample")["Condition"]

# Normal群のサンプル名リストを取得（Conditionが'Normal'の行をフィルタ）
normal_samples = sample_info.query("Condition == 'Normal'")["Sample"].tolist()
# Tumor群のサンプル名リストを取得（Conditionが'Tumor'の行をフィルタ）
tumor_samples  = sample_info.query("Condition == 'Tumor'")["Sample"].tolist()
# タンパク質数と各群のサンプル数を表示して、データの規模を確認する
print(f"Proteins: {len(df)}, Normal: {len(normal_samples)}, Tumor: {len(tumor_samples)}")
```

### Welch's t-test

```python
def welch_ttest(df, normal_samples, tumor_samples):
    """全タンパク質に Welch t-test を実行し DataFrame を返す。"""
    rows = []  # 各タンパク質の検定結果を一時的に格納するリスト
    # 全タンパク質を1つずつループして検定を行う
    for protein in df.index:
        # そのタンパク質のNormal群の発現量を取得し、欠損値（NaN）を除外する
        nv = df.loc[protein, normal_samples].dropna()
        # そのタンパク質のTumor群の発現量を取得し、欠損値（NaN）を除外する
        tv = df.loc[protein, tumor_samples].dropna()
        # 各群に2個以上のデータがないとt検定が計算できないためスキップする
        if len(nv) < 2 or len(tv) < 2:
            continue
        # Welchのt検定を実行（equal_var=Falseで等分散を仮定しない）
        # t_stat: t統計量（値が大きいほど2群の差が大きい）
        # p_val: p値（この差が偶然生じる確率。小さいほど有意）
        t_stat, p_val = stats.ttest_ind(tv, nv, equal_var=False)
        # log2 fold change = Tumor群の平均 - Normal群の平均（データが既にlog2スケールのため引き算で比を表す）
        # 正の値なら腫瘍で発現増加、負の値なら腫瘍で発現減少を意味する
        log2fc = tv.mean() - nv.mean()
        # 検定結果を辞書形式でリストに追加する
        rows.append({
            "Protein": protein,               # タンパク質名
            "Mean_Normal": nv.mean(),          # Normal群の平均発現量（log2スケール）
            "Mean_Tumor": tv.mean(),           # Tumor群の平均発現量（log2スケール）
            "Log2FC": log2fc,                  # 発現変化量（正=腫瘍で増加、負=腫瘍で減少）
            "T_statistic": t_stat,             # t統計量（2群の差の大きさを標準誤差で割った値）
            "P_value": p_val,                  # p値（帰無仮説「2群に差がない」が正しい確率）
            # -log10(p値): Volcanoプロットのy軸用。p値が小さいほど値が大きくなり上に表示される
            # max(p_val, 1e-300)でp値が0の場合のlog計算エラーを防ぐ
            "Neg_log10_P": -np.log10(max(p_val, 1e-300)),
        })
    # リストをDataFrameに変換する（1行=1タンパク質の検定結果）
    result = pd.DataFrame(rows)
    # 有意性の分類: まず全タンパク質を "NS"（Not Significant = 非有意）に設定する
    result["Significant"] = "NS"
    # p値が閾値未満 かつ log2FCが正の閾値を超える → "Up"（腫瘍で有意に増加）
    result.loc[(result["P_value"] < P_THRESHOLD) & (result["Log2FC"] >  LOG2FC_THRESHOLD), "Significant"] = "Up"
    # p値が閾値未満 かつ log2FCが負の閾値を下回る → "Down"（腫瘍で有意に減少）
    result.loc[(result["P_value"] < P_THRESHOLD) & (result["Log2FC"] < -LOG2FC_THRESHOLD), "Significant"] = "Down"
    return result  # 全タンパク質の検定結果DataFrameを返す
```

### t検定の実行

```python
# 上で定義したwelch_ttest関数を実行し、全タンパク質の検定結果を取得する
result_df = welch_ttest(df, normal_samples, tumor_samples)

# 腫瘍で有意に増加したタンパク質の数をカウントする
n_up   = (result_df["Significant"] == "Up").sum()
# 腫瘍で有意に減少したタンパク質の数をカウントする
n_down = (result_df["Significant"] == "Down").sum()
# 検定を実行できたタンパク質の総数を表示する
print(f"検定タンパク質数: {len(result_df)}")
# 有意差ありの内訳（Up: 腫瘍で増加、Down: 腫瘍で減少）を表示する
print(f"有意差あり: {n_up + n_down} (Up {n_up}, Down {n_down})")

# 検定結果をCSVファイルに保存する（後の解析やCOSMIC照合で使用する）
result_df.to_csv(f"{TABLE_DIR}/differential_proteins.csv", index=False)
```

### Volcanoプロット

```python
# 図と軸オブジェクトを作成する（figsize=(横8, 縦6)インチ）
fig, ax = plt.subplots(figsize=(8, 6))

# NS（非有意）→ Up（増加）→ Down（減少）の順にプロットする
# NSを先に描画することで、有意な点が上に表示されて見やすくなる
for sig, color, alpha in [("NS", NS, 0.3), ("Up", UP, 0.6), ("Down", DOWN, 0.6)]:
    # 各カテゴリに該当するタンパク質のブールマスクを作成する
    m = result_df["Significant"] == sig
    # 散布図を描画する（x=発現変化量, y=統計的有意性）
    # c: 点の色, s: 点のサイズ(10pt), alpha: 透明度（NSは薄く表示）
    ax.scatter(result_df.loc[m, "Log2FC"], result_df.loc[m, "Neg_log10_P"],
              c=color, s=10, alpha=alpha, label=f"{sig} ({m.sum()})")

# p値の閾値ライン（水平破線）: この線より上が統計的に有意
ax.axhline(-np.log10(P_THRESHOLD), color="gray", ls="--", lw=0.5)
# log2FCの正の閾値ライン（垂直破線）: この線より右が2倍以上の増加
ax.axvline( LOG2FC_THRESHOLD, color="gray", ls="--", lw=0.5)
# log2FCの負の閾値ライン（垂直破線）: この線より左が2倍以上の減少
ax.axvline(-LOG2FC_THRESHOLD, color="gray", ls="--", lw=0.5)

# x軸ラベル: 発現変化の方向と大きさを示す
ax.set_xlabel("Log2 Fold Change (Tumor / Normal)")
# y軸ラベル: 統計的有意性を示す（値が大きいほど有意）
ax.set_ylabel("-Log10(P-value)")
# グラフタイトル: 腫瘍 vs 非腫瘍の差分タンパク質量解析
ax.set_title("Differential Protein Abundance: Tumor vs Non-tumor")
# 凡例を右上に表示する（frameon=Falseで枠線を非表示にしてすっきりさせる）
ax.legend(frameon=False, loc="upper right")
# 上辺と右辺の枠線を非表示にする（学術論文スタイルのすっきりした見た目にする）
ax.spines[["top", "right"]].set_visible(False)

# 図をPNGファイルとして保存する（dpi=150で高解像度、bbox_inches="tight"で余白を最小化）
fig.savefig(f"{FIG_DIR}/fig_bonus_volcano.png", dpi=150, bbox_inches="tight")
# 画面に表示する
plt.show()
```

### 信頼楕円のヘルパー関数

```python
def confidence_ellipse(x, y, ax, n_std=2.0, **kwargs):
    """95% 信頼楕円を描画する。"""
    # データ点が2個未満では共分散を計算できないため、描画をスキップする
    if len(x) < 2:
        return
    # xとyの共分散行列を計算する（2x2行列: 分散と共分散を含む）
    cov = np.cov(x, y)
    # ピアソン相関係数を計算する（xとyの線形関係の強さ。-1〜+1の範囲）
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    # 楕円の横半径を計算する（相関が高いほど横長になる）
    ell_rx = np.sqrt(1 + pearson)
    # 楕円の縦半径を計算する（相関が高いほど縦短になる）
    ell_ry = np.sqrt(1 - pearson)
    # 原点中心の楕円オブジェクトを作成する（**kwargsで色や透明度を受け取る）
    ellipse = Ellipse((0, 0), width=ell_rx * 2, height=ell_ry * 2, **kwargs)
    # アフィン変換を組み合わせて楕円をデータの分布に合わせる
    transf = (transforms.Affine2D()
              .rotate_deg(45)  # 45度回転（共分散の主軸方向に合わせる）
              # n_std倍の標準偏差に合わせてスケールする（n_std=2で約95%信頼区間）
              .scale(np.sqrt(cov[0, 0]) * n_std, np.sqrt(cov[1, 1]) * n_std)
              # データの重心位置に移動する
              .translate(np.mean(x), np.mean(y)))
    # 変換を楕円に適用する（ax.transDataでデータ座標系に変換）
    ellipse.set_transform(transf + ax.transData)
    # 楕円をAxesに追加して描画する
    return ax.add_patch(ellipse)
```

### 全有意差タンパク質のヒートマップとPCA

```python
# --- 有意差タンパク質の抽出 ---
# 検定結果から有意差あり（UpまたはDown）のタンパク質だけを抽出する
sig_df = result_df[result_df["Significant"] != "NS"]
# 元の発現量データから有意差タンパク質の行だけを取り出す
sig_data = df.loc[df.index.isin(sig_df["Protein"])]
# 各タンパク質の発現変化方向（Log2FC）をSeriesとして取得する（色分けに使用）
direction = sig_df.set_index("Protein")["Log2FC"]

# 各サンプルの条件（Normal/Tumor）に応じた色をSeriesとして作成する（ヒートマップの行カラーバーに使用）
sample_colors  = conditions.reindex(df.columns).map(COND_MAP)
# 各タンパク質の発現変化方向に応じた色をSeriesとして作成する（ヒートマップの列カラーバーに使用）
# Log2FC > 0 なら赤（Up: 腫瘍で増加）、それ以外は青（Down: 腫瘍で減少）
protein_colors = pd.Series(
    sig_data.index.map(lambda p: UP if direction.get(p, 0) > 0 else DOWN),
    index=sig_data.index, name="Expression"
)

# --- Figure 2a: ヒートマップ ---
# 階層的クラスタリング付きヒートマップを作成する
g = sns.clustermap(
    sig_data.T,              # データを転置する（行=サンプル、列=タンパク質にする）
    method="ward",           # Ward法でクラスタリングする（分散を最小化する手法）
    cmap="RdBu_r",           # カラーマップ: 赤=高発現、青=低発現（反転版）
    z_score=1,               # 列方向（タンパク質ごと）に標準化する（平均0、標準偏差1に変換）
    row_colors=sample_colors,   # 行（サンプル）の横にNormal/Tumorの色バーを表示する
    col_colors=protein_colors,  # 列（タンパク質）の上にUp/Downの色バーを表示する
    xticklabels=False,       # x軸のタンパク質名ラベルを非表示にする（数が多すぎて読めないため）
    yticklabels=True,        # y軸のサンプル名ラベルを表示する
    figsize=(10, 8),         # 図のサイズを横10×縦8インチに設定する
    vmin=-3, vmax=3,         # 色の範囲を-3〜+3に制限する（外れ値による色の偏りを防ぐ）
)
# y軸のサンプル名のフォントサイズを小さくして見やすくする
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=7)
# ヒートマップをPNGファイルとして保存する
g.savefig(f"{FIG_DIR}/fig2a_heatmap_all.png", dpi=150, bbox_inches="tight")
plt.show()

# --- Figure 2b: PCA ---
# 主成分分析（PCA）オブジェクトを作成する（2次元に次元削減する）
pca = PCA(n_components=2)
# 有意差タンパク質の発現データでPCAを実行する（転置して行=サンプルにする）
# scores: 各サンプルのPC1, PC2座標を格納した配列
scores = pca.fit_transform(sig_data.T)

# PCAプロット用の図と軸を作成する
fig, ax = plt.subplots(figsize=(8, 6))
# Normal群とTumor群をそれぞれ異なる色でプロットする
for cond, color in [("Normal", NORMAL), ("Tumor", TUMOR)]:
    # 各サンプルが現在の条件に該当するかのブールマスクを作成する
    m = (conditions.reindex(df.columns) == cond).values
    # 散布図を描画する（s=80: 点のサイズ、alpha=0.8: やや透明、edgecolors="white": 白い縁取り）
    ax.scatter(scores[m, 0], scores[m, 1], c=color, s=80, alpha=0.8,
              label=cond, edgecolors="white")
    # 95%信頼楕円を描画する（群の分布範囲を視覚的に示す）
    confidence_ellipse(scores[m, 0], scores[m, 1], ax,
                       facecolor=color, alpha=0.15, edgecolor=color, lw=1.5)
# x軸ラベル: 第1主成分と寄与率（全分散のうちPC1が説明する割合）を表示する
ax.set_xlabel(f"Component 1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
# y軸ラベル: 第2主成分と寄与率を表示する
ax.set_ylabel(f"Component 2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
# グラフタイトル: 全有意差タンパク質を使ったPCA
ax.set_title("PCA using All Differentially Abundant Proteins")
# 凡例を表示する（frameon=Falseで枠線を非表示にする）
ax.legend(frameon=False)
# 上辺と右辺の枠線を非表示にする
ax.spines[["top", "right"]].set_visible(False)
# 薄い破線のグリッドを表示する（データの位置を読み取りやすくする）
ax.grid(True, alpha=0.3, ls="--")
# PCAプロットをPNGファイルとして保存する
fig.savefig(f"{FIG_DIR}/fig2b_pca_all.png", dpi=150, bbox_inches="tight")
plt.show()
```

### Top N 解析

```python
def plot_top_n(df, result_df, conditions):
    """Top 50/100/200 の clustering + PCA を一括描画する。"""
    # 有意差のあるタンパク質だけを抽出する（.copy()で元データを変更しないようにする）
    sig = result_df[result_df["Significant"] != "NS"].copy()
    # log2FCの絶対値を計算する（発現変化の大きさでソートするため。方向は問わない）
    sig["AbsLog2FC"] = sig["Log2FC"].abs()
    # 腫瘍で増加したタンパク質を発現変化の大きい順にソートする
    sig_up   = sig.query("Significant == 'Up'").sort_values("AbsLog2FC", ascending=False)
    # 腫瘍で減少したタンパク質を発現変化の大きい順にソートする
    sig_down = sig.query("Significant == 'Down'").sort_values("AbsLog2FC", ascending=False)

    # 各タンパク質の発現変化方向を取得する（色分けに使用する）
    direction = sig.set_index("Protein")["Log2FC"]
    # 各サンプルの条件に応じた色を作成する（ヒートマップの行カラーバー用）
    sample_colors = conditions.reindex(df.columns).map(COND_MAP)
    # 各Top Nに対応するFigure番号のラベル辞書（ヒートマップ用: Fig 2c, 2d, 2e）
    hm_label = {50: "c", 100: "d", 200: "e"}
    # 各Top Nに対応するFigure番号のラベル辞書（PCA用: Fig 2f, 2g, 2h）
    pc_label = {50: "f", 100: "g", 200: "h"}

    # Top 50, 100, 200 のそれぞれについて解析を実行する
    for n in [50, 100, 200]:
        # 上昇Top N個と低下Top N個のタンパク質名を結合する
        # min(n, len(...))で、データ数がNより少ない場合にエラーを防ぐ
        top_proteins = np.concatenate([
            sig_up.head(min(n, len(sig_up)))["Protein"].values,
            sig_down.head(min(n, len(sig_down)))["Protein"].values,
        ])
        # 元の発現量データからTop Nタンパク質の行だけを取り出す
        top_data = df.loc[df.index.isin(top_proteins)]
        # 各タンパク質の発現変化方向に応じた色をSeriesとして作成する（カラーバー用）
        prot_colors = pd.Series(
            top_data.index.map(lambda p: UP if direction.get(p, 0) > 0 else DOWN),
            index=top_data.index, name="Expression"
        )

        # 階層的クラスタリング付きヒートマップを作成する (Fig 2c-e)
        g = sns.clustermap(
            top_data.T,              # 転置して行=サンプル、列=タンパク質にする
            method="ward",           # Ward法でクラスタリングする
            cmap="RdBu_r",           # カラーマップ: 赤=高発現、青=低発現
            z_score=1,               # タンパク質ごとに標準化する
            row_colors=sample_colors,   # サンプルの条件を色バーで表示する
            col_colors=prot_colors,     # タンパク質のUp/Downを色バーで表示する
            xticklabels=False,       # タンパク質名を非表示にする
            yticklabels=True,        # サンプル名を表示する
            figsize=(8, 10),         # 図のサイズ（横8×縦10インチ）
            vmin=-3, vmax=3,         # 色の範囲を制限する
        )
        # サンプル名のフォントサイズを小さくして見やすくする
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=7)
        # ヒートマップをファイルに保存する（ファイル名にTop Nの数を含める）
        g.savefig(f"{FIG_DIR}/fig2{hm_label[n]}_clustering_top{n}.png",
                  dpi=150, bbox_inches="tight")
        plt.show()

        # Top NタンパクだけでPCAを実行する (Fig 2f-h)
        pca = PCA(n_components=2)  # 2次元に次元削減するPCAオブジェクトを作成する
        # PCAを実行してサンプルごとのPC1, PC2座標を得る
        scores = pca.fit_transform(top_data.T)
        # PCAプロット用の図と軸を作成する
        fig, ax = plt.subplots(figsize=(8, 6))
        # Normal群とTumor群をそれぞれプロットする
        for cond, color in [("Normal", NORMAL), ("Tumor", TUMOR)]:
            # 各サンプルが現在の条件に該当するかのブールマスクを作成する
            m = (conditions.reindex(df.columns) == cond).values
            # 散布図を描画する（s=80: 点のサイズ、edgecolors="white": 白い縁取り）
            ax.scatter(scores[m, 0], scores[m, 1], c=color, s=80, alpha=0.8,
                      label=cond, edgecolors="white")
            # 95%信頼楕円を描画して群の分布範囲を視覚的に示す
            confidence_ellipse(scores[m, 0], scores[m, 1], ax,
                               facecolor=color, alpha=0.15, edgecolor=color, lw=1.5)
        # x軸ラベル: 第1主成分と寄与率を表示する
        ax.set_xlabel(f"Component 1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        # y軸ラベル: 第2主成分と寄与率を表示する
        ax.set_ylabel(f"Component 2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        # グラフタイトル: 使用したTop Nの数を含める
        ax.set_title(f"PCA using Top {n} Differentially Abundant Proteins")
        # 凡例を枠なしで表示する
        ax.legend(frameon=False)
        # 上辺と右辺の枠線を非表示にする
        ax.spines[["top", "right"]].set_visible(False)
        # 薄い破線のグリッドを表示する
        ax.grid(True, alpha=0.3, ls="--")
        # PCAプロットをファイルに保存する
        fig.savefig(f"{FIG_DIR}/fig2{pc_label[n]}_pca_top{n}.png",
                    dpi=150, bbox_inches="tight")
        plt.show()
        # 実際に使用されたタンパク質数を表示する（Top N × 2方向の合計）
        print(f"Top {n}: {len(top_data)} proteins")
```

```python
# Top 50/100/200の解析を一括実行する（ヒートマップとPCAプロットが各Nごとに生成される）
plot_top_n(df, result_df, conditions)
```

---

## コード詳細

### stats.ttest_ind() のパラメータ

```python
# Welchのt検定を実行する: tvとnvの2群間で平均値に有意な差があるかを検定する
# equal_var=False で等分散を仮定しない（Welch版）。Studentのt検定より頑健
# t_stat: t統計量（2群の差の大きさ）、p_val: p値（差が偶然である確率）
t_stat, p_val = stats.ttest_ind(tv, nv, equal_var=False)
```

| パラメータ | 値 | 意味 |
|-----------|-----|------|
| 第1引数 | `tv` | 群1のデータ（腫瘍組織の発現量） |
| 第2引数 | `nv` | 群2のデータ（正常組織の発現量） |
| `equal_var` | `False` | **Welch版を使用**（等分散を仮定しない） |

- **t統計量（t_stat）**: 値が大きいほど2群の差が大きい
- **p値（p_val）**: この差が偶然生じる確率。0.05未満なら「統計的に有意」

### Log2 Fold Changeの解釈

| Log2FC | Fold Change | 意味 |
|--------|-------------|------|
| +3.0 | 8倍 | 腫瘍で8倍増加 |
| +1.0 | 2倍 | 腫瘍で2倍増加 |
| 0 | 1倍 | 変化なし |
| -1.0 | 0.5倍 | 腫瘍で半分に減少 |
| -3.0 | 0.125倍 | 腫瘍で1/8に減少 |

### Volcanoプロットの読み方

```
           高い有意性
              ↑
      ●       |       ●
     青(Down) |     赤(Up)
   腫瘍で減少  |   腫瘍で増加
              |
  ← 大きな減少 ── 0 ── 大きな増加 →
              |
       灰色（有意でない）
```

## 実行方法

**スクリプトで一括実行する場合:**

```bash
python scripts/step_07_differential_expression.py
```

**Notebook でセルごとに実行する場合:**

```bash
jupyter notebook notebooks/step_08.ipynb
```

Notebook版ではVolcanoプロットやヒートマップの出力をインラインで確認できます。

## 可視化結果

### Volcanoプロット

![Volcanoプロット](images/fig_bonus_volcano.png)

Volcanoプロットは差分発現解析の結果を一望するための散布図です。X軸にlog2 fold change（発現変化の大きさ）、Y軸に-log10(p値)（統計的有意性）をとっています。赤い点（Up: 867個）は腫瘍で有意に増加したタンパク質、青い点（Down: 188個）は有意に減少したタンパク質、灰色（NS: 1,026個）は有意でないタンパク質を示します。破線は閾値（p=0.05、|log2FC|=1.0）を示しています。

### 有意差タンパク質のヒートマップとPCA

![全有意差タンパク質ヒートマップ](images/fig2a_heatmap_all.png)

有意差のある1,055タンパク質全てについて、発現量のヒートマップを作成しました。行がタンパク質、列がサンプルで、階層的クラスタリングにより類似したパターンのタンパク質とサンプルがまとめられています。NormalとTumorで明瞭に異なる発現パターンが観察されます。

![全有意差タンパク質PCA](images/fig2b_pca_all.png)

有意差タンパク質のみを使ったPCA解析です。有意差タンパク質に絞ることで、Normal/Tumor間の分離がより鮮明になっています。

### Top Nタンパク質による詳細解析

有意差の大きい順にTop 50、Top 100、Top 200のタンパク質を抽出し、クラスタリングとPCAを実行しました。少数の強力なマーカーだけでも群分離が可能かを検証しています。

#### Top 50 タンパク質

![Top 50 クラスタリング](images/fig2c_clustering_top50.png)

最も有意差の大きい50タンパク質による階層的クラスタリングです。わずか50個のタンパク質でもNormal/Tumorの完全な分離が達成されており、バイオマーカーパネルとしての実用可能性を示唆しています。

![Top 50 PCA](images/fig2f_pca_top50.png)

Top 50タンパク質によるPCA。少数の特徴量でも明瞭な群分離が得られています。

#### Top 100 タンパク質

![Top 100 クラスタリング](images/fig2d_clustering_top100.png)

Top 100タンパク質によるクラスタリングです。タンパク質数を増やすことで、より詳細な発現パターンの構造が見えてきます。

![Top 100 PCA](images/fig2g_pca_top100.png)

Top 100タンパク質によるPCA。分離パターンが安定していることが確認できます。

#### Top 200 タンパク質

![Top 200 クラスタリング](images/fig2e_clustering_top200.png)

Top 200タンパク質によるクラスタリングです。より広範な発現変動パターンを捉えつつ、群分離は維持されています。

![Top 200 PCA](images/fig2h_pca_top200.png)

Top 200タンパク質によるPCA。タンパク質数を増やしても群構造は安定しており、解析の頑健性が示されています。

## 本書データでの実測値

| 指標 | 本書 (sage) | 論文 (DIA-NN) |
|------|-------------|--------------|
| 検定対象タンパク質 | 2,081 | 10,329 |
| 有意差 (p<0.05, FC>2) | **1,055** | **2,642** |
| ↑ Up-regulated | **867** | **1,475** |
| ↓ Down-regulated | **188** | **1,167** |

本書の Up/Down 比は 4.61 と、論文の 1.26 より大きくなっています。これは sage の理論スペクトルベース検索が **高発現の Up-regulated タンパク質を捉えやすく、低発現の Down-regulated は相対的に見落としやすい** 特性を反映しています。

## まとめ

Welch's t-testで有意に変動するタンパク質を特定し、Volcanoプロットで結果を俯瞰しました。本書データでは 1,055 個の有意差タンパク質（Up 867 + Down 188）を特定でき、論文の 2,642 個には及ばないものの、バイオロジカルに重要な主要ドライバー遺伝子（KRAS, CTNNB1, PIK3CA 等）は捕捉できています。

> 前回: [#7 全体像の可視化](article-07-visualization.md)
> 次回: [#9 COSMIC照合](article-09-cosmic.md) — がん関連タンパク質の同定

#バイオインフォマティクス #プロテオミクス #統計解析 #Volcanoプロット #labcode
