---
title: "深層学習DIA基本可視化解析【論文再現シリーズ #13b】"
emoji: "🎯"
type: "tech"
topics: ["proteomics", "python", "deeplearning", "visualization", "labcode"]
published: false
---

# 深層学習DIA基本可視化解析

## はじめに

前回（[#13a 差分発現統計解析](article-13a-openms-differential-stats.md)）で特定した有意差タンパク質について、ヒートマップとPCA解析による基本的な可視化を行います。深層学習による高精度検出の効果を確認し、全有意差タンパク質での群分離能力を評価します。

:::message
**この記事で行う処理**
前回のWelch's t検定で特定した有意差タンパク質を使用して、階層的クラスタリング付きヒートマップとPCA解析を実行します。深層学習による包括的検出の優位性を視覚的に確認します。
:::

## 前提

- [#13a 深層学習DIA差分発現統計解析](article-13a-openms-differential-stats.md) が完了していること
- **対応Notebook**: [`notebooks/step_13b_basic.ipynb`](../notebooks/step_13b_basic.ipynb) — この記事のコードをセルごとに実行できます

## コード全文（対応Notebook: step_13b_basic.ipynb）

### ライブラリと設定

```python
import numpy as np           # 数値計算ライブラリ（配列操作・数学関数に使用）
import pandas as pd          # データフレーム操作ライブラリ（表形式データの読み込み・加工に使用）
import matplotlib.pyplot as plt  # グラフ描画ライブラリ（PCAプロット等の作成に使用）
from matplotlib.patches import Ellipse  # 楕円パッチ（PCAプロットに95%信頼楕円を描画するため）
import matplotlib.transforms as transforms  # 座標変換（楕円の回転・拡縮・移動に使用）
import seaborn as sns        # 統計的可視化ライブラリ（ヒートマップ等の高水準プロットに使用）
from sklearn.decomposition import PCA  # 主成分分析（PCA）クラス（次元削減による群分離の可視化に使用）

# --- パス ---
RESULTS  = "../results"              # 解析結果の保存先ディレクトリへのパス
FIG_DIR  = f"{RESULTS}/figures"      # 図の保存先ディレクトリへのパス

# --- カラー ---
NORMAL, TUMOR = "#3498DB", "#E74C3C"  # Normal群を青、Tumor群を赤で表示する色コード
UP, DOWN, NS  = "#E74C3C", "#3498DB", "#CCCCCC"  # 上昇(赤)・低下(青)・非有意(灰)の色コード
# 条件名から色への対応辞書（グラフ描画時に群ごとの色を自動で割り当てるため）
COND_MAP = {"Normal": NORMAL, "Tumor": TUMOR}
```

### データ読み込み

```python
# 前処理済みデータと統計解析結果を読み込む
# OpenMS + AlphaPeptDeepで前処理済みタンパク質発現データを読み込む（行=タンパク質、列=サンプル、値=log2発現量）
df = pd.read_csv(f"{RESULTS}/preprocessed_data_openms.csv", index_col=0)
# サンプル情報（各サンプルがNormalかTumorか）を読み込む
sample_info = pd.read_csv(f"{RESULTS}/sample_info_openms.csv")
# サンプル名をインデックスにして条件列だけを取り出す（後でサンプルから条件を引くための辞書的Series）
conditions = sample_info.set_index("Sample")["Condition"]
# 前回のt検定結果を読み込む（有意差タンパク質の特定に使用）
result_df = pd.read_csv(f"{RESULTS}/differential_proteins_openms.csv")

print(f"Total proteins: {len(df)}, Differential proteins: {(result_df['Significant'] != 'NS').sum()}")
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
# グラフタイトル: 深層学習DIA解析の全有意差タンパク質を使ったPCA
ax.set_title("Deep Learning DIA: PCA using All Differentially Abundant Proteins")
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

## コード詳細

### ヒートマップ可視化パラメータ

| パラメータ | 設定値 | 意味 |
|-----------|--------|------|
| `method="ward"` | Ward法 | 分散最小化によるクラスタリング |
| `cmap="RdBu_r"` | 反転Red-Blue | 赤=高発現、青=低発現 |
| `z_score=1` | 列方向標準化 | タンパク質ごとに平均0、SD=1に正規化 |
| `vmin=-3, vmax=3` | 色範囲制限 | 外れ値による色の偏りを防止 |

### PCA解析設定

| パラメータ | 意味 |
|-----------|------|
| `n_components=2` | 2次元に次元削減 |
| `confidence_ellipse(n_std=2.0)` | 95%信頼楕円（約2標準偏差） |
| `explained_variance_ratio_` | 各主成分の寄与率（分散説明率） |

## 実行方法

**スクリプトで一括実行する場合:**
```bash
python scripts/step_13b_openms_basic_visualization.py
```

**Notebook でセルごとに実行する場合:**
```bash
jupyter notebook notebooks/step_13b_basic.ipynb
```

## 可視化結果

### 有意差タンパク質のヒートマップとPCA

![深層学習DIA 全有意差タンパク質ヒートマップ](images/fig2a_heatmap_all.png)

深層学習により検出された有意差タンパク質について、発現量のヒートマップを作成しました。**従来手法（Sage: 1,055個）を大幅に上回る検出数**により、より詳細な発現パターンの階層構造が観察されます。

![深層学習DIA 全有意差タンパク質PCA](images/fig2b_pca_all.png)

深層学習による高精度検出タンパク質のPCA解析です。豊富なマーカー候補により、Normal/Tumor間の分離がより鮮明になっています。

## まとめ

OpenMS + AlphaPeptDeepによる深層学習DIA解析で、全有意差タンパク質を用いた基本的な可視化解析を実行しました。**深層学習による高精度検出**により、従来の理論スペクトル検索では見落とされていた低発現タンパク質や微細な発現変動も捉えることができています。

次回では、発現変化の大きさに応じたTop N解析により、より詳細な群分離能力の評価を行います。

> 前回: [#13a 深層学習DIA差分発現統計解析](article-13a-openms-differential-stats.md)
> 次回: [#13c 深層学習DIAトップN解析](article-13c-openms-differential-topn.md) — Top50/100/200解析

#バイオインフォマティクス #プロテオミクス #深層学習 #可視化解析 #PCA #ヒートマップ #labcode