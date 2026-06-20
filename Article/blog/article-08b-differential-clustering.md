---
title: "差分発現タンパク質の基本可視化【論文再現シリーズ #8b】"
emoji: "🎯"
type: "tech"
topics: ["proteomics", "python", "clustering", "labcode"]
published: false
---

# 差分発現タンパク質の基本可視化

## はじめに

前回（[#8a Volcanoプロット](article-08a-differential-volcano.md)）で特定した1,055個の有意差タンパク質を使って、**階層的クラスタリング**と**主成分分析（PCA）**を実行します。有意差タンパク質だけに絞ることで、Normal/Tumor群の分離パターンがより明瞭になることを確認します。

:::message
**この記事で行う処理**
有意差タンパク質全体（1,055個）でヒートマップとPCAを作成し、Normal/Tumor群の分離パターンを可視化します。階層的クラスタリングにより類似発現パターンのタンパク質をグルーピングし、PCAで群分離の質を定量的に評価します。
:::

## 前提

- [#8a Volcanoプロット](article-08a-differential-volcano.md) が完了していること
- 差分発現タンパク質データ（`differential_proteins.csv`）が生成済み
- **対応Notebook**: [`notebooks/step_08.ipynb`](../notebooks/step_08.ipynb) — この記事のコードをセルごとに実行できます

## コード全文（対応Notebook: step_08.ipynb）

### ライブラリと設定（追加）

```python
import matplotlib.pyplot as plt  # グラフ描画ライブラリ（PCAプロット等の作成に使用）
from matplotlib.patches import Ellipse  # 楕円パッチ（PCAプロットに95%信頼楕円を描画するため）
import matplotlib.transforms as transforms  # 座標変換（楕円の回転・拡縮・移動に使用）
import seaborn as sns        # 統計的可視化ライブラリ（ヒートマップ等の高水準プロットに使用）
from sklearn.decomposition import PCA  # 主成分分析（PCA）クラス（次元削減による群分離の可視化に使用）
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


## 実行方法

**スクリプトで一括実行する場合:**

```bash
python scripts/step_07_differential_expression.py
```

**Notebook でセルごとに実行する場合:**

```bash
jupyter notebook notebooks/step_08.ipynb
```

## 可視化結果

### 有意差タンパク質のヒートマップとPCA

![全有意差タンパク質ヒートマップ](images/fig2a_heatmap_all.png)

有意差のある1,055タンパク質全てについて、発現量のヒートマップを作成しました。階層的クラスタリング（Ward法）により類似したパターンのタンパク質とサンプルが自動的にグルーピングされています。**Normal群とTumor群で明瞭に異なる発現パターン**が観察され、有意差タンパク質に絞った効果が確認できます。

![全有意差タンパク質PCA](images/fig2b_pca_all.png)

有意差タンパク質のみを使ったPCA解析です。**第1主成分と第2主成分で群が明確に分離**されており、Normal/Tumor間の分子的差異が主成分の方向として捉えられています。有意差タンパク質に絞ることで、ノイズが除去され群分離がより鮮明になっています。

## コード詳細

### sns.clustermap() のパラメータ

| パラメータ | 値 | 意味 |
|-----------|-----|------|
| `method` | `"ward"` | Ward法を使用（分散を最小化） |
| `cmap` | `"RdBu_r"` | 赤-青の逆順カラーマップ |
| `z_score` | `1` | 列方向（タンパク質ごと）に標準化 |
| `row_colors` | サンプル色 | Normal/Tumorを色バーで表示 |
| `col_colors` | タンパク質色 | Up/Downを色バーで表示 |

### PCAの解釈

- **第1主成分（PC1）**: 最大の分散方向。通常、群間の主要な差異を捉える
- **第2主成分（PC2）**: 第2の分散方向。群内のサブ構造を捉えることが多い
- **寄与率**: 各主成分が説明する全分散の割合（%表示）

## まとめ

有意差タンパク質1,055個について、**階層的クラスタリングとPCA解析**による可視化を実行しました。主要成果は以下の通りです：

### 可視化解析成果

1. **群分離の確認**: ヒートマップ・PCAともにNormal/Tumor群の明瞭な分離を確認
2. **発現パターンの構造化**: 階層的クラスタリングにより類似発現パターンの自動グルーピング
3. **主成分の解釈**: 第1・第2主成分が群間差異を効果的に捉えることを確認
4. **有意差絞り込みの効果**: 全タンパク質より群分離が鮮明になることを確認

これらの結果により、有意差タンパク質が大腸がんの分子的特徴を適切に反映していることが示されました。

次回では、これらの有意差タンパク質の中から最も変化の大きいタンパク質群を段階的に選択し、バイオマーカーパネルとしての実用性を評価します。

> 前回: [#8a Volcanoプロット](article-08a-differential-volcano.md)
> 次回: [#8c トップN解析実行](article-08c-differential-topn-analysis.md) — バイオマーカー段階的解析

#バイオインフォマティクス #プロテオミクス #クラスタリング #PCA #labcode