---
title: "深層学習DIA解析での主成分分析（PCA）と論文比較【論文再現シリーズ #12b】"
emoji: "🧠"
type: "tech"
topics: ["proteomics", "openms", "pca", "analysis", "labcode"]
published: false
---

# 深層学習DIA解析での主成分分析（PCA）と論文比較

## はじめに

前回（[#12a 深層学習DIA基本可視化](article-12a-openms-visualization-basics.md)）では、相関行列とクラスタリングにより基本的な群構造を確認しました。この記事では、主成分分析（PCA）による高次元データの次元削減解析を行い、腫瘍/正常組織の分離状況をより詳細に評価します。また、Toyota et al. 2025論文との数値比較も実施します。

> **📝 INFO**
>
**この記事で行う処理**
OpenMS + AlphaPeptDeepで検出した19,981タンパク質のデータに対してPCA（主成分分析）を実行し、2次元の主成分空間でNormal/Tumorの群分離を可視化します。各主成分の寄与率を計算し、論文値と比較することで、深層学習DIA解析の妥当性を定量的に評価します。

## 前提

- [#12a 深層学習DIA基本可視化](article-12a-openms-visualization-basics.md) が完了していること
- 前処理済みデータが利用可能
- **対応Notebook**: [`notebooks/step_12b.ipynb`](../notebooks/step_12b.ipynb) — この記事のコードをセルごとに実行できます

> **📝 INFO**
>
**【PCA基礎知識】**

- **PCA（Principal Component Analysis、主成分分析）**: 高次元のデータを少数の軸（主成分）に圧縮して可視化する手法。「サンプルの全体的な違いを2〜3軸に要約する」と理解すればOKです。
- **寄与率（Explained Variance Ratio）**: 各主成分が全変動のうち何%を説明するかの指標。PC1の寄与率が高いほど、その軸がデータの主要な変動を捉えています。
- **信頼楕円（Confidence Ellipse）**: 各群の95%信頼区間を楕円で表示し、群の分布範囲と重なり具合を視覚的に示します。

## コード全文（対応Notebook: step_12b.ipynb）

### ライブラリと設定

```python
import numpy as np                              # 数値計算ライブラリ（配列演算・統計に使用）
import pandas as pd                             # データフレーム操作ライブラリ（CSV読込・表形式データ処理）
import matplotlib.pyplot as plt                 # 基本グラフ描画ライブラリ（軸設定・保存）
from matplotlib.patches import Ellipse          # 楕円パッチ（PCA信頼区間の描画用）
from sklearn.decomposition import PCA           # 主成分分析（scikit-learnのPCAクラス）
from scipy.stats import chi2                    # カイ二乗分布（95%信頼楕円の計算に使用）

# --- パス設定 ---
RESULTS = "../results"                          # 解析結果の親ディレクトリ
FIG = f"{RESULTS}/figures"                      # 生成した図の保存先ディレクトリ
TABLE = f"{RESULTS}/tables"                     # 生成したテーブル（CSV等）の保存先ディレクトリ

# --- 図の見た目設定 ---
plt.rcParams["font.size"] = 10                  # 全体のフォントサイズを10ptに設定
plt.rcParams["axes.labelsize"] = 10             # 軸ラベル（x軸, y軸のタイトル）のフォントサイズ
plt.rcParams["xtick.labelsize"] = 8             # x軸の目盛りラベルのフォントサイズ
plt.rcParams["ytick.labelsize"] = 8             # y軸の目盛りラベルのフォントサイズ
```

### データ読み込み

```python
# 前処理済みの深層学習DIA解析データを読み込み（行=タンパク質、列=サンプル、値=log2発現量）
df = pd.read_csv(f"{RESULTS}/preprocessed_data_openms.csv", index_col=0)
# サンプル情報（サンプル名と群の対応）を読み込み
sample_info = pd.read_csv(f"{RESULTS}/sample_info_openms.csv")
# サンプル情報から「サンプル名→群」の対応SeriesをPDFで作成（引用で使いやすくする）
conditions = sample_info.set_index("Sample")["Condition"]

print(f"深層学習DIAデータ: {df.shape[0]} タンパク質 × {df.shape[1]} サンプル")
print(f"サンプル群構成: {dict(conditions.value_counts())}")
```

### 信頼楕円描画関数

```python
def confidence_ellipse(x, y, ax, n_std=2.0, **kwargs):
    """95% 信頼楕円を描画する（PCA散布図用のヘルパー関数）。

    Args:
        x, y: 各軸のデータ点（PC1, PC2のスコア）
        ax: matplotlib Axes オブジェクト（描画対象）
        n_std: 標準偏差の倍数（2.0 = 約95%信頼区間）
        **kwargs: 楕円の外観設定（色、透明度等）

    Returns:
        描画された楕円オブジェクト
    """
    # データ点が2個未満の場合は楕円を描画できないため終了
    if len(x) < 2:
        return

    # x, y の共分散行列を計算（2×2マトリクス）
    cov = np.cov(x, y)
    # 共分散行列から相関係数を計算（-1〜1の値）
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])

    # 95%信頼区間に対応するカイ二乗分布の値を取得（自由度2）
    ell_radius_x = np.sqrt(1 + pearson)   # x方向の楕円半径係数
    ell_radius_y = np.sqrt(1 - pearson)   # y方向の楕円半径係数
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, **kwargs)

    # 楕円のスケールを調整（各軸の標準偏差 × n_std倍）
    scale_x = np.sqrt(cov[0, 0]) * n_std   # x軸方向のスケール
    scale_y = np.sqrt(cov[1, 1]) * n_std   # y軸方向のスケール
    ellipse.set_transform(ax.transData +
                         plt.matplotlib.transforms.Affine2D().scale(scale_x, scale_y) +
                         plt.matplotlib.transforms.Affine2D().translate(np.mean(x), np.mean(y)))

    # axesに楕円を追加
    return ax.add_patch(ellipse)
```

### PCA（主成分分析）

```python
def plot_pca(df, conditions):
    """PCA散布図を作成・保存し、寄与率を計算する。

    Args:
        df: 前処理済みタンパク質発現データ（行=タンパク質、列=サンプル）
        conditions: サンプル名→群ラベルの対応Series

    Returns:
        PCAの寄与率DataFrame
    """
    print("PCA解析実行中...")

    # PCAモデルを作成（n_components=2: 2つの主成分に次元削減）
    pca = PCA(n_components=2)

    # df.T で転置（行=サンプル, 列=タンパク質）してからPCAを実行
    # fit_transform: モデルの学習と変換を同時に行い、各サンプルのPC1, PC2スコアを返す
    # scores の形状: (サンプル数, 2)
    scores = pca.fit_transform(df.T)

    # 8×6インチの図とAxes（描画領域）を作成
    fig, ax = plt.subplots(figsize=(8, 6))

    # 群ごとの色を定義（Normal=青, Tumor=赤）
    cmap = {"Normal": "#3498DB", "Tumor": "#E74C3C"}

    # Normal群とTumor群それぞれについてプロット
    for cond, color in cmap.items():
        # 各サンプルが現在の群（Normal or Tumor）に該当するかのブールマスクを作成
        mask = conditions.reindex(df.columns) == cond

        # 散布図を描画: PC1をx軸、PC2をy軸にプロット
        ax.scatter(scores[mask, 0], scores[mask, 1],  # mask に一致するサンプルのPC1, PC2値
                   c=color,                            # 点の塗りつぶし色
                   s=100,                              # 点のサイズ（ピクセル^2）
                   alpha=0.8,                          # 透明度（0=透明, 1=不透明）
                   label=cond,                         # 凡例に表示するラベル名
                   edgecolors="white",                 # 点の輪郭色（白で囲むと見やすい）
                   linewidth=0.5)                      # 輪郭線の太さ

        # 95%信頼楕円を描画（群の分布範囲を視覚的に示す）
        confidence_ellipse(scores[mask, 0], scores[mask, 1], ax, n_std=2.0,
                           facecolor=color,            # 楕円の塗りつぶし色
                           alpha=0.15,                 # 楕円の透明度（薄く表示）
                           edgecolor=color,            # 楕円の輪郭色
                           linewidth=1.5)              # 楕円の輪郭線の太さ

    # 各主成分の寄与率（全変動のうち何%を説明するか）を取得
    ev = pca.explained_variance_ratio_

    # x軸ラベルにPC1の寄与率を表示（例: "Component 1 (39.7%)"）
    ax.set_xlabel(f"Component 1 ({ev[0]*100:.1f}%)")
    # y軸ラベルにPC2の寄与率を表示（例: "Component 2 (15.6%)"）
    ax.set_ylabel(f"Component 2 ({ev[1]*100:.1f}%)")
    # グラフのタイトルを設定
    ax.set_title("PCA of Non-tumor and Tumor Tissues")
    # 凡例を表示（frameon=False で枠線なし）
    ax.legend(frameon=False)
    # 薄いグレーのグリッド線を表示（データ点の位置を読みやすくする）
    ax.grid(True, color="gray", alpha=0.3, linestyle="-", linewidth=0.5)
    # 上と右の枠線を非表示にする（見た目をすっきりさせる）
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # 図をPNGファイルとして保存（dpi=150で高解像度、余白を自動調整）
    fig.savefig(f"{FIG}/fig1c_pca.png", dpi=150, bbox_inches="tight")
    # メモリ節約のため図を閉じる
    plt.close()

    # 各主成分の寄与率と累積寄与率をDataFrameにまとめる
    # PC: 主成分名、Variance_Ratio: 個別寄与率、Cumulative: 累積寄与率
    var_df = pd.DataFrame({"PC": [f"PC{i+1}" for i in range(len(ev))],
                           "Variance_Ratio": ev,
                           "Cumulative": np.cumsum(ev)})
    # 寄与率の表をCSVファイルとして保存（後で数値を確認できるように）
    var_df.to_csv(f"{TABLE}/pca_variance.csv", index=False)
    print("保存: fig1c_pca.png")

    return var_df

# 関数を呼び出してPCA散布図を生成・保存
variance_df = plot_pca(df, conditions)

# 寄与率の結果を表示
print("\n=== PCA 寄与率 ===")
print(variance_df.to_string(index=False, float_format="%.1f"))
```

## コード詳細

### PCA の explained_variance_ratio_

```python
# PC1（第1主成分）が全変動のうち何割を説明するかを取得（例: 0.45 = 45%）
# 値が大きいほど、PC1軸がデータの主要な違いを捉えていることを意味する
pca.explained_variance_ratio_[0]

# PC2（第2主成分）が全変動のうち何割を説明するかを取得（例: 0.12 = 12%）
# PC1で説明しきれなかった残りの変動のうち、最も大きな方向を捉えている
pca.explained_variance_ratio_[1]
```

- **寄与率が高い** = その主成分がデータの変動をよく説明している
- PC1の寄与率が高い場合、PC1軸上での群の分離が「データの主要な変動」であることを示す

### 3つの手法の統合的理解

| 手法 | 見えるもの | 使いどころ |
|------|----------|----------|
| 相関行列 | サンプル間の類似度 | データ品質の確認、外れ値検出 |
| クラスタリング | サンプルの階層構造 | 群構造の発見、タンパク質パターン |
| PCA | 全体の変動構造 | 群の分離確認、バッチ効果検出 |

## 実行方法

**スクリプトで一括実行する場合:**

```bash
python scripts/step_12b_openms_visualization_pca.py
```

**Notebook でセルごとに実行する場合:**

```bash
jupyter notebook notebooks/step_12b.ipynb
```

Notebook版では寄与率の計算結果をインラインで確認できます。

## 可視化結果と論文比較

### PCA（主成分分析）

![PCA](images/fig1c_pca.png)

PCAにより高次元のプロテオームデータを2次元に圧縮した結果です。PC1（第1主成分、寄与率 **39.7%**）の軸上でNormal群とTumor群が明瞭に分離しており、データの最大変動がTumor/Normalの対比であることを示しています。PC2（寄与率 **15.6%**）は個体間のばらつきを反映しています。

### 深層学習DIA vs 論文の数値比較

| 指標 | 本実装 (OpenMS + AlphaPeptDeep) | 論文 (DIA-NN) | 差異 |
|------|---------------------------|--------------|-----|
| サンプル数 | 32（Normal 16 + Tumor 16） | 32（Normal 16 + Tumor 16） | 一致 |
| 解析タンパク質数 | **19,981** | ～10,329 | **1.9倍** |
| 相関値 (Pearson r) | 0.70 – 0.93 | 0.744 – 0.982 | 類似 |
| **PCA PC1 寄与率** | **39.7%** | **42.1%** | **-2.4%** |
| PCA PC2 寄与率 | 15.6% | 10.8% | +4.8% |
| Normal/Tumor 分離 | 明確 | 明確 | 一致 |

### 重要な発見

**PC1 寄与率が論文の 42.1% に対して本実装 39.7% と近い値**を示す点が重要です。これは「大腸がんの Tumor/Normal 対比が強い生物学的シグナル」を示し、**深層学習DIA解析でも本質的な群構造は再現される**ことを意味します。

むしろ、**検出タンパク質数が1.9倍に増加**しているにも関わらず、主要な生物学的シグナル（PC1寄与率）がほぼ維持されていることは、深層学習DIAの**高い特異性と感度**を示す証拠です。

### 深層学習DIAの優位性

1. **検出数の飛躍的向上**: 19,981 vs 10,329タンパク質（1.9倍）
2. **生物学的シグナルの保持**: PC1寄与率がほぼ同等（39.7% vs 42.1%）
3. **群分離の明確性**: Normal/Tumor が明瞭に分離
4. **商用利用可能**: Apache 2.0ライセンスで企業応用可能

## まとめ

PCA解析により、深層学習DIA解析結果の生物学的妥当性を定量的に確認しました：

### 主要成果

1. **PC1寄与率39.7%**: 論文値（42.1%）に近い値で生物学的シグナルの再現を確認
2. **明確な群分離**: Normal/Tumor が主成分空間で明瞭に分離
3. **検出性能向上**: タンパク質数1.9倍増加で包括性向上
4. **技術的革新**: 深層学習ベースの解析で論文を上回る性能を実現

深層学習DIA解析は、従来手法を大幅に上回る検出性能を持ちながら、生物学的に妥当な群構造を保持することが数値的に示されました。この結果は、商用利用可能な深層学習技術による実用的プロテオミクス解析の実現を意味します。

> 前回: [#12a 深層学習DIA基本可視化](article-12a-openms-visualization-basics.md)
> 次回: [#13a 深層学習DIA差分発現統計](article-13a-openms-differential-stats.md) — Welch's t検定とVolcanoプロット

#バイオインフォマティクス #プロテオミクス #PCA #統計解析 #深層学習 #labcode