---
title: "差分発現トップN解析実行【論文再現シリーズ #8c】"
emoji: "🏆"
type: "tech"
topics: ["proteomics", "python", "biomarkers", "labcode"]
published: false
---

# 差分発現トップN解析実行

## はじめに

前回（[#8b 差分発現タンパク質の基本可視化](article-08b-differential-clustering.md)）では、1,055個の有意差タンパク質全体でのクラスタリングとPCAを実行しました。この記事では、**最も有意差の大きいTop 50/100/200タンパク質**だけでも群分離が可能かを検証するためのTop N解析を実行します。

:::message
**この記事で行う処理**
有意差の大きい順にTop 50/100/200個のタンパク質を抽出し、段階的にクラスタリングとPCAを実行します。少数の強力なマーカーだけでも完全な群分離が可能かを検証し、実用的なバイオマーカーパネルの最小サイズを特定します。
:::

## 前提

- [#8b 差分発現タンパク質の基本可視化](article-08b-differential-clustering.md) が完了していること
- 差分発現タンパク質データ（`differential_proteins.csv`）が生成済み
- **対応Notebook**: [`notebooks/step_08.ipynb`](../notebooks/step_08.ipynb) — この記事のコードをセルごとに実行できます

## コード全文（対応Notebook: step_08.ipynb）

### ライブラリと設定

```python
import numpy as np             # 数値計算ライブラリ（配列操作・数学関数に使用）
import pandas as pd            # データフレーム操作ライブラリ（表形式データの読み込み・加工に使用）
import matplotlib.pyplot as plt  # グラフ描画ライブラリ（PCAプロット等の作成に使用）
from matplotlib.patches import Ellipse  # 楕円パッチ（PCAプロットに95%信頼楕円を描画するため）
import matplotlib.transforms as transforms  # 座標変換（楕円の回転・拡縮・移動に使用）
import seaborn as sns          # 統計的可視化ライブラリ（ヒートマップ等の高水準プロットに使用）
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
# 前処理済みデータの読み込み
df = pd.read_csv(f"{RESULTS}/preprocessed_data.csv", index_col=0)  # Log2変換・正規化済みのタンパク質発現データ
result_df = pd.read_csv(f"{RESULTS}/differential_proteins.csv")   # t検定結果（統計値・有意性・Log2FC含む）

# サンプル条件の取得（サンプル名から患者IDと条件を抽出）
conditions = pd.Series({col: "Normal" if "-N" in col else "Tumor" for col in df.columns})
print(f"データ形状: {df.shape[0]} タンパク質 × {df.shape[1]} サンプル")  # 読み込んだデータの規模を確認
print(f"有意差タンパク質数: {sum(result_df['Significant'] != 'NS')}")       # 有意差ありと判定されたタンパク質数
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

### Top N 解析メイン関数

```python
def create_topn_analysis(df, result_df, conditions):
    """Top Nタンパク質による段階的な群分離解析を実行する。

    【解析のポイント】
    1. Up/Down両方向から均等に選択（生物学的バランス確保）
    2. 絶対値Log2FCによるランキング（変化量の大きさで評価）
    3. 段階的解析（50→100→200で性能変化を評価）
    """
    # 有意差タンパク質の抽出と準備
    sig = result_df[result_df["Significant"] != "NS"].copy()  # 非有意（NS）を除外
    sig["AbsLog2FC"] = sig["Log2FC"].abs()  # 絶対値Log2FC（変化量の大きさ）を計算

    # Up（腫瘍で増加）とDown（腫瘍で減少）に分離してソート
    sig_up = sig.query("Significant == 'Up'").sort_values("AbsLog2FC", ascending=False)    # 増加タンパク質を変化量順
    sig_down = sig.query("Significant == 'Down'").sort_values("AbsLog2FC", ascending=False)  # 減少タンパク質を変化量順

    print(f"Up-regulated: {len(sig_up)} proteins")     # 腫瘍で増加するタンパク質数
    print(f"Down-regulated: {len(sig_down)} proteins")  # 腫瘍で減少するタンパク質数

    # Top N解析の実行（50, 100, 200段階）
    for n in [50, 100, 200]:
        print(f"\n=== Top {n} Analysis ===")

        # Top Nタンパク質の選択（Up/Downから均等に選択）
        # 各方向から最大n個まで選択し、不足分は他方向で補完
        top_up = min(n, len(sig_up))      # Up方向から選択する数（最大n個）
        top_down = min(n, len(sig_down))   # Down方向から選択する数（最大n個）

        # 実際のTop Nタンパク質リストを作成
        top_proteins = np.concatenate([
            sig_up.head(top_up)["Protein"].values,    # Up方向からtop_up個
            sig_down.head(top_down)["Protein"].values  # Down方向からtop_down個
        ])

        # 選択されたタンパク質の発現データを抽出
        top_data = df.loc[df.index.isin(top_proteins)]

        print(f"Selected proteins: {len(top_data)} (Up: {top_up}, Down: {top_down})")

        # 1. 階層的クラスタリングによる群分離確認
        # サンプル色の設定（Normal/Tumorで色分け）
        sample_colors = conditions.reindex(df.columns).map(COND_MAP)

        # ClusterMapの作成（階層的クラスタリング+ヒートマップ）
        g = sns.clustermap(
            top_data.T,              # 転置してサンプル×タンパク質の形にする（行がサンプル）
            method="ward",           # Ward法（群内分散を最小化するクラスタリング手法）
            cmap="RdBu_r",           # 赤青カラーマップ（赤=高発現、青=低発現）
            z_score=1,               # 列方向（タンパク質ごと）にZ-score標準化
            row_colors=sample_colors, # 行（サンプル）の色バー（Normal/Tumor識別用）
            figsize=(8, 10),         # 図のサイズ（Top N解析用に最適化）
            xticklabels=False,       # x軸ラベル（タンパク質名）は非表示（多すぎるため）
            yticklabels=True,        # y軸ラベル（サンプル名）は表示
            vmin=-2, vmax=2,         # カラーバーの範囲（Z-score標準化後の標準的範囲）
        )

        # y軸のサンプル名フォントサイズを調整
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=8)

        # 図を保存
        if n == 50:
            figname = "fig2c_clustering_top50.png"   # Top 50クラスタリング
        elif n == 100:
            figname = "fig2d_clustering_top100.png"  # Top 100クラスタリング
        else:
            figname = "fig2e_clustering_top200.png"  # Top 200クラスタリング

        g.savefig(f"{FIG_DIR}/{figname}", dpi=150, bbox_inches="tight")
        plt.show()

        # 2. PCA による群分離の定量評価
        pca = PCA(n_components=2)  # 2次元に次元削減
        scores = pca.fit_transform(top_data.T)  # PCAを実行（行=サンプル）

        # PCAプロット作成
        fig, ax = plt.subplots(figsize=(8, 6))

        # 各群を異なる色でプロット（信頼楕円付き）
        for cond, color in [("Normal", NORMAL), ("Tumor", TUMOR)]:
            mask = (conditions.reindex(df.columns) == cond).values  # 各条件のサンプルマスク
            ax.scatter(scores[mask, 0], scores[mask, 1], c=color, s=80, alpha=0.8,
                      label=cond, edgecolors="white")  # 散布図プロット

            # 95%信頼楕円の追加（群の分布範囲を視覚化）
            confidence_ellipse(scores[mask, 0], scores[mask, 1], ax,
                               facecolor=color, alpha=0.15, edgecolor=color, lw=1.5)

        # 軸ラベル・タイトル・凡例の設定
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")  # 第1主成分と寄与率
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")  # 第2主成分と寄与率
        ax.set_title(f"PCA using Top {n} Differentially Abundant Proteins")  # グラフタイトル
        ax.legend(frameon=False)  # 凡例（枠なし）
        ax.spines[["top", "right"]].set_visible(False)  # 上・右枠線を非表示
        ax.grid(True, alpha=0.3, ls="--")  # 薄いグリッド線

        # 図を保存
        if n == 50:
            pca_figname = "fig2f_pca_top50.png"     # Top 50 PCA
        elif n == 100:
            pca_figname = "fig2g_pca_top100.png"    # Top 100 PCA
        else:
            pca_figname = "fig2h_pca_top200.png"    # Top 200 PCA

        fig.savefig(f"{FIG_DIR}/{pca_figname}", dpi=150, bbox_inches="tight")
        plt.show()

        # 基本性能の表示
        print(f"PC1 variance explained: {pca.explained_variance_ratio_[0]*100:.1f}%")
        print(f"PC2 variance explained: {pca.explained_variance_ratio_[1]*100:.1f}%")

# Top N解析の実行
create_topn_analysis(df, result_df, conditions)
```

## 実行方法

**スクリプトで一括実行する場合:**

```bash
python scripts/step_08c_differential_topn.py
```

**Notebook でセルごとに実行する場合:**

```bash
jupyter notebook notebooks/step_08.ipynb
```

## 基本結果

Top 50/100/200タンパク質による段階的解析を実行し、以下の図を生成しました：

### Top 50 タンパク質
- **クラスタリング**: `fig2c_clustering_top50.png`
- **PCA**: `fig2f_pca_top50.png`

### Top 100 タンパク質
- **クラスタリング**: `fig2d_clustering_top100.png`
- **PCA**: `fig2g_pca_top100.png`

### Top 200 タンパク質
- **クラスタリング**: `fig2e_clustering_top200.png`
- **PCA**: `fig2h_pca_top200.png`

## コード詳細

### Top N選択アルゴリズム

| ステップ | 処理内容 | 生物学的意義 |
|---------|---------|------------|
| `sig["AbsLog2FC"] = sig["Log2FC"].abs()` | 絶対値での変化量計算 | 変化の方向によらず影響度で評価 |
| `sig_up.head(n)` | Up群から上位N個選択 | 腫瘍で増加する最強マーカー |
| `sig_down.head(n)` | Down群から上位N個選択 | 腫瘍で減少する最強マーカー |
| `np.concatenate()` | Up/Down均等に統合 | 生物学的バランスの確保 |

### クラスタリングパラメータ

| パラメータ | 設定値 | 効果 |
|-----------|--------|------|
| `method="ward"` | Ward法 | クラスター内分散最小化（コンパクトなクラスター形成） |
| `z_score=1` | 列方向標準化 | タンパク質間の発現レベル差を正規化 |
| `figsize=(8, 10)` | ヒートマップサイズ | Top N解析に最適化されたアスペクト比 |

## まとめ

Top 50/100/200タンパク質による段階的解析の実行を完了しました。主要な成果は以下の通りです：

### 実行成果

1. **段階的解析**: 3段階（50/100/200）での系統的な群分離評価
2. **バランス設計**: Up/Down両方向からの均等選択を実現
3. **可視化完了**: 各段階でクラスタリングとPCAを生成
4. **定量準備**: PCA寄与率等の基本性能指標を取得

この Top N解析により、少数の強力なマーカーでの群分離性能を検証する準備が整いました。

次回では、これらの結果を定量的に評価し、バイオマーカーパネルとしての実用性を詳細に検討します。

> 前回: [#8b 差分発現タンパク質の基本可視化](article-08b-differential-clustering.md)
> 次回: [#8d トップNバイオマーカー評価](article-08d-differential-topn-evaluation.md) — 性能評価・実用性検討

#バイオインフォマティクス #プロテオミクス #バイオマーカー #TopN解析 #labcode