---
title: "Pythonで論文のFigure 1を再現する：相関行列・クラスタリング・PCA【論文再現シリーズ #6】"
emoji: "📊"
type: "tech"
topics: ["proteomics", "python", "visualization", "labcode"]
published: false
---

# Pythonで論文のFigure 1を再現する

## はじめに

この記事では、Toyota et al. 2025 のFigure 1（プロテオームの全体像）をPythonで再現します。腫瘍組織と正常組織のタンパク質発現パターンが異なることを、3つの異なる可視化手法で確認します。

## 前提

- [#5 前処理](article-05-preprocess.md) が完了していること
- `results/preprocessed_data.csv` が存在

## 再現するFigure

| パネル | 内容 | 手法 |
|--------|------|------|
| Figure 1a | 相関行列ヒートマップ | ピアソン相関 + クラスターマップ |
| Figure 1b | 階層的クラスタリング | Ward法 + ヒートマップ |
| Figure 1c | PCA | 主成分分析 |

## スクリプト全文: step_03_overview_visualization.py

> **注意:** スクリプトが長いため、ここではコアとなる3つの関数を抜粋して掲載します。完全版は `scripts/step_03_overview_visualization.py` を参照してください。

### (a) 相関行列ヒートマップ

```python
def plot_correlation_matrix(df, sample_info):
    """
    Figure 1a: サンプル間の相関行列をヒートマップとして可視化する関数。

    【相関行列とは？】
      各サンプル同士の「発現パターンの類似度」をピアソン相関係数で数値化したもの。
      - 相関係数 = 1.0: 完全に同じパターン
      - 相関係数 = 0.0: 全く関係がない
      NormalとTumorが異なる発現パターンを持つなら、群間の相関が群内より低くなる。
    """
    print("  (a) 相関行列ヒートマップ")

    # ピアソン相関係数の計算
    # df.corr() はデフォルトでピアソン相関係数を計算する
    # 結果はサンプル数 × サンプル数の正方行列になる
    corr = df.corr(method="pearson")

    # サンプルの色分け（Normal → 青, Tumor → 赤）
    conditions = sample_info.set_index("Sample")["Condition"]
    colors = conditions.map({"Normal": "#4EAED1", "Tumor": "#E8524A"})
    row_colors = colors.reindex(corr.index)

    # sns.clustermap() は seaborn のヒートマップ＋クラスタリング関数
    g = sns.clustermap(
        corr,                   # 描画するデータ（相関行列）
        method="average",       # 【クラスタリング手法】UPGMA法（群平均法）
        metric="correlation",   # 【距離指標】1 - ピアソン相関係数
        cmap="RdBu_r",          # 【カラーマップ】赤=高い相関、青=低い相関
        vmin=0.85, vmax=1.0,    # 色の範囲（0.85〜1.0の微妙な差を可視化）
        figsize=(10, 10),       # 図のサイズ（横10インチ × 縦10インチ）
        row_colors=row_colors,  # 行の左に表示するカラーバー
        col_colors=row_colors,  # 列の上に表示するカラーバー
    )

    filepath = os.path.join(FIG_DIR, "fig1a_correlation.png")
    g.savefig(filepath, dpi=150, bbox_inches="tight")
    plt.close()
```

### (b) 階層的クラスタリング

```python
def plot_hierarchical_clustering(df, sample_info):
    """
    Figure 1b: 全タンパク質の発現量を使って、サンプルをクラスタリングする。

    【階層的クラスタリングとは？】
      似たサンプルを自動でツリー状に分類する手法。
      教師なし（ラベルなし）でデータの構造を発見できる。
    """
    print("  (b) 階層的クラスタリング")

    conditions = sample_info.set_index("Sample")["Condition"]
    sample_colors = conditions.map({"Normal": "#4EAED1", "Tumor": "#E8524A"})

    g = sns.clustermap(
        df.T,                   # 転置: サンプル（行）× タンパク質（列）
        method="ward",          # 【Ward法】分散最小化基準でクラスタを結合
                                #   最もよく使われるクラスタリング手法
        cmap="RdBu_r",          # カラーマップ
        center=0,               # 色の中心値を0にする（正=赤、負=青）
        figsize=(14, 8),
        row_colors=sample_colors.reindex(df.columns),
        z_score=0,              # 【Z-score正規化】行方向で標準化
                                #   各サンプルの平均を0、SDを1にする
        xticklabels=False,      # タンパク質名は多すぎるので非表示
        yticklabels=True,       # サンプル名は表示
    )

    filepath = os.path.join(FIG_DIR, "fig1b_clustering.png")
    g.savefig(filepath, dpi=150, bbox_inches="tight")
    plt.close()
```

### (c) PCA（主成分分析）

```python
def plot_pca(df, sample_info):
    """
    Figure 1c: PCA で高次元データを2次元に圧縮して可視化する。

    【PCA（主成分分析）とは？】
      数千次元のデータ（各タンパク質が1次元）を、情報の損失を最小限にしながら
      2〜3次元に圧縮する手法。PC1（第1主成分）がデータの変動を最も多く説明し、
      PC2（第2主成分）が次に多く説明する。
    """
    print("  (c) 主成分分析（PCA）")

    conditions = sample_info.set_index("Sample")["Condition"]

    # PCA の実行
    # n_components=2: 2次元に圧縮
    pca = PCA(n_components=2)

    # fit_transform: データに PCA を適用し、変換後の座標を返す
    # df.T: 転置して サンプル（行）× タンパク質（列）にする
    scores = pca.fit_transform(df.T)

    # プロット作成
    fig, ax = plt.subplots(figsize=(8, 6))

    for condition, color, marker in [
        ("Normal", "#4EAED1", "o"),   # 青い丸
        ("Tumor", "#E8524A", "s"),    # 赤い四角
    ]:
        mask = conditions.reindex(df.columns) == condition
        ax.scatter(
            scores[mask, 0],        # PC1のスコア（X軸）
            scores[mask, 1],        # PC2のスコア（Y軸）
            c=color, marker=marker,
            s=80, alpha=0.8,        # s=80: 点の大きさ、alpha=0.8: 透明度
            label=condition,
            edgecolors="white",     # 点の縁を白に（見やすくする）
        )

    # 軸ラベルに寄与率を表示
    # explained_variance_ratio_: 各主成分がデータ全体の変動の何%を説明するか
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title("PCA of Non-tumor and Tumor Tissues")
    ax.legend(frameon=False)

    filepath = os.path.join(FIG_DIR, "fig1c_pca.png")
    fig.savefig(filepath, dpi=150, bbox_inches="tight")
    plt.close()
```

---

## コード詳細

### sns.clustermap() のパラメータ解説

| パラメータ | 値 | 意味 |
|-----------|-----|------|
| `method` | `"ward"` | Ward法: 分散最小化でクラスタを結合。最もポピュラー |
| `metric` | `"correlation"` | 距離指標: 1 - ピアソン相関（相関が高いほど近い） |
| `cmap` | `"RdBu_r"` | Red-Blue reversed: 赤=高値、青=低値 |
| `vmin/vmax` | `0.85/1.0` | 色の範囲（この間の差を可視化） |
| `z_score` | `0` | 行方向でZ-score正規化（平均0、SD1に標準化） |
| `row_colors` | 色リスト | 行の横に表示するカラーバー（群の色分け） |

### PCA の explained_variance_ratio_

```python
pca.explained_variance_ratio_[0]  # PC1が説明する変動の割合（例: 0.45 = 45%）
pca.explained_variance_ratio_[1]  # PC2が説明する変動の割合（例: 0.12 = 12%）
```

- **寄与率が高い** = その主成分がデータの変動をよく説明している
- PC1の寄与率が高い場合、PC1軸上での群の分離が「データの主要な変動」であることを示す

### 3つの手法の使い分け

| 手法 | 見えるもの | 使いどころ |
|------|----------|----------|
| 相関行列 | サンプル間の類似度 | データ品質の確認、外れ値検出 |
| クラスタリング | サンプルの階層構造 | 群構造の発見、タンパク質パターン |
| PCA | 全体の変動構造 | 群の分離確認、バッチ効果検出 |

## 実行方法

```bash
python scripts/step_03_overview_visualization.py
```

## 可視化結果

以下に、本書データで生成した3つの可視化結果を示します。

### (a) 相関行列ヒートマップ

![相関行列ヒートマップ](images/fig1a_correlation.png)

サンプル間のピアソン相関係数を行列として可視化したものです。Normal群（青）とTumor群（赤）がそれぞれクラスターを形成し、群内の相関が高い（赤色）一方、群間の相関は相対的に低い（淡い色）ことが読み取れます。これは、腫瘍組織と正常組織が異なるタンパク質発現プロファイルを持つことを示しています。

### (b) 階層的クラスタリング

![階層的クラスタリング](images/fig1b_clustering.png)

全2,234タンパク質の発現量をもとに、Ward法による階層的クラスタリングを行いました。デンドログラム（樹形図）を見ると、NormalサンプルとTumorサンプルが明確に2つの枝に分離しています。教師なし手法（ラベルを使わない分類）でも、腫瘍/正常の群構造を正しく再現できることが確認されました。

### (c) PCA（主成分分析）

![PCA](images/fig1c_pca.png)

PCAにより高次元のプロテオームデータを2次元に圧縮した結果です。PC1（第1主成分、寄与率 **45.3%**）の軸上でNormal群とTumor群が明瞭に分離しており、データの最大変動がTumor/Normalの対比であることを示しています。PC2（寄与率 **14.5%**）は個体間のばらつきを反映しています。

## 本書データでの実測値（sage + 18 ファイル）

| 指標 | 本書 (sage) | 論文 (DIA-NN) |
|------|-------------|--------------|
| サンプル数 | 18（Normal 9 + Tumor 9） | 32（Normal 16 + Tumor 16） |
| 解析タンパク質数 | 2,234 | ～10,329 |
| 相関値 (Pearson r) | 0.70 – 0.93 | 0.744 – 0.982 |
| **PCA PC1 寄与率** | **45.3%** | **42.1%** |
| PCA PC2 寄与率 | 14.5% | 10.8% |
| Normal/Tumor 分離 | 明確 | 明確 |

PC1 寄与率が論文の **42.1%** に対して本書 **45.3%** と近い値を示す点が重要です。これは「大腸がんの Tumor/Normal 対比が強い生物学的シグナル」を示し、ツール（DIA-NN → sage）やサンプル数（32 → 18）を変えても **本質的な群構造は再現される** ことを意味します。

## まとめ

3つの可視化すべてで、腫瘍組織と正常組織が明確に分離されることを確認しました。本書のデータでは PC1 寄与率が論文に近い 45.3% となり、商用クリアな sage-proteomics でも論文と同等の生物学的シグナルを検出できることが示されました。

> 前回: [#5 前処理](article-05-preprocess.md)
> 次回: [#7 差分発現解析](article-07-differential.md) — Welch's t-testとVolcanoプロット

#バイオインフォマティクス #プロテオミクス #Python #PCA #labcode
