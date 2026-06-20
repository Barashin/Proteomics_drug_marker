# regenerate_all_figures.py の解説

## このコードの役割

Toyota et al. 2025論文の**全Figure（1a〜3）を論文スタイルで一括再生成**するマスタースクリプトです。sage→前処理パイプラインで得られた`preprocessed_data.csv`を入力として、相関解析・PCA・差分発現・クラスタリング・COSMIC解析・ステージ別解析の全てを実行し、論文準拠の色・レイアウトで可視化します。

## 全体フロー

```
preprocessed_data.csv → 相関解析 → PCA → 差分発現(t検定) → クラスタリング → COSMIC解析 → ステージ解析 → 全Figure出力
```

| | ファイル | 内容 |
|-|---------|------|
| 入力 | `results/preprocessed_data.csv` | sage由来の前処理済みタンパク質発現データ |
| 入力 | `results/sample_info.csv` | サンプル情報（Normal/Tumor分類） |
| 入力 | `results/clinical_info.csv` | 臨床情報（ステージ分類） |
| 出力 | `results/figures/*.png` | 論文Figure 1a〜3の全図 |
| 出力 | `blog/images/*.png` | 同じ図のブログ用コピー |
| 出力 | `results/tables/*.csv` | 解析結果テーブル（相関・差分発現・クラスター・ANOVA等） |

**所要時間の目安**: 約5-8分（全Figure生成・PCA・クラスタリング・統計解析込み）

## セクション別解説

### セクション1: データ読み込みと色設定

```python
# mzML → sage → preprocess パイプライン由来のデータを読み込み
df = pd.read_csv(os.path.join(RESULTS_DIR, "preprocessed_data.csv"), index_col=0)
sample_info = pd.read_csv(os.path.join(RESULTS_DIR, "sample_info.csv"))
clinical = pd.read_csv(os.path.join(RESULTS_DIR, "clinical_info.csv"))

# 論文準拠カラースキーム
HIST_NONTUMOR = "#3498DB"  # 青（正常組織）
HIST_TUMOR = "#E74C3C"     # 赤（腫瘍組織）
STAGE_COLORS = {
    "Normal": "#3498DB", "I": "#66BB6A", "II": "#FFD54F",
    "III": "#FFA726", "IV": "#AB47BC"
}
COLOR_UP_PROTEIN = "#E8A0BF"    # ピンク（腫瘍で上昇）
COLOR_DOWN_PROTEIN = "#A0C4E8"  # 水色（腫瘍で減少）
```

**目的**: sage-proteomics由来の前処理済みデータを読み込み、論文Figure準拠の色設定を定義する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `pd.read_csv(index_col=0)` | 1列目（タンパク質名）をインデックスとして読み込み |
| `stage_map.get()` | サンプル名→ステージのマッピング辞書 |
| `save_fig()` | results/figures/とblog/images/の両方に図を保存 |

**引数の意味**:
- `HIST_NONTUMOR/TUMOR`: 組織型色（Figure 1b/2a-hの行色バー用）
- `STAGE_COLORS`: ステージ別色（Normal→青からStage IV→紫のグラデーション）
- `COLOR_UP/DOWN_PROTEIN`: タンパク質発現方向色（Figure 1b/2a-hの列色バー用）

**つまずきやすいポイント**:
- 色コードは論文Figureに厳密準拠（独自色は論文との整合性を損なう）
- データパイプラインは mzML→sage→preprocess の順で処理済みが前提

### セクション2: Figure 1a（相関行列ヒートマップ）

```python
def plot_fig1a():
    corr = df.corr(method="pearson")  # 全サンプル間のピアソン相関係数行列

    # 行/列色バー = 組織型（Normal青/Tumor赤）
    histology_colors = conditions.map({
        "Normal": HIST_NONTUMOR, "Tumor": HIST_TUMOR
    }).reindex(corr.index)

    g = sns.clustermap(
        corr, method="average", metric="correlation",
        cmap="Reds", vmin=0.7, vmax=1.0,
        row_colors=histology_colors, col_colors=histology_colors,
        linewidths=0, xticklabels=True, yticklabels=True,
    )
```

**目的**: 全サンプル間のピアソン相関係数を計算し、階層的クラスタリングで並べ替えたヒートマップを作成する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `df.corr(method="pearson")` | 全サンプルペア間のピアソン相関係数行列を計算 |
| `sns.clustermap(method="average")` | Average linkageで階層クラスタリング＋ヒートマップ |
| `vmin=0.7, vmax=1.0` | 相関係数の表示範囲（0.7〜1.0の高相関域に限定） |

**引数の意味**:
- `metric="correlation"`: 1-ピアソン相関を距離として使用
- `cmap="Reds"`: 赤系カラーマップ（高相関=濃い赤）
- `row_colors/col_colors`: 組織型を示すカラーバー

**つまずきやすいポイント**:
- 相関行列は対称行列なので行・列のクラスタリング結果は同じになる
- vmin=0.7設定により低相関（<0.7）は同じ薄い色で表示される

### セクション3: Figure 1b（階層クラスタリング）

```python
def plot_fig1b():
    # 二重行色バー: 組織型 + 臨床ステージ
    row_colors_df = pd.DataFrame({
        "Histology": hist_colors,         # Normal青/Tumor赤
        "Clinical stage": stage_colors,    # Normal青→Stage IV紫
    })

    # 列色バー: タンパク質発現方向（腫瘍vs正常の平均差）
    normal_mean = df[normal_samples].mean(axis=1)
    tumor_mean = df[tumor_samples].mean(axis=1)
    fc = tumor_mean - normal_mean  # log2スケールでの差
    protein_colors = fc.apply(lambda x: COLOR_UP_PROTEIN if x > 0 else COLOR_DOWN_PROTEIN)

    g = sns.clustermap(
        df.T, method="ward", metric="euclidean",
        cmap="RdBu_r", center=0, vmin=-3, vmax=3,
        row_colors=row_colors_df, col_colors=col_colors_df,
        z_score=1,  # 列方向（タンパク質方向）にZ-score正規化
    )
```

**目的**: 全データをWard法階層クラスタリングで並べ替え、組織型・ステージ・タンパク質方向の3軸情報を色バーで表現したヒートマップを作成する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `df.T` | データを転置（サンプル×タンパク質→タンパク質×サンプル） |
| `method="ward"` | Ward連結法（クラスター内分散最小化） |
| `z_score=1` | 列方向（タンパク質ごと）にZ-score正規化 |
| `cmap="RdBu_r"` | 青-白-赤カラーマップ（reversed） |

**引数の意味**:
- `center=0`: Z-score=0を白色に設定
- `vmin=-3, vmax=3`: Z-scoreの標準的表示範囲
- `row_colors=row_colors_df`: 2つの色バー（組織型・ステージ）を同時表示

**つまずきやすいポイント**:
- Z-score=1（列方向）により、各タンパク質の発現パターンが平均0・標準偏差1に正規化される
- 二重色バーは論文Figure 1bの複雑な情報表現を再現

### セクション4: Figure 1c（PCA）

```python
def plot_fig1c():
    pca = PCA(n_components=2)
    scores = pca.fit_transform(df.T)  # サンプル×PC1/PC2のスコア行列

    # Normal組織: 全て緑色（論文スタイル）
    ax.scatter(scores[normal_idx, 0], scores[normal_idx, 1],
               c="#2E8B57", marker="o", s=100, alpha=0.85)

    # Tumor組織: 臨床ステージ別色分け
    for i in tumor_idx:
        stage = stage_map.get(sample, "III")  # デフォルトStage III
        color = STAGE_COLORS.get(stage, "#E74C3C")
        ax.scatter(scores[i, 0], scores[i, 1], c=color, marker="o", s=100)

    # 破線信頼楕円（論文スタイル）
    confidence_ellipse(scores[normal_idx, 0], scores[normal_idx, 1], ax, n_std=2.0,
                       linestyle="--", linewidth=1.5)
```

**目的**: 主成分分析による次元削減で、組織型・ステージの分離を2次元平面上で可視化する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `PCA(n_components=2)` | 2主成分への次元削減 |
| `pca.fit_transform(df.T)` | サンプル×遺伝子 → サンプル×PC1/PC2 |
| `confidence_ellipse()` | 95%信頼楕円の描画（正規分布仮定） |
| `explained_variance_ratio_` | 各主成分の寄与率（%表示用） |

**引数の意味**:
- `n_std=2.0`: 95%信頼区間（2σ）の楕円
- `linestyle="--"`: 破線スタイル（論文準拠）
- `scores[normal_idx, 0]`: Normal組織のPC1スコア

**つまずきやすいポイント**:
- Normal組織は単一色（緑）、Tumor組織はステージ別色分けの論文デザイン
- 寄与率を軸ラベルに含めて主成分の説明力を明示

### セクション5: Welch t検定（差分発現解析）

```python
def run_welch_ttest():
    results = []
    for protein in df.index:
        nv = df.loc[protein, normal_samples].dropna()  # Normal組織の発現値
        tv = df.loc[protein, tumor_samples].dropna()   # Tumor組織の発現値

        # Welch t検定（等分散性を仮定しない）
        t_stat, p_val = stats.ttest_ind(tv, nv, equal_var=False)
        log2fc = tv.mean() - nv.mean()  # log2 fold change

        results.append({
            "Protein": protein, "Log2FC": log2fc,
            "T_statistic": t_stat, "P_value": p_val,
        })

    # 有意性判定: p<0.05 かつ |log2FC|>1
    result_df["Significant"] = "NS"
    result_df.loc[(result_df["P_value"] < 0.05) & (result_df["Log2FC"] > 1.0), "Significant"] = "Up"
    result_df.loc[(result_df["P_value"] < 0.05) & (result_df["Log2FC"] < -1.0), "Significant"] = "Down"
```

**目的**: Normal vs Tumor間でタンパク質発現に有意差があるかをWelch t検定で判定し、ボルケーノプロット用データを作成する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `stats.ttest_ind(equal_var=False)` | Welch t検定（不等分散t検定） |
| `dropna()` | 欠損値を除外して純粋なデータで検定 |
| `-np.log10(max(p_val, 1e-300))` | ボルケーノプロット用の-log10(p値)変換 |

**引数の意味**:
- `equal_var=False`: 分散の等質性を仮定しないWelch t検定
- `|log2FC| > 1`: 2倍以上の変化（生物学的意義ありの閾値）
- `p < 0.05`: 統計的有意性の閾値

**つまずきやすいポイント**:
- log2FCが正ならTumorで上昇、負なら減少を意味する
- 統計的有意性（p<0.05）と生物学的意義（|FC|>2）を組み合わせた判定

### セクション6: Figure 2a-h（サブセットヒートマップ・PCA）

```python
def plot_fig2_topN(result_df):
    # Top N/2 上昇 + Top N/2 下降タンパク質を選択
    for n, fig_c, fig_f in [(50, "fig2c_clustering_top50.png", "fig2f_pca_top50.png"),
                             (100, "fig2d_clustering_top100.png", "fig2g_pca_top100.png"),
                             (200, "fig2e_clustering_top200.png", "fig2h_pca_top200.png")]:

        n_half = n // 2
        top_up = sig_up.head(min(n_half, len(sig_up)))      # 上昇Top25/50/100
        top_down = sig_down.head(min(n_half, len(sig_down)))  # 下降Top25/50/100
        top_combined = pd.concat([top_up, top_down])

        # クラスタリング（Figure 2c-e）
        g = sns.clustermap(top_data.T, method="ward", z_score=1, ...)

        # PCA（Figure 2f-h）
        pca = PCA(n_components=2)
        scores = pca.fit_transform(top_data.T)
```

**目的**: 差分発現タンパク質のうち変化量上位N個（50/100/200）を選抜し、サブセット解析用のヒートマップ・PCAを生成する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `sig_up.head(n_half)` | Log2FC降順でTop N/2の上昇タンパク質を選択 |
| `sig_down.head(n_half)` | Log2FC昇順でTop N/2の下降タンパク質を選択 |
| `pd.concat([top_up, top_down])` | 上昇・下降タンパク質を結合 |

**引数の意味**:
- `n_half = n // 2`: バランス良く上昇・下降タンパク質を選択
- `[(50, ...), (100, ...), (200, ...)]`: 3つのサブセットサイズ
- `z_score=1`: タンパク質ごとのZ-score正規化（Figure 1bと同様）

**つまずきやすいポイント**:
- 上昇・下降のバランスを保ちながらTop Nを選択する仕組み
- 同じ色設定・レイアウトでFigure 1との一貫性を保つ

### セクション7: COSMIC解析

```python
def plot_cosmic():
    COSMIC_CRC_GENES = [
        "APC", "TP53", "KRAS", "BRAF", "PIK3CA", "SMAD4", ...  # 65遺伝子
    ]
    COSMIC_ALL_GENES = list(set(COSMIC_CRC_GENES + [
        "ABL1", "AKT1", "ALK", "AR", "BRCA1", "BRCA2", ...    # 389遺伝子
    ]))

    identified = set(df.index)  # 本研究で同定されたタンパク質
    overlap_all = identified & set(COSMIC_ALL_GENES)    # 積集合
    overlap_crc = identified & set(COSMIC_CRC_GENES)    # 積集合

    coverage_all = len(overlap_all) / len(set(COSMIC_ALL_GENES)) * 100
    coverage_crc = len(overlap_crc) / len(set(COSMIC_CRC_GENES)) * 100
```

**目的**: COSMIC（がん関連遺伝子データベース）との照合により、本研究のプロテオミクス解析の網羅性をカバー率として定量評価する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `set(df.index)` | 同定タンパク質をset（集合）に変換 |
| `identified & set(COSMIC_ALL_GENES)` | 積集合演算（共通要素の抽出） |
| `len(overlap) / len(cosmic) * 100` | カバー率（%）の計算 |

**引数の意味**:
- `COSMIC_CRC_GENES`: 大腸がん特異的な65遺伝子
- `COSMIC_ALL_GENES`: 全がん種の389遺伝子
- `overlap_all/crc`: 本研究で検出できたCOSMIC遺伝子

**つまずきやすいポイント**:
- COSMIC完全リストは有料のため論文記載の代表的遺伝子のサブセットを使用
- カバー率70%以上は非常に優秀な成績（技術的制約を考慮）

### セクション8: Figure 3（ステージ別解析）

```python
def plot_fig3():
    # One-way ANOVA (ステージ間差の検定)
    for protein in df.index:
        groups = []  # 各ステージの発現値リスト
        for stage in available_stages:
            vals = df.loc[protein, stage_samples].dropna()
            if len(vals) >= 2:
                groups.append(vals.values)
        f_stat, p_val = stats.f_oneway(*groups)  # F検定実行

    # FDR補正
    reject, fdr, _, _ = multipletests(anova_df["P_value"], method="fdr_bh")

    # 階層クラスタリング（30クラスター）
    Z = linkage(z_data.values, method="ward")
    clusters = fcluster(Z, t=n_clust, criterion="maxclust")

    # 2×2パネル：各パネル=ヒートマップ+プロファイルプロット
    fig = plt.figure(figsize=(14, 12))
    outer_gs = GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)
```

**目的**: 臨床ステージ（Normal/I/II/III/IV）間でのタンパク質発現変化を統計的に解析し、論文Figure 3形式の4パネル複合図を作成する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `stats.f_oneway(*groups)` | One-way ANOVA（一元配置分散分析） |
| `multipletests(method="fdr_bh")` | Benjamini-Hochberg法FDR補正 |
| `GridSpec(2, 2)` | 2×2パネルレイアウト |
| `subgridspec(2, 1, height_ratios=[3, 2])` | 各パネル内の上下分割 |

**引数の意味**:
- `available_stages`: データに実際に存在するステージ
- `t=n_clust`: 30クラスターに分割（論文準拠）
- `height_ratios=[3, 2]`: ヒートマップ:プロファイル=3:2の高さ比

**つまずきやすいポイント**:
- One-way ANOVAは3群以上の平均値差を一度に検定（多重比較回避）
- 4パネル構成で増加・減少パターンのクラスターを視覚的に整理

## 【深掘り】 一括再生成スクリプトとは？

- **ひとことで**: 個別解析スクリプトを統合し、論文全Figureを一回の実行で完全再現するマスタースクリプト
- **利点**: 個別スクリプト実行の手間を削減、色・レイアウトの一貫性保証、ブログ用画像も同時生成
- **論文再現の要点**: 元論文と同じ色設定・統計手法・可視化パラメータを厳密に再現
- **商用利用対応**: sage-proteomics等の商用利用可能ツール由来データを使用
- **品質管理**: results/figuresとblog/imagesの二重保存でバックアップ効果

## 【深掘り】 信頼楕円とは？

- **ひとことで**: 2次元散布図上でデータ点の95%が含まれる範囲を楕円で表現する統計的可視化手法
- **計算原理**: データの共分散行列から楕円の長軸・短軸・回転角を算出
- **PCAでの意味**: 主成分空間での各群（Normal/Tumor）の分布範囲を視覚化
- **パラメータ**: n_std=2.0で95%信頼区間（2σ範囲）
- **論文での使用**: 群間分離の程度を直感的に示す（楕円が重なれば分離不良）

## 【深掘り】 Ward連結法とは？

- **ひとことで**: 階層クラスタリングでクラスター内の分散（ばらつき）を最小化する連結方法
- **他の手法との違い**: Single linkage（最近傍）、Average linkage（平均距離）より均一なクラスターを形成
- **計算方法**: 各結合ステップで「結合後のクラスター内分散の増加量」を最小化
- **利点**: コンパクトで球状のクラスターが得られ、生物学的解釈しやすい
- **プロテオミクスでの意義**: 発現パターンが似たタンパク質群を明確に分離

## よくある質問

**Q: なぜ個別スクリプトがあるのに一括再生成スクリプトが必要なのですか？**
A: ①色・レイアウトの一貫性確保 ②論文Figure番号との正確な対応 ③ブログ記事用画像の同時生成 ④全解析の一括実行による効率化のためです。

**Q: sage由来データとDIA-NN由来データで結果は変わりますか？**
A: 検出タンパク質数・定量精度に違いが出る可能性があります。sageはオープンソースでDIA-NNより検出数がやや少ないですが、商用利用制限がありません。

**Q: COSMIC遺伝子リストは論文完全版と同じですか？**
A: 論文補足データの代表的遺伝子のサブセットを使用しています。完全なCOSMICリストは有料ライセンスが必要なためです。

**Q: Figure 3の4パネルはどう選択されますか？**
A: Normal→Stage IVの変化量（trend）が最大の2クラスター（増加）と最小の2クラスター（減少）を自動選択します。論文補足データがある場合はそれを優先使用します。

**Q: 信頼楕円が重なっている場合の解釈は？**
A: Normal/Tumor群の分離が不十分であることを示します。より多くの差分発現タンパク質や異なる次元削減手法が必要かもしれません。

**Q: Z-score正規化のaxis=1（列方向）の意味は？**
A: 各タンパク質について全サンプルにわたって平均0・標準偏差1に正規化することです。これにより発現量の絶対値ではなくパターンの形状で比較できます。

## 出力されるファイルの内容

### figures/とblog/images/の全Figure
```
fig1a_correlation.png       # Figure 1a: サンプル間相関ヒートマップ
fig1b_clustering.png        # Figure 1b: 全データ階層クラスタリング
fig1c_pca.png              # Figure 1c: 全データPCA
fig_bonus_volcano.png      # ボーナス: ボルケーノプロット
fig2a_heatmap_all.png      # Figure 2a: 全差分発現タンパク質ヒートマップ
fig2b_pca_all.png          # Figure 2b: 全差分発現タンパク質PCA
fig2c_clustering_top50.png  # Figure 2c: Top50クラスタリング
fig2d_clustering_top100.png # Figure 2d: Top100クラスタリング
fig2e_clustering_top200.png # Figure 2e: Top200クラスタリング
fig2f_pca_top50.png        # Figure 2f: Top50 PCA
fig2g_pca_top100.png       # Figure 2g: Top100 PCA
fig2h_pca_top200.png       # Figure 2h: Top200 PCA
fig_cosmic_coverage.png    # COSMIC遺伝子カバー率
fig3_stage_heatmap.png     # Figure 3: ステージ別4パネル複合図
```

### tables/の解析結果CSV
```
correlation_matrix.csv     # 全サンプル間相関係数行列
pca_variance.csv          # PCA各成分の寄与率
differential_proteins.csv # Welch t検定結果（全タンパク質）
cosmic_overlap.csv       # COSMIC遺伝子との重複リスト
anova_results.csv        # One-way ANOVA結果
cluster_assignments.csv  # 階層クラスタリング割り当て
```

## 重要な注意点

- **データパイプライン依存**: mzML→sage→前処理の順で処理済みデータが前提
- **実行時間**: 全Figure生成に5-8分必要（大量のPCA・クラスタリング計算のため）
- **メモリ使用量**: 数千タンパク質×数十サンプルの行列計算で数GB必要
- **色の一貫性**: 論文Figure準拠色を変更すると論文との比較が困難になる
- **ファイル上書き**: 既存のfigures/とblog/images/の同名ファイルは上書きされる

## 次のステップ

この再生成スクリプトの出力は以下で利用されます：
- ブログ記事への図の埋め込み（blog/images/から参照）
- 論文との視覚的比較・検証
- 追加解析用のベースライン図として参照
- 技術書・プレゼンテーション用の高品質Figure素材