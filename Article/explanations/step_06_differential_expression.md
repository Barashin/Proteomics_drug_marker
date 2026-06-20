# step_06_differential_expression.py の解説

## このコードの役割

前処理済みプロテオミクスデータを用いて差分発現解析（Differential Expression Analysis）を実行し、腫瘍（Tumor）と正常（Normal）組織間で発現量が統計的に有意に異なるタンパク質を特定するスクリプトです。論文 Figure 2 に対応する可視化も含みます。

## 全体フロー

```
前処理済みデータ → Welch's t検定 → Volcano プロット → 有意差タンパク質抽出 → ヒートマップ・PCA
```

| | ファイル | 内容 |
|-|---------|------|
| 入力 | `results/preprocessed_data.csv` | Log2変換・前処理済みタンパク質定量データ |
| 入力 | `results/sample_info.csv` | サンプル群情報（Normal/Tumor） |
| 出力 | `results/tables/differential_proteins.csv` | 全タンパク質の統計検定結果 |
| 出力 | `results/figures/fig_bonus_volcano.png` | Volcano プロット（差分発現概観） |
| 出力 | `results/figures/fig2a_heatmap_all.png` | 全有意差タンパク質ヒートマップ |
| 出力 | `results/figures/fig2b_pca_all.png` | 全有意差タンパク質PCA |
| 出力 | `results/figures/fig2c-h_*.png` | Top N タンパク質解析（N=50,100,200） |

**所要時間の目安**: 約3-5分（タンパク質数による）

## セクション別解説

### セクション1: 統計的閾値と色設定

```python
# p値の閾値: 統計的に有意とみなす上限
P_VALUE_THRESHOLD = 0.05

# Log2 Fold Change の閾値: log2(2) = 1.0 → fold change > 2
LOG2FC_THRESHOLD = 1.0

# カラーパレット定義
COLOR_UP = "#E74C3C"        # 腫瘍で発現上昇
COLOR_DOWN = "#3498DB"      # 腫瘍で発現低下
COLOR_NS = "#CCCCCC"        # 有意差なし
```

**目的**: 生物学的に意味のある変化と統計的有意性の基準を設定し、一貫した色分けで可視化する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `P_VALUE_THRESHOLD = 0.05` | 5%の偽陽性率でタンパク質を有意と判定 |
| `LOG2FC_THRESHOLD = 1.0` | 2倍以上の発現変化を生物学的に意味ありと判定 |

**引数の意味**:
- `p < 0.05`: 統計的有意性（「偶然ではなく真の差がある」の信頼度95%）
- `log2FC > 1.0`: fold change > 2倍（腫瘍で発現2倍以上増加）
- `log2FC < -1.0`: fold change < 0.5倍（腫瘍で発現1/2以下に減少）

**つまずきやすいポイント**:
- Log2スケールでは差の引き算がfold changeになる
- 閾値は研究分野により異なる（プロテオミクスでは1.5-2倍が一般的）

### セクション2: Welch's t検定実行

```python
def welch_ttest(df, normal_samples, tumor_samples):
    for protein in df.index:
        normal_vals = df.loc[protein, normal_samples].dropna()
        tumor_vals = df.loc[protein, tumor_samples].dropna()

        # 最低2サンプル必要（t検定の前提条件）
        if len(normal_vals) < 2 or len(tumor_vals) < 2:
            continue

        # Welch's t検定を実行
        t_stat, p_val = stats.ttest_ind(tumor_vals, normal_vals, equal_var=False)

        # Log2 Fold Change計算
        log2fc = tumor_vals.mean() - normal_vals.mean()
```

**目的**: 各タンパク質について Normal群 vs Tumor群 の平均値差を統計的に検定し、偶然ではない真の差があるかを判定する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `stats.ttest_ind()` | 独立2標本のt検定（scipy統計ライブラリ） |
| `equal_var=False` | Welch補正適用（分散が等しくない前提） |
| `dropna()` | 欠損値を除外してから検定実行 |

**引数の意味**:
- `equal_var=False`: 2群の分散が等しいと仮定しない（現実的）
- `t_stat`: t統計量（正なら腫瘍群が高い、負なら低い）
- `p_val`: 有意確率（0に近いほど有意差大）
- `log2fc`: Log2 fold change（データが既にlog2なので平均の差）

**つまずきやすいポイント**:
- 最低サンプル数（各群2以上）が必要
- 欠損値があると検定できない
- Welchの補正により自由度が非整数になる

### セクション3: 有意性分類

```python
# 全タンパク質を初期値"NS"で初期化
result_df["Significant"] = "NS"

# Up-regulated条件マスク
up_mask = (result_df["P_value"] < P_VALUE_THRESHOLD) & (result_df["Log2FC"] > LOG2FC_THRESHOLD)

# Down-regulated条件マスク
down_mask = (result_df["P_value"] < P_VALUE_THRESHOLD) & (result_df["Log2FC"] < -LOG2FC_THRESHOLD)

# 条件に合致する行を"Up"または"Down"に分類
result_df.loc[up_mask, "Significant"] = "Up"
result_df.loc[down_mask, "Significant"] = "Down"
```

**目的**: p値とfold changeの両方の条件を満たすタンパク質を Up/Down/NS の3カテゴリに分類する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `&` | ブールマスクのAND演算（両方の条件を満たす） |
| `result_df.loc[mask, column]` | 条件に合う行の特定列を一括更新 |

**引数の意味**:
- `"Up"`: p<0.05 かつ log2FC>1.0（統計的に有意な増加）
- `"Down"`: p<0.05 かつ log2FC<-1.0（統計的に有意な減少）
- `"NS"`: 上記以外（Not Significant、有意差なし）

### セクション4: Volcano プロット描画

```python
def plot_volcano(result_df):
    fig, ax = plt.subplots(figsize=(8, 6))

    # 3カテゴリごとに色分けして描画
    for sig, color, alpha in [
        ("NS", COLOR_NS, 0.3),     # 灰色・半透明
        ("Up", COLOR_UP, 0.6),     # 赤色
        ("Down", COLOR_DOWN, 0.6), # 青色
    ]:
        mask = result_df["Significant"] == sig
        ax.scatter(
            result_df.loc[mask, "Log2FC"],       # X軸: fold change
            result_df.loc[mask, "Neg_log10_P"],  # Y軸: -log10(p値)
            c=color, s=10, alpha=alpha,
            label=f"{sig} ({mask.sum()})",
        )
```

**目的**: 全タンパク質の統計検定結果を1枚の散布図で俯瞰し、有意差のあるタンパク質の分布パターンを可視化する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `ax.scatter()` | 散布図描画（matplotlib） |
| `ax.axhline()`/`ax.axvline()` | 閾値線の描画（水平線・垂直線） |
| `-np.log10(P_VALUE_THRESHOLD)` | p値を-log10変換（小さいp値ほど上に位置） |

**引数の意味**:
- X軸 `Log2FC`: 発現変化量（右=増加、左=減少）
- Y軸 `Neg_log10_P`: 統計的有意性（上=有意、下=非有意）
- 右上象限: 有意に増加、左上象限: 有意に減少

**つまずきやすいポイント**:
- p=0のとき log10(0)=-∞ になるため、`max(p_val, 1e-300)`で下限設定
- 火山状分布になるのは、大部分が中央（変化なし）に集中するため

### セクション5: 階層的クラスタリングヒートマップ

```python
def plot_heatmap_all(df, result_df, sample_info):
    # 有意差タンパク質のみ抽出
    sig_df = result_df[result_df["Significant"] != "NS"].copy()
    sig_data = df.loc[df.index.isin(sig_proteins)]

    g = sns.clustermap(
        sig_data.T,                # 転置: 行=サンプル、列=タンパク質
        method="ward",             # ウォード法
        cmap="RdBu_r",            # 赤-青カラーマップ
        z_score=1,                 # 列方向（タンパク質）でZ-score標準化
        row_colors=sample_row_colors,   # 左カラーバー: Normal/Tumor
        col_colors=protein_colors,      # 上カラーバー: Up/Down
    )
```

**目的**: 有意差のある全タンパク質について、サンプルとタンパク質を同時に階層的クラスタリングし、発現パターンの全体構造を可視化する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `sns.clustermap()` | 階層クラスタリング付きヒートマップ（seaborn） |
| `method="ward"` | ウォード法（クラスタ内分散最小化） |
| `z_score=1` | 列方向（タンパク質ごと）でZ-score標準化 |
| `row_colors`/`col_colors` | サイドバーの色分け |

**引数の意味**:
- `z_score=1`: 各タンパク質を平均0・SD1に標準化（相対パターン比較）
- `cmap="RdBu_r"`: 赤（高発現）-白（平均）-青（低発現）
- `row_colors`: 左側にNormal（青）・Tumor（赤）のカラーバー
- `col_colors`: 上側にUp（赤）・Down（青）のタンパク質方向カラーバー

### セクション6: PCA散布図

```python
def plot_pca_all(df, result_df, sample_info):
    # 有意差タンパク質でPCA実行
    pca = PCA(n_components=2)
    scores = pca.fit_transform(sig_data.T)

    # Normal/Tumor別に散布図描画
    for condition, color in [("Normal", COLOR_NORMAL), ("Tumor", COLOR_TUMOR)]:
        mask = conditions.reindex(df.columns) == condition
        ax.scatter(scores[mask, 0], scores[mask, 1], c=color, label=condition)
        # 95%信頼楕円も描画
        confidence_ellipse(scores[mask, 0], scores[mask, 1], ax, n_std=2.0)
```

**目的**: 有意差タンパク質の高次元発現パターンを2次元に圧縮し、Normal群・Tumor群の分離度を定量的に評価する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `PCA(n_components=2)` | 主成分分析で2次元に圧縮（scikit-learn） |
| `fit_transform()` | PCAの学習と座標変換を同時実行 |
| `confidence_ellipse()` | 95%信頼楕円描画（カスタム関数） |
| `explained_variance_ratio_` | 各主成分の寄与率 |

**引数の意味**:
- `sig_data.T`: サンプル（行）×タンパク質（列）に転置
- `scores`: 各サンプルのPC1・PC2座標
- `n_std=2.0`: 2標準偏差（約95%信頼領域）の楕円

### セクション7: Top Nタンパク質解析

```python
def plot_top_n_analysis(df, result_df, sample_info, n_top=50):
    # Up/Down別に|fold change|順でソート
    sig_up = sig_df[sig_df["Significant"] == "Up"].copy()
    sig_down = sig_df[sig_df["Significant"] == "Down"].copy()
    sig_up = sig_up.sort_values("AbsLog2FC", ascending=False)
    sig_down = sig_down.sort_values("AbsLog2FC", ascending=False)

    for n in [50, 100, 200]:
        # 上昇Top N + 低下Top N を選択
        top_up = sig_up.head(n)["Protein"].values
        top_down = sig_down.head(n)["Protein"].values
        top_proteins = np.concatenate([top_up, top_down])
```

**目的**: 最も発現変化の大きいタンパク質に絞って解析することで、ノイズを除去し、より明確な群分離パターンを観察する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `sort_values("AbsLog2FC", ascending=False)` | 絶対値fold changeの降順ソート |
| `np.concatenate()` | 上昇・下降タンパク質リストの結合 |
| `head(n)` | 上位N個を選択 |

**引数の意味**:
- `AbsLog2FC`: |log2FC|（変化の大きさ、方向無視）
- `n in [50, 100, 200]`: 論文準拠の3段階タンパク質数
- 上昇Top50 + 低下Top50 = 合計100タンパク質を使用

## 【深掘り】 Welch's t検定とは？

- **ひとことで**: 2群の平均値差が統計的に有意かを判定するが、2群の分散が異なっても正確に検定できる改良版t検定
- **定義**: Student's t検定の「等分散仮定」を緩和し、Welch-Satterthwaiteの式で自由度を補正
- **どんなとき使う**: 生物学データで2群の分散が異なることが多いため、プロテオミクス解析では標準的
- **落とし穴**: 各群に最低2サンプル必要、欠損値があると検定不可
- **参考**: scipy.stats.ttest_ind(equal_var=False)

## 【深掘り】 Volcano プロットとは？

- **ひとことで**: 差分発現解析の結果を「統計的有意性」と「生物学的意義」の2軸で俯瞰できる散布図
- **定義**: X軸=fold change（変化量）、Y軸=-log10(p値)（有意性）のプロット
- **名前の由来**: データ点の分布が火山の噴火のような形になることから
- **読み方**: 右上=有意に増加、左上=有意に減少、下部=有意差なし
- **閾値線**: 垂直線（fold change閾値）と水平線（p値閾値）で4象限に分割
- **参考**: 論文Figure「Volcano plot showing differential protein abundance」

## 【深掘り】 Z-score標準化とは？

- **ひとことで**: 各タンパク質の発現値を「平均0、標準偏差1」に変換し、異なるタンパク質を同じスケールで比較可能にする変換
- **定義**: z = (値 - 平均) / 標準偏差。サンプル間の相対的パターンを強調
- **計算例**: タンパク質Aで平均=100、SD=20なら、値120 → z=(120-100)/20=1.0
- **ヒートマップでの効果**: 高発現・低発現タンパク質を同じ色範囲で可視化可能
- **注意**: 絶対発現量情報は失われ、相対パターンのみ表現
- **参考**: seaborn clustermap z_score=1（列方向標準化）

## 【深掘り】 95%信頼楕円とは？

- **ひとことで**: 2次元散布図上で、データの散らばりを表す楕円で、約95%のデータ点がこの楕円内に含まれる
- **定義**: 2変量正規分布を仮定し、共分散行列から楕円の形状・向き・大きさを決定
- **PCAでの意味**: 各群（Normal/Tumor）のPC1-PC2空間での散らばり範囲を可視化
- **楕円の重なり**: 重なりが小さいほど2群の分離が良好
- **パラメータ**: n_std=2.0で約95%信頼領域（2標準偏差）
- **参考**: matplotlib.patches.Ellipse + Affine変換

## よくある質問

**Q: p<0.05 と fold change>2 の両方が必要なのはなぜ？**
A: p値は統計的有意性（偶然ではない）、fold changeは生物学的意義（変化が大きい）を示します。両方の条件で、「統計的に信頼でき、かつ生物学的に意味のある変化」を特定できます。

**Q: Welch's t検定と通常のt検定の違いは？**
A: 通常のt検定は2群の分散が等しいと仮定しますが、Welch's t検定はこの仮定を不要にします。生物学データでは群間で分散が異なることが多いため、Welch版がより適切です。

**Q: Top Nタンパク質解析で50/100/200を選ぶ理由は？**
A: 変化量の大きいタンパク質ほど生物学的に重要な可能性が高いためです。段階的に解析することで、最適なタンパク質数を探索し、ノイズの影響を評価できます。

**Q: Z-score標準化はいつ使う？**
A: ヒートマップで異なる発現レベルのタンパク質を同じ色範囲で比較したい場合です。標準化により、高発現タンパク質が低発現タンパク質のパターンを隠すことを防げます。

**Q: 自分のデータで閾値を変えたい場合は？**
A: `P_VALUE_THRESHOLD`と`LOG2FC_THRESHOLD`を変更してください。より厳しい基準なら0.01と1.5、より緩い基準なら0.1と0.5などが考えられます。

## 出力される統計結果の形式

### differential_proteins.csv
```
Protein,Mean_Normal,Mean_Tumor,Log2FC,T_statistic,P_value,Neg_log10_P,Significant
TP53,12.3,14.8,2.5,3.2,0.003,2.5,Up
EGFR,15.1,13.2,-1.9,-2.8,0.012,1.9,Down
MYC,10.5,10.7,0.2,0.3,0.76,0.1,NS
```

- **Protein**: タンパク質名（遺伝子シンボル）
- **Mean_Normal/Tumor**: 各群の平均Log2発現値
- **Log2FC**: Tumor - Normal の差（Log2 fold change）
- **T_statistic**: Welch's t検定の統計量
- **P_value**: 有意確率（小さいほど有意）
- **Significant**: Up（増加）/Down（減少）/NS（有意差なし）

## 重要な注意点

- **多重検定補正**: この実装では未補正。必要に応じてBenjamini-Hochberg法等を適用
- **サンプル数**: 各群最低3-5サンプル推奨（統計検出力確保）
- **欠損値の扱い**: 欠損があるタンパク質は検定から除外される
- **メモリ使用量**: 大規模データ（1万タンパク質超）では数GB必要

## 次のステップ

この解析結果は以下のステップで利用されます：
- Step 7: `step_07_cosmic_analysis.py`での既知がん関連遺伝子との照合
- Step 8: `step_08_stage_analysis.py`でのステージ別解析
- バイオマーカー候補タンパク質の機能解析・パスウェイ解析