# step_05_overview_visualization.py の解説

## このコードの役割

前処理済みプロテオミクスデータを用いて、データ全体の傾向を可視化するスクリプトです。論文 Figure 1 に相当する3種類の図（相関行列、階層的クラスタリング、PCA）を作成し、Normal群とTumor群の分離を確認します。

## 全体フロー

```
前処理済みデータ → 相関行列ヒートマップ → 階層的クラスタリング → PCA散布図 → 図の保存
```

| | ファイル | 内容 |
|-|---------|------|
| 入力 | `results/preprocessed_data.csv` | Log2変換・欠損値補完済みタンパク質マトリクス |
| 入力 | `results/sample_info.csv` | サンプル群情報（Normal/Tumor） |
| 出力 | `results/figures/fig1a_correlation.png` | サンプル間相関行列ヒートマップ |
| 出力 | `results/figures/fig1b_clustering.png` | 教師なし階層的クラスタリング |
| 出力 | `results/figures/fig1c_pca.png` | PCA散布図（ステージ別色分け） |

**所要時間の目安**: 約2-3分（サンプル・タンパク質数による）

## セクション別解説

### セクション1: データ読み込みと設定

```python
COLOR_NORMAL = "#4EAED1"    # 青色: Non-tumor（正常組織）を示す
COLOR_TUMOR = "#E8524A"     # 赤色: Tumor（腫瘍組織）を示す

def load_data():
    df = pd.read_csv(os.path.join(RESULTS_DIR, "preprocessed_data.csv"), index_col=0)
    sample_info = pd.read_csv(os.path.join(RESULTS_DIR, "sample_info.csv"))
    return df, sample_info
```

**目的**: 前処理済みデータと臨床情報を読み込み、論文準拠の色設定で一貫した可視化を準備する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `pd.read_csv(index_col=0)` | CSV読み込み時に1列目（タンパク質名）をインデックスに指定 |
| `os.path.join()` | OS非依存のファイルパス結合 |

**引数の意味**:
- `index_col=0`: 1列目（Protein列）を行インデックスとして扱う
- `COLOR_NORMAL/TUMOR`: 論文原図と同じ色コード（視覚的一貫性）

**つまずきやすいポイント**:
- CSVの1行目がヘッダー、1列目がタンパク質名の形式が前提
- sample_info.csvにSample列とCondition列が必須

### セクション2: 相関行列ヒートマップ（Figure 1a）

```python
# ピアソン相関係数の計算
corr = df.corr(method="pearson")

# サンプルの色分け（行・列に沿ったカラーバー）
conditions = sample_info.set_index("Sample")["Condition"]
colors = conditions.map({"Normal": COLOR_NORMAL, "Tumor": COLOR_TUMOR})

g = sns.clustermap(
    corr,                   # 描画するデータ（相関行列）
    method="average",       # クラスタリング手法: UPGMA法（群平均法）
    metric="correlation",   # 距離指標: 1 - ピアソン相関係数
    cmap="Reds",            # カラーマップ: 赤系グラデーション
    vmin=0.7, vmax=1.0,     # 色の範囲固定
    row_colors=row_colors,  # 行の左にカラーバー
    col_colors=row_colors,  # 列の上にカラーバー
)
```

**目的**: サンプル間の発現パターン類似度をピアソン相関係数で定量化し、Normal・Tumor群の自然な分離を確認する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `df.corr(method="pearson")` | 全サンプル間のピアソン相関係数行列を計算 |
| `sns.clustermap()` | ヒートマップ+デンドログラム+クラスタリングを統合実行 |
| `method="average"` | UPGMA法による階層クラスタリング |
| `metric="correlation"` | 距離=1-相関係数（高相関=近距離） |

**引数の意味**:
- `vmin=0.7, vmax=1.0`: 相関係数0.7-1.0を色の全範囲にマップ
- `cmap="Reds"`: 高い相関（類似）ほど濃い赤で表示
- `row_colors/col_colors`: Normal（青）・Tumor（赤）のカラーバー

**つまずきやすいポイント**:
- 相関行列は正方行列（サンプル数×サンプル数）になる
- クラスタリングにより行・列の順序が自動変更される

### セクション3: 階層的クラスタリング（Figure 1b）

```python
# タンパク質ごとのfold change方向を計算（列カラーバー用）
normal_samples_list = sample_info[sample_info["Condition"] == "Normal"]["Sample"].tolist()
tumor_samples_list = sample_info[sample_info["Condition"] == "Tumor"]["Sample"].tolist()

normal_mean = df[[s for s in normal_samples_list if s in df.columns]].mean(axis=1)
tumor_mean = df[[s for s in tumor_samples_list if s in df.columns]].mean(axis=1)
fc = tumor_mean - normal_mean

protein_direction_colors = fc.apply(
    lambda x: "#E74C3C" if x > 0 else "#3498DB"  # 赤=Up, 青=Down
)

g = sns.clustermap(
    df.T,                   # 転置: サンプル（行）× タンパク質（列）
    method="ward",          # ウォード法
    metric="euclidean",     # ユークリッド距離
    cmap="RdBu_r",          # Red-Blue (reversed): 赤=高発現、青=低発現
    z_score=1,              # 列方向（タンパク質方向）で標準化
    vmin=-3, vmax=3,        # Z-scoreの±3範囲
)
```

**目的**: 教師なし（群ラベルを使わない）でサンプル・タンパク質を同時クラスタリングし、生物学的な発現パターンを発見する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `df.T` | データフレーム転置（行列入れ替え） |
| `method="ward"` | ウォード法（クラスタ内分散最小化） |
| `z_score=1` | 列方向（タンパク質方向）でZ-score標準化 |
| `fc.apply(lambda x:...)` | 各タンパク質にfold change方向判定を適用 |

**引数の意味**:
- `z_score=1`: 各タンパク質を平均0・標準偏差1に標準化
- `cmap="RdBu_r"`: 赤（高発現）-白（平均）-青（低発現）
- `vmin=-3, vmax=3`: Z-score±3範囲をフル色域に割り当て

**つまずきやすいポイント**:
- 転置（df.T）により行・列の意味が変わる
- Z-score標準化で絶対的発現量は失われ、相対パターンのみ表示

### セクション4: 95%信頼楕円ヘルパー関数

```python
def confidence_ellipse(x, y, ax, n_std=2.0, **kwargs):
    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, **kwargs)
    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    transf = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)
```

**目的**: 2変量データ（PCAスコア等）の散らばりを楕円で表現し、群の広がりと重なりを視覚化する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `np.cov(x, y)` | 2変量の共分散行列を計算 |
| `matplotlib.patches.Ellipse` | 楕円パッチオブジェクト |
| `transforms.Affine2D()` | 2次元アフィン変換（回転・拡大・平行移動） |

**引数の意味**:
- `n_std=2.0`: 標準偏差の2倍（約95%信頼領域）
- `pearson`: ピアソン相関係数（楕円の傾き）
- `scale_x/scale_y`: 各軸方向の楕円半径

### セクション5: PCA散布図（Figure 1c）

```python
# PCA の実行
pca = PCA(n_components=2)
scores = pca.fit_transform(df.T)

# ステージ別カラー（論文準拠）
stage_colors = {
    "Normal": "#3498DB",   # 青（Non-tumor tissue）
    "I":      "#66BB6A",   # 緑（Stage I）
    "II":     "#FFD54F",   # 黄（Stage II / IIA）
    "III":    "#FFA726",   # オレンジ（Stage III / IIIB）
    "IV":     "#AB47BC",   # 紫（Stage IV / IVC）
}

ax.scatter(
    scores[i, 0],              # X座標: PC1スコア
    scores[i, 1],              # Y座標: PC2スコア
    c=color,                   # 点の色（ステージ別）
    s=120,                     # 点のサイズ
    alpha=0.85,                # 透過度
)

ax.set_xlabel(f"Component 1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"Component 2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
```

**目的**: 高次元プロテオームデータを2次元に圧縮し、Normal・Tumor・ステージ別の分離パターンを可視化する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `PCA(n_components=2)` | 主成分分析オブジェクト（上位2成分） |
| `fit_transform(df.T)` | PCAの学習と変換を同時実行 |
| `explained_variance_ratio_` | 各主成分の寄与率（情報保持率） |
| `ax.scatter()` | 散布図プロット |

**引数の意味**:
- `n_components=2`: PC1・PC2の2次元で可視化
- `df.T`: サンプル（行）×タンパク質（列）に転置
- `s=120`: マーカーサイズ（論文原図に合わせ大きめ）
- `alpha=0.85`: 透明度（重複時の視認性向上）

**つまずきやすいポイント**:
- PCAは入力データの転置（df.T）が必要
- 寄与率が低いと2D可視化の信頼性が下がる

## 【深掘り】 ピアソン相関係数とは？

- **ひとことで**: 2つのサンプルの発現パターンが「どれくらい似ているか」を-1から1の数値で表す指標
- **定義**: 2つの変数間の線形関係の強さを測る統計量。0に近いほど無関係、1に近いほど正の相関、-1に近いほど負の相関
- **プロテオミクスでの意味**: 相関係数が高い（0.8以上）なら、2つのサンプルで似たタンパク質が似たレベルで発現している
- **どんなとき使う**: サンプル品質確認、バッチ効果検出、群間差の事前評価
- **参考**: 論文Methods「Pearson correlation analysis」

## 【深掘り】 階層的クラスタリングとは？

- **ひとことで**: データを「似ている度合い」で自動的にグループ分けし、ツリー状の関係図（デンドログラム）を作る手法
- **定義**: 距離行列に基づいて、近い（類似する）サンプル同士から順番に結合していく教師なし分類
- **ウォード法**: クラスタ結合時にクラスタ内分散増加が最小になるペアを優先選択
- **「教師なし」の意味**: Normal/Tumorのラベル情報を使わずに、発現データのみでクラスタリング
- **期待**: 生物学的に意味のある群が自然に分離すれば、データに真の群間差が存在する証拠
- **参考**: 論文Methods「Unsupervised hierarchical clustering」

## 【深掘り】 Z-score標準化とは？

- **ひとことで**: 各タンパク質の発現値を「平均0、標準偏差1」に変換し、相対的な高低パターンで比較可能にする変換
- **定義**: z = (値 - 平均) / 標準偏差。タンパク質間の絶対発現量の違いを除去
- **計算例**: タンパク質Aが平均100、SD=20なら、値120はz=(120-100)/20=1.0になる
- **どんなとき使う**: 異なる発現レベルのタンパク質をヒートマップで同等に視覚化したい場合
- **落とし穴**: 絶対的な発現量情報は失われ、相対パターンのみ表現される
- **参考**: 論文Methods「Z-score normalization」

## 【深掘り】 主成分分析（PCA）とは？

- **ひとことで**: 数千次元の高次元データを情報をなるべく保ったまま2-3次元に圧縮し、散布図で可視化可能にする手法
- **定義**: データの分散が最大になる方向（主成分）を順番に見つけ、その軸で新しい座標系を構築
- **PC1/PC2**: 第1主成分（最大分散方向）・第2主成分（PC1と直交する次の最大分散方向）
- **寄与率**: 各主成分が元データのどれだけの情報を保持しているかの割合
- **どんなとき使う**: 高次元データの全体構造把握、群の分離確認、外れ値検出
- **期待**: PC1・PC2だけで元データの60%以上を説明できれば、2D可視化が元データをよく反映
- **参考**: 論文Methods「Principal component analysis」

## よくある質問

**Q: 相関行列で何が分かるの？**
A: サンプル間の発現パターンの類似度です。Normal同士・Tumor同士の相関が高く、Normal-Tumor間の相関が低ければ、2群に異なる発現パターンがあることの証拠になります。

**Q: クラスタリングで「教師なし」の意味は？**
A: Normal/Tumorのラベル情報を一切使わず、発現データのみでサンプルを分類することです。それでも群が自然に分離すれば、ラベルに依存しない真の生物学的差があることを示します。

**Q: PCAの寄与率はどれくらいあれば良い？**
A: PC1+PC2で50%以上なら可視化に信頼性があります。70%以上なら非常に良好です。30%未満だと2D可視化では元データの情報が大幅に失われています。

**Q: Z-score標準化は必要？**
A: ヒートマップでは必要です。タンパク質ごとに発現レベルが桁違いに異なるため、標準化しないと高発現タンパク質のパターンのみが支配的になり、低発現タンパク質の重要なパターンが見えなくなります。

**Q: 色分けは論文通りでないとダメ？**
A: 再現性の観点から論文と同じ色を使うことを推奨します。読者が論文と比較しやすく、一貫性のある視覚表現になります。

## 出力される図の内容

### Figure 1a: 相関行列ヒートマップ
- **色**: 薄い（低相関）～濃い赤（高相関）
- **カラーバー**: 行・列の端にNormal（青）・Tumor（赤）
- **デンドログラム**: 似たサンプル同士が自動的に近くに配置

### Figure 1b: 階層的クラスタリング
- **色**: 赤（高発現）-白（平均）-青（低発現）
- **行**: サンプル（Normal・Tumor別にクラスター化）
- **列**: タンパク質（発現パターン別にクラスター化）
- **カラーバー**: サンプル群（行）・発現方向（列）

### Figure 1c: PCA散布図
- **軸**: PC1（最大分散方向）・PC2（次の最大分散方向）
- **色**: ステージ別（Normal=青、Stage I=緑、II=黄、III=橙、IV=紫）
- **寄与率**: 軸ラベルに各主成分の情報保持率を表示

## 重要な注意点

- **データ品質**: 前処理が不十分だとクラスタリング結果が不安定になる
- **サンプル数**: 少なすぎる（n<10）とPCA・クラスタリングの信頼性が低下
- **メモリ使用量**: 大規模データ（1万タンパク質超）では数GB必要
- **図の解像度**: dpi=150は仮設定、論文投稿時は300-600に上げる

## 次のステップ

この可視化結果は以下の解析で参照されます：
- 外れ値サンプルの特定と除外判断
- Step 6: `step_06_differential_expression.py`での群設定妥当性確認
- 論文Figure 1との視覚的比較検証