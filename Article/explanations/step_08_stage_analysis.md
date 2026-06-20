# step_08_stage_analysis.py の解説

## このコードの役割

がんの**臨床ステージ（I〜IV期）**ごとにタンパク質発現パターンの変化を統計的に解析し、論文Figure 3を再現するスクリプトです。One-way ANOVA（一元配置分散分析）によるステージ間差の検定、階層的クラスタリング、複合パネル可視化を行います。

## 全体フロー

```
前処理データ＋臨床情報 → ステージ別中央値算出 → One-way ANOVA → FDR補正 → 階層クラスタリング → Figure3可視化
```

| | ファイル | 内容 |
|-|---------|------|
| 入力 | `results/preprocessed_data.csv` | 前処理済みタンパク質発現データ |
| 入力 | `results/sample_info.csv` | サンプル情報（条件・患者情報含む） |
| 入力 | `results/clinical_info.csv` | 臨床情報（ステージ分類、存在する場合） |
| 出力 | `results/tables/anova_results.csv` | ANOVA検定結果（F統計量・p値・FDR値） |
| 出力 | `results/tables/cluster_assignments.csv` | 階層クラスタリング結果 |
| 出力 | `results/figures/fig3_stage_heatmap.png` | Figure 3: ヒートマップ+ラインプロット複合パネル |

**所要時間の目安**: 約2-3分（ANOVA計算・クラスタリング・可視化）

## セクション別解説

### セクション1: データ読み込みとステージ情報の付与

```python
def load_data():
    df = pd.read_csv(os.path.join(RESULTS_DIR, "preprocessed_data.csv"), index_col=0)
    sample_info = pd.read_csv(os.path.join(RESULTS_DIR, "sample_info.csv"))
    return df, sample_info

def assign_stages(sample_info):
    clinical_path = os.path.join(RESULTS_DIR, "clinical_info.csv")
    if os.path.exists(clinical_path):
        clinical = pd.read_csv(clinical_path)
        # サンプル名 → ステージのマッピング
        stage_map = {}
        for _, row in clinical.iterrows():
            stage_map[row["Sample_N"]] = "Normal"
            stage_map[row["Sample_T"]] = row["Stage"]
        sample_info["Stage"] = sample_info["Sample"].apply(
            lambda x: stage_map.get(re.sub(r"_dup\d+$", "", x))
        )
```

**目的**: 前処理済みデータを読み込み、各サンプルにがんの臨床ステージ（Normal/I/II/III/IV）を割り当てる。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `pd.read_csv(index_col=0)` | 1列目（タンパク質名）をインデックスとして読み込み |
| `os.path.exists()` | 臨床情報ファイルの存在確認 |
| `re.sub(r"_dup\d+$", "", x)` | 重複Gene Symbol接尾語（_dup1等）を除去 |

**引数の意味**:
- `clinical_info.csv`: 各患者のステージ分類を記録するファイル
- `Sample_N`: 正常組織サンプル名
- `Sample_T`: 腫瘍組織サンプル名
- `Stage`: I, II, III, IV の臨床ステージ分類

**つまずきやすいポイント**:
- 臨床情報がない場合はCondition列（Normal/Tumor）で代替される
- 重複Gene Symbol処理で付与された`_dupN`接尾語を除去してマッピングする必要がある

### セクション2: One-way ANOVA（一元配置分散分析）

```python
def run_anova(df, sample_info):
    stage_order = ["Normal", "I", "II", "III", "IV"]
    available = [s for s in sample_info["Stage"].unique() if pd.notna(s)]
    stages_with_normal = [s for s in stage_order if s in available]

    results = []
    for protein in df.index:
        groups = []
        for stage in stages_with_normal:
            samples = sample_info[sample_info["Stage"] == stage]["Sample"].tolist()
            vals = df.loc[protein, [s for s in samples if s in df.columns]].dropna()
            if len(vals) >= 2:
                groups.append(vals.values)

        if len(groups) < 2:
            continue

        # One-way ANOVAの実行
        f_stat, p_val = stats.f_oneway(*groups)
        results.append({
            "Protein": protein,
            "F_statistic": f_stat,
            "P_value": p_val,
        })
```

**目的**: 各タンパク質について、5つのステージ（Normal, I, II, III, IV）間で発現量に有意差があるかをOne-way ANOVAで検定する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `stats.f_oneway(*groups)` | 3群以上の平均値の差を一度に検定（F統計量とp値を返す） |
| `dropna()` | 欠損値を除外してクリーンなデータで検定 |
| `len(vals) >= 2` | 分散計算に必要な最低サンプル数の確認 |

**引数の意味**:
- `*groups`: リストを展開して各ステージのデータを引数として渡す
- `F_statistic`: グループ間ばらつき÷グループ内ばらつき（大きいほど有意差大）
- `P_value`: 帰無仮説（全群の平均が等しい）が正しい確率

**つまずきやすいポイント**:
- ANOVAはt検定と違い3群以上を一度に比較できる（多重比較問題の回避）
- 各グループに最低2サンプルが必要（分散計算のため）
- ステージの生物学的順序（Normal→IV）を保持することが重要

### セクション3: FDR補正（多重検定補正）

```python
# FDR補正（Benjamini-Hochberg法）
if len(result_df) > 0:
    reject, fdr, _, _ = multipletests(result_df["P_value"], method="fdr_bh")
    result_df["FDR"] = fdr
    result_df["Significant"] = fdr < FDR_THRESHOLD  # 0.01
```

**目的**: 数千個のタンパク質を同時検定することによる多重比較問題を解決し、偽陽性率を制御する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `multipletests(method="fdr_bh")` | Benjamini-Hochberg法でFDR調整済みp値を計算 |
| `FDR_THRESHOLD = 0.01` | 偽発見率1%以下の厳しい閾値設定 |

**引数の意味**:
- `fdr_bh`: Benjamini-Hochberg法（最も広く使われるFDR補正法）
- `reject`: 帰無仮説を棄却するかどうかのブール配列
- `fdr`: FDR補正後のp値（調整済みp値）

**つまずきやすいポイント**:
- 生のp値ではなく必ずFDR補正後の値で判定する
- 0.01という閾値は一般的な0.05より厳しい設定（高精度を優先）

### セクション4: 階層的クラスタリング

```python
def cluster_and_plot(median_df, anova_df, n_clusters=N_CLUSTERS):
    # 有意なタンパク質の抽出
    sig_proteins = anova_df[anova_df["Significant"]]["Protein"].tolist()
    sig_data = median_df.loc[median_df.index.isin(sig_proteins)]

    # Z-score正規化（行方向 = 各タンパク質ごと）
    z_data = sig_data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    z_data = z_data.dropna()

    # 階層的クラスタリングの実行
    Z = linkage(z_data.values, method="ward")
    clusters = fcluster(Z, t=actual_clusters, criterion="maxclust")
```

**目的**: ANOVAで有意と判定されたタンパク質を、発現パターンの類似性に基づいて30個のクラスターに分類する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `apply(lambda x: (x - x.mean()) / x.std(), axis=1)` | 各タンパク質を行方向にZ-score正規化 |
| `linkage(method="ward")` | Ward法で階層的クラスタリングの連結行列を計算 |
| `fcluster(criterion="maxclust")` | 指定数のクラスターに分割 |

**引数の意味**:
- `axis=1`: 行方向（各タンパク質の全ステージ値）でZ-score計算
- `method="ward"`: クラスター内分散を最小化する連結法
- `t=actual_clusters`: 目標クラスター数（30、データ数で制限）

**つまずきやすいポイント**:
- Z-score正規化により「発現量の絶対値」ではなく「パターンの形」で分類される
- Ward法は均一で解釈しやすいクラスターを作りやすい

### セクション5: Figure 3複合パネルの作成

```python
def plot_figure3(df, sample_info, z_data, clusters):
    # 論文の補足データ（Table S9-S12）からクラスター情報を読み込み
    paper_clusters = load_paper_clusters(RAW_DIR)

    # 2x2の4パネル構成
    fig = plt.figure(figsize=(16, 20))
    outer_gs = GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.3)

    for idx, info in enumerate(selected):
        # 各パネルの内部レイアウト: 上段ヒートマップ、下段ラインプロット
        inner_gs = GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer_gs[idx],
            height_ratios=[3, 1], hspace=0.15
        )

        # ヒートマップ: ステージ別中央値をZ-score化
        cluster_z = median_matrix.apply(
            lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x * 0, axis=1
        )

        sns.heatmap(
            cluster_z,
            cmap="RdYlGn_r",          # 論文準拠カラーマップ
            center=0, vmin=-2, vmax=2,
            xticklabels=True, yticklabels=False,
            cbar_kws={"shrink": 0.5, "label": "Z-score"},
        )
```

**目的**: 論文Figure 3を再現する2×2の複合パネル（各パネル＝ヒートマップ+ラインプロット）を作成する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `GridSpec(2, 2)` | 2×2の外側レイアウト（4パネル）を作成 |
| `GridSpecFromSubplotSpec(2, 1, height_ratios=[3, 1])` | 各パネル内の上下分割（ヒートマップ:ラインプロット=3:1） |
| `sns.heatmap(cmap="RdYlGn_r")` | 赤=高発現、緑=低発現の論文準拠カラーマップ |
| `load_paper_clusters()` | 論文補足データ（Table S9-S12）のクラスター情報を読み込み |

**引数の意味**:
- `height_ratios=[3, 1]`: ヒートマップを3倍の高さ、ラインプロットを1倍の高さに設定
- `center=0, vmin=-2, vmax=2`: Z-scoreの標準的な表示範囲
- `"RdYlGn_r"`: Red-Yellow-Green reversed（赤=高、黄=中、緑=低）

**つまずきやすいポイント**:
- 完全再現のため論文の補足データ（Table S9-S12）を読み込む仕組みが必要
- 4つのクラスター（増加2個、減少2個）を適切に選択・配置する必要がある

## 【深掘り】 One-way ANOVA とは？

- **ひとことで**: 3群以上のグループ間で平均値に差があるかを1回の検定で判定する統計手法
- **定義**: F統計量（グループ間ばらつき÷グループ内ばらつき）を計算し、F分布と比較して有意性を判定
- **t検定との違い**: t検定は2群比較のみ、ANOVAは3群以上を一度に比較（多重比較問題を回避）
- **F統計量の意味**: 大きいほど「グループ間の差が個体差を上回る」ことを示す
- **がん研究での使いどころ**: Normal vs Stage I vs Stage II vs Stage III vs Stage IV の5群比較
- **参考**: Toyota et al. 2025 Methods "One-way ANOVA followed by Benjamini-Hochberg correction"

## 【深掘り】 FDR補正（偽発見率制御）とは？

- **ひとことで**: 大量の検定を行う際に、偽陽性の割合を統計的に制御する手法
- **多重比較問題**: 10,000個のタンパク質でp<0.05検定すると、差がなくても約500個が偶然有意になる
- **FDRの定義**: 有意と判定したもののうち、本当は差がないものの割合（False Discovery Rate）
- **Benjamini-Hochberg法**: p値を小さい順に並べ、順位に応じた閾値で段階的に判定
- **FDR<0.01の意味**: 有意判定した100個のうち1個以下が偽陽性という厳しい基準
- **なぜ必要**: プロテオミクスでは数千〜万のタンパク質を同時解析するため必須
- **参考**: Benjamini & Hochberg (1995) "Controlling the False Discovery Rate"

## 【深掘り】 Z-score正規化とは？

- **ひとことで**: 各データを「平均0、標準偏差1」に変換して異なる尺度のデータを統一する手法
- **計算式**: Z = (値 - 平均) / 標準偏差
- **目的**: 発現量の絶対値ではなく「発現パターンの形状」でタンパク質を分類・比較
- **例**: 100→200（2倍増加）と1→2（2倍増加）は同じZ-scoreパターンになる
- **行方向vs列方向**: 行方向=各タンパク質の全ステージ値で正規化、列方向=各ステージの全タンパク質で正規化
- **ヒートマップでの効果**: 赤/青のコントラストが発現変化の方向性を直感的に表現
- **参考**: プロテオミクス解析における標準的な前処理手法

## 【深掘り】 階層的クラスタリング（Ward法）とは？

- **ひとことで**: 似たパターンを持つデータを段階的にグループ化し、樹形図（デンドログラム）で表現する手法
- **ボトムアップ方式**: 各データが独立クラスター → 最も似たペアを結合 → 指定数まで繰り返し
- **Ward法の特徴**: クラスター内の分散（ばらつき）を最小化するように結合する方法
- **距離の計算**: ユークリッド距離でタンパク質間の「発現パターンの近さ」を測定
- **linkage関数**: (n-1)×4の連結行列を返し、各行が1回の結合操作を記録
- **fcluster関数**: 連結行列から指定数のクラスターに「切断」して分割
- **プロテオミクスでの意義**: 同じ生物学的機能・経路のタンパク質が同じクラスターに集まりやすい

## 【深掘り】 臨床ステージとは？

- **ひとことで**: がんの進行度を表す国際的な分類システム（TNM分類に基づく）
- **Stage I**: 早期がん（腫瘍が小さく、リンパ節転移なし）
- **Stage II**: 中期がん（腫瘍がやや大きい、または限定的なリンパ節転移）
- **Stage III**: 進行がん（周囲組織への浸潤、広範囲なリンパ節転移）
- **Stage IV**: 末期がん（遠隔転移あり、他臓器への転移）
- **なぜ重要**: ステージが上がるほど予後が悪く、治療法も異なる
- **プロテオミクス解析での意義**: ステージ進行に伴う分子変化を追跡し、バイオマーカーを発見
- **本論文での知見**: Cluster 3/14（増加）とCluster 20/25（減少）でステージ特異的変化を確認

## よくある質問

**Q: なぜt検定ではなくOne-way ANOVAを使うのですか？**
A: 5つのステージ（Normal, I, II, III, IV）を比較するため。t検定で全ペア比較すると10回（5C2）の検定が必要で多重比較問題が深刻化しますが、ANOVAなら1回の検定で済みます。

**Q: FDR 0.01は厳しすぎませんか？**
A: プロテオミクスでは数千のタンパク質を同時解析するため、厳しい基準が必要です。一般的な0.05では偽陽性が多くなりすぎます。高精度な候補選定を優先した設定です。

**Q: Z-score正規化はなぜ行方向（axis=1）なのですか？**
A: 各タンパク質の「ステージ間での変化パターン」を抽出するため。列方向だと各ステージ内でのタンパク質間比較になり、ステージ進行の情報が失われます。

**Q: 論文のクラスターデータ（Table S9-S12）がない場合はどうなりますか？**
A: 自前のクラスタリング結果から「Normal→Stage IV」の変化量（trend）が最大・最小の2クラスターずつを選んでフォールバック表示します。

**Q: 臨床情報（clinical_info.csv）がない場合はどうなりますか？**
A: sample_info.csvのCondition列（Normal/Tumor）をステージ代替として使用し、2群比較の解析になります。詳細なステージ別解析は行えません。

**Q: なぜクラスター数を30に設定しているのですか？**
A: 論文の解析パラメータに準拠した設定です。生物学的に意味のある分類の粒度と計算効率のバランスを考慮した値です。

## 出力されるファイルの内容

### anova_results.csv
```
Protein,F_statistic,P_value,FDR,Significant
TP53,15.234,0.00001,0.0023,True
EGFR,8.452,0.00089,0.0156,False
MYC,22.156,0.00000,0.0001,True
```

- **F_statistic**: F統計量（大きいほどステージ間差が顕著）
- **P_value**: 生のp値（FDR補正前）
- **FDR**: FDR補正後のp値（多重比較を考慮）
- **Significant**: FDR < 0.01 かどうかのブール値

### cluster_assignments.csv
```
Protein,Cluster
TP53,14
EGFR,3
MYC,20
```

- **Protein**: タンパク質名（Gene Symbol）
- **Cluster**: 階層的クラスタリングで割り当てられたクラスター番号（1〜30）

### fig3_stage_heatmap.png
- **レイアウト**: 2×2の4パネル構成
- **各パネル**: 上段ヒートマップ + 下段ラインプロット
- **パネル分類**: (a)(b)増加クラスター、(c)(d)減少クラスター
- **カラーマップ**: RdYlGn_r（赤=高発現、緑=低発現）
- **ステージ色**: Normal→青、I→緑、II→黄、III→オレンジ、IV→赤

## 重要な注意点

- **論文補足データ依存**: Table S9-S12がある場合とない場合で結果表示が変わる
- **ステージ情報の必要性**: clinical_info.csvがないと詳細なステージ別解析ができない
- **FDR閾値の厳しさ**: 0.01は高精度だが、有意タンパク質数が少なくなる可能性がある
- **計算コスト**: 数千タンパク質のANOVA + クラスタリングで2-3分の処理時間
- **メモリ使用量**: Z-score計算・クラスタリング・可視化で数GBのメモリが必要

## 次のステップ

この解析結果は以下で利用されます：
- 論文Discussion: ステージ特異的タンパク質変化の生物学的解釈
- バイオマーカー候補の臨床的意義評価
- 創薬標的探索への応用（ステージ早期で変化するタンパク質の同定）
- 他のがん種との比較解析