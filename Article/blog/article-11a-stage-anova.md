---
title: "大腸がんステージ別タンパク質変動をANOVAで検出する【論文再現シリーズ #11a】"
emoji: "📊"
type: "tech"
topics: ["proteomics", "statistics", "anova", "labcode"]
published: false
---

# 大腸がんステージ別タンパク質変動をANOVAで検出する

## はじめに

大腸がんはStage I〜IVに分類されます。この記事では、**One-way ANOVA**を使用して各ステージ間でタンパク質発現量に統計的な差があるかを検定し、疾患進行に関連するタンパク質を特定します。

> **📊 INFO**
>
**この記事で行う処理**
大腸がんのステージ（I〜IV）と正常組織の5群間で、各タンパク質の発現量に差があるかをOne-way ANOVAで検定し、FDR補正で偽陽性を制御します。有意差が検出されたタンパク質は次の記事でクラスタリング解析を行い、疾患進行パターンを可視化します。

> **📝 INFO**
>
**【用語メモ】**

- **ANOVA（Analysis of Variance、分散分析）**: 3群以上の平均を比較する統計手法。ここでは Normal / Stage I / II / III / IV の5群間で発現量が有意に異なるタンパク質を検出します。
- **FDR（False Discovery Rate、偽発見率）**: 有意と判定されたタンパク質のうち、実は偽陽性である割合の期待値。`FDR < 0.01` は「有意と判定された中で偽陽性率は1%以下」を意味し、Benjamini-Hochberg 法で補正します。

## 前提

- [#9 COSMIC照合](article-09-cosmic.md) が完了していること
- ステージ情報（Table S1）が利用可能
- **対応Notebook**: [`notebooks/step_10.ipynb`](../notebooks/step_10.ipynb) — この記事のコードをセルごとに実行できます

## コード全文（対応Notebook: step_10.ipynb）

### ライブラリと設定

```python
import re                      # 正規表現モジュール: サンプル名から重複ランサフィックスを除去するために使用
import numpy as np             # 数値計算ライブラリ: 配列操作やクラスター番号の一意抽出に使用
import pandas as pd            # データ分析ライブラリ: CSV読み込み、DataFrame操作の中心
import matplotlib.pyplot as plt  # グラフ描画ライブラリ: ヒートマップやラインプロットの描画に使用
from scipy import stats        # 科学計算ライブラリの統計モジュール: ANOVAのf_oneway関数を使用
from statsmodels.stats import multitest  # 多重検定補正ライブラリ: FDR補正を行う

# パス設定
RESULTS_DIR = "results"
FIGURES_DIR = f"{RESULTS_DIR}/figures"

# パラメータ設定
FDR_THRESHOLD = 0.01  # FDR閾値 1%（Toyota et al.と同設定）
```

### データ読み込みとステージ付与

```python
# 前処理済みデータを読み込む（#6 前処理の出力）
df = pd.read_csv(f"{RESULTS_DIR}/preprocessed_data.csv", index_col=0)
print(f"データサイズ: {df.shape[0]} proteins × {df.shape[1]} samples")

# サンプル情報を読み込む
sample_info = pd.read_csv(f"{RESULTS_DIR}/sample_info.csv", index_col=0)

# Toyota et al. 論文のTable S1に基づいてステージ情報を付与
# 実際の臨床情報からステージを取得（簡略版マッピング）
stage_mapping = {
    'CRC01-N': 'Normal', 'CRC01-T': 'Stage_II',
    'CRC02-N': 'Normal', 'CRC02-T': 'Stage_I',
    'CRC03-N': 'Normal', 'CRC03-T': 'Stage_II',
    'CRC04-N': 'Normal', 'CRC04-T': 'Stage_III',
    'CRC05-N': 'Normal', 'CRC05-T': 'Stage_II',
    'CRC06-N': 'Normal', 'CRC06-T': 'Stage_I',
    'CRC07-N': 'Normal', 'CRC07-T': 'Stage_III',
    'CRC08-N': 'Normal', 'CRC08-T': 'Stage_IV',
    'CRC09-N': 'Normal', 'CRC09-T': 'Stage_II',
    'CRC10-N': 'Normal', 'CRC10-T': 'Stage_III',
    'CRC11-N': 'Normal', 'CRC11-T': 'Stage_I',
    'CRC12-N': 'Normal', 'CRC12-T': 'Stage_IV',
    'CRC13-N': 'Normal', 'CRC13-T': 'Stage_II',
    'CRC14-N': 'Normal', 'CRC14-T': 'Stage_III',
    'CRC15-N': 'Normal', 'CRC15-T': 'Stage_II',
    'CRC16-N': 'Normal', 'CRC16-T': 'Stage_I'
}

# サンプル名からサフィックス（_rep1など）を除去してクリーンな名前を取得
def clean_sample_name(name):
    """サンプル名から重複実行サフィックスを除去する"""
    return re.sub(r'_rep\d+$', '', name)  # _rep1, _rep2などを除去

# ステージ情報を付与
sample_stage = []
for sample in df.columns:
    clean_name = clean_sample_name(sample)
    stage = stage_mapping.get(clean_name, 'Unknown')
    sample_stage.append(stage)

# ステージ別サンプル数の確認
stage_counts = pd.Series(sample_stage).value_counts()
print("ステージ別サンプル数:")
print(stage_counts)
```

### One-way ANOVA

```python
def perform_anova_analysis(data, stages):
    """
    各タンパク質について5群間（Normal, Stage I-IV）でOne-way ANOVAを実行

    Args:
        data: タンパク質×サンプルの発現量DataFrame
        stages: 各サンプルのステージ情報のlist

    Returns:
        結果DataFrame（タンパク質名、F統計量、p値、補正p値）
    """
    results = []

    # ステージごとにサンプルインデックスを分類
    stage_groups = {}
    for i, stage in enumerate(stages):
        if stage not in stage_groups:
            stage_groups[stage] = []
        stage_groups[stage].append(i)

    print(f"解析対象群: {list(stage_groups.keys())}")

    # 各タンパク質について統計検定を実行
    for protein in data.index:
        try:
            # 各ステージのデータを取得
            group_data = []
            for stage in ['Normal', 'Stage_I', 'Stage_II', 'Stage_III', 'Stage_IV']:
                if stage in stage_groups:
                    indices = stage_groups[stage]
                    values = data.loc[protein, data.columns[indices]].dropna()  # 欠損値除去
                    if len(values) > 1:  # 最低2サンプル以上必要
                        group_data.append(values)

            # ANOVAを実行（3群以上必要）
            if len(group_data) >= 3:
                f_stat, p_val = stats.f_oneway(*group_data)
                results.append({
                    'protein': protein,
                    'f_statistic': f_stat,
                    'p_value': p_val
                })
            else:
                # データ不足の場合はNaN
                results.append({
                    'protein': protein,
                    'f_statistic': np.nan,
                    'p_value': np.nan
                })

        except Exception as e:
            # エラーが発生した場合もNaN
            results.append({
                'protein': protein,
                'f_statistic': np.nan,
                'p_value': np.nan
            })

    return pd.DataFrame(results).set_index('protein')

# ANOVA解析の実行
print("One-way ANOVA解析を実行中...")
anova_results = perform_anova_analysis(df, sample_stage)

# 有効な結果のみを抽出（NaNを除去）
valid_results = anova_results.dropna()
print(f"統計検定完了: {len(valid_results)} タンパク質")
print(f"有効p値の範囲: {valid_results['p_value'].min():.2e} - {valid_results['p_value'].max():.2e}")
```

### FDR補正とフィルタリング

```python
# Benjamini-Hochberg法によるFDR補正
if len(valid_results) > 0:
    rejected, p_adjusted, _, _ = multitest.multipletests(
        valid_results['p_value'],
        alpha=FDR_THRESHOLD,
        method='fdr_bh'
    )

    # 結果をデータフレームに追加
    valid_results['p_adjusted'] = p_adjusted
    valid_results['significant'] = rejected

    # 有意なタンパク質を抽出
    significant_proteins = valid_results[valid_results['significant']]

    print(f"FDR < {FDR_THRESHOLD} の有意タンパク質: {len(significant_proteins)} 個")
    print(f"全体に占める割合: {len(significant_proteins)/len(valid_results)*100:.1f}%")

    # 結果を保存
    valid_results.to_csv(f"{RESULTS_DIR}/anova_stage_results.csv")
    significant_proteins.to_csv(f"{RESULTS_DIR}/anova_stage_significant.csv")

    print(f"結果を保存: anova_stage_results.csv ({len(valid_results)} proteins)")
    print(f"有意差タンパク質を保存: anova_stage_significant.csv ({len(significant_proteins)} proteins)")
else:
    print("エラー: 有効な統計検定結果がありません")
```

## コード詳細

### t検定 vs ANOVA の使い分け

| 比較 | t検定 | ANOVA |
|------|-------|-------|
| **群数** | 2群のみ | 3群以上 |
| **本シリーズでの用途** | Normal vs Tumor（#8差分発現） | Normal vs Stage I-IV（本記事） |
| **帰無仮説** | 2群の平均が等しい | 全群の平均が等しい |
| **検定統計量** | t統計量 | F統計量 |

### FDR補正の重要性

多重検定補正なしでα=0.05を使用すると、1,000タンパク質で約50個の偽陽性が発生します。FDR補正により偽陽性率を1%以下に制御し、信頼性の高い発見を保証します。

## 実行方法

**スクリプトで一括実行する場合:**
```bash
python scripts/step_10_stage_analysis.py
```

**Notebook でセルごとに実行する場合:**
```bash
jupyter notebook notebooks/step_10.ipynb
```

## まとめ

One-way ANOVAにより、大腸がんのステージ進行に伴って有意に変動するタンパク質を統計的に特定しました。この結果は次の記事でクラスタリング解析を行い、ステージ進行のパターンを可視化します。

> 前回: [#9 COSMIC照合](article-09-cosmic.md)
> 次回: [#11b ステージ別クラスタリング可視化](article-11b-stage-clustering.md) — 有意差タンパク質の発現パターン解析

#バイオインフォマティクス #プロテオミクス #統計解析 #ANOVA #labcode