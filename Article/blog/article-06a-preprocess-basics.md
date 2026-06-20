---
title: "Pythonでプロテオミクスデータの基本前処理【論文再現シリーズ #6a】"
emoji: "🧹"
type: "tech"
topics: ["proteomics", "python", "preprocessing", "labcode"]
published: false
---

# Pythonでプロテオミクスデータの基本前処理

## はじめに

プロテオミクスデータの前処理は、多くの論文でPerseus（MaxQuant付属ソフト）が使われます。しかしPerseusは**商用利用に制限がある**ため、本シリーズではPythonで完全代替します。この記事では、Log2変換・欠損値フィルタリングまでの基本前処理をPythonで実装します。

> **📝 INFO**
>
**この記事で行う処理**
質量分析の生の強度値は桁が大きく分布が歪んでいるため、統計解析に適した形に整えます。具体的には、①Log2変換で正規分布に近づけ、②有効値が少なすぎるタンパク質を除去（70%ルール）します。これにより、統計検定や可視化に適したクリーンなデータが準備できます。


## 前提

- [#1 環境構築](article-01-setup.md) が完了していること
- [#4b sage実行](article-04b-sage-execution.md) で `results/protein_matrix_from_sage.csv` が生成されていること
- **対応Notebook**: [`notebooks/step_06.ipynb`](../notebooks/step_06.ipynb) — この記事のコードをセルごとに実行できます

## 本書で扱うデータ規模

本シリーズでは論文データの **32 ファイル（16 患者分: CRC01-CRC16、各 Normal/Tumor）** を扱います。sage による解析の結果として以下のデータが入力になります：

| 項目 | 値 |
|------|-----|
| タンパク質数 | **2,110** |
| サンプル数 | **32** (Normal 16 + Tumor 16) |
| 論文全体のタンパク質数 | 10,329 (参考) |

前処理後：

| 段階 | タンパク質数 |
|------|-------------|
| sage 出力 | 2,110 |
| Log2 変換 | 2,110 |
| 70% 有効値フィルタ | **2,081**（29 除去） |

## 前処理の全体像

```
DIA解析の出力（タンパク質定量マトリクス）
  → Log2変換（正規分布に近づける）
  → 欠損値フィルタリング（70%ルール）
  → 欠損値補完（Perseus互換: downshift法）
  → 前処理済みデータ
```

## コード全文（対応Notebook: step_06.ipynb）

以下が基本前処理コードです。Perseusの主要機能をPythonで再現しています。Notebook版（`notebooks/step_06.ipynb`）ではセルごとに実行できます。

```python
import os       # ファイルパスの結合・操作に使う標準ライブラリ
import numpy as np   # 数値計算ライブラリ（log2変換・乱数生成などに使用）
import pandas as pd  # データフレーム操作ライブラリ（CSV読み書き・欠損値処理などに使用）

# --- パス ---
# 結果ファイルを格納するディレクトリへの相対パス（notebookから実行する場合の基準）
RESULTS_DIR = "../results"
# 前処理済みデータのCSVファイルパス（step_04のsageで作成したタンパク質定量マトリクス）
INPUT_CSV = os.path.join(RESULTS_DIR, "protein_matrix_from_sage.csv")

# --- Perseus互換パラメータ ---
# 有効値の最低割合（論文の設定: 70%）
# この割合未満の有効値しかないタンパク質は信頼性が低いため除去する
VALID_RATIO = 0.70
```

### 1. データ読み込み

```python
# CSVファイルを読み込み、1列目（タンパク質名）をインデックスに設定する
df = pd.read_csv(INPUT_CSV, index_col=0)
# インデックスの名前を"Protein"に設定（後の処理で参照しやすくするため）
df.index.name = "Protein"

# サンプル名の末尾 -N / -T で Normal（正常組織） / Tumor（腫瘍組織） を判別
# リスト内包表記で、列名に"-N"を含むサンプル名だけを抽出する
normal_samples = [c for c in df.columns if "-N" in c]
# 列名に"-T"を含むサンプル名だけを抽出する
tumor_samples  = [c for c in df.columns if "-T" in c]
# Normal群とTumor群をまとめた辞書を作成（後のフィルタリングで群ごとに処理するため）
groups = {"Normal": normal_samples, "Tumor": tumor_samples}

# 読み込んだデータの概要を表示（行数=タンパク質数、列数=サンプル数）
print(f"読み込み: {df.shape[0]} タンパク質 × {df.shape[1]} サンプル")
# 各群のサンプル数を表示して、正しく分類されたか確認する
print(f"  Normal: {len(normal_samples)}, Tumor: {len(tumor_samples)}")
```

### 2. Log2変換

```python
# 質量分析の強度値は 10^6〜10^9 と桁が大きく分布が偏る。
# Log2変換で正規分布に近づけ、差が「倍率」に対応するようになる（差1 = 2倍変化）。

# 全サンプルの中央値の中央値を計算し、データが既にlog2スケールか判定する基準にする
# df.median()で各列（サンプル）の中央値を取得し、さらにその中央値を取る
median_val = df.median().median()

# 中央値が100より大きい場合は生の強度値と判断してLog2変換を実行する
if median_val > 100:
    # 値0をNaN（欠損値）に置換してからLog2変換する（log2(0)=-∞を避けるため）
    df = np.log2(df.replace(0, np.nan))
    # 変換前の中央値を表示して、変換が正しく行われたことを確認する
    print(f"Log2変換実行（変換前の中央値: {median_val:.1f}）")
else:
    # 中央値が100以下の場合、既にlog2スケールに変換済みと判断してスキップする
    print(f"既にlog2スケール（中央値: {median_val:.1f}）→ スキップ")
```

### 3. 有効値フィルタリング（70%ルール）

```python
def filter_by_valid_ratio(df, groups, ratio=VALID_RATIO):
    """いずれかの群で有効値割合 >= ratio を満たすタンパク質を残す。

    【なぜフィルタリングが必要か？】
      - ほとんどのサンプルで検出されないタンパク質は統計的に信頼できない
      - 欠損値が多すぎると補完の精度も低下する
    """
    # 全タンパク質をFalse（除去対象）で初期化する。条件を満たしたものだけTrueにする
    keep = pd.Series(False, index=df.index)
    # 各群（Normal, Tumor）について有効値割合を確認するループ
    for samples in groups.values():
        # notna()で欠損でないセルをTrue/Falseに変換し、行方向に合計して有効値数を求める
        # それをサンプル数で割って有効値の「割合」を計算する（0.0〜1.0の範囲）
        valid = df[samples].notna().sum(axis=1) / len(samples)
        # OR演算（|=）で「いずれかの群で基準を満たせばTrue」にする
        # これにより、Normal群またはTumor群のどちらかで70%以上有効なら残す
        keep |= (valid >= ratio)
    # keepがTrueのタンパク質だけを残したDataFrameを返す
    return df[keep]

# フィルタリング前のタンパク質数を記録しておく（除去数の表示に使う）
n_before = len(df)
# 70%ルールでフィルタリングを実行し、結果をdfに上書きする
df = filter_by_valid_ratio(df, groups)
# フィルタリング結果を表示: 何個のタンパク質が残り、何個除去されたかを確認する
print(f"フィルタリング: {n_before} → {len(df)} タンパク質（{n_before - len(df)} 除去）")
```

## コード詳細

### Log2変換

質量分析データの強度値は以下の特徴があります：

| 項目 | 変換前 | Log2変換後 |
|------|---------|------------|
| **値の範囲** | 10⁶ 〜 10⁹ | 20 〜 30 |
| **分布** | 偏った分布（右に長い裾） | 正規分布に近い |
| **解釈** | 絶対強度 | **倍率**（差1=2倍変化） |

```python
# 例: 強度が 1,000,000 → log2(1,000,000) ≈ 19.9
# 倍率の解釈: log2値が 20 → 21 なら2倍増加
```

### 欠損値フィルタリング（70%ルール）

プロテオミクスデータでは、検出限界以下のタンパク質は「欠損値（NaN）」として記録されます。

| フィルタリング基準 | 説明 |
|------------------|------|
| **70%ルール** | Normal群またはTumor群のどちらかで70%以上検出されたタンパク質を残す |
| **なぜ70%？** | 統計検定に十分なサンプルサイズを確保する慣例的な基準 |
| **OR条件** | Normal群で30%しか検出されなくても、Tumor群で80%検出されれば残す |

## 実行方法

**Notebook でセルごとに実行する場合:**

```bash
jupyter notebook notebooks/step_06.ipynb
```

Notebook版では各処理の途中経過を確認しながら進められます。

## まとめ

Log2変換と有効値フィルタリングの基本前処理により、統計解析に適したデータ形式が準備できました。本書のデータでは、sage由来の2,110タンパク質から70%フィルタで2,081タンパク質が残りました。

> 前回: [#4b sage実行](article-04b-sage-execution.md)
> 次回: [#6b 欠損値補完](article-06b-preprocess-imputation.md) — Perseus互換downshift法による欠損値補完

#バイオインフォマティクス #プロテオミクス #Python #前処理 #labcode