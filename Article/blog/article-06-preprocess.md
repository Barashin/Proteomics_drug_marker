---
title: "Perseusなしで！Pythonでプロテオミクスデータを前処理する【論文再現シリーズ #6】"
emoji: "🧹"
type: "tech"
topics: ["proteomics", "python", "bioinformatics", "labcode"]
published: false
---

# Perseusなしで！Pythonでプロテオミクスデータを前処理する

## はじめに

プロテオミクスデータの前処理は、多くの論文でPerseus（MaxQuant付属ソフト）が使われます。しかしPerseusは**商用利用に制限がある**ため、本シリーズではPythonで完全代替します。この記事では、Log2変換・欠損値フィルタリング・欠損値補完をPythonで実装します。

:::message
**この記事で行う処理**
質量分析の生の強度値は桁が大きく分布が歪んでいるため、統計解析に適した形に整えます。具体的には、①Log2変換で正規分布に近づけ、②有効値が少なすぎるタンパク質を除去し（70%ルール）、③検出限界以下で欠損になった値をPerseus互換のdownshift法で補完します。これにより、欠損のないクリーンなlog2スケールのマトリクスが得られ、以降の統計検定や可視化が正しく実行できるようになります。
:::

## 前提

- [#1 環境構築](article-01-setup.md) が完了していること
- [#4 sageによるDIA解析](article-04-sage.md) または [#5 OpenMS + AlphaPeptDeep](article-05-openms.md) で `results/protein_matrix_from_*.csv` が生成されていること
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
| 欠損値補完後 | 2,081（922 個の値を補完） |

## 前処理の全体像

```
DIA解析の出力（タンパク質定量マトリクス）
  → Log2変換（正規分布に近づける）
  → 欠損値フィルタリング（70%ルール）
  → 欠損値補完（Perseus互換: downshift法）
  → 前処理済みデータ
```

## コード全文（対応Notebook: step_06.ipynb）

以下が前処理コードの全文です。Perseusの3つの主要機能をPythonで再現しています。Notebook版（`notebooks/step_06.ipynb`）ではセルごとに実行できます。

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
# downshift: 有効値の平均から何SD下にシフトするか
#   → 2.4 = 「検出されなかったタンパク質は有効値の平均より2.4SD低い」と仮定
# width: 補完値の分布幅（SDの倍率）
#   → 0.3 = 元のSDの30%の幅で補完値をばらつかせる

# Perseus互換のdownshiftパラメータ: 有効値の平均から2.4SD下にシフトして補完値の中心を決める
DOWNSHIFT = 2.4
# Perseus互換のwidthパラメータ: 補完値のばらつきを元のSDの30%に設定する
WIDTH = 0.3

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

### 4. 欠損値補完（Perseus互換 downshift法）

```python
def impute_downshift(df, downshift=DOWNSHIFT, width=WIDTH):
    """Perseus互換: 各サンプルごとに低値側から欠損を補完する。

    【補完の考え方】
      「検出されなかった」= 存在しないのではなく「低すぎて検出限界以下」
      → 検出された値の分布から推測して低い値を割り当てる

    【手順（各サンプルごと）】
      1. 有効値の平均(μ)と標準偏差(σ)を計算
      2. 補完値の中心 = μ - downshift × σ
      3. 補完値の幅 = width × σ
      4. この分布からランダムにサンプリング
    """
    # 元のDataFrameを変更しないようにコピーを作成する
    df = df.copy()
    # 補完した値の総数をカウントする変数（最後に報告用に使う）
    total = 0
    # 各サンプル（列）を順番に処理するループ
    for col in df.columns:
        # そのサンプルの有効値（NaNでない値）だけを取り出す
        valid = df[col].dropna()
        # 有効値が0個の場合は補完できないのでスキップする
        if len(valid) == 0:
            continue
        # 補完値の中心: 有効値の平均からdownshift×SD分だけ低い値にする
        # 例: 平均20, SD=2, downshift=2.4 → 中心は 20 - 2.4×2 = 15.2
        imp_mean = valid.mean() - downshift * valid.std()
        # 補完値のばらつき: 元のSDのwidth倍（0.3なら30%）に設定する
        imp_std  = width * valid.std()
        # そのサンプルでNaN（欠損値）の位置をTrue/Falseのマスクとして取得する
        mask = df[col].isna()
        # 欠損値の個数を数える
        n = mask.sum()
        # 欠損値が1個以上ある場合のみ補完を実行する
        if n > 0:
            # 正規分布N(imp_mean, imp_std)からn個の乱数を生成し、欠損位置に代入する
            df.loc[mask, col] = np.random.normal(imp_mean, imp_std, n)
            # 補完した個数を累計に加算する
            total += n
    # 補完済みDataFrameと補完値の総数を返す
    return df, total

# 再現性のため乱数シードを固定する（42は慣例的によく使われる値）
np.random.seed(42)
# downshift法で欠損値を補完し、補完後のデータと補完数を受け取る
df, n_imputed = impute_downshift(df)
# 補完した値の総数とパラメータを表示する
print(f"欠損値補完: {n_imputed} 値（downshift={DOWNSHIFT}, width={WIDTH}）")
# 補完後に残っている欠損値の数を確認する（0であれば全て補完完了）
print(f"残り欠損: {df.isna().sum().sum()}")
```

### 5. 保存

```python
# 前処理済みのタンパク質定量マトリクスをCSVファイルとして保存する
df.to_csv(os.path.join(RESULTS_DIR, "preprocessed_data.csv"))

# サンプル情報（各サンプルがNormalかTumorか）をDataFrameにまとめる
# 後続のステップ（差次解析・可視化）で群の情報が必要になるため、別ファイルとして保存する
sample_info = pd.DataFrame({
    # サンプル名のリスト: Normal群とTumor群を結合する
    "Sample": normal_samples + tumor_samples,
    # 条件ラベル: Normal群のサンプル数分"Normal"、Tumor群のサンプル数分"Tumor"を生成する
    "Condition": ["Normal"] * len(normal_samples) + ["Tumor"] * len(tumor_samples)
})
# サンプル情報をCSVファイルとして保存する（index=Falseで行番号を出力しない）
sample_info.to_csv(os.path.join(RESULTS_DIR, "sample_info.csv"), index=False)

# 保存完了のメッセージを表示して、最終的なデータサイズを確認する
print(f"保存完了: {df.shape[0]} タンパク質 × {df.shape[1]} サンプル")
print(f"  → preprocessed_data.csv")
print(f"  → sample_info.csv")
```

---

## コード詳細

### Log2変換

```python
# 値0をNaN（欠損値）に置換した後、全セルをLog2変換する（log2(0)=-∞を避けるため先にNaN化）
df = np.log2(df.replace(0, np.nan))
```

- **`df.replace(0, np.nan)`**: 値が0のセルをNaN（欠損値）に置換します。Log2(0)は-∞になってしまうため、先にNaNにしておきます。
- **`np.log2()`**: NumPyのlog2関数。各セルの値をLog2変換します。NaNはそのままNaNとして保持されます。
- **変換の意味**: 例えば値が1,000,000の場合、Log2(1,000,000) ≈ 20 になります。Log2スケールでは差が「倍率」に直結します（差1.0 = 2倍の変化）。

---

### 欠損値フィルタリング（70%ルール）

```python
# 各タンパク質について、指定群のサンプル中の有効値（非欠損値）の割合を計算する
valid = df[samples].notna().sum(axis=1) / len(samples)
# OR代入演算子で「いずれかの群で基準割合以上ならTrue」にする（どちらかの群で70%以上あれば残す）
keep |= (valid >= ratio)
```

- **`notna()`**: 各セルがNaN（欠損）でなければTrue、NaNならFalseを返します。
- **`sum(axis=1)`**: 行方向（各タンパク質）のTrueの数を合計します。axis=0なら列方向。
- **`/ len(samples)`**: サンプル数で割って「有効値の割合」を計算します。
- **`|=`**: OR演算の代入。いずれかの群で基準を満たせば残します。

---

### Perseus互換 downshift 欠損値補完

```python
# 補完値の中心: 有効値の平均からdownshift×SD分だけ低い値（検出限界以下を模擬する）
imp_mean = valid.mean() - downshift * valid.std()
# 補完値のばらつき: 元のSDのwidth倍（0.3なら30%）に絞って自然なばらつきを与える
imp_std  = width * valid.std()
# 正規分布N(imp_mean, imp_std)からn個の乱数を生成し、欠損位置（mask=True）に代入する
df.loc[mask, col] = np.random.normal(imp_mean, imp_std, n)
```

- **`valid.mean() - downshift * valid.std()`**: 有効値の平均から `downshift` × 標準偏差 分だけ低い値を補完の中心にします。例えば平均20、SD=2、downshift=2.4の場合、補完の中心は 20 - 2.4×2 = **15.2** です。
- **`width * valid.std()`**: 補完値のばらつきを元のSDの `width` 倍に設定します。0.3なら元のSDの30%の幅です。
- **`np.random.normal()`**: 正規分布から指定個数の乱数を生成します。

**なぜdownshiftするのか？**: 質量分析で「検出されなかった」ということは、そのタンパク質の量が「検出限界より低かった」ことを意味します。そのため、検出された値よりも低い値で補完するのが合理的です。

---

## 実行方法

**スクリプトで一括実行する場合:**

```bash
micromamba activate crc-proteomics
python scripts/step_04_preprocess.py
```

**Notebook でセルごとに実行する場合:**

```bash
jupyter notebook notebooks/step_06.ipynb
```

Notebook版では各処理の途中経過を確認しながら進められます。

## まとめ

Perseusの3つの主要機能をPython約80行で再現しました。パラメータ（downshift=2.4, width=0.3）を論文と合わせれば、同等の結果が得られます。本書のデータでは、sage 由来の 2,110 タンパク質 → 70% フィルタで 2,081 → 欠損値補完で 922 個の値を補完、という流れでクリーンな log2 マトリクスが得られました。

> 前回: [#5 OpenMS + AlphaPeptDeep](article-05-openms.md)
> 次回: [#7 全体像の可視化](article-07-visualization.md) — 相関行列・クラスタリング・PCA

#バイオインフォマティクス #プロテオミクス #Python #Perseus #labcode
