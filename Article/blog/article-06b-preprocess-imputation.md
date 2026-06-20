---
title: "Perseus互換欠損値補完でプロテオミクスデータを仕上げる【論文再現シリーズ #6b】"
emoji: "🔮"
type: "tech"
topics: ["proteomics", "python", "imputation", "labcode"]
published: false
---

# Perseus互換欠損値補完でプロテオミクスデータを仕上げる

## はじめに

前回（[#6a 基本前処理](article-06a-preprocess-basics.md)）でLog2変換と有効値フィルタリングを完了しました。この記事では、残った欠損値をPerseus互換の**downshift法**で補完し、統計解析に適したクリーンなデータを完成させます。

> **📝 INFO**
>
**この記事で行う処理**
「検出されなかった」タンパク質の値をdownshift法で補完します。この方法では「検出限界以下の低い値が存在した」と仮定し、検出されたタンパク質の分布から推測した低い値をランダムサンプリングで割り当てます。Perseus互換パラメータ（downshift=2.4, width=0.3）を使用して論文と同等の前処理を実現します。

## 前提

- [#6a 基本前処理](article-06a-preprocess-basics.md) が完了していること
- 有効値フィルタリング済みのデータが準備されていること
- **対応Notebook**: [`notebooks/step_06.ipynb`](../notebooks/step_06.ipynb) — この記事のコードをセルごとに実行できます

## 欠損値補完の理論

### なぜ補完が必要か？

プロテオミクスデータの欠損値は「存在しない」のではなく「検出限界以下で検出できなかった」ことを意味します。統計解析では欠損値があると多くの手法が適用できないため、合理的な推定値で補完します。

```
検出されたタンパク質の分布:
     │ ●●●●●●●●●
強度 │ ●●●●●●●●●●● ← この辺りが検出の中心
     │ ●●●●●●●●
     │ ●●●●●
     │ ●●○○○○○○ ← 検出限界以下（○が補完対象）
     └─────────────
```

## コード全文（対応Notebook: step_06.ipynb）

### Perseus互換パラメータ設定

```python
import numpy as np   # 数値計算ライブラリ（乱数生成・統計計算に使用）
import pandas as pd  # データフレーム操作ライブラリ（CSV読み書き・データ操作に使用）
import os           # ファイルパス操作ライブラリ

# --- Perseus互換パラメータ ---
# downshift: 有効値の平均から何SD下にシフトするか
#   → 2.4 = 「検出されなかったタンパク質は有効値の平均より2.4SD低い」と仮定
# width: 補完値の分布幅（SDの倍率）
#   → 0.3 = 元のSDの30%の幅で補完値をばらつかせる

# Perseus互換のdownshiftパラメータ: 有効値の平均から2.4SD下にシフトして補完値の中心を決める
DOWNSHIFT = 2.4
# Perseus互換のwidthパラメータ: 補完値のばらつきを元のSDの30%に設定する
WIDTH = 0.3

# --- パス設定 ---
RESULTS_DIR = "../results"  # 結果ファイルを格納するディレクトリへの相対パス
```

### 欠損値補完（Perseus互換 downshift法）

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

### データ保存

```python
# 前処理済みのタンパク質定量マトリクスをCSVファイルとして保存する
df.to_csv(os.path.join(RESULTS_DIR, "preprocessed_data.csv"))

# サンプル情報（各サンプルがNormalかTumorか）をDataFrameにまとめる
# 後続のステップ（差次解析・可視化）で群の情報が必要になるため、別ファイルとして保存する

# サンプル名の末尾 -N / -T で Normal（正常組織） / Tumor（腫瘍組織） を判別
normal_samples = [c for c in df.columns if "-N" in c]
tumor_samples  = [c for c in df.columns if "-T" in c]

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

## コード詳細

### Perseus互換 downshift 欠損値補完

```python
# 補完値の中心: 有効値の平均からdownshift×SD分だけ低い値（検出限界以下を模擬する）
imp_mean = valid.mean() - downshift * valid.std()
# 補完値のばらつき: 元のSDのwidth倍（0.3なら30%）に絞って自然なばらつきを与える
imp_std  = width * valid.std()
# 正規分布N(imp_mean, imp_std)からn個の乱数を生成し、欠損位置（mask=True）に代入する
df.loc[mask, col] = np.random.normal(imp_mean, imp_std, n)
```

| パラメータ | 設定値 | 意味 |
|-----------|--------|------|
| **downshift** | 2.4 | 補完値の中心を有効値平均から2.4SD下に設定 |
| **width** | 0.3 | 補完値のばらつきを元のSDの30%に設定 |

### 補完の具体例

有効値の分布が平均=20、標準偏差=2の場合：

```python
# 補完値の分布パラメータ
imp_mean = 20 - 2.4 * 2 = 15.2  # 中心値
imp_std = 0.3 * 2 = 0.6          # ばらつき幅

# 補完値は N(15.2, 0.6) からサンプリング
# → 約14.0〜16.4の範囲で補完される
```

### なぜdownshiftするのか？

質量分析で「検出されなかった」ということは、そのタンパク質の量が「検出限界より低かった」ことを意味します。そのため、検出された値よりも低い値で補完するのが合理的です。

## 本書データでの実測値

| 項目 | 値 |
|------|-----|
| 補完前の欠損値数 | 922個 |
| 補完後の欠損値数 | **0個**（完全補完） |
| 最終タンパク質数 | **2,081個** |
| 最終サンプル数 | **32個** |

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

## 出力ファイル

| ファイル | 内容 |
|---------|------|
| `preprocessed_data.csv` | 前処理済みタンパク質定量マトリクス（2,081×32） |
| `sample_info.csv` | サンプル情報（Sample, Conditionカラム） |

これらのファイルが以降のすべての統計解析・可視化の入力データになります。

## まとめ

Perseus互換のdownshift法により、欠損値のないクリーンなlog2スケールマトリクスが完成しました。パラメータ（downshift=2.4, width=0.3）を論文と合わせることで、同等の結果が得られます。本書のデータでは922個の欠損値を補完し、統計解析に適した形に整えました。

> 前回: [#6a 基本前処理](article-06a-preprocess-basics.md)
> 次回: [#7a Sage基本可視化](article-07a-visualization-basics.md) — 相関行列・クラスタリング

#バイオインフォマティクス #プロテオミクス #Python #Perseus #欠損値補完 #labcode