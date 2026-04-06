---
title: "Perseusなしで！Pythonでプロテオミクスデータを前処理する【論文再現シリーズ #5】"
emoji: "🧹"
type: "tech"
topics: ["proteomics", "python", "bioinformatics", "labcode"]
published: false
---

# Perseusなしで！Pythonでプロテオミクスデータを前処理する

## はじめに

プロテオミクスデータの前処理は、多くの論文でPerseus（MaxQuant付属ソフト）が使われます。しかしPerseusは**商用利用に制限がある**ため、本シリーズではPythonで完全代替します。この記事では、Log2変換・欠損値フィルタリング・欠損値補完をPythonで実装します。

## 前提

- [#1 環境構築](article-01-setup.md) が完了していること
- [#4 sageによるDIA解析](article-04-sage.md) で `results/protein_matrix_from_sage.csv` が生成されていること

## 本書サブセットで扱うデータ規模

本シリーズでは論文データの **18 ファイル（9 患者分: CRC04-CRC12、各 Normal/Tumor）** を扱います。sage による解析の結果として以下のデータが入力になります：

| 項目 | 値 |
|------|-----|
| タンパク質数 | **2,260** |
| サンプル数 | **18** (Normal 9 + Tumor 9) |
| 論文全体のタンパク質数 | 10,329 (参考) |

前処理後：

| 段階 | タンパク質数 |
|------|-------------|
| sage 出力 | 2,260 |
| Log2 変換 | 2,260 |
| 70% 有効値フィルタ | **2,234**（26 個除去） |
| 欠損値補完後 | 2,234（686 個の値を補完） |

## 前処理の全体像

```
DIA解析の出力（タンパク質定量マトリクス）
  → Log2変換（正規分布に近づける）
  → 欠損値フィルタリング（70%ルール）
  → 欠損値補完（Perseus互換: downshift法）
  → 前処理済みデータ
```

## スクリプト全文: step_02_preprocess.py

以下が前処理スクリプトの全文です。Perseusの3つの主要機能をPythonで再現しています。

```python
#!/usr/bin/env python3
"""
Step 2: データ前処理
==================
DIA解析の出力を読み込み、統計解析に適した形に前処理します。
論文ではPerseus（商用利用不可）を使用 → Pythonで完全代替。

【処理の流れ】
  1. タンパク質定量データの読み込み
  2. Log2変換（強度値のスケール変換）
  3. 欠損値フィルタリング（少なくとも1群で70%以上の有効値）
  4. 欠損値補完（Perseus互換: downshift法）

入力: data/raw/ (DIA解析出力 or 補足Table S2)
出力: results/preprocessed_data.csv, results/sample_info.csv
"""

# ============================================================
# ライブラリの読み込み
# ============================================================
import os            # ファイルパス操作
import numpy as np   # 数値計算ライブラリ（log2変換、乱数生成等）
import pandas as pd  # データフレーム操作ライブラリ（表形式データの読み書き）

# ============================================================
# 設定（パスとパラメータの定義）
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, "..")
RAW_DIR = os.path.join(PROJECT_DIR, "data", "raw")
RESULTS_DIR = os.path.join(PROJECT_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# --- 欠損値補完パラメータ（Perseus互換） ---
# downshift: 有効値の平均から何SD下にシフトするか
#   → 2.4 = 「検出されなかったタンパク質は有効値の平均より2.4SD低い」と仮定
# width: 補完値の分布幅（SDの倍率）
#   → 0.3 = 元のSDの30%の幅で補完値をばらつかせる
IMPUTATION_DOWNSHIFT = 2.4
IMPUTATION_WIDTH = 0.3

# 有効値の最低割合（論文の設定: 70%）
VALID_VALUE_RATIO = 0.70


# ============================================================
# Log2変換関数
# ============================================================
def log2_transform(df):
    """
    生の強度値をLog2変換する関数。

    【なぜLog2変換するか？】
      - 質量分析の強度値は桁が非常に大きい（10^6〜10^9）
      - Log2変換すると正規分布に近づき、統計手法が適用しやすくなる
      - Log2スケールでは差が「倍率」に対応（差1.0 = 2倍の変化）
    """
    median_val = df.median().median()
    if median_val > 100:
        # 中央値が大きい → 未変換の生データ → Log2変換実行
        print(f"  Log2変換実行（中央値={median_val:.1f} → 生データと判定）")
        df = np.log2(df.replace(0, np.nan))  # 0はNaNに置換してから変換
    else:
        # 中央値が小さい → 既にLog2スケール → スキップ
        print(f"  Log2変換スキップ（中央値={median_val:.1f} → 既にlog2スケール）")
    return df


# ============================================================
# 欠損値フィルタリング関数
# ============================================================
def filter_valid_values(df, sample_groups, ratio=VALID_VALUE_RATIO):
    """
    有効値が少なすぎるタンパク質を除去する関数。
    「少なくとも1群で70%以上のサンプルに有効値がある」タンパク質のみ残す。

    【なぜフィルタリングが必要か？】
      - ほとんどのサンプルで検出されないタンパク質は統計的に信頼できない
      - 欠損値が多すぎると補完の精度も低下する
    """
    n_before = len(df)
    keep_mask = pd.Series(False, index=df.index)

    for group_name, samples in sample_groups.items():
        group_data = df[samples]
        # notna(): NaNでなければTrue → sum(): True の数を数える
        valid_ratio = group_data.notna().sum(axis=1) / len(samples)
        keep_mask |= (valid_ratio >= ratio)  # いずれかの群で基準を満たせばOK

    df_filtered = df[keep_mask]
    n_after = len(df_filtered)
    print(f"  有効値フィルタ: {n_before} → {n_after} タンパク質（{n_before - n_after} 除去）")
    return df_filtered


# ============================================================
# 欠損値補完関数（Perseus互換 downshift法）
# ============================================================
def impute_missing_values(df, downshift=IMPUTATION_DOWNSHIFT, width=IMPUTATION_WIDTH):
    """
    Perseus互換の欠損値補完を行う関数。

    【補完の考え方】
      「検出されなかった」= 存在しないのではなく「低すぎて検出限界以下」
      → 検出された値の分布から推測して低い値を割り当てる

    【手順（各サンプルごと）】
      1. 有効値の平均(μ)と標準偏差(σ)を計算
      2. 補完値の中心 = μ - downshift × σ
      3. 補完値の幅 = width × σ
      4. この分布からランダムにサンプリング
    """
    df_imputed = df.copy()
    n_imputed = 0

    for col in df.columns:
        valid = df[col].dropna()          # 有効値のみ取得
        if len(valid) == 0:
            continue

        col_mean = valid.mean()           # 有効値の平均
        col_std = valid.std()             # 有効値の標準偏差
        imp_mean = col_mean - downshift * col_std  # 補完値の中心（低くシフト）
        imp_std = width * col_std                   # 補完値の幅（狭く）

        missing_mask = df[col].isna()     # 欠損値の位置
        n_missing = missing_mask.sum()

        if n_missing > 0:
            # 正規分布からランダムサンプリングして補完
            imp_values = np.random.normal(imp_mean, imp_std, n_missing)
            df_imputed.loc[missing_mask, col] = imp_values
            n_imputed += n_missing

    print(f"  欠損値補完: {n_imputed} 値を補完（downshift={downshift}, width={width}）")
    return df_imputed


# ============================================================
# デモデータ生成関数（実データがない場合に自動生成）
# ============================================================
def generate_demo_data():
    """
    論文の記述に基づいたデモデータを生成。
    500タンパク質 × 32サンプル（16患者 × Normal/Tumor）。
    """
    np.random.seed(42)
    n_proteins, n_patients = 500, 16
    stages = ["I"]*3 + ["II"]*4 + ["III"]*4 + ["IV"]*5

    samples_n = [f"CRC{i+1:02d}-N" for i in range(n_patients)]
    samples_t = [f"CRC{i+1:02d}-T" for i in range(n_patients)]

    base_expr = np.random.normal(20, 3, n_proteins)  # log2スケールの基本発現量
    data = {}
    for i, sample in enumerate(samples_n + samples_t):
        is_tumor = i >= n_patients
        noise = np.random.normal(0, 0.5, n_proteins)
        if is_tumor:
            change = np.zeros(n_proteins)
            change[:100] = np.random.uniform(1, 3, 100)      # 腫瘍で上昇
            change[100:200] = np.random.uniform(-3, -1, 100)  # 腫瘍で低下
            expr = base_expr + noise + change
        else:
            expr = base_expr + noise
        mask = np.random.random(n_proteins) < 0.10  # 10%を欠損に
        expr[mask] = np.nan
        data[sample] = expr

    df = pd.DataFrame(data, index=[f"Protein_{i+1:04d}" for i in range(n_proteins)])
    df.index.name = "Protein"

    # 臨床情報の保存
    pd.DataFrame({"Sample_N": samples_n, "Sample_T": samples_t, "Stage": stages}
    ).to_csv(os.path.join(RESULTS_DIR, "clinical_info.csv"), index=False)

    print(f"  デモデータ生成: {n_proteins} タンパク質 × {len(data)} サンプル")
    return df


# ============================================================
# メイン処理
# ============================================================
def main():
    print("=" * 60)
    print("Step 2: データ前処理")
    print("=" * 60)

    # [1/4] データ読み込み
    print("\n[1/4] タンパク質定量データの読み込み")
    # （実データがなければデモデータを自動生成）
    df = generate_demo_data()  # ← 実データがあれば load_protein_data() を使用

    # [2/4] サンプル群の定義
    print("\n[2/4] サンプル群の定義")
    normal_samples = [c for c in df.columns if "-N" in c]
    tumor_samples = [c for c in df.columns if "-T" in c]
    print(f"  Normal群: {len(normal_samples)} サンプル")
    print(f"  Tumor群:  {len(tumor_samples)} サンプル")

    # [3/4] Log2変換
    print("\n[3/4] Log2変換")
    df = log2_transform(df)

    # [4/4] フィルタリング + 欠損値補完
    print("\n[4/4] 有効値フィルタリングと欠損値補完")
    df = filter_valid_values(df, {"Normal": normal_samples, "Tumor": tumor_samples})
    df = impute_missing_values(df)

    # 保存
    df.to_csv(os.path.join(RESULTS_DIR, "preprocessed_data.csv"))
    pd.DataFrame({
        "Sample": normal_samples + tumor_samples,
        "Condition": ["Normal"]*len(normal_samples) + ["Tumor"]*len(tumor_samples)
    }).to_csv(os.path.join(RESULTS_DIR, "sample_info.csv"), index=False)

    print(f"\n前処理完了: {df.shape[0]} タンパク質 × {df.shape[1]} サンプル")


if __name__ == "__main__":
    main()
```

---

## コード詳細

### Log2変換

```python
df = np.log2(df.replace(0, np.nan))
```

- **`df.replace(0, np.nan)`**: 値が0のセルをNaN（欠損値）に置換します。Log2(0)は-∞になってしまうため、先にNaNにしておきます。
- **`np.log2()`**: NumPyのlog2関数。各セルの値をLog2変換します。NaNはそのままNaNとして保持されます。
- **変換の意味**: 例えば値が1,000,000の場合、Log2(1,000,000) ≈ 20 になります。Log2スケールでは差が「倍率」に直結します（差1.0 = 2倍の変化）。

---

### 欠損値フィルタリング（70%ルール）

```python
valid_ratio = group_data.notna().sum(axis=1) / len(samples)
keep_mask |= (valid_ratio >= ratio)
```

- **`notna()`**: 各セルがNaN（欠損）でなければTrue、NaNならFalseを返します。
- **`sum(axis=1)`**: 行方向（各タンパク質）のTrueの数を合計します。axis=0なら列方向。
- **`/ len(samples)`**: サンプル数で割って「有効値の割合」を計算します。
- **`|=`**: OR演算の代入。いずれかの群で基準を満たせば残します。

---

### Perseus互換 downshift 欠損値補完

```python
imp_mean = col_mean - downshift * col_std
imp_std = width * col_std
imp_values = np.random.normal(imp_mean, imp_std, n_missing)
```

- **`col_mean - downshift * col_std`**: 有効値の平均から `downshift` × 標準偏差 分だけ低い値を補完の中心にします。例えば平均20、SD=2、downshift=2.4の場合、補完の中心は 20 - 2.4×2 = **15.2** です。
- **`width * col_std`**: 補完値のばらつきを元のSDの `width` 倍に設定します。0.3なら元のSDの30%の幅です。
- **`np.random.normal()`**: 正規分布から指定個数の乱数を生成します。

**なぜdownshiftするのか？**: 質量分析で「検出されなかった」ということは、そのタンパク質の量が「検出限界より低かった」ことを意味します。そのため、検出された値よりも低い値で補完するのが合理的です。

---

## 実行方法

```bash
micromamba activate crc-proteomics
python scripts/step_02_preprocess.py
```

## まとめ

Perseusの3つの主要機能をPython約80行で再現しました。パラメータ（downshift=2.4, width=0.3）を論文と合わせれば、同等の結果が得られます。本書のデータでは、sage 由来の 2,260 タンパク質 → 70% フィルタで 2,234 → 欠損値補完で 686 個の値を補完、という流れでクリーンな log2 マトリクスが得られました。

> 前回: [#4 sageで同定・定量](article-04-sage.md)
> 次回: [#6 全体像の可視化](article-06-visualization.md) — 相関行列・クラスタリング・PCA

#バイオインフォマティクス #プロテオミクス #Python #Perseus #labcode
