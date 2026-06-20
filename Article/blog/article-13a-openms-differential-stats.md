---
title: "深層学習DIA解析での差分発現統計解析【論文再現シリーズ #13a】"
emoji: "📊"
type: "tech"
topics: ["proteomics", "python", "deeplearning", "statistics", "labcode"]
published: false
---

# 深層学習DIA解析での差分発現統計解析

## はじめに

この記事では、OpenMS + AlphaPeptDeepによる深層学習DIA解析結果から、腫瘍 vs 正常組織で有意に異なるタンパク質をWelch's t検定で統計的に特定します。Volcanoプロットで結果を可視化し、深層学習による高精度検出の効果を確認します。

> **📝 INFO**
>
**この記事で行う処理**
OpenMS + AlphaPeptDeepで前処理済みの19,981タンパク質について、各タンパク質のNormal群とTumor群をWelchのt検定で比較し、統計的に有意（p<0.05）かつ発現変化が大きい（fold change > 2倍）タンパク質を特定します。結果をVolcanoプロットで可視化し、腫瘍で増加・減少したタンパク質を論文レベルの精度で把握します。

## 前提

- [#12b 深層学習DIA PCA解析](article-12b-openms-visualization-pca.md) が完了していること
- **対応Notebook**: [`notebooks/step_13.ipynb`](../notebooks/step_13.ipynb) — この記事のコードをセルごとに実行できます

## コード全文（対応Notebook: step_13.ipynb）

### ライブラリと設定

```python
import numpy as np           # 数値計算ライブラリ（配列操作・数学関数に使用）
import pandas as pd          # データフレーム操作ライブラリ（表形式データの読み込み・加工に使用）
import matplotlib.pyplot as plt  # グラフ描画ライブラリ（Volcanoプロット等の作成に使用）
from scipy import stats      # 統計検定ライブラリ（Welch's t-testの実行に使用）

# --- パス ---
RESULTS  = "../results"              # 解析結果の保存先ディレクトリへのパス
FIG_DIR  = f"{RESULTS}/figures"      # 図の保存先ディレクトリへのパス
TABLE_DIR = f"{RESULTS}/tables"      # テーブル（CSV等）の保存先ディレクトリへのパス

# --- 閾値 ---
# p値の閾値: 0.05未満なら「統計的に有意」と判定する（慣例的に広く使われる基準）
P_THRESHOLD    = 0.05
# log2 fold changeの閾値: 1.0以上なら2倍以上の発現変化があることを意味する
LOG2FC_THRESHOLD = 1.0

# --- カラー ---
UP, DOWN, NS  = "#E74C3C", "#3498DB", "#CCCCCC"  # 上昇(赤)・低下(青)・非有意(灰)の色コード
```

### データ読み込み

```python
# データ読み込み & Normal/Tumor サンプル分離
# OpenMS + AlphaPeptDeepで前処理済みタンパク質発現データを読み込む（行=タンパク質、列=サンプル、値=log2発現量）
df = pd.read_csv(f"{RESULTS}/preprocessed_data_openms.csv", index_col=0)
# サンプル情報（各サンプルがNormalかTumorか）を読み込む
sample_info = pd.read_csv(f"{RESULTS}/sample_info_openms.csv")
# サンプル名をインデックスにして条件列だけを取り出す（後でサンプルから条件を引くための辞書的Series）
conditions = sample_info.set_index("Sample")["Condition"]

# Normal群のサンプル名リストを取得（Conditionが'Normal'の行をフィルタ）
normal_samples = sample_info.query("Condition == 'Normal'")["Sample"].tolist()
# Tumor群のサンプル名リストを取得（Conditionが'Tumor'の行をフィルタ）
tumor_samples  = sample_info.query("Condition == 'Tumor'")["Sample"].tolist()
# タンパク質数と各群のサンプル数を表示して、データの規模を確認する
print(f"Proteins: {len(df)}, Normal: {len(normal_samples)}, Tumor: {len(tumor_samples)}")
```

### Welch's t-test

```python
def welch_ttest(df, normal_samples, tumor_samples):
    """全タンパク質に Welch t-test を実行し DataFrame を返す。"""
    rows = []  # 各タンパク質の検定結果を一時的に格納するリスト
    # 全タンパク質を1つずつループして検定を行う
    for protein in df.index:
        # そのタンパク質のNormal群の発現量を取得し、欠損値（NaN）を除外する
        nv = df.loc[protein, normal_samples].dropna()
        # そのタンパク質のTumor群の発現量を取得し、欠損値（NaN）を除外する
        tv = df.loc[protein, tumor_samples].dropna()
        # 各群に2個以上のデータがないとt検定が計算できないためスキップする
        if len(nv) < 2 or len(tv) < 2:
            continue
        # Welchのt検定を実行（equal_var=Falseで等分散を仮定しない）
        # t_stat: t統計量（値が大きいほど2群の差が大きい）
        # p_val: p値（この差が偶然生じる確率。小さいほど有意）
        t_stat, p_val = stats.ttest_ind(tv, nv, equal_var=False)
        # log2 fold change = Tumor群の平均 - Normal群の平均（データが既にlog2スケールのため引き算で比を表す）
        # 正の値なら腫瘍で発現増加、負の値なら腫瘍で発現減少を意味する
        log2fc = tv.mean() - nv.mean()
        # 検定結果を辞書形式でリストに追加する
        rows.append({
            "Protein": protein,               # タンパク質名
            "Mean_Normal": nv.mean(),          # Normal群の平均発現量（log2スケール）
            "Mean_Tumor": tv.mean(),           # Tumor群の平均発現量（log2スケール）
            "Log2FC": log2fc,                  # 発現変化量（正=腫瘍で増加、負=腫瘍で減少）
            "T_statistic": t_stat,             # t統計量（2群の差の大きさを標準誤差で割った値）
            "P_value": p_val,                  # p値（帰無仮説「2群に差がない」が正しい確率）
            # -log10(p値): Volcanoプロットのy軸用。p値が小さいほど値が大きくなり上に表示される
            # max(p_val, 1e-300)でp値が0の場合のlog計算エラーを防ぐ
            "Neg_log10_P": -np.log10(max(p_val, 1e-300)),
        })
    # リストをDataFrameに変換する（1行=1タンパク質の検定結果）
    result = pd.DataFrame(rows)
    # 有意性の分類: まず全タンパク質を "NS"（Not Significant = 非有意）に設定する
    result["Significant"] = "NS"
    # p値が閾値未満 かつ log2FCが正の閾値を超える → "Up"（腫瘍で有意に増加）
    result.loc[(result["P_value"] < P_THRESHOLD) & (result["Log2FC"] >  LOG2FC_THRESHOLD), "Significant"] = "Up"
    # p値が閾値未満 かつ log2FCが負の閾値を下回る → "Down"（腫瘍で有意に減少）
    result.loc[(result["P_value"] < P_THRESHOLD) & (result["Log2FC"] < -LOG2FC_THRESHOLD), "Significant"] = "Down"
    return result  # 全タンパク質の検定結果DataFrameを返す
```

### t検定の実行

```python
# 上で定義したwelch_ttest関数を実行し、全タンパク質の検定結果を取得する
result_df = welch_ttest(df, normal_samples, tumor_samples)

# 腫瘍で有意に増加したタンパク質の数をカウントする
n_up   = (result_df["Significant"] == "Up").sum()
# 腫瘍で有意に減少したタンパク質の数をカウントする
n_down = (result_df["Significant"] == "Down").sum()
# 検定を実行できたタンパク質の総数を表示する
print(f"検定タンパク質数: {len(result_df)}")
# 有意差ありの内訳（Up: 腫瘍で増加、Down: 腫瘍で減少）を表示する
print(f"有意差あり: {n_up + n_down} (Up {n_up}, Down {n_down})")

# 検定結果をCSVファイルに保存する（後の解析やCOSMIC照合で使用する）
result_df.to_csv(f"{TABLE_DIR}/differential_proteins_openms.csv", index=False)
```

### Volcanoプロット

```python
# 図と軸オブジェクトを作成する（figsize=(横8, 縦6)インチ）
fig, ax = plt.subplots(figsize=(8, 6))

# NS（非有意）→ Up（増加）→ Down（減少）の順にプロットする
# NSを先に描画することで、有意な点が上に表示されて見やすくなる
for sig, color, alpha in [("NS", NS, 0.3), ("Up", UP, 0.6), ("Down", DOWN, 0.6)]:
    # 各カテゴリに該当するタンパク質のブールマスクを作成する
    m = result_df["Significant"] == sig
    # 散布図を描画する（x=発現変化量, y=統計的有意性）
    # c: 点の色, s: 点のサイズ(10pt), alpha: 透明度（NSは薄く表示）
    ax.scatter(result_df.loc[m, "Log2FC"], result_df.loc[m, "Neg_log10_P"],
              c=color, s=10, alpha=alpha, label=f"{sig} ({m.sum()})")

# p値の閾値ライン（水平破線）: この線より上が統計的に有意
ax.axhline(-np.log10(P_THRESHOLD), color="gray", ls="--", lw=0.5)
# log2FCの正の閾値ライン（垂直破線）: この線より右が2倍以上の増加
ax.axvline( LOG2FC_THRESHOLD, color="gray", ls="--", lw=0.5)
# log2FCの負の閾値ライン（垂直破線）: この線より左が2倍以上の減少
ax.axvline(-LOG2FC_THRESHOLD, color="gray", ls="--", lw=0.5)

# x軸ラベル: 発現変化の方向と大きさを示す
ax.set_xlabel("Log2 Fold Change (Tumor / Normal)")
# y軸ラベル: 統計的有意性を示す（値が大きいほど有意）
ax.set_ylabel("-Log10(P-value)")
# グラフタイトル: 深層学習DIA解析での差分タンパク質量解析
ax.set_title("Deep Learning DIA: Differential Protein Abundance (Tumor vs Non-tumor)")
# 凡例を右上に表示する（frameon=Falseで枠線を非表示にしてすっきりさせる）
ax.legend(frameon=False, loc="upper right")
# 上辺と右辺の枠線を非表示にする（学術論文スタイルのすっきりした見た目にする）
ax.spines[["top", "right"]].set_visible(False)

# 図をPNGファイルとして保存する（dpi=150で高解像度、bbox_inches="tight"で余白を最小化）
fig.savefig(f"{FIG_DIR}/fig_volcano_openms.png", dpi=150, bbox_inches="tight")
# 画面に表示する
plt.show()
```

## コード詳細

### Welch's t-testの特徴

| パラメータ | 設定値 | 意味 |
|-----------|--------|------|
| `equal_var=False` | False | 2群の分散が等しくないと仮定（Welchのt検定） |
| `P_THRESHOLD` | 0.05 | 統計的有意性の閾値（5%有意水準） |
| `LOG2FC_THRESHOLD` | 1.0 | 発現変化の閾値（2倍変化 = log2 1.0） |

### Log2 Fold Changeの解釈

| 値 | 意味 | 発現変化倍率 |
|----|------|-------------|
| +1.0 | 腫瘍で2倍増加 | 2.0倍 |
| +2.0 | 腫瘍で4倍増加 | 4.0倍 |
| -1.0 | 腫瘍で半減 | 0.5倍 |
| -2.0 | 腫瘍で1/4に減少 | 0.25倍 |

### Volcanoプロットの読み方

- **x軸**: Log2 Fold Change（発現変化の方向と大きさ）
- **y軸**: -Log10(P-value)（統計的有意性、値が大きいほど有意）
- **右上（赤）**: 腫瘍で有意に増加したタンパク質
- **左上（青）**: 腫瘍で有意に減少したタンパク質
- **下部（灰）**: 統計的に有意でないタンパク質

## 実行方法

**スクリプトで一括実行する場合:**
```bash
python scripts/step_13a_openms_differential_stats.py
```

**Notebook でセルごとに実行する場合:**
```bash
jupyter notebook notebooks/step_13.ipynb
```

## まとめ

OpenMS + AlphaPeptDeepによる深層学習DIA解析で、腫瘍 vs 正常組織の差分発現タンパク質を統計的に特定しました。Welch's t検定により有意差ありタンパク質を抽出し、Volcanoプロットで可視化しました。

高度な可視化解析（ヒートマップ、PCA、Top N解析）については次章で詳しく解説します。

> 前回: [#12b 深層学習DIA PCA解析](article-12b-openms-visualization-pca.md)
> 次回: [#13b 深層学習DIA高度可視化解析](article-13b-openms-differential-visualization.md) — ヒートマップとPCAによる詳細解析

#バイオインフォマティクス #プロテオミクス #深層学習 #統計解析 #labcode