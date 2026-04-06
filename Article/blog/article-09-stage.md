---
title: "大腸がんのステージ進行に伴うタンパク質変動をクラスター分析で可視化する【論文再現シリーズ #9】"
emoji: "📈"
type: "tech"
topics: ["proteomics", "python", "clustering", "labcode"]
published: false
---

# 大腸がんのステージ進行に伴うタンパク質変動をクラスター分析で可視化する

## はじめに

大腸がんはStage I〜IVに分類されます。この記事では、各ステージでのタンパク質発現パターンを解析し、**疾患進行に伴って一貫して変動するタンパク質群**を特定します。Toyota et al. 2025のFigure 3の再現です。

## 前提

- [#8 COSMIC照合](article-08-cosmic.md) が完了していること
- ステージ情報（Table S1）が利用可能

## スクリプト全文: step_06_stage_analysis.py

> コアとなるANOVA関数とクラスター分析関数を抜粋します。完全版は `scripts/step_06_stage_analysis.py` を参照。

### コアロジック: One-way ANOVA

```python
def run_anova(df, sample_info):
    """
    各タンパク質について One-way ANOVA を実行する関数。

    【One-way ANOVAとは？】
      3群以上の平均値に差があるかを検定する手法。
      t検定は2群専用だが、ANOVAは3群以上に対応。
      本解析では Normal, Stage I, II, III, IV の5群を比較する。

    【FDR補正（Benjamini-Hochberg法）とは？】
      数千回の検定を同時に行うと、偶然だけで有意になるものが出てくる。
      例: 10,000回検定 × p < 0.05 → 偶然500個が有意に
      FDR補正はこの「偽陽性の割合」を制御する手法。
    """
    # ステージのリスト（Normalを含む5群）
    stages_list = ["Normal", "I", "II", "III", "IV"]

    results = []
    for protein in df.index:
        # 各ステージのサンプルの発現値を取得
        groups = []
        for stage in stages_list:
            samples = sample_info[sample_info["Stage"] == stage]["Sample"].tolist()
            vals = df.loc[protein, [s for s in samples if s in df.columns]].dropna()
            if len(vals) >= 2:  # 各群最低2サンプル必要
                groups.append(vals.values)

        if len(groups) < 2:
            continue

        # --- One-way ANOVA の実行 ---
        # stats.f_oneway(*groups): 可変長引数で複数の群を渡す
        # 戻り値: F統計量（群間差/群内差の比）とp値
        f_stat, p_val = stats.f_oneway(*groups)

        results.append({
            "Protein": protein,
            "F_statistic": f_stat,
            "P_value": p_val,
        })

    result_df = pd.DataFrame(results)

    # --- FDR補正（Benjamini-Hochberg法） ---
    # multipletests() で多重検定補正を実行
    # method="fdr_bh": Benjamini-Hochberg法（最もポピュラーなFDR補正法）
    from statsmodels.stats.multitest import multipletests
    reject, fdr, _, _ = multipletests(result_df["P_value"], method="fdr_bh")
    result_df["FDR"] = fdr
    result_df["Significant"] = fdr < 0.01  # FDR < 0.01 を有意とする

    return result_df
```

### コアロジック: クラスター分析

```python
def cluster_and_plot(median_df, anova_df, n_clusters=30):
    """
    有意なタンパク質を階層的クラスタリングし、30クラスターに分割する関数。

    【階層的クラスタリングの手順】
      1. 各タンパク質のステージ別中央値をZ-score正規化
         → 全タンパク質を「同じスケール」に揃える
      2. Ward法で階層的にクラスタを結合
         → 最も分散が小さくなるようにペアを結合していく
      3. 30クラスターに分割（論文のパラメータ）
         → 似た発現パターンのタンパク質をグループ化

    【Z-score正規化とは？】
      各タンパク質の発現値を「平均0、標準偏差1」に変換する処理。
      これにより、発現量の絶対値ではなく「変動パターン」で比較できるようになる。
      計算式: z = (x - 平均) / 標準偏差
    """
    # 有意なタンパク質のみ抽出
    sig_proteins = anova_df[anova_df["Significant"]]["Protein"].tolist()
    sig_data = median_df.loc[median_df.index.isin(sig_proteins)]

    # --- Z-score正規化（行方向 = 各タンパク質ごと） ---
    # apply(lambda x: ..., axis=1) は各行に対して関数を適用
    z_data = sig_data.apply(
        lambda x: (x - x.mean()) / x.std(),  # Z-score の計算式
        axis=1   # axis=1 = 行方向（各タンパク質）
    ).dropna()

    # --- 階層的クラスタリング ---
    from scipy.cluster.hierarchy import linkage, fcluster

    # linkage(): 階層的クラスタリングの連結行列を計算
    # method="ward": Ward法（分散最小化基準）
    Z = linkage(z_data.values, method="ward")

    # fcluster(): 連結行列から指定数のクラスターに分割
    # t=30: 30クラスターに分割（論文のパラメータ）
    # criterion="maxclust": 最大クラスター数で分割
    clusters = fcluster(Z, t=n_clusters, criterion="maxclust")

    return z_data, clusters
```

### コアロジック: プロファイルプロット

```python
def plot_profile_plots(z_data, clusters):
    """
    Figure 3 のプロファイルプロットを作成する関数。

    【プロファイルプロットとは？】
      X軸にステージ（Normal, I, II, III, IV）、Y軸にZ-score をとり、
      各タンパク質の発現変動パターンを線グラフで描画したもの。
      灰色の細い線 = 個々のタンパク質
      黒い太い線 = クラスター全体の平均パターン

    【上昇クラスター】: Normalより腫瘍で高く、ステージとともに上昇
      → バイオマーカー候補（腫瘍促進因子）
    【下降クラスター】: Normalより腫瘍で低く、ステージとともに下降
      → 腫瘍抑制因子の候補
    """
    # ステージのカラー（論文準拠）
    # 青 → 緑 → 黄 → オレンジ → 赤 で進行を直感的に表現
    stage_colors = {
        "Normal": "#4EAED1",  # 青
        "I": "#66BB6A",       # 緑
        "II": "#FFD54F",      # 黄
        "III": "#FFA726",     # オレンジ
        "IV": "#E8524A",      # 赤
    }

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    stages = z_data.columns.tolist()

    for idx, cl in enumerate(selected_clusters):
        ax = axes.flatten()[idx]
        mask = clusters == cl
        cluster_data = z_data.iloc[mask]

        # 個別タンパク質のプロファイル（薄い灰色の線）
        for _, row in cluster_data.iterrows():
            ax.plot(stages, row.values, color="gray", alpha=0.1, linewidth=0.5)

        # 平均プロファイル（太い黒線）
        mean_profile = cluster_data.mean(axis=0)
        ax.plot(stages, mean_profile.values, color="#333333", linewidth=2.5, marker="o")

        ax.set_title(f"Cluster {cl} ({mask.sum()} proteins)")
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.5)

    plt.tight_layout()
    plt.savefig("results/figures/fig3_profile_plots.png", dpi=150, bbox_inches="tight")
```

---

## コード詳細

### t検定 vs ANOVA の使い分け

| 手法 | 群数 | 使いどころ |
|------|------|----------|
| **t検定** | 2群 | Normal vs Tumor の比較（Step 4） |
| **ANOVA** | 3群以上 | Normal vs Stage I vs II vs III vs IV の比較（Step 6） |

### FDR補正の重要性

```python
# FDR補正なし: 10,000タンパク質 × p < 0.05 → 偶然500個が有意に！
# FDR補正あり: 偽陽性の割合を1%以下に制御
reject, fdr, _, _ = multipletests(p_values, method="fdr_bh")
```

- `fdr_bh` = Benjamini-Hochberg法（FDRを制御する最もポピュラーな方法）
- FDR < 0.01: 「有意」と判定されたタンパク質のうち、偽陽性は1%以下

### 論文の主な発見

| クラスター | タンパク質数 | パターン | 意味 |
|-----------|-----------|---------|------|
| Cluster 3 | 34 | Normal→IV で一貫上昇 | 腫瘍促進因子の候補 |
| Cluster 14 | 1,324 | Normal→IV で上昇 | 広範な腫瘍関連変動 |
| Cluster 20 | 16 | Normal→IV で一貫下降 | 腫瘍抑制因子の候補 |
| Cluster 25 | 1,062 | Normal→IV で下降 | 正常組織機能の喪失 |

### 本書データでの再現結果

> **⚠ ステージ構成の制約**
> 本書のサブセット（CRC04-CRC12）には **Stage I と Stage IV が含まれず、Normal (9人) + Stage II (4人) + Stage III (5人) の 3 段階** のみです。論文の 5 段階（Normal → I → II → III → IV）とは異なる構成で解析します。

| 指標 | 本書 (sage, 2234タンパク質) | 論文 (DIA-NN, 10329タンパク質) |
|------|----------------------------|-------------------------------|
| ANOVA 有意 (FDR<0.01) | **151** | 記載なし（論文はクラスター先行） |
| クラスター数 | **30** | 30 |

### ステージ別ヒートマップとクラスタープロファイル

![ステージ別ヒートマップ](images/fig3_stage_heatmap.png)

ANOVA で有意（FDR<0.01）と判定された151タンパク質について、ステージ別の発現パターンをヒートマップとプロファイルプロットで可視化しています。30クラスターに分割した結果から代表的なクラスターを抽出し、Normal → Stage II → Stage III にかけて段階的に変動するパターンが観察されます。灰色の細い線は個々のタンパク質、黒い太い線はクラスター平均を示しています。

特筆すべきは、サンプル数を増やし9患者18ファイルで再解析したことで、**ANOVAで151個のタンパク質が有意**（FDR<0.01）と判定され、30クラスターへの分割が可能になった点です。Stage I/IVが含まれないサブセットであるにもかかわらず、疾患進行に伴うプロテオーム変動パターンという生物学的シグナルは商用クリアな sage パイプラインでも検出できることを示しています。

## 実行方法

```bash
python scripts/step_08_stage_analysis.py
```

## まとめ

ステージ別解析により、疾患進行に伴って一貫して変動するタンパク質クラスターを同定しました。本書データでは ANOVA 有意な 151 タンパク質を 30 クラスターに分割し、Normal → Stage II → Stage III にかけて段階的に変動する候補マーカー群を特定しました。

> 前回: [#8 COSMIC照合](article-08-cosmic.md)
> 次回: [#10 まとめと次のステップ](article-10-conclusion.md) — 論文化への道筋

#バイオインフォマティクス #プロテオミクス #クラスター分析 #大腸がん #labcode
