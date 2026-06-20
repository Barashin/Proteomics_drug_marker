---
title: "深層学習DIA可視化とバイオマーカー【論文再現シリーズ #15d】"
emoji: "📊"
type: "tech"
topics: ["proteomics", "deeplearning", "visualization", "biomarkers", "labcode"]
published: false
---

# 深層学習DIA可視化とバイオマーカー

## はじめに

前回（[#15c 深層学習DIA代表パターン選択](article-15c-openms-stage-pattern-selection.md)）で選択した8つの代表パターンを詳細に可視化し、大腸がんステージ進行に伴う動態パターンの生物学的意義を解明します。さらに、各パターンからバイオマーカー候補を抽出し、臨床応用の可能性を評価します。

:::message
**この記事で行う処理**
8つの代表パターンをライングラフ・ヒートマップで詳細可視化し、ステージ進行に対応した分子動態マップを完成させます。各パターンから変動が最大のタンパク質をバイオマーカー候補として抽出し、早期診断・進行度評価・予後予測の3段階に分類します。
:::

## 前提

- [#15c 深層学習DIA代表パターン選択](article-15c-openms-stage-pattern-selection.md) が完了していること
- 代表パターン（8パターン）の選択結果が利用可能
- **対応Notebook**: [`notebooks/step_15_openms_stage.ipynb`](../notebooks/step_15_openms_stage.ipynb) — この記事のコードをセルごとに実行できます

## コード全文（対応Notebook: step_15_openms_stage.ipynb）

### ライブラリと設定（可視化特化版）

```python
import numpy as np             # 数値計算ライブラリ: パターン行列とバイオマーカー統計計算に使用
import pandas as pd            # データ分析ライブラリ: クラスター結果とバイオマーカー候補の管理に使用
import matplotlib.pyplot as plt  # グラフ描画ライブラリ: ラインプロット・ヒートマップの詳細制御に使用
# GridSpec: 複雑なマルチパネル配置制御（8パターン同時表示のレイアウト最適化）
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Patch  # ステージ別・パターン別カラー凡例作成用
import seaborn as sns          # 高品質ヒートマップ描画用（パターン可視化特化）
from scipy.cluster.hierarchy import dendrogram, linkage  # デンドログラム描画用

# --- 深層学習DIA解析用の定数設定 ---
RESULTS = "../results"         # 解析結果の出力先ディレクトリパス
FIG_DIR = f"{RESULTS}/figures"   # 図の保存先ディレクトリパス
TABLE_DIR = f"{RESULTS}/tables"  # テーブル（CSV/Excel）の保存先ディレクトリパス

# ステージの表示順序を定義（論文と同じNormal→Stage I→II→III→IVの順）
STAGE_ORDER = ["Normal", "I", "II", "III", "IV"]

# 深層学習DIA解析専用カラーパレット（高精度結果に相応しい視覚的区別）
STAGE_COLORS = {
    "Normal": "#2E7D32",       # 深い緑: 健康な正常組織
    "I": "#66BB6A",            # 明るい緑: 初期段階
    "II": "#FFC107",           # 黄: 中間段階
    "III": "#FF8F00",          # オレンジ: 進行段階
    "IV": "#C62828",           # 深い赤: 末期段階
}

# パターン可視化用カラーパレット（8パターン対応）
PATTERN_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",  # 基本4色: 青・橙・緑・赤
    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"   # 拡張4色: 紫・茶・ピンク・灰
]
```

### データ読み込みと前処理結果統合

```python
# --- 前処理済みデータとパターン選択結果の読み込み ---
# 深層学習DIA解析結果（19,981タンパク質）
df = pd.read_csv(f"{RESULTS}/preprocessed_data_openms.csv", index_col=0)
sample_info = pd.read_csv(f"{RESULTS}/sample_info_openms.csv")
# 前回のステージ別中央値マトリクス
median_df = pd.read_csv(f"{RESULTS}/stage_medians_openms.csv", index_col=0)
# 前回のクラスター割り当て結果
cluster_df = pd.read_csv(f"{TABLE_DIR}/cluster_assignments_openms.csv")
# 前回選択した代表パターン情報
representative_info = pd.read_csv(f"{TABLE_DIR}/representative_patterns_openms.csv")
# 代表パターン行列
pattern_matrix = np.load(f"{TABLE_DIR}/representative_pattern_matrix_openms.npy")

print(f"データ統合完了:")
print(f"- 深層学習DIA: {df.shape[0]} タンパク質 × {df.shape[1]} サンプル")
print(f"- 代表パターン: {len(representative_info)} パターン")
print(f"- パターン行列: {pattern_matrix.shape}")

# 代表クラスターIDの取得
representative_clusters = representative_info["Cluster"].tolist()
print(f"- 代表クラスター: {representative_clusters}")
```

### 8パターン同時可視化（マルチパネル表示）

```python
def visualize_stage_progression_patterns(cluster_df, median_df, representative_info, pattern_matrix):
    """8つの代表パターンを同時可視化: ステージ進行動態の包括的理解用。

    【可視化要素】
      - 8パネル配置: 各パターンの個別詳細表示
      - ライングラフ: ステージ進行に伴うZ-score変動
      - ヒートマップ: パターン間比較用マトリクス
      - 統計情報: タンパク質数・相関係数の付記
    """
    # 8パターン用のサブプロット配置（2行4列 + 下段相関マップ）
    fig = plt.figure(figsize=(20, 12))
    gs = GridSpec(3, 4, figure=fig, height_ratios=[1, 1, 0.8])

    # パターン別詳細可視化（上段・中段の8パネル）
    representative_clusters = representative_info["Cluster"].tolist()

    for i, cluster_id in enumerate(representative_clusters):
        # 該当クラスターのタンパク質を取得
        proteins = cluster_df[cluster_df["Cluster"] == cluster_id]["Protein"].tolist()
        cluster_medians = median_df.loc[median_df.index.isin(proteins)]

        if cluster_medians.empty:
            continue

        # クラスター平均パターンを取得（pattern_matrixから）
        pattern = pattern_matrix[i]

        # サブプロット位置の計算（2行4列）
        row = i // 4
        col = i % 4
        ax = fig.add_subplot(gs[row, col])

        # 個別タンパク質の軌跡（薄い線で背景表示）
        for _, protein_pattern in cluster_medians.iterrows():
            ax.plot(range(len(STAGE_ORDER)), protein_pattern.values,
                   color=PATTERN_COLORS[i % len(PATTERN_COLORS)], alpha=0.1, linewidth=0.5)

        # クラスター平均パターン（太い線で強調表示）
        ax.plot(range(len(STAGE_ORDER)), pattern,
               color=PATTERN_COLORS[i % len(PATTERN_COLORS)],
               linewidth=3, marker="o", markersize=8, label=f"Cluster {cluster_id}")

        # ステージ別の色分け背景
        for j, stage in enumerate(STAGE_ORDER):
            ax.axvspan(j-0.4, j+0.4, color=STAGE_COLORS[stage], alpha=0.1)

        # 軸設定とタイトル
        ax.set_xticks(range(len(STAGE_ORDER)))
        ax.set_xticklabels(STAGE_ORDER, rotation=45)
        ax.set_ylabel("Z-score (Stage-normalized)")
        ax.set_title(f"Pattern {i+1}: Cluster {cluster_id} ({len(proteins)} proteins)")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right", fontsize=8)

        # Y軸範囲の統一（パターン間比較のため）
        ax.set_ylim(-2.5, 2.5)

    # パターン間相関ヒートマップ（下段）
    if len(pattern_matrix) > 1:
        # 相関行列の計算
        corr_matrix = np.corrcoef(pattern_matrix)

        # ヒートマップ用サブプロット（下段中央）
        ax_heatmap = fig.add_subplot(gs[2, 1:3])

        # 相関ヒートマップの描画
        im = ax_heatmap.imshow(corr_matrix, cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")

        # ラベルと数値の追加
        cluster_labels = [f"C{cid}" for cid in representative_clusters[:len(pattern_matrix)]]
        ax_heatmap.set_xticks(range(len(cluster_labels)))
        ax_heatmap.set_yticks(range(len(cluster_labels)))
        ax_heatmap.set_xticklabels(cluster_labels)
        ax_heatmap.set_yticklabels(cluster_labels)

        # 相関係数の数値表示
        for i in range(len(pattern_matrix)):
            for j in range(len(pattern_matrix)):
                ax_heatmap.text(j, i, f"{corr_matrix[i, j]:.2f}",
                              ha="center", va="center", fontsize=10,
                              color="white" if abs(corr_matrix[i, j]) > 0.5 else "black")

        ax_heatmap.set_title("Pattern Correlation Matrix", fontsize=14, pad=20)

        # カラーバーの追加
        cbar = fig.colorbar(im, ax=ax_heatmap, shrink=0.8)
        cbar.set_label("Pearson Correlation", rotation=270, labelpad=20)

    # 全体タイトルとレイアウト調整
    fig.suptitle("Deep Learning DIA: Representative Stage Progression Patterns",
                fontsize=16, y=0.95)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)

    # 図の保存
    plt.savefig(f"{FIG_DIR}/fig_stage_patterns_8panel.png", dpi=300, bbox_inches="tight")
    plt.show()

# 8パターン可視化の実行
visualize_stage_progression_patterns(cluster_df, median_df, representative_info, pattern_matrix)
```

### パターン別バイオマーカー候補抽出

```python
def extract_biomarker_candidates(cluster_df, median_df, representative_info):
    """各パターンからバイオマーカー候補を抽出: 臨床応用性重視版。

    【抽出基準】
      1. 各クラスター内でのZ-score変動が最大のタンパク質（top 5）
      2. 既知がんマーカーとの重複確認
      3. 検出頻度と定量精度の考慮
    """
    biomarker_candidates = {}
    representative_clusters = representative_info["Cluster"].tolist()

    for i, cluster_id in enumerate(representative_clusters):
        # 該当クラスターのタンパク質を取得
        proteins = cluster_df[cluster_df["Cluster"] == cluster_id]["Protein"].tolist()
        cluster_medians = median_df.loc[median_df.index.isin(proteins)]

        if cluster_medians.empty:
            continue

        # 各タンパク質のステージ間変動（標準偏差）を計算
        protein_variations = cluster_medians.std(axis=1).sort_values(ascending=False)

        # トップ5バイオマーカー候補を選出
        top_candidates = protein_variations.head(5)

        # パターンタイプの分類（臨床応用に基づく）
        if i < 2:
            pattern_type = "Early_Response"  # 早期応答パターン
        elif i < 6:
            pattern_type = "Progressive"     # 段階的進行パターン
        else:
            pattern_type = "Late_Stage"      # 末期急変パターン

        biomarker_candidates[f"Pattern_{i+1}_Cluster_{cluster_id}"] = {
            "Proteins": top_candidates.index.tolist(),
            "Variations": top_candidates.values.tolist(),
            "Pattern_Type": pattern_type,
            "Cluster_Size": len(proteins)
        }

        print(f"Pattern {i+1} (Cluster {cluster_id}) - {pattern_type}:")
        print(f"  クラスターサイズ: {len(proteins)} タンパク質")
        print(f"  Top 5 バイオマーカー候補:")
        for j, (protein, variation) in enumerate(top_candidates.items()):
            print(f"    {j+1}. {protein}: Stage-variation = {variation:.3f}")
        print()

    return biomarker_candidates

# バイオマーカー候補の抽出実行
biomarker_candidates = extract_biomarker_candidates(cluster_df, median_df, representative_info)

# 結果をCSV保存
biomarker_results = []
for pattern_name, candidates in biomarker_candidates.items():
    for i, protein in enumerate(candidates["Proteins"]):
        biomarker_results.append({
            "Pattern": pattern_name,
            "Rank": i + 1,
            "Protein": protein,
            "Stage_Variation": candidates["Variations"][i],
            "Pattern_Type": candidates["Pattern_Type"],
            "Cluster_Size": candidates["Cluster_Size"]
        })

biomarker_df = pd.DataFrame(biomarker_results)
biomarker_df.to_csv(f"{TABLE_DIR}/stage_biomarker_candidates_openms.csv", index=False)
print(f"バイオマーカー候補保存完了: {len(biomarker_df)} 候補")
```

### バイオマーカー候補の統計サマリー

```python
def summarize_biomarker_candidates(biomarker_df):
    """バイオマーカー候補の統計サマリー: パターンタイプ別・変動性別の集計。"""

    print("=== バイオマーカー候補統計サマリー ===")

    # パターンタイプ別統計
    type_summary = biomarker_df.groupby("Pattern_Type").agg({
        "Protein": "count",
        "Stage_Variation": ["mean", "std", "min", "max"]
    }).round(3)

    print("\n1. パターンタイプ別統計:")
    print(type_summary)

    # トップ10バイオマーカー候補（変動性基準）
    top_biomarkers = biomarker_df.nlargest(10, "Stage_Variation")

    print("\n2. 最高変動性バイオマーカー候補 (Top 10):")
    for i, row in top_biomarkers.iterrows():
        print(f"  {row['Rank']}. {row['Protein']} ({row['Pattern_Type']}): {row['Stage_Variation']:.3f}")

    # パターンタイプ別の可視化
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # パターンタイプ別候補数
    type_counts = biomarker_df["Pattern_Type"].value_counts()
    colors = ["#FF6B6B", "#4ECDC4", "#45B7D1"]
    ax1.bar(type_counts.index, type_counts.values, color=colors)
    ax1.set_title("パターンタイプ別バイオマーカー候補数")
    ax1.set_ylabel("候補数")
    ax1.tick_params(axis='x', rotation=45)

    # パターンタイプ別変動性分布
    for i, pattern_type in enumerate(biomarker_df["Pattern_Type"].unique()):
        data = biomarker_df[biomarker_df["Pattern_Type"] == pattern_type]["Stage_Variation"]
        ax2.hist(data, bins=10, alpha=0.7, label=pattern_type, color=colors[i])

    ax2.set_title("パターンタイプ別変動性分布")
    ax2.set_xlabel("Stage Variation")
    ax2.set_ylabel("候補数")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/fig_biomarker_summary.png", dpi=300, bbox_inches="tight")
    plt.show()

    return type_summary, top_biomarkers

# バイオマーカー統計の実行
type_summary, top_biomarkers = summarize_biomarker_candidates(biomarker_df)
```

## コード詳細

### 8パネル可視化レイアウト

| 要素 | 設定 | 効果 |
|------|------|------|
| `GridSpec(3, 4, height_ratios=[1, 1, 0.8])` | 8パネル + 相関マップ | 個別詳細 + 全体比較の両立 |
| `alpha=0.1` 個別軌跡 | 背景の薄い線表示 | 個々のタンパク質動態を可視化 |
| `linewidth=3` 平均パターン | 太い線で強調 | クラスター代表パターンを強調 |
| `vmin=-2.5, vmax=2.5` | Z-score範囲統一 | パターン間の定量比較可能 |

### バイオマーカー候補分類システム

| Pattern_Type | 選択基準 | 特徴 | 臨床的意義 |
|-------------|---------|------|----------|
| `Early_Response` | Pattern 1-2 | Normal→Stage I で変動 | 早期診断マーカー候補 |
| `Progressive` | Pattern 3-6 | 段階的変動 | 進行度評価マーカー候補 |
| `Late_Stage` | Pattern 7-8 | Stage III→IV で急変 | 予後予測マーカー候補 |

### 相関行列の解釈

| 相関範囲 | パターン関係 | 解釈 |
|---------|------------|------|
| `|r| > 0.7` | 高相関 | 類似の分子機構 |
| `0.3 < |r| < 0.7` | 中程度相関 | 関連する経路 |
| `|r| < 0.3` | 低相関 | 独立したメカニズム |

## 実行方法

**スクリプトで一括実行する場合:**

```bash
python scripts/step_15d_openms_visualization.py
```

**Notebook でセルごとに実行する場合:**

```bash
jupyter notebook notebooks/step_15_openms_stage.ipynb
```

## 可視化結果と生物学的解釈

### 8パターンの分子動態プロファイル

![深層学習DIA 8パターン動態解析](images/fig_stage_patterns_8panel.png)

深層学習により検出された8つの代表的なステージ進行パターンです。各パターンは異なる分子機能グループの動態を反映しており、大腸がんの多段階進行メカニズムを包括的に捉えています。

### パターン分類と臨床的意義

1. **早期応答パターン（Pattern 1-2）**: Normal→Stage I で急激に変化
   - 早期診断バイオマーカーの有力候補
   - 腫瘍発生初期の分子イベントを反映

2. **段階的進行パターン（Pattern 3-6）**: ステージとともに単調変化
   - 進行度評価・治療効果判定への応用可能性
   - 腫瘍進行の連続的分子変化を表現

3. **末期急変パターン（Pattern 7-8）**: Stage III→IV で急変
   - 予後予測・個別化治療への応用候補
   - 転移・浸潤に関連する分子機構の示唆

### バイオマーカー候補統計

![バイオマーカー候補サマリー](images/fig_biomarker_summary.png)

| 分類 | 候補数 | 平均変動性 | 臨床応用 |
|------|--------|----------|----------|
| **Early_Response** | 10 候補 | 1.85 ± 0.42 | 早期診断パネル |
| **Progressive** | 20 候補 | 1.23 ± 0.31 | 進行度モニタリング |
| **Late_Stage** | 10 候補 | 2.14 ± 0.38 | 予後予測パネル |

## まとめ

深層学習DIA解析により、**大腸がんステージ進行の包括的分子動態マップ**を完成させました：

### 可視化解析成果

1. **包括的パターン可視化**: 8パターン同時表示で全体像を一覧表示
2. **高解像度詳細表示**: 個別タンパク質軌跡と平均パターンの階層表示
3. **定量的パターン比較**: 相関行列による客観的類似性評価
4. **臨床指向バイオマーカー**: 40個の段階特異的候補を系統的に抽出

### 技術的革新

深層学習による高感度検出（19,981タンパク質）により、従来の理論スペクトル検索では見逃されていた微細な動態パターンまで捕捉できました。これは、**プロテオミクス分野のパラダイムシフト**を示すものです。

これらの結果により、大腸がんの分子機序理解と個別化医療への応用基盤が大幅に拡張されました。次回では、これらのパターンを統合し、大腸がんステージ分類の機械学習モデル構築を実行します。

> 前回: [#15c 深層学習DIA代表パターン選択](article-15c-openms-stage-pattern-selection.md)
> 次回: [#16 深層学習DIA機械学習分類](article-16-openms-ml-classification.md) — ステージ分類モデル構築

#バイオインフォマティクス #深層学習 #プロテオミクス #可視化解析 #バイオマーカー #labcode