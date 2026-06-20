---
title: "ステージ別タンパク質発現パターンをクラスタリングで可視化【論文再現シリーズ #11b】"
emoji: "📈"
type: "tech"
topics: ["proteomics", "clustering", "visualization", "labcode"]
published: false
---

# ステージ別タンパク質発現パターンをクラスタリングで可視化

## はじめに

前回（[#11a ANOVA解析](article-11a-stage-anova.md)）で特定した有意差タンパク質について、**階層的クラスタリング**を実行してステージ進行に伴う発現パターンを可視化します。Toyota et al. 2025のFigure 3を再現します。

> **📈 INFO**
>
**この記事で行う処理**
ANOVA有意差タンパク質を階層的クラスタリングで30グループに分割し、各グループのステージ別発現プロファイルを作成します。最終的にヒートマップとプロファイルプロットでステージ進行パターンを可視化し、疾患進行のバイオマーカー候補を特定します。

## 前提

- [#11a ANOVA解析](article-11a-stage-anova.md) が完了していること
- `anova_stage_significant.csv` が生成されていること
- **対応Notebook**: [`notebooks/step_10.ipynb`](../notebooks/step_10.ipynb) の後半部分

## コード全文（クラスタリング・可視化部分）

### ライブラリの追加導入

```python
import seaborn as sns          # 統計データ可視化ライブラリ: ヒートマップ描画に使用
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # 複雑なレイアウト制御
from matplotlib.patches import Patch  # 凡例用のカラーパッチを作成
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram  # 階層的クラスタリング
```

### 有意差タンパク質データの準備

```python
# 前回の結果を読み込み
significant_proteins_info = pd.read_csv(f"{RESULTS_DIR}/anova_stage_significant.csv", index_col=0)
significant_protein_names = significant_proteins_info.index.tolist()

print(f"ANOVA有意差タンパク質数: {len(significant_protein_names)}")

# 有意差タンパク質の発現データを抽出
significant_data = df.loc[significant_protein_names].copy()

# ステージ別中央値を計算（クラスタリング用）
stage_groups = {
    'Normal': [i for i, s in enumerate(sample_stage) if s == 'Normal'],
    'Stage_I': [i for i, s in enumerate(sample_stage) if s == 'Stage_I'],
    'Stage_II': [i for i, s in enumerate(sample_stage) if s == 'Stage_II'],
    'Stage_III': [i for i, s in enumerate(sample_stage) if s == 'Stage_III'],
    'Stage_IV': [i for i, s in enumerate(sample_stage) if s == 'Stage_IV']
}

# 各ステージの中央値を計算
median_by_stage = pd.DataFrame()
for stage, indices in stage_groups.items():
    if len(indices) > 0:
        stage_data = significant_data.iloc[:, indices].median(axis=1)
        median_by_stage[stage] = stage_data

print(f"ステージ別中央値マトリクス: {median_by_stage.shape[0]} proteins × {median_by_stage.shape[1]} stages")
median_by_stage.head()
```

### 階層的クラスタリング

```python
# Z-score正規化（行方向: 各タンパク質の発現パターンを標準化）
from scipy.stats import zscore
median_zscore = median_by_stage.apply(zscore, axis=1)

# 欠損値やinf値をチェック・除去
median_zscore_clean = median_zscore.replace([np.inf, -np.inf], np.nan).dropna()
print(f"クリーニング後: {median_zscore_clean.shape[0]} proteins")

# 階層的クラスタリング（Ward法、ユークリッド距離）
linkage_matrix = linkage(median_zscore_clean.values, method='ward', metric='euclidean')

# クラスター数を指定（Toyota論文に基づき30クラスター）
n_clusters = 30
cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')

# クラスター情報をデータフレームに追加
median_zscore_clean['cluster'] = cluster_labels
cluster_counts = pd.Series(cluster_labels).value_counts().sort_index()

print(f"クラスタリング完了: {n_clusters} clusters")
print(f"クラスター別タンパク質数: {cluster_counts.head()}")
```

### Figure 3: ステージ別ヒートマップとプロファイルプロット

```python
# 論文のFigure 3を再現: ヒートマップ + プロファイルプロット
def create_stage_heatmap_with_profiles(data_zscore, cluster_labels, n_clusters=30):
    """
    ステージ別ヒートマップとクラスタープロファイルプロットを作成
    Toyota et al. 2025 Figure 3の再現
    """
    # クラスター順にデータを並び替え
    data_with_clusters = data_zscore.copy()
    data_with_clusters['cluster'] = cluster_labels
    data_sorted = data_with_clusters.sort_values('cluster')

    # プロット用データ（クラスター列を除去）
    plot_data = data_sorted.drop('cluster', axis=1)

    # 複合レイアウトの設定
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 3],
                  hspace=0.3, wspace=0.3)

    # === メインヒートマップ（右下）===
    ax_heatmap = fig.add_subplot(gs[1, 0])

    # ヒートマップの描画（クラスター境界付きで）
    sns.heatmap(plot_data,
                cmap='RdBu_r',           # 赤（高発現）-白-青（低発現）
                center=0,                # Z-scoreの中心を0に設定
                vmin=-3, vmax=3,         # 色の範囲を±3に制限
                cbar_kws={'label': 'Z-score'},
                xticklabels=True,
                yticklabels=False,       # タンパク質名は非表示（多すぎるため）
                ax=ax_heatmap)

    ax_heatmap.set_title('Protein Expression by Cancer Stage\n(Clustered by Expression Pattern)',
                        fontsize=14, fontweight='bold')
    ax_heatmap.set_xlabel('Cancer Stage', fontsize=12)
    ax_heatmap.set_ylabel(f'{len(plot_data)} Proteins (ANOVA significant)', fontsize=12)

    # クラスター境界線を追加
    cluster_boundaries = []
    current_pos = 0
    for cluster_id in sorted(data_sorted['cluster'].unique()):
        cluster_size = (data_sorted['cluster'] == cluster_id).sum()
        current_pos += cluster_size
        cluster_boundaries.append(current_pos)

    for boundary in cluster_boundaries[:-1]:  # 最後の境界は不要
        ax_heatmap.axhline(y=boundary, color='white', linewidth=0.5, alpha=0.7)

    # === プロファイルプロット（左下）===
    ax_profiles = fig.add_subplot(gs[1, 1])

    # 各クラスターの平均プロファイルを計算・描画
    stage_order = ['Normal', 'Stage_I', 'Stage_II', 'Stage_III', 'Stage_IV']
    colors = plt.cm.tab20(np.linspace(0, 1, min(n_clusters, 20)))  # 最大20色

    for i, cluster_id in enumerate(sorted(data_sorted['cluster'].unique())):
        if i >= 20:  # 色の制限
            break
        cluster_data = data_sorted[data_sorted['cluster'] == cluster_id]
        cluster_profile = cluster_data[stage_order].mean()

        ax_profiles.plot(range(len(stage_order)), cluster_profile.values,
                        color=colors[i], alpha=0.7, linewidth=1.5,
                        label=f'Cluster {cluster_id} (n={len(cluster_data)})')

    ax_profiles.set_xticks(range(len(stage_order)))
    ax_profiles.set_xticklabels(stage_order, rotation=45)
    ax_profiles.set_ylabel('Mean Z-score')
    ax_profiles.set_title('Cluster Expression Profiles')
    ax_profiles.grid(True, alpha=0.3)
    ax_profiles.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

    # === 上部: ステージ別サンプル数バー（上左）===
    ax_bar = fig.add_subplot(gs[0, 0])
    stage_counts = pd.Series(sample_stage).value_counts()[stage_order]
    bars = ax_bar.bar(range(len(stage_order)), stage_counts.values,
                      color=['lightgray', 'lightblue', 'lightgreen', 'orange', 'red'],
                      alpha=0.7)

    ax_bar.set_xticks(range(len(stage_order)))
    ax_bar.set_xticklabels(stage_order)
    ax_bar.set_ylabel('Sample Count')
    ax_bar.set_title('Sample Distribution by Stage')

    # バーの上に数値を表示
    for i, (bar, count) in enumerate(zip(bars, stage_counts.values)):
        ax_bar.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                   str(count), ha='center', va='bottom', fontweight='bold')

    plt.suptitle('Figure 3: Stage-Associated Protein Expression Patterns\n' +
                 f'({len(plot_data)} ANOVA-significant proteins, {n_clusters} clusters)',
                 fontsize=16, fontweight='bold', y=0.95)

    # 保存
    plt.tight_layout()
    plt.savefig(f"{FIGURES_DIR}/fig3_stage_heatmap.png", dpi=300, bbox_inches='tight')
    plt.show()

# Figure 3の作成実行
create_stage_heatmap_with_profiles(median_zscore_clean, cluster_labels, n_clusters)
```

### クラスタープロファイルの詳細解析

```python
# 代表的なクラスター（進行パターン別）の抽出
def analyze_progression_patterns(data_zscore, cluster_labels):
    """
    ステージ進行パターンの分類と代表クラスターの同定
    """
    stage_order = ['Normal', 'Stage_I', 'Stage_II', 'Stage_III', 'Stage_IV']
    cluster_patterns = {}

    for cluster_id in sorted(np.unique(cluster_labels)):
        # 該当クラスターのタンパク質を抽出
        cluster_mask = cluster_labels == cluster_id
        cluster_data = data_zscore[cluster_mask]

        # クラスター平均プロファイル
        profile = cluster_data[stage_order].mean()
        cluster_size = len(cluster_data)

        # 進行パターンの分類（簡略版）
        normal_val = profile['Normal']
        stage4_val = profile['Stage_IV']

        if stage4_val - normal_val > 1:
            pattern_type = "Progressive Increase"
        elif normal_val - stage4_val > 1:
            pattern_type = "Progressive Decrease"
        else:
            pattern_type = "Complex/Stable"

        cluster_patterns[cluster_id] = {
            'pattern_type': pattern_type,
            'size': cluster_size,
            'profile': profile.values,
            'normal_to_stage4_change': stage4_val - normal_val
        }

    return cluster_patterns

# パターン解析の実行
progression_patterns = analyze_progression_patterns(median_zscore_clean, cluster_labels)

# パターン別のクラスター数を集計
pattern_summary = {}
for cluster_id, info in progression_patterns.items():
    pattern = info['pattern_type']
    if pattern not in pattern_summary:
        pattern_summary[pattern] = {'count': 0, 'total_proteins': 0}
    pattern_summary[pattern]['count'] += 1
    pattern_summary[pattern]['total_proteins'] += info['size']

print("ステージ進行パターン別の分類:")
for pattern, summary in pattern_summary.items():
    print(f"- {pattern}: {summary['count']} clusters, {summary['total_proteins']} proteins")
```

## コード詳細

### 階層的クラスタリングの選択理由

- **Ward法**: クラスター内分散を最小化する手法。バランスの取れたクラスターを生成
- **ユークリッド距離**: 発現パターンの類似性を効果的に測定
- **Z-score正規化**: タンパク質間の発現レベル差を除去し、パターンに着目

### 可視化の工夫点

| 要素 | 設定 | 理由 |
|------|------|------|
| **カラーマップ** | RdBu_r | 赤（高発現）-青（低発現）で直感的 |
| **スケール範囲** | ±3 Z-score | 極値の影響を抑制、パターンを強調 |
| **クラスター境界** | 白線 | 視覚的にグループを分離 |

## 実行方法

**Notebook でセルごとに実行する場合:**
```bash
jupyter notebook notebooks/step_10.ipynb
```
（クラスタリング・可視化セクションを実行）

## まとめ

ANOVA有意差タンパク質の階層的クラスタリングにより、大腸がんのステージ進行に伴う発現パターンを可視化しました。この解析により、疾患進行のバイオマーカー候補タンパク質群が明確になります。

**主な発見**:
- ステージ進行に伴って一貫して増加するタンパク質群
- 段階的に減少するタンパク質群
- 複雑な変動パターンを示すタンパク質群

> 前回: [#11a ANOVA解析](article-11a-stage-anova.md)
> 次回: [#11 深層学習DIA前処理](article-11-openms-preprocess.md) — OpenMS結果の前処理

#バイオインフォマティクス #プロテオミクス #クラスタリング #可視化 #labcode