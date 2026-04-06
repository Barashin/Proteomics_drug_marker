---
title: "DIA-MSプロテオミクス論文再現のまとめと論文化への道筋【論文再現シリーズ #10】"
emoji: "🎓"
type: "tech"
topics: ["proteomics", "bioinformatics", "research", "labcode"]
published: false
---

# DIA-MSプロテオミクス論文再現のまとめと論文化への道筋

## はじめに

シリーズ最終回です。これまでの記事でToyota et al. 2025のdry解析パイプラインを無料ツールのみで再現してきました。この記事では、全体のまとめと、**このパイプラインを使って自分の論文を書く方法**を紹介します。

## シリーズの振り返り

### 再現した解析パイプライン

```
DIA-MS データ (mzML, 18ファイル, 9患者分)
  → RAW→mzML変換（msconvert / ThermoRawFileParser）... #3
  → sage-proteomics（タンパク質同定・定量、6.5分で完了）... #4
  → 前処理（Log2, フィルタ, 補完）        ... #5
  → 全体像可視化（相関, クラスタリング, PCA） ... #6
  → 差分発現解析（t-test, Volcano）       ... #7
  → COSMIC照合（がん関連タンパク質）       ... #8
  → ステージ別解析（ANOVA, クラスター）    ... #9
```

### 本書で得た実測値（sage + 18ファイル）

| 段階 | 本書 | 論文 (DIA-NN + 32ファイル) |
|------|------|--------------------------|
| 同定タンパク質 | **2,260** | 10,329 |
| 前処理後 | **2,234** | - |
| 有意差 (p<0.05, FC>2) | **1,015** (↑818 / ↓197) | 2,642 (↑1,475 / ↓1,167) |
| PCA PC1 寄与率 | **45.3%** | **42.1%** |
| COSMIC 全がん関連 | 21 / 200 (10.5%) | 531 / 748 (71%) |
| COSMIC 主要ドライバー | **KRAS, CTNNB1, PIK3CA 等 21個** | 同等 |
| ステージ ANOVA有意 (FDR<0.01) | **151 タンパク質、30クラスター** | - |

### 使用ツール（すべて無料・商用利用可能）

| ツール | 用途 | ライセンス |
|--------|------|-----------|
| **sage-proteomics** | タンパク質同定・定量（DIA） | **MIT**（商用完全OK、Rust製） |
| pymzml | mzMLメタ情報検査 | MIT |
| Python + numpy/pandas | データ操作 | BSD |
| scipy | 統計検定（t-test, ANOVA） | BSD |
| scikit-learn | PCA | BSD |
| seaborn/matplotlib | ヒートマップ・プロット | BSD |
| statsmodels | FDR補正 | BSD |

**ポイント**: すべて **商用利用に一切の制約がない** ライセンスで、LabCode 書籍・企業の研究・受託解析のいずれでも安心して使用できます。

### 再現した論文のFigure

| Figure | 内容 | 記事 |
|--------|------|------|
| Figure 1a | 相関行列ヒートマップ | #6 |
| Figure 1b | 階層的クラスタリング | #6 |
| Figure 1c | PCA | #6 |
| Figure 2 | Volcanoプロット・Top N解析 | #7 |
| Figure 3 | ステージ別プロファイルプロット | #9 |

## 論文化への道筋

このパイプラインを使って論文を書くための3つのアプローチを紹介します。

### アプローチ1: 別の疾患データに適用（最もおすすめ）

同じDIA-MSパイプラインを、**別の疾患の公開データ**に適用します。

```
1. ProteomeXchange で別の疾患のDIA-MSデータを探す
2. 本パイプラインで解析
3. 疾患特異的なバイオマーカー候補を報告
```

**候補となるデータ例:**
- 肝がん（HCC）、膵がん、肺がん等のDIA-MSデータ
- 炎症性腸疾患（IBD）のプロテオミクスデータ
- 自己免疫疾患のプロテオミクスデータ

**想定投稿先:** Proteomes (IF 3.0), Proteomics Clinical Applications (IF 3.6), Journal of Proteome Research (IF 4.0)

### アプローチ2: 機械学習を追加

本論文は統計解析のみですが、**機械学習（Random Forest, XGBoost, LASSO等）**を追加してバイオマーカーパネルを構築すれば、新規性が高まります。

```
本パイプラインの差分発現結果
  → 特徴量選択（LASSO, Boruta等）
  → 機械学習モデル構築（RF, XGBoost）
  → ROC-AUC で診断性能を評価
  → 独立コホートで検証
```

**想定投稿先:** Advanced Science (IF 14.3), Nature Communications (IF 14.7)

### アプローチ3: 方法論の比較

Perseus vs Pythonの解析結果を比較し、**オープンソースパイプラインの妥当性**を示す方法論論文を書くことも可能です。

```
同一データに対して:
  ① Perseus で解析
  ② Python（本パイプライン）で解析
  → 結果を比較（相関、一致率、再現性）
  → オープンソース代替の有用性を主張
```

**想定投稿先:** Journal of Proteomics (IF 3.5), Bioinformatics (IF 5.8)

## パイプラインの拡張案

| 拡張 | 追加ツール | 内容 |
|------|----------|------|
| パスウェイ解析 | gseapy (Python) | GO/KEGG/GSEAエンリッチメント |
| ネットワーク解析 | networkx (Python) | PPI（タンパク質間相互作用）ネットワーク |
| 生存解析 | lifelines (Python) | 臨床データとの統合 |
| マルチオミクス | - | トランスクリプトームとの統合 |

## 全シリーズのコード

全スクリプトは一括実行可能です：

```bash
micromamba activate crc-proteomics
bash scripts/run_all.sh
```

## シリーズ全記事リンク

| # | タイトル |
|---|---------|
| #0 | [論文紹介と全体像](article-00-introduction.md) |
| #1 | [環境構築](article-01-setup.md) |
| #2 | [データ取得](article-02-data.md) |
| #3 | [RAW→mzML変換](article-03-convert.md) |
| #4 | [sageで同定・定量](article-04-sage.md) |
| #5 | [前処理](article-05-preprocess.md) |
| #6 | [全体像の可視化](article-06-visualization.md) |
| #7 | [差分発現解析](article-07-differential.md) |
| #8 | [COSMIC照合](article-08-cosmic.md) |
| #9 | [ステージ別解析](article-09-stage.md) |
| #10 | [まとめと次のステップ](article-10-conclusion.md) |

## おわりに

DIA-MSプロテオミクスの解析パイプラインを、無料ツールのみで構築しました。このパイプラインは大腸がんに限らず、**あらゆるDIA-MSデータに適用可能**です。

プロテオミクスは今後ますます重要になる分野です。本シリーズが、皆さんの研究の一助になれば幸いです。

> 前回: [#9 ステージ別解析](article-09-stage.md)

#バイオインフォマティクス #プロテオミクス #論文執筆 #labcode
