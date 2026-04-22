---
title: "DIA-MSプロテオミクスで大腸がんバイオマーカーを探る【論文再現シリーズ #0】"
emoji: "🔬"
type: "tech"
topics: ["proteomics", "bioinformatics", "python", "labcode"]
published: false
---

# DIA-MSプロテオミクスで大腸がんバイオマーカーを探る

:::message
**LC-MSを使った解析が初めての方へ**
本シリーズに入る前に、AJACSの解説動画をご覧いただくことをおすすめします。LC-MS（液体クロマトグラフィー質量分析）の基礎原理からデータ解析の考え方まで、わかりやすく解説されています。

- [AJACS LC-MS解説動画（前編）](https://youtu.be/I9cArPIAkrw?si=zurCmWwMXxLKukdm) — LC-MSの基本原理とプロテオミクスにおけるデータ取得の流れを解説
- [AJACS LC-MS解説動画（後編）](https://youtu.be/YSJ0BhvWWFw?si=ZQS4xxrVtdnpZnxH) — LC-MSデータの解析手法やバイオインフォマティクスへの応用を解説

AJACS（あじゃっくす）は、バイオサイエンスデータベースセンター（NBDC）が主催するバイオインフォマティクスのトレーニングプログラムです。初学者向けの講義動画が多数公開されており、プロテオミクスに限らず幅広い分野の基礎を学ぶことができます。
:::

## はじめに

この記事シリーズでは、DIA-MSプロテオミクスの論文（Toyota et al., Proteomes 2025）のdry解析パイプラインを**無料ツールだけで完全再現**します。

プロテオミクス（タンパク質の網羅的解析）は近年急速に発展している分野ですが、解析パイプラインの構築は初心者にとって大きなハードルです。この記事では、**コピペで動くコード**と**丁寧な解説**で、誰でもDIA-MSデータの解析ができるようになることを目指します。

## フォローする論文

**Toyota N, Konno R, et al.** "Identification of Cancer-Associated Proteins in Colorectal Cancer Using Mass Spectrometry." *Proteomes* 2025; 13(3):38.
DOI: https://doi.org/10.3390/proteomes13030038

:::message alert
**本シリーズに掲載する図について**

本シリーズ中の図（ヒートマップ・PCA・Volcanoプロット等）は、すべて **本書のパイプラインで独自に生成した再現図** です。論文の Figure を直接転載したものではありません。使用ツールの違い（DIA-NN → sage）や同定タンパク質数の差により、論文オリジナルの Figure とは**細部が異なる場合があります**。論文と対比する際は本シリーズの Figure 番号と論文の Figure 番号が1対1対応するとは限らない点にご注意ください。
:::

### 論文の概要

- 大腸がん患者16人の腫瘍/正常組織をDIA-MSで解析
- **10,329タンパク質**を同定
- COSMICデータベースと照合し**531のがん関連タンパク質**を特定
- ステージ進行に伴う**タンパク質変動パターン**を発見

## 解析パイプライン

```
DIA-MS (Orbitrap Exploris 480)
  → sage-proteomics (タンパク質同定・定量、MITライセンス)
  → 前処理 (Log2変換, 欠損値補完)
  → 可視化 (相関行列, クラスタリング, PCA)
  → 差分発現解析 (Welch's t-test)
  → COSMIC照合 (がん関連タンパク質)
  → ステージ別解析 (ANOVA, クラスター)
```

## 無料ツールへの置き換え

論文では Perseus（統計解析）と DIA-NN（DIA解析）を使用していますが、どちらも**商用利用に制限がある**ため、本シリーズでは **MIT/BSD ライセンスのツール** で完全代替します。

| 論文のツール | 本シリーズ | ライセンス |
|------------|----------|----------|
| DIA-NN | **sage-proteomics** | **MIT**（商用完全OK、Rust製） |
| Perseus v1.6.15.0 | Python (scipy, sklearn, statsmodels) | BSD |
| Seaborn / Matplotlib | Seaborn / Matplotlib | BSD |

本書では論文全32ファイル（16患者分）を扱い、sage で同定した **2,110 タンパク質**、有意差 **1,055 個**、COSMIC 主要ドライバー（KRAS, CTNNB1, PIK3CA 等）**22個** を検出しました。

## シリーズ構成

| # | タイトル | 内容 |
|---|---------|------|
| **#0** | **本記事** | 論文紹介と全体像 |
| #1 | 環境構築 | micromamba + Python環境 |
| #2 | データ取得 | ProteomeXchangeからダウンロード |
| #3 | RAW→mzML変換 | Thermo RAWファイルをmzMLに変換 |
| #4 | sageで同定・定量 | タンパク質の同定と定量（商用OK） |
| #5 | OpenMS + AlphaPeptDeep | 深層学習ベースの代替パイプライン |
| #6 | 前処理 | Log2変換・欠損値補完 |
| #7 | 全体像の可視化 | 相関行列・クラスタリング・PCA |
| #8 | 差分発現解析 | Welch's t-test・Volcanoプロット |
| #9 | COSMIC照合 | がん関連タンパク質の同定 |
| #10 | ステージ別解析 | ANOVA・クラスター分析 |
| #11 | まとめと次のステップ | 論文化への道筋 |

> 次回: [#1 環境構築](article-01-setup.md) — micromambaとPython環境のセットアップ

#バイオインフォマティクス #プロテオミクス #DIA-MS #labcode
