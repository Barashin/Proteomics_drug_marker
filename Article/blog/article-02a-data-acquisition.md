---
title: "ProteomeXchangeから大腸がんプロテオミクスデータを取得する【論文再現シリーズ #2a】"
emoji: "📥"
type: "tech"
topics: ["proteomics", "bioinformatics", "labcode"]
published: false
---

# ProteomeXchangeから大腸がんプロテオミクスデータを取得する

## はじめに

この記事では、Toyota et al. 2025 の公開データをProteomeXchangeからダウンロードします。プロテオミクス論文のデータは通常、ProteomeXchangeコンソーシアムを通じて公開されており、誰でも無料でアクセスできます。

> **📝 INFO**
>
**この記事で行う処理**
ProteomeXchange（国際プロテオミクスデータベース）から Toyota et al. 2025 論文のDIA-MSデータ（PXD058672）をブラウザ経由でダウンロードします。16患者分（Normal/Tumor各1ペア）の計32ファイル（約40GB）を取得し、以降の解析ステップで使用するmzMLデータを準備します。

## 前提

- [#1 環境構築](article-01-setup.md) が完了していること
- 約40GBのストレージ容量があること
- **対応Notebook**: [`notebooks/step_02.ipynb`](../notebooks/step_02.ipynb) — この記事のコードをセルごとに実行できます

## ProteomeXchangeとは

ProteomeXchange（https://www.proteomexchange.org）は、質量分析ベースのプロテオミクスデータを共有する国際コンソーシアムです。論文で使用されたRAWデータや解析結果がここに登録されています。

本論文のデータは以下に登録されています：

| リポジトリ | ID | 内容 |
|----------|-----|------|
| ProteomeXchange | **PXD058672** | RAWファイル, mzML |
| jPOST | **JPST003422** | DIA解析出力, 定量データ |

## ブラウザでGUIからダウンロードする

ブラウザからGUI操作でデータをダウンロードする方法です。

### Step 2-1: ProteomeXchangeでデータセットを検索する

ProteomeCentral（https://proteomecentral.proteomexchange.org）にアクセスし、トップページの **「Datasets」** をクリックします。

![ProteomeCentral トップページ — 「Datasets」をクリック](images/Datasets_choose.png)

Datasetsページに移動したら、左側の **Filter** 欄にデータセットID `PXD058672` を入力して検索します。

![検索バーに PXD058672 を入力](images/検索.png)

検索結果に該当するデータセット（PXD058672）が1件表示されます。**Toyota et al. (2025)** の論文であることを確認し、赤枠で囲まれたデータセット行をクリックして詳細ページに移動します。

![検索結果 — PXD058672 のデータセットをクリック](images/Dataset_click.png)

データセットの詳細ページが表示されたら、ページ下部の **「jPOST dataset」** のリンクをクリックして、jPOSTリポジトリに移動します。

https://repository.jpostdb.org/entry/JPST003422

![jPOST dataset のリンクをクリック](images/JPOST_dataset.png)

### Step 2-2: jPOSTリポジトリのデータセット概要を確認する

本論文のデータは **jPOSTrepo**（ProteomeXchangeのパートナーリポジトリ）にホストされています。データセット概要ページ（下図）で以下の情報を確認できます：

- **Project title**: 論文のタイトル
- **Keywords**: colorectal cancer, surgical specimens, tumor tissues 等
- **Dataset ID**: JPST003422
- **Publication(s)**: 論文のDOIとリンク

> **ポイント:** ページ下部の **"Files"** セクションに、ダウンロード可能なファイル一覧があります。

### Step 2-3: 「Download all」でデータをダウンロードする

ページ下部の **"Files"** セクションにある **「Download all」** ボタン（青色、約38.74 GB）をクリックして、全データを一括ダウンロードします。

![jPOSTデータセット概要ページ — Filesセクションの「Download all」ボタン](images/Download_all.png)

> **注意:** 全ファイルで約40GBあるため、ダウンロードには時間がかかります。安定したネットワーク環境での実行を推奨します。
>
> **補足:** 本シリーズではRAWデータをさらに処理したmzMLファイルを用いて解析を行います。そのため、mzMLファイルのみのダウンロードでも問題ありません。jPOSTの "Files" セクションからmzMLファイルだけを選択してダウンロードすれば、容量を大幅に節約できます。

## RAWデータの構成

jPOSTからダウンロードしたファイルは以下の通りです：

| ファイル名パターン | 形式 | 内容 |
|----------------|------|------|
| `CRC01-N.raw` | RAW | 患者1の**正常組織**のDIA-MSデータ |
| `CRC01-T.raw` | RAW | 患者1の**腫瘍組織**のDIA-MSデータ |
| `CRC02-N.raw` 〜 `CRC16-T.raw` | RAW | 患者2〜16の正常/腫瘍組織 |

- 全16患者 × 2条件（Normal/Tumor）= **32ファイル**
- 各ファイル約1〜2GB、合計約40GB
- 本シリーズでは、sage-proteomicsを使ってmzMLからタンパク質を同定・定量します

> **補足:** 本シリーズでは全16患者分（CRC01-CRC16）の32ファイルを使用して解析を行います。RAWファイルからmzMLへの変換手順については、[#3a RAW→mzML変換](article-03a-convert-basics.md)を参照してください。

### 実験プリセット情報

jPOSTの各ファイルには以下の実験条件が記録されています：

| タグ | 項目 | 値 | 意味 |
|:---:|------|-----|------|
| **S** | Sample | Human sample | ヒト由来サンプル |
| **F** | Fractionation | SP3 Digestion | SP3法（磁気ビーズ）でタンパク質を消化 |
| **E** | Enzyme/Mod. | Human-Car | ヒトプロテオーム + カルバミドメチル化（Cys修飾） |
| **M** | MS mode | Orbitrap Exploris 480-DIA | Orbitrap Exploris 480でDIAモード測定 |

## まとめ

ProteomeXchange / jPOST からDIA-MSのRAWデータとmzMLファイルをダウンロードしました。DIAデータの理論的背景とファイル形式の詳細については次章で解説します。

この手順は**他の論文のDIA-MSデータ**にもそのまま適用できます。ProteomeXchangeにデータが公開されている論文であれば、データセットIDを変えるだけで同じ流れでデータを取得できます。

> 前回: [#1 環境構築](article-01-setup.md)
> 次回: [#2b DIA-MSデータの理解](article-02b-data-formats.md) — DIA-MSの理論とファイル形式

#バイオインフォマティクス #プロテオミクス #ProteomeXchange #labcode