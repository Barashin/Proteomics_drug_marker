---
title: "Thermo RAWファイルをmzMLに変換する【論文再現シリーズ #3】"
emoji: "🔄"
type: "tech"
topics: ["proteomics", "bioinformatics", "msconvert", "labcode"]
published: false
---

# Thermo RAWファイルをmzMLに変換する

## はじめに

ProteomeXchange / jPOST で公開されているプロテオミクスデータの多くは、Thermo Fisher Scientific の質量分析計で取得された **.raw 形式** のみで提供されています。.rawはThermo独自のバイナリ形式なので、そのままでは多くの解析ツールで読み込めません。

この記事では、RAWファイルをオープンフォーマットである **mzML** に変換する方法を、**無料かつ商用利用可能**なツールで解説します。

## 前提

- [#2 データ取得](article-02-data.md) が完了し、RAWファイルが手元にあること
- **Windows PC**（または Windows 仮想マシン）を推奨

## プロテオミクス解析の全体像

まず、プロテオミクス解析の全体像を整理しておきましょう。質量分析による定量プロテオミクスは、大きく以下の流れで進みます：

```
  装置（LC-MS）
    ↓
  Rawファイル（.raw / .d / .wiff …）※ベンダー独自形式
    ↓  ← msconvert（ProteoWizard）でここを変換
  mzML（オープン標準形式）
    ↓  ← sage / MaxQuant / FragPipe など
  ペプチド・タンパク質 同定 + 定量テーブル
    ↓  ← Perseus / R / Python
  統計解析・可視化（差次発現、パスウェイ解析 …）
```

本記事では、この中の **Raw → mzML** の変換ステップを扱います。

## なぜ mzML に変換するのか

### ベンダーごとに Raw ファイル形式が異なる

質量分析計のメーカー（ベンダー）ごとに、Rawファイルの形式はバラバラです。**同じ「生データ」でも、ベンダーが違えば拡張子も内部構造も別物**になります。

| ベンダー | 拡張子 | 構造 | 代表的な装置 |
|---------|-------|------|-------------|
| **Thermo Fisher Scientific** | `.raw` | 単一バイナリファイル | Orbitrap シリーズ, Q Exactive, Exploris |
| **Bruker** | `.d` | ディレクトリ（フォルダ） | timsTOF, maXis |
| **SCIEX** | `.wiff` + `.wiff.scan` | 複数ファイル | TripleTOF, ZenoTOF |
| **Waters** | `.raw` | ディレクトリ（フォルダ） | Synapt, Xevo |
| **Agilent** | `.d` | ディレクトリ（フォルダ） | Q-TOF 6500シリーズ |

> **注意:** Thermoの `.raw` と Watersの `.raw` は**名前は同じでも全く別の形式**です。Bruker/Agilent/Watersは「ディレクトリ形式」なので、エクスプローラーでは単一ファイルに見えても実体はフォルダです。

### mzML による統一

これらのベンダー固有形式のままでは、解析ツールごとに対応可否が異なり、**同じパイプラインで異なるベンダーのデータを扱えません**。そこで登場するのが **mzML** です。

**mzML は HUPO-PSI（Human Proteome Organisation - Proteomics Standards Initiative）が定める国際標準のオープンフォーマット**で、どのベンダーのRawデータも同じ構造のXMLに変換できます。

```
Thermo .raw  ─┐
Bruker .d    ─┤
SCIEX .wiff  ─┼─→ msconvert ─→ .mzML ─→ 統一された解析パイプライン
Waters .raw  ─┤
Agilent .d   ─┘
```

| | ベンダー独自Raw | .mzML |
|---|---|---|
| **形式** | ベンダーごとにバラバラ | オープンXML（統一） |
| **読み込めるツール** | ベンダー対応ツールに限定 | **ほぼすべての解析ツール** |
| **OS依存** | Windows向けライブラリが必要なことが多い | OS非依存 |
| **標準規格** | ベンダー独自 | HUPO-PSI国際標準 |
| **再現性** | ベンダーツールに依存 | 誰でも同じ環境で再解析可能 |

**mzMLに変換しておくメリット:**

1. **ツール選択の自由度** — sage・OpenMS・MaxQuant・FragPipe・Skyline など、あらゆるツールで利用可能
2. **ベンダー横断的な解析** — 異なる装置で取得したデータを同じパイプラインで処理できる
3. **再現性・長期保存** — オープン標準なので、10年後も確実に読める
4. **OS非依存** — macOS / Linux でも読み込み可能

> **補足:** 一部のツールはRAW（Thermo形式）を直接読めますが、**sage・OpenMS・MaxQuant・多くのRパッケージはmzMLを要求**します。本シリーズのように複数ツールを組み合わせる場合、mzMLに統一しておくのが無難です。

---

## 変換ツールの選択

RAWからmzMLへの変換に使える主なツールは以下の2つです。どちらも無料・商用利用可能です。

| ツール | ライセンス | 推奨環境 | 特徴 |
|--------|-----------|---------|------|
| **msconvert**（ProteoWizard） | Apache 2.0 | Windows | 事実上の標準、GUI版あり、多機能 |
| **ThermoRawFileParser** | Apache 2.0 | Windows / Linux | 軽量・高速、CLI特化 |

> **重要:** いずれも**Thermo社の RawFileReader ライブラリ**（無料だが独自ライセンス、Windows専用DLL）に依存しています。このためmacOS（特にApple Silicon）では動作に制約があります（詳細は後述）。

本記事では **msconvert（ProteoWizard）** をメインに解説します。事実上の業界標準で、Windows環境であれば最も確実に動作します。

---

## 方法1: msconvert（ProteoWizard）【推奨】

### Step 1: ProteoWizard のインストール

1. [ProteoWizard 公式サイト](https://proteowizard.sourceforge.io/) にアクセス
2. **Download** ページから **Windows版** の最新インストーラをダウンロード
   - 64-bit版（`pwiz-setup-x.x.x.x-x86_64.msi`）を推奨
3. インストーラを実行してインストール

> **ライセンス同意:** インストール時にThermoベンダーライブラリのライセンスへの同意が求められます。

<!-- TODO(Windows): インストール画面のスクリーンショットを追加 -->

### Step 2: 単一ファイルの変換

PowerShell またはコマンドプロンプトで、以下のコマンドを実行します：

```powershell
msconvert CRC06-T.raw --mzML --filter "peakPicking vendor msLevel=1-" -o .
```

主要なオプション：

| オプション | 意味 |
|-----------|------|
| `--mzML` | 出力形式をmzMLに指定 |
| `--filter "peakPicking vendor msLevel=1-"` | ベンダーアルゴリズムでピーク検出（centroid化）— MS1〜MSn全レベル |
| `-o .` | 出力先ディレクトリ（`.` はカレント） |

<!-- TODO(Windows): 実際の実行ログをここに貼る -->

### Step 3: 全RAWファイルの一括変換

ディレクトリ内の全RAWファイルを一括変換するには、ワイルドカードを使います：

```powershell
msconvert *.raw --mzML --filter "peakPicking vendor msLevel=1-" -o mzml
```

`mzml/` サブディレクトリに、すべての変換結果が出力されます。

### Step 4: GUI版を使う場合

コマンドラインが苦手な場合は、**ProteoWizard GUI（MSConvertGUI.exe）** を使うと簡単です：

1. スタートメニューから **MSConvertGUI** を起動
2. **Browse** ボタンでRAWファイルを選択
3. **Output format** を **mzML** に設定
4. **Filters** セクションで **Peak Picking** を追加し、**Vendor** を有効化
5. **Start** ボタンをクリック

<!-- TODO(Windows): MSConvertGUI のスクリーンショットを追加 -->

### Step 5: 変換結果の確認

<!-- TODO(Windows): 変換前後のファイルサイズ、スキャン数などを記載 -->

```
変換前: CRC06-T.raw     (X.X GB)
変換後: CRC06-T.mzML    (X.X GB)
スキャン数: XX,XXX
```

---

## 方法2: ThermoRawFileParser（代替手段）

軽量なCLI特化ツールです。Windows / Linux（x86_64）で動作します。

### Windows の場合

1. [GitHub Releases](https://github.com/compomics/ThermoRawFileParser/releases) から **Windows版** をダウンロード
2. ZIP を解凍
3. PowerShell で実行：

```powershell
ThermoRawFileParser.exe -i CRC06-T.raw -o . -f 2
```

| オプション | 意味 |
|-----------|------|
| `-i` | 入力RAWファイル |
| `-o` | 出力ディレクトリ |
| `-f 2` | 出力形式（**2 = mzML**, 1 = indexed mzML, 0 = MGF） |

### Linux（x86_64）の場合

Docker経由が最も簡単です：

```bash
docker pull quay.io/biocontainers/thermorawfileparser:1.4.5--ha8f3691_0

docker run --rm \
  -v /path/to/data:/data \
  quay.io/biocontainers/thermorawfileparser:1.4.5--ha8f3691_0 \
  ThermoRawFileParser \
  -i /data/CRC06-T.raw \
  -o /data \
  -f 2
```

---

## ⚠️ macOS Apple Silicon での制約

**macOS（特にM1/M2/M3/M4 の Apple Silicon）では、RAW → mzML 変換は現時点で困難**です。以下の制約があります：

| 試行した方法 | 結果 |
|-------------|------|
| ThermoRawFileParser Docker v1.4.x（Mono経由） | Mono JITが x86→ARM64 エミュレーション環境でクラッシュ |
| msconvert Docker（Wine経由） | Wine が Rosetta/QEMU環境で assertion failure |
| ThermoRawFileParser v2.0.0 osx-arm64 ネイティブビルド | Thermo RawFileReader DLL がx86依存のため読み込み失敗 |

**根本原因:** Thermo社の RawFileReader ライブラリが Windows x86_64 向け専用で、ARM64ネイティブ版が存在しないためです。

**Apple Silicon Mac ユーザーへの推奨:**

- **Windows PC または Windows 仮想マシン**（UTM, Parallels, VMware Fusion等）で msconvert を実行する
- **Linux x86_64 サーバー**（クラウドVMなど）で Docker 経由の ThermoRawFileParser を使う

> **シリーズの方針:** 本シリーズではsage-proteomics（mzML入力）を使用するため、mzMLへの変換が必要です。macOSユーザーは上記いずれかの方法で Windows / Linux 環境での変換を行ってください。

---

## まとめ

Thermo RAW ファイルを mzML に変換する方法を解説しました。

| 環境 | 推奨方法 |
|------|---------|
| **Windows** | msconvert（ProteoWizard GUI / CLI） |
| **Linux x86_64** | ThermoRawFileParser（Docker） |
| **macOS Apple Silicon** | Windows環境またはLinuxサーバーでの変換を推奨 |

次の記事では、変換したmzMLファイルを使って、sage-proteomicsでタンパク質を同定・定量します。

> 前回: [#2 データ取得](article-02-data.md)
> 次回: [#4 sageで同定・定量](article-04-sage.md) — sage-proteomicsによるタンパク質同定・定量

#バイオインフォマティクス #プロテオミクス #mzML #msconvert #labcode
