---
title: "RAWファイルをmzMLに変換実行する【論文再現シリーズ #3b】"
emoji: "🔄"
type: "tech"
topics: ["proteomics", "msconvert", "conversion", "labcode"]
published: false
---

# RAWファイルをmzMLに変換実行する

## はじめに

前回（[#3a ProteoWizardセットアップ](article-03a-convert-basics.md)）でProteoWizardの環境設定が完了しました。この記事では、実際のデータファイル **CRC01-N.raw**（大腸がん患者の正常組織）を使って、RAW→mzML変換を実践します。

> **📝 INFO**
>
**この記事で行う処理**
Toyota et al. 2025 論文のThermo RAWファイル（CRC01-N.raw, 1.07GB）をmzML形式に変換します。コマンドライン版とGUI版の2つの方法で変換を行い、単一ファイルと複数ファイル一括変換の手順を習得します。変換されたmzMLファイルは次章のsage解析で使用します。

## 前提

- [#3a ProteoWizardセットアップ](article-03a-convert-basics.md) が完了していること
- **CRC01-N.raw** が手元にあること（1.07GB）
- **対応Notebook**: [`notebooks/step_03.ipynb`](../notebooks/step_03.ipynb) — この記事のコードをセルごとに実行できます

## 今回変換するファイル

```bash
📁 /home/shizuku/labcode/article/Proteomics_drug_marker/Inputs/
└── CRC01-N.raw (1.07 GB) ← 今回変換するファイル
```

**CRC01-N.raw の詳細:**
- **患者**: CRC01（大腸がん患者1番）
- **組織**: N（Normal、正常組織）
- **装置**: Thermo Orbitrap Exploris 480
- **測定法**: DIA-MS（Data-Independent Acquisition）

## コマンドライン版での変換

PATHが通ったので、簡単なコマンドで変換できます：

```cmd
# RAWファイルがあるディレクトリに移動
cd /d D:\proteomics_data\Inputs\

# CRC01-N.rawをmzMLに変換（PATHが通っているため短縮可能）
msconvert CRC01-N.raw --mzML --filter "peakPicking vendor msLevel=1-" -o .
```

**パスが通っていない場合は**、フルパスで実行：
```cmd
"C:\Users\Koki\AppData\Local\Apps\ProteoWizard 3.0.26095.be91649 64-bit\msconvert.exe" CRC01-N.raw --mzML --filter "peakPicking vendor msLevel=1-" -o .
```

### コマンドの詳細解説

```cmd
msconvert CRC01-N.raw --mzML --filter "peakPicking vendor msLevel=1-" -o .
```

| 部分 | 説明 |
|------|------|
| `CRC01-N.raw` | 入力ファイル（1.07GB） |
| `--mzML` | 出力形式をmzMLに指定 |
| `--filter "peakPicking vendor msLevel=1-"` | Thermoの標準アルゴリズムでピーク検出（MS1+MS2レベル） |
| `-o .` | 現在のディレクトリに出力 |

### 実行前の最終確認

変換を実行する前に、必要なファイルがあることを確認：

```cmd
# RAWファイルの存在確認
dir CRC01-N.raw

# 出力先の空き容量確認（約1.5GB必要）
dir
```

### 変換の実行と結果

変換には5-10分程度かかります：

```
Processing: CRC01-N.raw
Writing: CRC01-N.mzML
conversion completed: CRC01-N.mzML (1.2GB)
```

**変換結果:**
```
📁 Inputs/
├── CRC01-N.raw     (1.07 GB) ← 元ファイル
└── CRC01-N.mzML    (1.2 GB)  ← 変換後ファイル
```

## GUI版での変換

**MSConvertGUI**を使った変換手順を詳しく解説します。スタートメニューから **MSConvert** を起動してください。

![MSConvertGUI設定画面](images/MSCOnvertGUI.jpg)

**① ファイル選択（左上赤枠）**：

**単一ファイルの場合**：
1. **File:** フィールドの **Browse** ボタンをクリック
2. CRC01-N.rawファイルがある場所を選択
3. **CRC01-N.raw** を選択して **開く** をクリック
4. **🚨重要🚨 Add ボタンを押すことを忘れずに！**

**複数ファイルの場合**：
1. **Browse** ボタンをクリック
2. **Ctrl** キーを押しながら複数のRAWファイルを選択
   - CRC01-N.raw
   - CRC01-T.raw
   - CRC02-N.raw
   - ...
3. **開く** をクリック
4. **Add** ボタンをクリック
5. 左側のリストに全ファイルが表示されることを確認

> **注意**: ファイルを選択しただけでは変換リストに追加されません。必ず **Add** ボタンをクリックして、左側の変換ファイルリストにすべてのRAWファイルが表示されることを確認してください。

**② 出力形式の設定（左下赤枠）**：
- **Output format:** が **mzML** になっていることを確認
- デフォルトでmzMLが選択されているはずです

**③ Peak Picking設定（中央赤枠）**：
1. **Filters** セクションで **Peak Picking** が選択されていることを確認
2. **Vendor (does not work for UNIFI, and it MUST be the first filter!)** が表示されていることを確認

> **重要**: "Vendor (does not work for UNIFI, and it MUST be the first filter!)" の表示は正常です。これはThermoの標準アルゴリズムでピーク検出を行うという意味です。

**④ 変換実行（右下赤枠）**：
- すべての設定が完了したら **Start** ボタンをクリック

### 変換の実行と確認

**変換中の表示**：
- **Start** ボタンクリック後、変換が開始されます
- **Conversion Progress** ウィンドウが表示されます
- 進行状況がプログレスバーで表示されます
- **CRC01-N.raw (1.07GB)** の変換には **5-10分程度** かかります

![変換進行画面](images/run画面.png)

**変換プロセス**：
1. **Starting...** - 変換準備開始
2. **Opening file** - RAWファイルを読み込み中
3. **Calculating SHA1 checksum...** - ファイル整合性確認
4. **Processing...** - データ処理中
5. **Writing spectra: 2404/51585** - スペクトラム書き込み進行状況
6. **Writing mzML** - mzMLファイル生成中

**変換完了後**：

**単一ファイルの場合**：
- **51,585スペクトラム** の処理が完了
- **CRC01-N.mzML (約1.2GB)** が同じフォルダに生成されます

```
📁 出力結果（単一ファイル）
├── CRC01-N.raw     (1.07 GB) ← 元ファイル
└── CRC01-N.mzML    (1.2 GB)  ← 変換後ファイル
```

## 複数ファイルの一括変換

### コマンドライン版：ワイルドカード使用
```cmd
# 同一フォルダ内の全RAWファイルを一括変換
msconvert *.raw --mzML --filter "peakPicking vendor msLevel=1-" -o .

# 特定パターンのファイルのみ変換（例：CRCで始まるファイル）
msconvert CRC*.raw --mzML --filter "peakPicking vendor msLevel=1-" -o .

# 変換結果確認
dir *.mzML
```

**例：5つのRAWファイルを一括変換**
```cmd
# 変換前の確認
dir *.raw
CRC01-N.raw    (1.07 GB)
CRC01-T.raw    (1.15 GB)
CRC02-N.raw    (1.02 GB)
CRC02-T.raw    (1.08 GB)
CRC03-N.raw    (1.11 GB)

# 一括変換実行
msconvert *.raw --mzML --filter "peakPicking vendor msLevel=1-" -o .

# 変換後の確認
dir *.mzML
CRC01-N.mzML   (1.2 GB)
CRC01-T.mzML   (1.3 GB)
CRC02-N.mzML   (1.1 GB)
CRC02-T.mzML   (1.2 GB)
CRC03-N.mzML   (1.2 GB)
```

### GUI版での一括変換

**複数ファイルの場合**：
- 各ファイルが順次処理されます
- すべての変換が完了するとウィンドウが閉じます

```
📁 出力結果（複数ファイル）
├── CRC01-N.raw     (1.07 GB) ← 元ファイル
├── CRC01-N.mzML    (1.2 GB)  ← 変換後ファイル
├── CRC01-T.raw     (1.15 GB) ← 元ファイル
├── CRC01-T.mzML    (1.3 GB)  ← 変換後ファイル
├── CRC02-N.raw     (1.02 GB) ← 元ファイル
└── CRC02-N.mzML    (1.1 GB)  ← 変換後ファイル
```

> **一括変換の利点**: 複数ファイルの場合、設定を一度行うだけで全ファイルに同じ処理が適用されるため、一貫性のある変換結果が得られます。

> **GUI版のメリット**: コマンドラインが苦手でも直感的に設定でき、設定ミスを防げます。特に **Peak Picking** の **Vendor** 設定が視覚的に確認できるのが利点です。

## DIA-MSデータの特徴確認

変換されたmzMLファイルの中身を確認してみましょう：

```xml
<!-- CRC01-N.mzMLの一部 -->
<mzML>
  <cvList>...</cvList>
  <fileDescription>
    <fileContent>
      <cvParam name="DIA" />  <!-- ← DIA測定であることを確認 -->
    </fileContent>
  </fileDescription>
  <run id="CRC01-N">
    <spectrumList count="XXXX">  <!-- スペクトラム数 -->
```

**確認ポイント:**
- **DIA**: Data-Independent Acquisitionで測定されていること
- **スペクトラム数**: 通常数万〜数十万スペクトラム
- **ファイルサイズ**: RAWとほぼ同等（圧縮効果により若干増減）

> **スペクトラム数の確認**: 進行画面の "writing spectra: 2404/51585" から、このDIA-MSデータには **51,585個のスペクトラム** が含まれていることがわかります。これは典型的なDIA測定のスペクトラム数です。

## まとめ

**CRC01-N.raw** をProteoWizardで **CRC01-N.mzML** に変換しました。

### 変換前後の比較
```
変換前: CRC01-N.raw  (1.07 GB) - Thermo独自形式
変換後: CRC01-N.mzML (1.2 GB)  - 標準オープン形式
```

### 次のステップ

この mzML ファイルは sage-proteomics でタンパク質同定・定量に使用できます。Toyota et al. 2025 の再現解析において、このファイルから**約2,000個のタンパク質**を同定できる予定です。

> 前回: [#3a ProteoWizardセットアップ](article-03a-convert-basics.md)
> 次回: [#4a sageの基礎](article-04a-sage-fundamentals.md) — sage-proteomicsの原理とセットアップ

#バイオインフォマティクス #プロテオミクス #ProteoWizard #msconvert #labcode