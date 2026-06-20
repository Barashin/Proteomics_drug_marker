---
title: "ProteoWizardセットアップとインストール【論文再現シリーズ #3a】"
emoji: "⚙️"
type: "tech"
topics: ["proteomics", "msconvert", "proteowizard", "setup", "labcode"]
published: false
---

# ProteoWizardセットアップとインストール

## はじめに

この記事では、RAW→mzML変換に必要な **ProteoWizard** のインストールと環境設定を行います。現在のProteoWizardは従来とインストール先が変わっているため、最新の手順を詳しく解説します。

> **📝 INFO**
>
**この記事で行う処理**
ProteoWizard（オープンソースの質量分析データ変換ツール）をWindowsにインストールし、コマンドプロンプトから`msconvert`コマンドが使える状態にします。最新バージョンではAppDataフォルダにインストールされるため、正確な場所の特定と環境変数PATHの設定まで行います。

## 前提

- [#2b DIA-MSの理解](article-02b-data-formats.md) が完了していること
- **Windows PC**（推奨）または Windows仮想マシン
- **対応Notebook**: [`notebooks/step_03.ipynb`](../notebooks/step_03.ipynb) — この記事のコードをセルごとに実行できます

> **注意**: macOS Apple Siliconでは制約があります（記事末尾で説明）

## なぜProteoWizardが必要なのか

**Thermo独自のRAWファイル**は以下の問題があります：

1. **独自バイナリ形式**: Thermoでしか読み込めない
2. **解析ツールの制限**: sage、MaxQuant、OpenMSなどは直接読み込み不可
3. **OS依存**: Windows専用ライブラリが必要

**mzMLは国際標準のオープンフォーマット**で：

- どの解析ツールでも読み込み可能
- OS非依存（Windows/Linux/macOS）
- 長期保存に適している（10年後も確実に読める）

## ProteoWizardとは

**ProteoWizard** は、質量分析データの変換・処理を行うオープンソースツール群です。

- **開発**: Vanderbilt大学 + Pacific Northwest National Laboratory
- **ライセンス**: Apache 2.0（商用利用可能）
- **msconvert**: RAW→mzML変換の事実上の標準ツール
- **対応形式**: Thermo、Bruker、SCIEX、Waters、Agilentなど全メーカー

## ProteoWizardのインストール

### Step 1: ダウンロードとインストール

1. [ProteoWizard 公式サイト](https://proteowizard.sourceforge.io/) にアクセス
2. **Download** → **Windows** から最新版をダウンロード
   ```
   pwiz-setup-3.0.XXXX-x86_64.msi (約500MB)
   ```
3. MSIファイルを実行してインストール

**インストール先について**:
- **記載上**: `C:\Program Files\ProteoWizard`
- **実際のデフォルト**: `C:\Users\[ユーザー名]\AppData\Local\Apps\ProteoWizard 3.0.XXXX...`

> **重要**: 最近のバージョンでは、管理者権限を必要としない **AppDataフォルダ** がデフォルトのインストール先になっています。このため、多くの解説サイトで示されている `Program Files` とは異なる場所にインストールされます。

### Step 2: インストール場所の確認と検索方法

**現在のProteoWizardは、デフォルトでAppDataフォルダにインストール**されます：

```cmd
# 最も可能性が高い場所（現在のデフォルト）
dir "C:\Users\%USERNAME%\AppData\Local\Apps\ProteoWizard*"

# 従来の場所（古いバージョンや手動指定した場合）
dir "C:\Program Files\ProteoWizard\"
dir "C:\Program Files (x86)\ProteoWizard\"
```

#### ProteoWizardが見つからない場合の検索方法

**1. コマンドプロンプトでシステム全体検索**
```cmd
# Cドライブ全体からmsconvert.exeを検索（時間がかかります）
dir /s C:\*msconvert.exe*

# より高速：特定のディレクトリから検索
dir /s "C:\Users\%USERNAME%\*msconvert.exe*"
dir /s "C:\Program Files\*msconvert.exe*"
```

**2. PowerShellでの高速検索**
```powershell
# PowerShellを開いて実行（より高速）
Get-ChildItem -Path C:\ -Name "msconvert.exe" -Recurse -ErrorAction SilentlyContinue

# ProteoWizardフォルダを検索
Get-ChildItem -Path C:\ -Name "*ProteoWizard*" -Recurse -ErrorAction SilentlyContinue
```

**3. Windowsの検索機能を使用**
- **Windowsキー** を押してスタートメニューを開く
- **msconvert** または **ProteoWizard** で検索
- 見つかったアプリを右クリック → **ファイルの場所を開く**

**4. エクスプローラーでの検索**
- **エクスプローラー** を開く
- **Cドライブ** を選択
- 検索ボックスに **msconvert.exe** と入力
- 検索結果からパスを確認

**5. whereコマンド（PATHが通っている場合のみ）**
```cmd
# パスが通っていれば場所を表示
where msconvert.exe
```

**典型的なインストール先**:
```
C:\Users\Koki\AppData\Local\Apps\ProteoWizard 3.0.26095.be91649 64-bit\
```

> **検索のコツ**: AppDataフォルダは隠しフォルダのため、エクスプローラーで手動で探す場合は、**表示** → **隠しファイル** を有効にしてください。

> **なぜAppDataなのか**: 管理者権限を不要とし、ユーザー個別のインストールを可能にするため、最近のバージョンではAppDataフォルダがデフォルトになっています。

## 環境変数PATHの設定

### Step 3: 環境変数PATHに追加

毎回長いパスを入力するのを避けるため、PATHに追加します：

#### 方法1: GUIで追加（推奨）

1. **環境変数設定画面を開く**
   ```cmd
   # 直接開く場合
   rundll32 sysdm.cpl,EditEnvironmentVariables
   ```
   または **Windowsキー + R** → `sysdm.cpl` → **詳細設定** → **環境変数**

2. **ユーザー環境変数の編集**
   - **ユーザー環境変数** セクションで **Path** を選択
   - **編集** をクリック
   - **新規** をクリック
   - 以下のパスを追加：
     ```
     C:\Users\Koki\AppData\Local\Apps\ProteoWizard 3.0.26095.be91649 64-bit
     ```
   - **OK** → **OK** → **OK** で保存

> **なぜGUIが推奨なのか**: コマンド（setx）にはPATH長の制限（1024文字）があり、既存のPATHが長い場合に新しいパスが切り捨てられてしまいます。GUIではこの制限がありません。

#### 方法2: PowerShellで追加（コマンド派向け）
```powershell
# PowerShellを管理者として実行
$oldPath = [Environment]::GetEnvironmentVariable("Path", "User")
$newPath = $oldPath + ";C:\Users\Koki\AppData\Local\Apps\ProteoWizard 3.0.26095.be91649 64-bit"
[Environment]::SetEnvironmentVariable("Path", $newPath, "User")
```

> **注意**: `setx PATH "%PATH%;..."`は文字数制限で失敗することが多いため推奨しません。

### Step 4: 動作確認

**重要**: コマンドプロンプトを再起動してから実行

```cmd
# 新しいコマンドプロンプトを開いて確認
msconvert --help
```

**成功時の出力例**：
```
Usage: msconvert [options] [filemasks]
Convert mass spec data file formats.

Return value: # of failed files.

Options:
  -f [ --filelist ] arg              : specify text file containing filenames
  -o [ --outdir ] arg (=.)           : set output directory ('-' for stdout, '.' for current directory) [.]
  -c [ --config ] arg                : configuration file (optionName=value)
  --outfile arg                      : Override the name of output file.
  -e [ --ext ] arg                   : set extension for output files [mzML|mzXML|mgf|txt|mz5|mzMLb]
  --mzML                             : write mzML format [default]
  --mzXML                            : write mzXML format
  --mz5                              : write mz5 format
  --mzMLb                            : write mzMLb format
  --mzMLbChunkSize arg (=1048576)    : mzMLb dataset chunk size in bytes
  --mzMLbCompressionLevel arg (=4)   : mzMLb GZIP compression level (0-9)
  --mgf                              : write Mascot generic format
  --text                             : write ProteoWizard internal text format
  --ms1                              : write MS1 format
  --cms1                             : write CMS1 format
  --ms2                              : write MS2 format
  --cms2                             : write CMS2 format
  ...
```

✅ **このようなヘルプ画面が表示されれば、ProteoWizardの設定は完了です！**

**確認ポイント**：
- `msconvert` コマンドが認識されている
- 各種出力フォーマット（mzML、mzXML、mgf等）が利用可能
- `--mzML` がデフォルト形式として設定されている

**エラーが出る場合**：
```cmd
'msconvert' は、内部コマンドまたは外部コマンド、
操作可能なプログラムまたはバッチ ファイルとして認識されていません。
```
→ PATH設定を再確認し、コマンドプロンプトを再起動してください

## macOS Apple Siliconでの制約

- **ProteoWizard** はWindows/Linux向けでmacOSネイティブ版はありません
- **Wine経由** での実行は可能ですが、Thermoライブラリの制約で不安定
- **仮想マシン** （Parallels等）でWindows環境を構築することを推奨します

## まとめ

ProteoWizardのインストールと環境設定が完了しました。最新バージョンのAppDataフォルダへのインストールと正確なPATH設定により、`msconvert`コマンドが使用可能になりました。

実際のRAWファイル変換手順については次章で詳しく解説します。

> 前回: [#2b DIA-MSの理解](article-02b-data-formats.md)
> 次回: [#3b RAW→mzML変換実行](article-03b-convert-execution.md) — 実際のファイル変換とトラブルシューティング

#バイオインフォマティクス #プロテオミクス #ProteoWizard #環境設定 #labcode