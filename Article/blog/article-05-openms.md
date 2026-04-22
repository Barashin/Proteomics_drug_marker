---
title: "OpenMS + AlphaPeptDeepで深層学習DIA解析する【論文再現シリーズ #5】"
emoji: "🧠"
type: "tech"
topics: ["proteomics", "openms", "alphapeptdeep", "bioinformatics", "labcode"]
published: false
---

# OpenMS + AlphaPeptDeepで深層学習DIA解析する

## はじめに

この記事では、#4 で紹介した sage-proteomics の**代替パイプライン**として、**OpenMS + AlphaPeptDeep** による DIA 解析を解説します。sage が理論スペクトルベース（in silico 消化）でペプチドを同定するのに対し、AlphaPeptDeep は**深層学習モデル**で保持時間・フラグメント強度を予測したスペクトルライブラリを生成し、OpenSWATH がそのライブラリを使って DIA データを検索します。

深層学習ベースのライブラリ予測は DIA-NN が先駆けですが、DIA-NN は商用利用に有料ライセンスが必要です。AlphaPeptDeep（Apache 2.0）と OpenMS（BSD）の組み合わせなら、**完全無料・商用利用可能**で同等のアプローチが取れます。

:::message
**この記事で行う処理**
mzMLファイル（32サンプル）とヒトプロテオームFASTAを入力として、AlphaPeptDeepで深層学習による予測スペクトルライブラリを生成し、OpenSWATHでDIAデータを検索し、PyProphetでFDR推定とスコアリングを行い、最終的にタンパク質×サンプルの定量マトリクス（protein_matrix_from_openms.csv）を出力します。sageの出力と同じ形式なので、以降の前処理・可視化ステップにそのまま接続できます。
:::

## 処理の全体像

この記事で行う処理を、「入力データ」「処理の中身」「出力データ」の3段階に分けて説明します。

### 入力：mzMLファイルとFASTA

sage と同じ入力データを使います。

- **mzMLファイル**: 32サンプル（16患者 × Normal/Tumor）のDIA-MS生データ
- **FASTAファイル**: ヒト全タンパク質のアミノ酸配列（UniProt reviewed、約20,000エントリ）

### 処理：4ステップのパイプライン

```
FASTA（ヒト全タンパク質配列）
  → [Step 1] AlphaPeptDeep: 深層学習でスペクトルライブラリ予測
      → 各ペプチドの保持時間・フラグメント強度を予測
  → [Step 2] OpenSWATH: 予測ライブラリ vs DIAデータの照合
      → 各サンプルでペプチドを同定・定量
  → [Step 3] PyProphet: FDR推定とスコアリング
      → 統計的に信頼できる同定結果だけを選別
  → [Step 4] タンパク質マトリクス構築
      → ペプチドレベル → タンパク質レベルに集約
```

**sage との最大の違い**はStep 1です。sage は理論的なフラグメントイオンのm/z値だけで照合しますが、AlphaPeptDeep は深層学習で「このペプチドはRT何分に溶出し、どのフラグメントがどの相対強度で出現するか」まで予測します。この追加情報により、同定の精度と感度が向上する可能性があります。

### 出力：タンパク質×サンプルの定量マトリクス

sage と同じ形式の CSV ファイルが得られます。

```
               CRC01-N    CRC01-T    CRC02-N    CRC02-T   ...  CRC16-T
タンパク質A     18.5       22.3       19.1       21.8      ...  23.0
タンパク質B     15.2       14.8       15.9       15.1      ...  14.5
...
```

---

## 前提

- [#1 環境構築](article-01-setup.md) が完了していること
- [#2 データ取得](article-02-data.md) で mzML ファイルが手元にあること
- [#4 sage](article-04-sage.md) で FASTA ファイルが取得済みであること
- **対応Notebook**: [`notebooks/step_05.ipynb`](../notebooks/step_05.ipynb) — この記事のコードをセルごとに実行できます

:::message alert
**sage と OpenMS + AlphaPeptDeep は択一です**。どちらか一方のパイプラインで生成した `protein_matrix_from_*.csv` を #6 の前処理に渡せば、以降のステップは共通です。両方試して結果を比較することもできます。
:::

## OpenMS + AlphaPeptDeep とは

### ツール比較

| ツール | ライセンス | 商用利用 | 本書採用 |
|-------|-----------|---------|:---:|
| DIA-NN | 学術無料 / 商用有料 | ❌ | |
| FragPipe (MSFragger) | 学術無料 / 商用有料 | ❌ | |
| Spectronaut | 商用有料 | ⚠ | |
| sage-proteomics | MIT | ✅ | ✓ (#4) |
| **OpenMS + AlphaPeptDeep** | **BSD + Apache 2.0** | **✅** | **✓ (#5)** |
| MaxQuant / MaxDIA | 非商用無料 | ⚠ | |

### 各ツールの役割

| ツール | 役割 | ライセンス | 論文 |
|-------|------|-----------|------|
| **AlphaPeptDeep** | 深層学習でスペクトルライブラリを予測 | Apache 2.0 | Wen et al. *Nat. Commun.* 2024 |
| **OpenMS (OpenSWATH)** | 予測ライブラリを使ったDIA検索 | BSD-3-Clause | Rost et al. *Nat. Biotechnol.* 2014 |
| **PyProphet** | 統計的スコアリングとFDR推定 | BSD-3-Clause | Teleman et al. *Bioinformatics* 2015 |

:::message
**【深層学習ライブラリ予測とは】**

従来のDIA解析では、事前にDDA実験で取得したスペクトルライブラリが必要でした（= library-based）。sage のような library-free アプローチは理論m/z値だけで照合するため、保持時間やフラグメント強度パターンの情報を使いません。

AlphaPeptDeep は、大規模なDDA実験データで訓練した深層学習モデル（Transformer）を使い、FASTAから直接「仮想的なDDAライブラリ」を生成します。各ペプチドについて以下を予測します：

- **保持時間（iRT）**: そのペプチドがLC-MSの何分に溶出するか
- **フラグメント強度パターン**: 各b/yイオンの相対的な強度

この予測ライブラリを使うことで、実際のDDAライブラリを取得せずに、library-based と同等の精度でDIA解析ができます。DIA-NN のコアアイデアと同じアプローチを、オープンソースで実現したものです。
:::

## インストール

```bash
# crc-proteomics 環境をアクティベート
micromamba activate crc-proteomics

# OpenMS のインストール（OpenSWATHWorkflow を含む）
micromamba install -c bioconda -c conda-forge openms -y

# AlphaPeptDeep のインストール（peptdeep CLI を含む）
pip install alphapeptdeep

# PyProphet のインストール（FDR推定ツール）
pip install pyprophet

# バージョン確認
OpenSwathWorkflow --help 2>&1 | head -1
peptdeep --help 2>&1 | head -1
pyprophet --version
```

:::message alert
**インストールの注意点**: AlphaPeptDeep は PyTorch に依存しています。GPU がある環境では CUDA 版 PyTorch が自動で使われ、ライブラリ予測が大幅に高速化されます。CPU のみでも動作しますが、Step 1 の処理時間が長くなります（CPU: 約30分、GPU: 約5分）。
:::

## コード全文（対応Notebook: step_05.ipynb）

### ライブラリと設定

```python
# 標準ライブラリ: ファイル操作(os)、パターン検索(glob)、正規表現(re)、外部コマンド実行(subprocess)、システム操作(sys)、時間計測(time)
import os, glob, re, subprocess, sys, time
import pandas as pd   # データフレーム操作ライブラリ（表形式データの読み書き・加工に使う）
import numpy as np    # 数値計算ライブラリ（配列演算・統計に使用）
```

```python
# --- パス設定（プロジェクト内の各ディレクトリ・ファイルへの相対パスを定数として定義） ---
MZML_DIR    = "../data/raw/raw_mzML"                # mzML生データが格納されているディレクトリ
FASTA_PATH  = "../data/raw/human_proteome.fasta"     # ヒトプロテオームのFASTA配列ファイルのパス
RESULTS_DIR = "../results"                           # 解析結果の出力先ルートディレクトリ

# OpenMS + AlphaPeptDeep パイプライン専用の出力ディレクトリ
OPENMS_OUT  = os.path.join(RESULTS_DIR, "openms_output")  # OpenMS関連の全出力をまとめるディレクトリ
LIBRARY_DIR = os.path.join(OPENMS_OUT, "library")          # AlphaPeptDeepが生成するスペクトルライブラリの保存先
OSWR_DIR    = os.path.join(OPENMS_OUT, "openswath")        # OpenSWATHの検索結果（.osw）の保存先
PYPROPHET_DIR = os.path.join(OPENMS_OUT, "pyprophet")      # PyProphetのスコアリング結果の保存先

# 各ディレクトリが存在しなければ作成する（exist_ok=Trueで既存でもエラーにならない）
for d in [OPENMS_OUT, LIBRARY_DIR, OSWR_DIR, PYPROPHET_DIR]:
    os.makedirs(d, exist_ok=True)

# mzMLファイルの一覧を取得してファイル名順にソート
mzml_files = sorted(glob.glob(os.path.join(MZML_DIR, "*.mzML")))
# 見つかったファイル数を表示（期待値: 32）
print(f"{len(mzml_files)} mzML files found")
```

```python
# --- パイプライン共通パラメータ ---
# ペプチド FDR（偽発見率）の閾値（論文と同じ 1% = 0.01）
Q_THRESHOLD = 0.01
# トリプシン消化のミスクリーベージ許容回数（論文と同じ: 最大1回の切断漏れを許容）
MISSED_CLEAVAGES = 1
# ペプチド長の下限と上限（論文のDIA-NNパラメータに合わせる）
PEPTIDE_MIN_LEN = 7
PEPTIDE_MAX_LEN = 45
# 前駆体イオンの電荷状態の範囲（2価〜4価のみ対象）
PRECURSOR_CHARGE_MIN = 2
PRECURSOR_CHARGE_MAX = 4
# フラグメントイオンのm/z範囲（質量分析計の測定範囲に合わせる）
FRAGMENT_MIN_MZ = 200.0
FRAGMENT_MAX_MZ = 1800.0
# 質量精度（ppm単位）: プリカーサーとフラグメントの質量誤差許容範囲
MS1_PPM = 10.0
MS2_PPM = 10.0
```

### Step 1: AlphaPeptDeepでスペクトルライブラリ生成

AlphaPeptDeep の `peptdeep` CLI を使い、FASTA からペプチド配列を生成し、深層学習モデルで保持時間（iRT）・フラグメント強度を予測します。出力は OpenSWATH が読み込める TSV 形式のスペクトルライブラリです。

```python
# --- Step 1: AlphaPeptDeepによる予測スペクトルライブラリの生成 ---

# AlphaPeptDeepの設定をYAML形式で定義する
# peptdeep CLIは設定ファイル経由でパラメータを受け取る
PEPTDEEP_SETTINGS = os.path.join(OPENMS_OUT, "peptdeep_settings.yaml")

# 設定ファイルの内容を文字列として構築（YAMLフォーマット）
settings_content = f"""
# AlphaPeptDeep ライブラリ予測設定ファイル
model_mgr:
  # 使用するモデル: 汎用事前学習モデル（generic）
  external_ms2_model: ''
  external_rt_model: ''
  external_ccs_model: ''

library:
  # FASTAファイルのパス: ヒトプロテオーム配列
  fasta_files:
    - {os.path.abspath(FASTA_PATH)}
  # 出力ファイルのパス: OpenSWATH用のTSV形式
  output_tsv: {os.path.abspath(os.path.join(LIBRARY_DIR, "predicted_library.tsv"))}
  # 消化酵素の設定（トリプシン: K/Rの後で切断、Pの前は切断しない）
  enzyme: trypsin
  # ミスクリーベージの最大回数
  max_missed_cleavages: {MISSED_CLEAVAGES}
  # ペプチド長の範囲
  min_peptide_length: {PEPTIDE_MIN_LEN}
  max_peptide_length: {PEPTIDE_MAX_LEN}
  # 前駆体電荷の範囲
  min_precursor_charge: {PRECURSOR_CHARGE_MIN}
  max_precursor_charge: {PRECURSOR_CHARGE_MAX}
  # フラグメントイオンのm/z範囲
  min_fragment_mz: {FRAGMENT_MIN_MZ}
  max_fragment_mz: {FRAGMENT_MAX_MZ}
  # 固定修飾: システインのカルバミドメチル化（IAA処理による標準修飾）
  fix_modifications:
    - Carbamidomethyl@C
  # 可変修飾: なし（シンプルな設定で実行）
  var_modifications: []
  max_var_mod_num: 0
  # フラグメントイオンの種類: b/yイオン
  fragment_types:
    - b
    - y
  # フラグメントの最大電荷
  max_fragment_charge: 2
"""

# 設定ファイルをディスクに書き出す
with open(PEPTDEEP_SETTINGS, "w") as f:
    f.write(settings_content)
print(f"AlphaPeptDeep設定ファイルを保存: {PEPTDEEP_SETTINGS}")
```

```python
# peptdeep CLIコマンドを構築: FASTAからスペクトルライブラリを予測する
cmd_peptdeep = [
    "peptdeep",         # AlphaPeptDeepのCLIコマンド
    "library",          # ライブラリ予測モード
    "--settings",       # 設定ファイルを指定するフラグ
    PEPTDEEP_SETTINGS   # 設定ファイルのパス
]

print("Running AlphaPeptDeep library prediction...")
print(f"  FASTA: {FASTA_PATH}")
print(f"  Output: {LIBRARY_DIR}")
t0 = time.time()  # 実行開始時刻を記録（所要時間の計算用）

# peptdeep を子プロセスとして起動し、ログをリアルタイム出力する
# stdout=PIPEで標準出力をキャプチャ、stderr=STDOUTでエラー出力も統合、text=Trueで文字列として扱う
proc = subprocess.Popen(cmd_peptdeep, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
# peptdeepの出力を1行ずつ読み取り、Notebook上にリアルタイム表示するループ
for line in proc.stdout:
    sys.stdout.write(line)  # 各ログ行をそのまま標準出力に書き出す
proc.wait()  # プロセスの終了を待機（終了コードが確定する）

# 所要時間を分単位で表示し、終了コードも出力（0なら正常終了）
elapsed = (time.time() - t0) / 60
print(f"\nAlphaPeptDeep完了: {elapsed:.1f} min (exit code {proc.returncode})")
```

```python
# 生成されたスペクトルライブラリの確認
PREDICTED_LIB = os.path.join(LIBRARY_DIR, "predicted_library.tsv")
# ライブラリファイルが存在するか確認する
if os.path.exists(PREDICTED_LIB):
    # ライブラリをDataFrameとして読み込み（大きいファイルなのでnrows=5で先頭5行だけ確認）
    df_lib = pd.read_csv(PREDICTED_LIB, sep="\t", nrows=5)
    # ライブラリのファイルサイズをMB単位で表示する
    lib_size = os.path.getsize(PREDICTED_LIB) / 1024**2
    print(f"ライブラリ生成成功: {lib_size:.1f} MB")
    # ライブラリの列名を表示して、正しいフォーマットで生成されたか確認する
    print(f"列: {list(df_lib.columns)}")
    # 先頭5行を表示する
    df_lib.head()
else:
    # ファイルが存在しない場合はエラーメッセージを表示する
    print(f"ERROR: ライブラリファイルが見つかりません: {PREDICTED_LIB}")
```

### Step 2: OpenSWATHでDIA解析

OpenSWATH（OpenSwathWorkflow）は、予測スペクトルライブラリと DIA データを照合し、各サンプルでペプチドを同定・定量します。出力は SQLite 形式の `.osw` ファイルです。

:::message
**【OpenSWATH の動作原理】**

OpenSWATH は以下の手順で DIA データを検索します：

1. **抽出イオンクロマトグラム（XIC）の構築**: ライブラリの各ペプチドについて、理論m/z ± 質量精度の範囲で DIA データからシグナルを抽出
2. **ピーク検出**: 保持時間（iRT）付近でピークグループを検出
3. **スコアリング**: ピークグループの品質を複数の指標（RT一致度、フラグメント強度相関、同位体パターンなど）でスコアリング
4. **出力**: スコア付きのピークグループを `.osw`（SQLite）ファイルに保存

AlphaPeptDeep が予測した保持時間とフラグメント強度パターンがあるため、ピーク検出とスコアリングの精度が向上します。
:::

```python
# --- Step 2: OpenSWATH による DIA 検索 ---

# ライブラリをOpenSWATH形式（.pqp）に変換する
# OpenSWATHはTSVライブラリも読めるが、PQP（SQLite）形式のほうが高速
PQP_LIB = os.path.join(LIBRARY_DIR, "predicted_library.pqp")

# OpenMS の TargetedFileConverter でTSV→PQP変換を行う
cmd_convert = [
    "TargetedFileConverter",  # OpenMSのファイル変換ツール
    "-in", PREDICTED_LIB,    # 入力: AlphaPeptDeepが生成したTSVライブラリ
    "-out", PQP_LIB          # 出力: OpenSWATH用のPQPライブラリ（SQLite形式）
]

print("TSVライブラリをPQP形式に変換中...")
t0 = time.time()  # 変換開始時刻を記録

# TargetedFileConverterを子プロセスとして実行する
proc = subprocess.Popen(cmd_convert, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
# 出力を1行ずつ表示する
for line in proc.stdout:
    sys.stdout.write(line)
proc.wait()  # プロセスの終了を待機する

# 変換の所要時間と終了コードを表示する
print(f"PQP変換完了: {(time.time() - t0):.1f}s (exit code {proc.returncode})")
```

```python
# 各mzMLファイルに対してOpenSWATHを実行する関数
def run_openswath(mzml_path, pqp_lib, output_dir, ms1_ppm, ms2_ppm):
    """1つのmzMLファイルに対してOpenSwathWorkflowを実行する。

    Args:
        mzml_path: 入力mzMLファイルのパス
        pqp_lib: PQP形式のスペクトルライブラリのパス
        output_dir: 出力ディレクトリ（.oswファイルの保存先）
        ms1_ppm: MS1の質量精度（ppm）
        ms2_ppm: MS2の質量精度（ppm）

    Returns:
        出力.oswファイルのパス
    """
    # 入力ファイル名から拡張子を除去してサンプル名を取得する（例: CRC01-N）
    sample_name = os.path.splitext(os.path.basename(mzml_path))[0]
    # 出力.oswファイルのパスを構築する（各サンプルごとに個別のファイル）
    osw_out = os.path.join(output_dir, f"{sample_name}.osw")

    # OpenSwathWorkflowのコマンドライン引数を構築する
    cmd = [
        "OpenSwathWorkflow",          # OpenSWATHのメインコマンド
        "-in", mzml_path,             # 入力: DIAデータのmzMLファイル
        "-tr", pqp_lib,               # ライブラリ: PQP形式の予測スペクトルライブラリ
        "-out_osw", osw_out,          # 出力: OSW形式（SQLite）の結果ファイル
        "-min_upper_edge_dist", "1",  # DIA窓境界からの最小距離（品質フィルタ）
        "-mz_extraction_window", str(ms2_ppm),       # MS2のm/z抽出窓（ppm単位）
        "-mz_extraction_window_unit", "ppm",          # m/z抽出窓の単位をppmに指定
        "-mz_extraction_window_ms1", str(ms1_ppm),   # MS1のm/z抽出窓（ppm単位）
        "-mz_extraction_window_ms1_unit", "ppm",      # MS1のm/z抽出窓の単位をppmに指定
        "-use_ms1_traces",            # MS1レベルのXIC（抽出イオンクロマトグラム）も使用する
        "-Scoring:stop_report_after_feature", "5",    # 各遷移グループで報告する特徴量の最大数
        "-Scoring:TransitionGroupPicker:min_peak_width", "10",  # ピークの最小幅（秒）
        "-threads", "4"               # 並列処理のスレッド数（環境に応じて調整可能）
    ]

    # OpenSwathWorkflowを子プロセスとして実行する
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    # 出力を1行ずつリアルタイム表示する
    for line in proc.stdout:
        sys.stdout.write(line)
    proc.wait()  # プロセスの終了を待機する

    return osw_out  # 出力ファイルのパスを返す
```

```python
# 全32サンプルに対してOpenSWATHを順次実行する
print(f"OpenSWATH実行開始: {len(mzml_files)} ファイル")
t0_all = time.time()  # 全体の開始時刻を記録する

# 各サンプルの出力ファイルパスを格納するリスト
osw_files = []

# 全mzMLファイルに対してOpenSWATHを順番に実行するループ
for i, mzml_path in enumerate(mzml_files):
    # 現在の処理ファイル番号とファイル名を表示する
    sample_name = os.path.basename(mzml_path)
    print(f"\n[{i+1}/{len(mzml_files)}] Processing: {sample_name}")
    t0 = time.time()  # 個別ファイルの開始時刻を記録する

    # OpenSWATHを実行し、出力ファイルパスを取得する
    osw_path = run_openswath(mzml_path, PQP_LIB, OSWR_DIR, MS1_PPM, MS2_PPM)
    # 出力ファイルパスをリストに追加する
    osw_files.append(osw_path)

    # 個別ファイルの処理時間を分単位で表示する
    elapsed = (time.time() - t0) / 60
    print(f"  完了: {elapsed:.1f} min")

# 全体の所要時間を表示する
total_elapsed = (time.time() - t0_all) / 60
print(f"\nOpenSWATH全完了: {total_elapsed:.1f} min ({len(osw_files)} ファイル)")
```

### Step 3: PyProphetでFDR推定

PyProphet は OpenSWATH の出力（`.osw` ファイル）に対して統計的スコアリングを行い、**偽発見率（FDR）** を推定します。これにより、信頼性の高い同定結果だけを抽出できます。

:::message
**【PyProphet のスコアリング】**

OpenSWATH が出力するスコアは「ピークグループの品質を示す生スコア」であり、統計的なFDRではありません。PyProphet は以下の手順で統計的な信頼度を付与します：

1. **semi-supervised learning**: target（本物のペプチド候補）と decoy（偽のペプチド候補）のスコア分布を学習
2. **d-score計算**: 複数のサブスコアを1つの統合スコア（d-score）に変換
3. **FDR推定**: target/decoy分析に基づき、各同定のq-valueを計算
4. **階層的FDR**: MS2レベル → ペプチドレベル → タンパク質レベルの順にFDRを制御
:::

```python
# --- Step 3: PyProphet によるスコアリングとFDR推定 ---

# 全サンプルの.oswファイルを1つに統合する（PyProphetはマルチファイル入力に対応）
# 統合することで、サンプル間の情報を利用したglobal FDR推定が可能になる
MERGED_OSW = os.path.join(PYPROPHET_DIR, "merged.osw")

# pyprophet merge コマンドで全.oswファイルを1つに統合する
cmd_merge = [
    "pyprophet", "merge",       # PyProphetのマージモード
    "--out", MERGED_OSW,        # 統合後のファイルパス
] + osw_files                   # 入力: 全32サンプルの.oswファイル

print("PyProphet: .oswファイルを統合中...")
t0 = time.time()  # 統合開始時刻を記録する

# pyprophet merge を子プロセスとして実行する
proc = subprocess.Popen(cmd_merge, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
# 出力を1行ずつ表示する
for line in proc.stdout:
    sys.stdout.write(line)
proc.wait()  # 統合の完了を待機する

print(f"統合完了: {(time.time() - t0):.1f}s (exit code {proc.returncode})")
```

```python
# MS2レベルのスコアリング: 各ピークグループにd-scoreとq-valueを付与する
cmd_score = [
    "pyprophet", "score",        # PyProphetのスコアリングモード
    "--in", MERGED_OSW,          # 入力: 統合された.oswファイル
    "--level", "ms2",            # スコアリングレベル: MS2（フラグメントイオンレベル）
    "--ss_initial_fdr", "0.15",  # semi-supervised learningの初期FDR閾値
    "--ss_iteration_fdr", "0.05" # 繰り返し学習時のFDR閾値
]

print("PyProphet: MS2レベルスコアリング中...")
t0 = time.time()  # スコアリング開始時刻を記録する

# pyprophet score を子プロセスとして実行する
proc = subprocess.Popen(cmd_score, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
# 出力を1行ずつ表示する
for line in proc.stdout:
    sys.stdout.write(line)
proc.wait()  # スコアリングの完了を待機する

print(f"MS2スコアリング完了: {(time.time() - t0) / 60:.1f} min (exit code {proc.returncode})")
```

```python
# ペプチドレベルのFDR推定: MS2スコアをペプチドレベルに統合する
cmd_peptide = [
    "pyprophet", "peptide",      # ペプチドレベルのFDR推定モード
    "--in", MERGED_OSW,          # 入力: スコアリング済みの.oswファイル
    "--context", "global"        # FDRのコンテキスト: 全サンプルを通じた global FDR
]

print("PyProphet: ペプチドレベルFDR推定中...")
t0 = time.time()  # 開始時刻を記録する

# pyprophet peptide を子プロセスとして実行する
proc = subprocess.Popen(cmd_peptide, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
for line in proc.stdout:
    sys.stdout.write(line)
proc.wait()  # 完了を待機する

print(f"ペプチドFDR完了: {(time.time() - t0):.1f}s (exit code {proc.returncode})")
```

```python
# タンパク質レベルのFDR推定: ペプチドスコアをタンパク質レベルに統合する
cmd_protein = [
    "pyprophet", "protein",      # タンパク質レベルのFDR推定モード
    "--in", MERGED_OSW,          # 入力: ペプチドFDR推定済みの.oswファイル
    "--context", "global"        # FDRのコンテキスト: global FDR
]

print("PyProphet: タンパク質レベルFDR推定中...")
t0 = time.time()  # 開始時刻を記録する

# pyprophet protein を子プロセスとして実行する
proc = subprocess.Popen(cmd_protein, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
for line in proc.stdout:
    sys.stdout.write(line)
proc.wait()  # 完了を待機する

print(f"タンパク質FDR完了: {(time.time() - t0):.1f}s (exit code {proc.returncode})")
```

```python
# スコアリング結果をTSV形式でエクスポートする（後段のPython処理で読み込むため）
EXPORT_TSV = os.path.join(PYPROPHET_DIR, "pyprophet_export.tsv")

cmd_export = [
    "pyprophet", "export",         # エクスポートモード
    "--in", MERGED_OSW,            # 入力: 全レベルのFDR推定が完了した.oswファイル
    "--out", EXPORT_TSV,           # 出力: TSV形式のエクスポートファイル
    "--max_global_peptide_qvalue", str(Q_THRESHOLD),   # ペプチドFDRフィルタ（1%以下）
    "--max_global_protein_qvalue", str(Q_THRESHOLD)    # タンパク質FDRフィルタ（1%以下）
]

print("PyProphet: 結果をTSVにエクスポート中...")
t0 = time.time()  # 開始時刻を記録する

# pyprophet export を子プロセスとして実行する
proc = subprocess.Popen(cmd_export, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
for line in proc.stdout:
    sys.stdout.write(line)
proc.wait()  # 完了を待機する

# エクスポートファイルのサイズと所要時間を表示する
if os.path.exists(EXPORT_TSV):
    export_size = os.path.getsize(EXPORT_TSV) / 1024**2
    print(f"エクスポート完了: {export_size:.1f} MB, {(time.time() - t0):.1f}s")
else:
    print("ERROR: エクスポートファイルが生成されませんでした")
```

### Step 4: タンパク質マトリクス構築

PyProphet がエクスポートした TSV ファイルから、sage と同じ形式の「タンパク質 × サンプル」定量マトリクスを構築します。

```python
# --- Step 4: タンパク質マトリクスの構築 ---

# PyProphetのエクスポートTSVを読み込む
df = pd.read_csv(EXPORT_TSV, sep="\t")
print(f"PyProphet export: {len(df)} rows, {len(df.columns)} columns")
# 列名を表示して、必要なカラムが含まれているか確認する
print(f"列名: {list(df.columns[:10])} ...")
```

```python
# PyProphetのエクスポートには以下の主要列が含まれる:
#   - ProteinName: タンパク質名（UniProtアクセッション）
#   - PeptideSequence: ペプチド配列（修飾含む）
#   - filename: 元のmzMLファイル名
#   - Intensity: ペプチドの定量強度
#   - m_score: PyProphetのq-value（FDR）

# filenameカラムからサンプル名を抽出する（ディレクトリパスと拡張子を除去）
df["sample"] = df["filename"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
# サンプル名の一覧を表示して正しく抽出できているか確認する
print(f"サンプル数: {df['sample'].nunique()}")
print(f"サンプル一覧: {sorted(df['sample'].unique())[:5]} ...")
```

```python
# FASTAからGene Symbolマッピング辞書を構築する（sage記事と同じ関数）
def parse_fasta_gene_map(fasta_path):
    """UniProt FASTA ヘッダから accession -> Gene Symbol の辞書を作る。

    ヘッダ例: >sp|P04637|P53_HUMAN ... GN=TP53 ...
    GN= がなければ Entry Name の先頭語 (例: P53) を使う。
    """
    id_to_gene = {}  # アクセッション番号→Gene Symbolの辞書を初期化する
    # 「GN=」の後に続く空白以外の文字列（=Gene Symbol）を抽出する正規表現パターン
    gn_re = re.compile(r"\bGN=(\S+)")
    # FASTAファイルを開いて1行ずつ読み込む
    with open(fasta_path) as f:
        for line in f:  # FASTAファイルの全行を走査する
            if not line.startswith(">"):
                continue  # アミノ酸配列行はスキップ（ヘッダ行のみ処理する）
            # ヘッダ行を「|」で最大3つに分割する（例: "sp", "P04637", "P53_HUMAN ..."）
            parts = line[1:].split("|", 2)
            if len(parts) >= 3:
                acc = parts[1]                     # アクセッション番号（例: P04637）
                entry_name = parts[2].split()[0]   # エントリ名（例: P53_HUMAN）
            else:
                acc = line[1:].split()[0]           # 非標準ヘッダの場合のフォールバック
                entry_name = acc
            m = gn_re.search(line)  # 「GN=遺伝子名」パターンを検索する
            # GN=が見つかればその遺伝子名、なければエントリ名の「_」前を使用する
            id_to_gene[acc] = m.group(1) if m else entry_name.split("_")[0]
    return id_to_gene  # {アクセッション番号: Gene Symbol} の辞書を返す

# FASTAを解析してGene Symbolマッピングを構築する
id_to_gene = parse_fasta_gene_map(FASTA_PATH)
print(f"FASTA entries: {len(id_to_gene)}")
```

```python
# ProteinNameからアクセッション番号を抽出し、Gene Symbolに変換する
# PyProphetの出力では ProteinName が "1/sp|P04637|P53_HUMAN" のような形式になることがある
df["acc"] = df["ProteinName"].str.extract(r"sp\|(\w+)\|", expand=False)
# 抽出できなかった場合はProteinNameをそのまま使う
df["acc"] = df["acc"].fillna(df["ProteinName"])
# アクセッション番号をGene Symbolに変換する（マッピングできない場合はアクセッション番号を維持）
df["gene"] = df["acc"].map(id_to_gene).fillna(df["acc"])

# 変換結果を確認する
print(f"Gene Symbol変換: {df['gene'].nunique()} ユニーク遺伝子")
```

```python
# ペプチドレベルの強度をタンパク質レベルに集約する
# 各サンプル×各タンパク質について、ペプチド強度の合計を計算する

# ピボットテーブルを構築: 行=gene（タンパク質）、列=sample、値=Intensityの合計
matrix = df.pivot_table(
    index="gene",          # 行: Gene Symbol（タンパク質名）
    columns="sample",      # 列: サンプル名（CRC01-N, CRC01-Tなど）
    values="Intensity",    # 値: ペプチドの定量強度
    aggfunc="sum"          # 集約方法: 同一タンパク質の全ペプチド強度を合計
)

# 値が0の箇所をNaN（未検出）に置換する
matrix = matrix.replace(0, np.nan)
# 全サンプルでNaN（全く検出されなかった）のタンパク質行を削除する
matrix = matrix.dropna(how="all")
# インデックス名を設定する（CSV出力時のヘッダになる）
matrix.index.name = "Protein"

# マトリクスのサイズを表示する（行数=タンパク質数、列数=サンプル数）
print(f"タンパク質マトリクス: {matrix.shape[0]} proteins x {matrix.shape[1]} samples")
# 先頭5行を表示して中身を確認する
matrix.head()
```

```python
# タンパク質マトリクスをCSVファイルとして保存する
out_csv = os.path.join(RESULTS_DIR, "protein_matrix_from_openms.csv")
# DataFrameをCSV形式で保存する（インデックス=タンパク質名も含める）
matrix.to_csv(out_csv)
# 保存先パスとマトリクスサイズを表示して確認する
print(f"保存完了: {out_csv}")
print(f"  {matrix.shape[0]} タンパク質 × {matrix.shape[1]} サンプル")
```

## sage との比較

| 指標 | sage (#4) | OpenMS + AlphaPeptDeep (#5) |
|------|-----------|---------------------------|
| **アプローチ** | 理論スペクトル（in silico 消化） | 深層学習ライブラリ予測 |
| **ライブラリ** | 不要（library-free） | AlphaPeptDeep が予測生成 |
| **FDR推定** | sage 内蔵 | PyProphet（外部ツール） |
| **RT予測** | なし | あり（深層学習モデル） |
| **フラグメント強度予測** | なし | あり（深層学習モデル） |
| **商用ライセンス** | MIT | BSD + Apache 2.0 |
| **インストール** | 単一バイナリ | 3ツール必要 |
| **処理時間（32ファイル）** | 約10分 | 約60-90分（環境依存） |
| **セットアップの簡単さ** | ★★★★★ | ★★★ |
| **同定タンパク質数（期待）** | 2,110 | 2,000-3,000（環境・パラメータ依存） |

:::message
**どちらを選ぶべきか？**

- **手軽さ・速さ重視** → sage（#4）がおすすめ。単一バイナリ、10分で完了
- **深層学習ベースを試したい** → OpenMS + AlphaPeptDeep（#5）。DIA-NNと同じアプローチをオープンソースで体験
- **両方試して比較** → 両方の `protein_matrix_from_*.csv` を #6 の前処理に渡して、下流の結果を比較するのも勉強になります

どちらの出力も同じ形式（タンパク質 × サンプルのCSV）なので、以降の解析ステップに影響はありません。
:::

## 実行方法

**スクリプトで一括実行する場合:**

```bash
micromamba activate crc-proteomics
python scripts/step_06_openms.py
```

**Notebook でセルごとに実行する場合:**

```bash
jupyter notebook notebooks/step_05.ipynb
```

:::message alert
**処理時間について**: OpenMS + AlphaPeptDeep パイプラインは sage より時間がかかります（合計60-90分程度）。特に Step 1（ライブラリ予測）と Step 2（OpenSWATH検索）が時間を要します。GPU がある環境では Step 1 が大幅に高速化されます。Notebook で実行する場合は、各ステップの出力ファイルが既に存在するかチェックし、スキップするロジックを入れると再実行が楽になります。
:::

## まとめ

OpenMS + AlphaPeptDeep は、DIA-NN と同じ「深層学習ライブラリ予測 → DIA 検索」アプローチを**完全オープンソース（BSD + Apache 2.0）**で実現するパイプラインです。sage の library-free アプローチと比べて手順は多いですが、保持時間とフラグメント強度の予測を活用することで、同定精度の向上が期待できます。

本シリーズでは sage（#4）をメインパイプラインとして採用していますが、この記事のパイプラインで生成した `protein_matrix_from_openms.csv` でも、以降の前処理・可視化・統計解析ステップはそのまま実行できます。

> 前回: [#4 sageで同定・定量](article-04-sage.md)
> 次回: [#6 前処理](article-06-preprocess.md) — Log2変換と欠損値補完

Sources:
- [OpenMS Documentation](https://openms.de/)
- [AlphaPeptDeep GitHub](https://github.com/MannLabs/alphapeptdeep)
- [PyProphet GitHub](https://github.com/PyProphet/pyprophet)
- [Wen et al. 2024 Nat. Commun.](https://doi.org/10.1038/s41467-024-46062-3)
- [Rost et al. 2014 Nat. Biotechnol.](https://doi.org/10.1038/nbt.2841)

#バイオインフォマティクス #プロテオミクス #OpenMS #AlphaPeptDeep #labcode
