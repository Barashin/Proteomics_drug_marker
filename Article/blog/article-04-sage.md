---
title: "商用利用OKのsage-proteomicsでDIA-MSデータを解析する【論文再現シリーズ #4】"
emoji: "🔮"
type: "tech"
topics: ["proteomics", "sage", "bioinformatics", "labcode"]
published: false
---

# 商用利用OKのsage-proteomicsでDIA-MSデータを解析する

## はじめに

この記事では、DIA-MSの生データからタンパク質を同定・定量する方法を解説します。論文（Toyota et al. 2025）では **DIA-NN v1.8.1** を使用していますが、DIA-NN は**商用利用に有料ライセンスが必要**です。

本シリーズでは **sage-proteomics（MIT ライセンス、完全無料・商用可）** で代替します。

:::message
**この記事で行う処理**
DIA-MSの生データ（mzMLファイル）をsage-proteomicsに入力し、ヒトプロテオームのFASTAデータベースとの照合によってペプチドを同定（どのタンパク質由来か特定）し、各サンプルにおける強度（量）を定量します。最終的にペプチドレベルの定量値をタンパク質レベルに集約し、「タンパク質×サンプル」の定量マトリクスを作成します。このマトリクスが以降のすべての統計解析・可視化の入力データになります。
:::

## 生データから定量マトリクスまで：処理の全体像

この記事で行う処理を、「入力データ」「処理の中身」「出力データ」の3段階に分けて説明します。

### 入力：mzMLファイル（質量分析の生データ）とは

質量分析計に組織サンプルを入れると、サンプル中のペプチド（タンパク質の断片）が **質量/電荷比（m/z）** と **シグナル強度** のペアとして記録されます。これを「マススペクトル」と呼びます。

```
1つのmzMLファイルの中身（イメージ）:

スペクトル#1  → [(m/z=300.15, 強度=1200), (m/z=450.23, 強度=8500), ...]
スペクトル#2  → [(m/z=301.02, 強度=900),  (m/z=512.41, 強度=15000), ...]
...
スペクトル#3500 → [...]

→ 1ファイルあたり数千枚のスペクトル。
  「どのm/zにどれくらいのシグナルがあったか」の生記録。
```

この段階では **「どのタンパク質がどれくらいあるか」は全くわかりません**。m/zの数値の羅列があるだけです。DIA（Data Independent Acquisition）方式では、質量範囲を一定幅のウィンドウで網羅的にスキャンするため、サンプル中のほぼすべてのペプチドを漏れなく記録できるのが特徴です。

### 処理：sageが行う「同定」と「定量」

sage は以下の手順で、m/zの数値データを「タンパク質名 × 強度」に変換します。

**ステップ1：理論スペクトルの生成（in silico消化）**

ヒトの全タンパク質のアミノ酸配列が書かれたFASTAファイルを読み込み、トリプシン（タンパク質を切る酵素）で消化した場合に生じるペプチド断片を計算機上で生成します。さらに、各ペプチドが質量分析計の中で壊れたときに生じるフラグメントイオンの理論的なm/z値を計算します。

```
FASTAファイル（ヒト全タンパク質 約2万個の配列）
  ↓ トリプシン消化をシミュレーション
約300万個のペプチド断片（理論値）
  ↓ フラグメントイオンのm/zを計算
約6,300万個の理論フラグメント
```

**ステップ2：照合（同定）**

実測スペクトル（mzMLの中身）と理論スペクトルを比較し、「このm/zのピークはこのペプチドのフラグメントだ」とマッチングします。十分な数のフラグメントが一致すれば、そのペプチドが「同定」されます。

**ステップ3：定量（LFQ）**

同定されたペプチドについて、mzMLファイル内のシグナル強度（ピーク面積）を計測し、各サンプルにおける量を数値化します。これが **Label-Free Quantification（LFQ）** です。

### 出力：タンパク質×サンプルの定量マトリクス

sageの直接の出力はペプチドレベルの定量値（`lfq.tsv`）です。これを後段のスクリプトでタンパク質レベルに集約し、最終的に以下の形のCSVファイルが得られます。

```
               CRC01-N    CRC01-T    CRC02-N    CRC02-T   ...  CRC16-T
タンパク質A     18.5       22.3       19.1       21.8      ...  23.0
タンパク質B     15.2       14.8       15.9       15.1      ...  14.5
タンパク質C      NaN       12.3       11.8       13.0      ...  12.7
...
（2,110行 × 32列）
```

- **行**：同定されたタンパク質（2,110個）
- **列**：各サンプル（16患者 × Normal/Tumor = 32サンプル）
- **値**：そのサンプルにおけるタンパク質の強度（量）。NaNは検出限界以下で検出されなかったことを意味します

このマトリクスが以降すべてのステップ（前処理→可視化→統計検定→データベース照合）の入力データになります。

---

| ツール | ライセンス | 商用利用 | 本書採用 |
|-------|-----------|---------|:---:|
| DIA-NN | 学術無料 / 商用有料 | ❌ | |
| FragPipe (MSFragger) | 学術無料 / 商用有料 | ❌ | |
| Spectronaut | 商用有料 | ⚠ | |
| **sage-proteomics** | **MIT** | **✅** | **✓** |
| OpenMS + AlphaPeptDeep | BSD + Apache 2.0 | ✅ | |
| MaxQuant / MaxDIA | 非商用無料 | ⚠ | |

## 前提

- [#1 環境構築](article-01-setup.md) が完了していること
- [#2 データ取得](article-02-data.md) で mzML ファイルが手元にあること
- ヒトプロテオームのFASTAファイル（下記「FASTAファイルの準備」で取得）
- **対応Notebook**: [`notebooks/step_04.ipynb`](../notebooks/step_04.ipynb) — この記事のコードをセルごとに実行できます

## FASTAファイルの準備（ヒトプロテオーム配列）

sage はmzMLファイルだけでは動きません。「どのタンパク質が存在しうるか」の候補リストとして、**ヒト全タンパク質のアミノ酸配列を記録したFASTAファイル**が必要です。

### FASTAファイルとは

FASTAファイルは、タンパク質（または遺伝子）の配列をテキストで記述するフォーマットです。`>` で始まるヘッダ行と、それに続くアミノ酸配列の組で構成されます。

```
>sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens GN=TP53
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYPQGLNGTVNLPGRNSFEV
...
```

- `sp|P04637|P53_HUMAN`: UniProt ID とエントリ名
- `GN=TP53`: 遺伝子シンボル（Gene Symbol）。sageの出力をタンパク質名に変換する際に使います
- 以降の行: アミノ酸の1文字表記による配列

ヒトの場合、約 **20,000 個** のタンパク質エントリが含まれます。

### UniProtからダウンロードする

FASTAファイルは **UniProt**（https://www.uniprot.org）からダウンロードできます。ヒトの参照プロテオーム（reviewed、Swiss-Prot）を取得します。

**ブラウザの場合:**

1. https://www.uniprot.org/uniprotkb?query=reviewed:true+AND+organism_id:9606 にアクセス
2. ページ上部の **「Download」** ボタンをクリック
3. Format で **「FASTA (canonical)」** を選択してダウンロード

**コマンドラインの場合:**

```bash
# UniProt からヒト reviewed プロテオームをダウンロード
curl -o data/raw/human_proteome.fasta \
  "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=reviewed:true+AND+organism_id:9606"

# 確認: エントリ数を表示
grep -c "^>" data/raw/human_proteome.fasta
# → 約 20,000 エントリ
```

:::message alert
**FASTAファイルのバージョンについて**: UniProtは定期的に更新されるため、ダウンロード時期によってエントリ数が若干異なります。論文と完全に同じ結果を再現したい場合は、論文で使用されたバージョンのFASTAを使うのが理想ですが、本シリーズでは最新版で問題ありません。
:::

## sage-proteomics とは

| 項目 | 内容 |
|------|------|
| 開発者 | Michael Lazear |
| 初版 | 2022年 |
| 言語 | Rust（単一バイナリ、依存なし） |
| ライセンス | **MIT**（商用完全OK） |
| DIA サポート | ✅ library-free DIA ネイティブ対応 |
| 論文 | Lazear MR. *J. Proteome Res.* 2023. DOI: 10.1021/acs.jproteome.3c00486 |
| GitHub | https://github.com/lazear/sage |

**sage の強み**:

- **単一バイナリ**でインストールが超簡単
- **Rust 製**で DIA-NN 並みに高速（本書のデータ 32ファイル を **10.7分** で解析）
- **MIT ライセンス**なので企業・出版物でも自由に使える
- **library-free** で FASTA を渡すだけで in silico 消化から定量まで完結

## sage のインストール

```bash
micromamba activate crc-proteomics
micromamba install -c bioconda -c conda-forge sage-proteomics -y
sage --version
# → sage 0.14.6
```

## sage 設定ファイル (sage_config.json)

論文の DIA-NN パラメータに合わせた設定ファイルを作ります。

```json
{
  "database": {
    "bucket_size": 32768,
    "fragment_min_mz": 200.0,
    "fragment_max_mz": 1800.0,
    "peptide_min_mass": 600.0,
    "peptide_max_mass": 4000.0,
    "enzyme": {
      "missed_cleavages": 1,
      "min_len": 7, "max_len": 45,
      "cleave_at": "KR", "restrict": "P", "c_terminal": true
    },
    "static_mods": { "C": 57.02146 },
    "variable_mods": {},
    "max_variable_mods": 1,
    "ion_kinds": ["b", "y"],
    "decoy_tag": "rev_", "generate_decoys": true,
    "fasta": "data/raw/human_proteome.fasta"
  },
  "quant": {
    "lfq": true,
    "lfq_settings": {
      "peak_scoring": "Hybrid",
      "integration": "Sum",
      "spectral_angle": 0.7,
      "ppm_tolerance": 10.0,
      "combine_charge_states": true
    }
  },
  "precursor_tol": { "ppm": [-10, 10] },
  "fragment_tol":  { "ppm": [-10, 10] },
  "precursor_charge": [2, 4],
  "isotope_errors": [0, 0],
  "deisotope": false,
  "chimera": true,
  "wide_window": true,
  "predict_rt": false,
  "min_peaks": 15, "max_peaks": 150,
  "min_matched_peaks": 4,
  "max_fragment_charge": 2,
  "report_psms": 1,
  "output_directory": "results/sage_output"
}
```

### 論文の DIA-NN パラメータとの対応

| 論文（DIA-NN） | sage の対応 |
|---------------|------------|
| `--fasta-search`（library-free） | `wide_window: true` + FASTA指定 |
| `--cut K*,R*`（トリプシン） | `enzyme.cleave_at: "KR"` |
| `--missed-cleavages 1` | `enzyme.missed_cleavages: 1` |
| `--min-pep-len 7` | `enzyme.min_len: 7` |
| `--max-pep-len 45` | `enzyme.max_len: 45` |
| `--pr-charges 2-4` | `precursor_charge: [2, 4]` |
| `--min-fr-mz 200` / `--max-fr-mz 1800` | `fragment_min_mz: 200` / `fragment_max_mz: 1800` |
| MS1/MS2 accuracy 10 ppm | `precursor_tol.ppm: [-10,10]` / `fragment_tol.ppm` |
| Cys carbamidomethylation | `static_mods: { "C": 57.02146 }` |
| FDR <1% | 出力の `peptide_q < 0.01` で後段フィルタ |
| MBR (Match Between Runs) | sage 0.14 では実装なし（将来実装予定） |

## sage の実行

```bash
sage \
  --fasta data/raw/human_proteome.fasta \
  --output_directory results/sage_output \
  --disable-telemetry-i-dont-want-to-improve-sage \
  scripts/sage_config.json \
  data/raw/raw_mzML/*.mzML
```

### 実行ログの要点

```
[INFO sage] generated 62769656 fragments, 3066221 peptides in 6029ms
[INFO sage] processing files 0 .. 5
[INFO sage] - search:     72700 ms (3507 spectra/s)
...
[INFO sage] finished in 262s
```

**約 10.7 分で 32 ファイルの DIA 解析が完了**します。DIA-NN と同等以上の速度です。

## 出力ファイル

```
results/sage_output/
├── results.json       (パラメータ記録)
├── results.sage.tsv   (PSM テーブル)
└── lfq.tsv            (ペプチド × サンプル定量マトリクス、27,586 行)
```

### lfq.tsv の中身

```
peptide   charge  proteins              q_value  score  CRC01-N.mzML   ...   CRC16-T.mzML
GLGTDEDTIIDIITHR  -1  sp|P08133|ANXA6_HUMAN  0.00122  0.815  94488098.3  ...   ...
DVVIC[+57.02146]PDASLEDAK  -1  sp|Q99497|PARK7_HUMAN  0.00207  0.742  28403038.6  ...   ...
```

- `charge = -1`: charge state は `combine_charge_states: true` により統合済み
- `q_value`: ペプチドレベル FDR。`< 0.01` でフィルタして論文と同じ 1% 基準に

## タンパク質マトリクスへの集約

sage の lfq.tsv はペプチドレベルなので、ダウンストリーム解析（第5章以降）で使うために **タンパク質レベルに集約** します。ここでは Notebook（`step_04.ipynb`）の全コードを 3 パートに分けて掲載します。

### Part 1: mzML ファイルの検査

sage を実行する前に、入力となる mzML ファイルが正しく揃っているかを **pymzml** で検査します。各ファイルの MS1/MS2 スペクトル数や保持時間（RT）範囲を確認し、欠損や破損がないことを確かめます。

```python
# 標準ライブラリ: ファイル操作(os)、パターン検索(glob)、正規表現(re)、時間計測(time)、外部コマンド実行(subprocess)、システム操作(sys)
import os, glob, re, time, subprocess, sys
import pandas as pd   # データフレーム操作ライブラリ（表形式データの読み書き・加工に使う）
import pymzml         # mzMLファイル（質量分析データ）をPythonで読み込むためのライブラリ
```

```python
# --- パス設定（プロジェクト内の各ディレクトリ・ファイルへの相対パスを定数として定義） ---
MZML_DIR    = "../data/raw/raw_mzML"                # mzML生データが格納されているディレクトリ
FASTA_PATH  = "../data/raw/human_proteome.fasta"     # ヒトプロテオームのFASTA配列ファイルのパス
RESULTS_DIR = "../results"                           # 解析結果の出力先ルートディレクトリ
SAGE_CONFIG = "../scripts/sage_config.json"          # sageの設定ファイル（検索パラメータを記述）
SAGE_OUT    = os.path.join(RESULTS_DIR, "sage_output")  # sageの出力先ディレクトリ（lfq.tsvなどが生成される）
TABLES_DIR  = os.path.join(RESULTS_DIR, "tables")       # 検査結果テーブルの保存先ディレクトリ

# 出力先ディレクトリが存在しなければ作成する（exist_ok=Trueで既存でもエラーにならない）
os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(SAGE_OUT, exist_ok=True)
```

```python
def inspect_mzml(path):
    """mzML を走査して MS1/MS2 数と RT 範囲を返す。"""
    # pymzmlでmzMLファイルを開き、スペクトルを1枚ずつ読み込めるReaderオブジェクトを作成
    reader = pymzml.run.Reader(path)
    ms1, ms2 = 0, 0  # MS1（前駆体スキャン）とMS2（フラグメントスキャン）のカウンターを初期化
    # 保持時間（RT）の最小値・最大値を追跡するための初期値（infは無限大）
    rt_min, rt_max = float("inf"), float("-inf")

    # mzMLファイル内の全スペクトルを1枚ずつ走査するループ
    for spec in reader:
        rt = spec.scan_time_in_minutes()  # そのスペクトルの保持時間（分単位）を取得
        if rt is not None:
            # 保持時間が取得できた場合、最小値・最大値を更新
            rt_min, rt_max = min(rt_min, rt), max(rt_max, rt)
        if spec.ms_level == 1:    # MS1スペクトル（前駆体イオンの全体像）の場合
            ms1 += 1              # MS1カウンターを加算
        elif spec.ms_level == 2:  # MS2スペクトル（フラグメントイオン）の場合
            ms2 += 1              # MS2カウンターを加算

    # 検査結果を辞書にまとめて返す（後でDataFrameの1行になる）
    return {
        "file": os.path.basename(path),                   # ファイル名のみ（ディレクトリパスを除去）
        "size_MB": round(os.path.getsize(path) / 1024**2, 1),  # ファイルサイズをMB単位に変換（小数1桁）
        "ms1": ms1,                    # MS1スペクトルの総数
        "ms2": ms2,                    # MS2スペクトルの総数
        "rt_min": round(rt_min, 2),    # 保持時間の最小値（分、小数2桁）
        "rt_max": round(rt_max, 2),    # 保持時間の最大値（分、小数2桁）
    }
```

```python
# mzML ファイル一覧を取得して検査を実行
# MZML_DIRディレクトリ内の全.mzMLファイルをglobで検索し、ファイル名順にソート
mzml_files = sorted(glob.glob(os.path.join(MZML_DIR, "*.mzML")))
print(f"{len(mzml_files)} mzML files found")  # 見つかったファイル数を表示（期待値: 32）

# リスト内包表記で全mzMLファイルに対してinspect_mzml関数を適用し、検査結果のリストを作成
records = [inspect_mzml(f) for f in mzml_files]
# 検査結果のリスト（辞書のリスト）をpandas DataFrameに変換（表形式にする）
df_inv = pd.DataFrame(records)
# 検査結果をCSVファイルとして保存（index=Falseで行番号は含めない）
df_inv.to_csv(os.path.join(TABLES_DIR, "mzml_inventory.csv"), index=False)
df_inv  # Jupyter Notebook上でDataFrameの内容を表示
```

### Part 2: sage の Python 実行

sage をコマンドラインではなく Python の `subprocess` 経由で実行します。こうすると実行時間の計測やログのリアルタイム表示を Notebook 内で完結できます。

```python
# sageコマンドの引数をリストとして構築（subprocessに渡す形式）
cmd = [
    "sage",                          # 実行するコマンド名（sage本体）
    "--fasta", FASTA_PATH,           # FASTAファイルのパスを指定（タンパク質配列データベース）
    "--output_directory", SAGE_OUT,  # 出力先ディレクトリを指定
    "--disable-telemetry-i-dont-want-to-improve-sage",  # テレメトリ（使用統計の送信）を無効化
    SAGE_CONFIG,                     # 設定ファイル（検索パラメータ）のパス
] + mzml_files  # 全mzMLファイルのパスをコマンド引数に追加

print(f"Running sage on {len(mzml_files)} files...")  # 解析対象ファイル数を表示
t0 = time.time()  # 実行開始時刻を記録（後で所要時間を計算するため）

# sage を子プロセスとして起動し、ログをリアルタイム出力
# stdout=PIPEで標準出力をキャプチャ、stderr=STDOUTでエラー出力も標準出力に統合、text=Trueで文字列として扱う
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
# sageの出力を1行ずつ読み取り、Notebook上にリアルタイム表示するループ
for line in proc.stdout:
    sys.stdout.write(line)  # 各ログ行をそのまま標準出力に書き出す
proc.wait()  # sageプロセスの終了を待機（終了コードが確定する）

# 所要時間を分単位で表示し、終了コードも出力（0なら正常終了）
print(f"\nDone in {(time.time() - t0) / 60:.1f} min (exit code {proc.returncode})")
```

### Part 3: タンパク質マトリクス構築

sage の `lfq.tsv`（ペプチドレベル LFQ）を FASTA の `GN=` フィールドで Gene Symbol にマップし、Gene × サンプルの定量マトリクスに集約します。

**ステップ 3-1: FASTA から Gene Symbol マッピング関数を定義**

```python
def parse_fasta_gene_map(fasta_path):
    """UniProt FASTA ヘッダから accession -> Gene Symbol の辞書を作る。

    ヘッダ例: >sp|P04637|P53_HUMAN ... GN=TP53 ...
    GN= がなければ Entry Name の先頭語 (例: P53) を使う。
    """
    id_to_gene = {}  # UniProtアクセッション番号をキー、Gene Symbolを値とする辞書を初期化
    # 「GN=」の後に続く空白以外の文字列（=Gene Symbol）を抽出する正規表現パターン
    gn_re = re.compile(r"\bGN=(\S+)")
    # FASTAファイルを開いて1行ずつ読み込む
    with open(fasta_path) as f:
        for line in f:  # FASTAファイルの全行を走査
            if not line.startswith(">"):
                continue  # 「>」で始まらない行はアミノ酸配列なのでスキップ（ヘッダ行のみ処理）
            # ヘッダ行を「|」で最大3つに分割（例: "sp", "P04637", "P53_HUMAN ..."）
            parts = line[1:].split("|", 2)
            if len(parts) >= 3:
                # UniProt標準形式（sp|アクセッション|エントリ名）の場合
                acc = parts[1]                     # アクセッション番号（例: P04637）
                entry_name = parts[2].split()[0]   # エントリ名（例: P53_HUMAN）
            else:
                # 標準形式でないヘッダの場合、先頭の単語をアクセッション番号として使う
                acc = line[1:].split()[0]
                entry_name = acc
            m = gn_re.search(line)  # ヘッダ行から「GN=遺伝子名」パターンを検索
            # GN=が見つかればその遺伝子名を、なければエントリ名の「_」前の部分を代替として使用
            id_to_gene[acc] = m.group(1) if m else entry_name.split("_")[0]
    return id_to_gene  # {アクセッション番号: Gene Symbol} の辞書を返す
```

**ステップ 3-2: lfq.tsv の読み込みと FDR フィルタリング**

```python
# ペプチド FDR（偽発見率）の閾値を定義（論文と同じ 1% = 0.01）
# FDRが1%未満のペプチドのみを使うことで、誤同定を統計的に制御する
Q_THRESHOLD = 0.01

# sageが出力したlfq.tsv（ペプチドレベルのLFQ定量値テーブル）のパスを構築
LFQ_TSV = os.path.join(SAGE_OUT, "lfq.tsv")
# タブ区切り（TSV）ファイルをpandas DataFrameとして読み込む
df = pd.read_csv(LFQ_TSV, sep="\t")
print(f"lfq.tsv: {len(df)} peptides, {len(df.columns)} columns")  # 行数（ペプチド数）と列数を表示

# FDR フィルタ: q_valueが閾値未満のペプチドだけを残す（信頼性の低い同定を除外）
# .copy()で新しいDataFrameを作成し、元データへの参照を切る（SettingWithCopyWarning防止）
df = df[df["q_value"] < Q_THRESHOLD].copy()
print(f"After q < {Q_THRESHOLD}: {len(df)} peptides")  # フィルタ後のペプチド数を表示

# メタデータ列（ペプチド配列、電荷、タンパク質名など）の名前を集合として定義
meta_cols = {"peptide", "charge", "proteins", "q_value", "score", "spectral_angle"}
# メタデータ列以外の列名を抽出 → これらがサンプルごとの強度値列（CRC01-N.mzML, CRC01-T.mzMLなど）
sample_cols = [c for c in df.columns if c not in meta_cols]

# 列名から「.mzML」サフィックスを除去して、サンプル名をクリーンにする辞書を作成
rename = {c: c.replace(".mzML", "") for c in sample_cols}
df = df.rename(columns=rename)       # DataFrameの列名を一括変更
sample_cols = list(rename.values())   # 変更後のサンプル名リストを更新
print(f"Samples: {len(sample_cols)}")  # サンプル数を表示（期待値: 32）
```

**ステップ 3-3: Gene Symbol マッピングとマトリクス集約**

```python
# FASTAファイルを解析して、アクセッション番号→Gene Symbolの変換辞書を構築
id_to_gene = parse_fasta_gene_map(FASTA_PATH)
print(f"FASTA entries: {len(id_to_gene)}")  # FASTAに含まれるエントリ数を表示（約20,000）

# proteins列（例: "sp|Q99497|PARK7_HUMAN;sp|..."）から代表アクセッション番号を抽出
# まず「;」で分割して最初のタンパク質を取り、次に「|」で分割して2番目の要素（=アクセッション番号）を取得
df["acc"] = df["proteins"].str.split(";").str[0].str.split("|").str[1]
# アクセッション番号をGene Symbolに変換（マッピングできなかった場合はアクセッション番号をそのまま使用）
df["gene"] = df["acc"].map(id_to_gene).fillna(df["acc"])

# 同じGene Symbolを持つペプチドの強度を合計して、タンパク質レベルに集約
# min_count=1: 全てNaNの場合はNaNを返す（0にしない）。生の強度値なのでsumが適切
matrix = df.groupby("gene")[sample_cols].sum(min_count=1)
# 値が0の箇所をNaN（欠損値=未検出）に置き換え、全サンプルでNaNの行（全く検出されなかったタンパク質）を削除
matrix = matrix.replace(0, pd.NA).dropna(how="all")
matrix.index.name = "Protein"  # インデックス名を「Protein」に設定（CSV出力時のヘッダになる）
# マトリクスのサイズを表示（行数=タンパク質数、列数=サンプル数）
print(f"Protein matrix: {matrix.shape[0]} proteins x {matrix.shape[1]} samples")
matrix.head()  # 先頭5行を表示して中身を確認
```

**ステップ 3-4: CSV 保存**

```python
# 出力CSVファイルのパスを構築（タンパク質×サンプルの定量マトリクス）
out_csv = os.path.join(RESULTS_DIR, "protein_matrix_from_sage.csv")
# DataFrameをCSVファイルとして保存（インデックス=タンパク質名も含める）
# このCSVが以降の全ステップ（前処理・可視化・統計検定・データベース照合）の入力データになる
matrix.to_csv(out_csv)
print(f"Saved: {out_csv}")  # 保存先パスを表示して確認
```

## 結果サマリー（本書での実測値）

| 指標 | 本書 (sage, 32ファイル) | 論文 (DIA-NN, 32ファイル) |
|------|----------------------|------------------------|
| ペプチド数（生） | **27,586** | 記載なし |
| ペプチド数 (FDR<1%) | **8,391** | 記載なし |
| **タンパク質数** | **2,110** | **10,329** |
| 解析時間 | **10.7分** (Apple Silicon) | 記載なし |

タンパク質数は論文 (10,329) より少ない (2,110) ですが、これは以下の理由によります：

1. sage は理論スペクトルベース（DIA-NN は深層学習ライブラリ予測）
2. sage は MBR 未実装

**商用完全クリアで論文と同等のパイプラインを学べる** のが本書の立ち位置です。

## まとめ

sage-proteomics は MIT ライセンスで、商用利用にも制約がない DIA 解析エンジンです。単一バイナリで導入が簡単、Rust 製で高速、library-free DIA をネイティブにサポートし、DIA-NN の商用代替として最適です。

> 前回: [#3 RAW→mzML変換](article-03-convert.md)
> 次回: [#5 OpenMS + AlphaPeptDeepによるDIA解析](article-05-openms.md) — 深層学習ベースの代替パイプライン

Sources:
- [sage GitHub](https://github.com/lazear/sage)
- [Lazear 2023 J. Proteome Res.](https://doi.org/10.1021/acs.jproteome.3c00486)
- [DIA-NN License](https://github.com/vdemichev/DiaNN)

#バイオインフォマティクス #プロテオミクス #sage #DIA #labcode
