# step_03_diann.py の解説

## このコードの役割

OpenMSのOpenSWATHワークフローを使用してDIA-MSデータからタンパク質を同定・定量するスクリプトです。論文で使用されたDIA-NNの商用利用可能な代替実装です。

## 全体フロー

```
FASTAファイル → in silico消化 → スペクトルライブラリ → OpenSWATH → FDR制御 → タンパク質マトリクス
```

| | ファイル | 内容 |
|-|---------|------|
| 入力 | `data/raw/*.mzML` | DIA-MSの生データファイル |
| 入力 | `data/raw/human_proteome.fasta.gz` | ヒトプロテオーム配列データベース |
| 出力 | `results/openms_output/digested_peptides.fasta` | in silico消化で生成されたペプチドリスト |
| 出力 | `results/openms_output/*_openswath.tsv` | 各mzMLファイルの解析結果 |

**所要時間の目安**: 数時間（データサイズとCPU性能による）

## セクション別解説

### セクション1: FASTAファイルの準備

```python
def prepare_fasta():
    """gzip圧縮FASTAを展開する。"""
    with gzip.open(FASTA_GZ, "rb") as f_in:
        with open(FASTA, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
```

**目的**: gzip圧縮されたFASTAファイルを展開し、OpenMSツールで読み込み可能な形式にする。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `gzip.open()` | gzip圧縮ファイルを読み込み専用で開く |
| `shutil.copyfileobj()` | ファイルオブジェクト間で効率的にデータをコピーする |

**つまずきやすいポイント**:
- FASTAファイルが存在しない場合はStep 02を先に実行する必要がある
- バイナリモード（"rb", "wb"）でファイルを開く必要がある

### セクション2: In silico消化（Digestor）

```python
def run_digestor(fasta_path):
    cmd = [
        "Digestor",
        "-in", fasta_path,
        "-out", output_path,
        "-enzyme", "Trypsin",           # トリプシン消化
        "-missed_cleavages", "1",       # ミスクリーベージ1回許容
        "-min_length", "7",             # ペプチド最小長7
        "-max_length", "45",            # ペプチド最大長45
    ]
```

**目的**: タンパク質配列を酵素（トリプシン）で仮想的に消化し、生成されるペプチドのリストを作成する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `subprocess.run()` | 外部コマンド（Digestor）を実行し、結果を取得 |
| `timeout=600` | 10分でタイムアウトして処理を中断 |

**引数の意味**:
- `-enzyme Trypsin`: トリプシン酵素（K、Rの後ろで切断）を使用
- `-missed_cleavages 1`: 切断されずに残る箇所を最大1箇所まで許容
- `-min_length 7`: 7アミノ酸未満の短いペプチドは除外
- `-max_length 45`: 45アミノ酸超の長いペプチドは除外（MS検出困難）

**つまずきやすいポイント**:
- OpenMSがインストールされていない場合、Digestorコマンドが見つからない
- メモリ不足でプロセスが終了する場合がある（大きなFASTAファイル）

### セクション3: OpenSWATH によるDIA解析

```python
def run_openswath(mzml_file, library_path):
    cmd = [
        "OpenSwathWorkflow",
        "-in", mzml_file,
        "-tr", library_path,            # スペクトルライブラリ
        "-out_tsv", output_path,
        # DIA-NN準拠のパラメータ
        "-rt_extraction_window", "600", # 保持時間ウィンドウ
        "-mz_extraction_window", "50",   # m/z抽出ウィンドウ (ppm)
    ]
```

**目的**: DIAデータ（mzML）をスペクトルライブラリと照合して、ペプチドとタンパク質を同定・定量する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `OpenSwathWorkflow` | OpenMSのDIA解析コマンドライン ツール |
| `-rt_extraction_window` | 保持時間の許容誤差（秒単位） |
| `-mz_extraction_window` | 質量数の許容誤差（ppm単位） |

**引数の意味**:
- `-in`: 入力mzMLファイル（DIA-MSデータ）
- `-tr`: トランジション（スペクトルライブラリ）ファイル
- `-rt_extraction_window 600`: 保持時間±600秒（10分）の範囲で検索
- `-mz_extraction_window 50`: 質量数±50ppmの範囲で検索

## 【深掘り】 DIA-MSとSWATHとは？

- **ひとことで**: 全ペプチドを同時に解析する質量分析手法
- **定義**: Data-Independent Acquisition（DIA）は、特定のペプチドを選ばずに全質量範囲を一括解析する手法。SWATHはその実装の一つ
- **どんなとき使う**: 定量精度を重視し、見落としを最小化したいプロテオミクス解析
- **落とし穴**: データ解析が複雑で、適切なスペクトルライブラリが必要
- **参考**: 論文 Methods セクションの「DIA-MS data acquisition」

## 【深掘り】 In silico消化とは？

- **ひとことで**: コンピュータ上でタンパク質を仮想的に切断すること
- **定義**: FASTAデータベース中のタンパク質配列を、実際の酵素の切断規則に従ってペプチド断片に分割する計算処理
- **どんなとき使う**: 理論ペプチドライブラリを作成し、実測データと照合するため
- **落とし穴**: 実際の切断は完全でないため、ミスクリーベージも考慮が必要
- **参考**: 論文 Methods セクションの「Protein database searching」

## よくある質問

**Q: なぜDIA-NNではなくOpenMSを使うの？**
A: DIA-NNは商用利用に有料ライセンスが必要ですが、OpenMSはBSD 3-Clauseライセンスで完全に無料。機能的にはほぼ同等の解析が可能です。

**Q: 処理に数時間かかるのはなぜ？**
A: 32個のmzMLファイルそれぞれに対してスペクトルライブラリとの照合計算を実行するため。CPUコア数を増やすとある程度高速化できます。

**Q: 自分のデータで試したい。何を変えればいい？**
A: `RAW_DIR`のパスを変更し、自分のmzMLファイルとFASTAファイルを配置してください。酵素がトリプシン以外の場合は`-enzyme`パラメータを変更が必要です。

**Q: OpenMSがインストールされていない場合は？**
A: 以下のコマンドでインストールできます：
```bash
micromamba install -c bioconda -c conda-forge openms -y
```

## 重要な注意点

- **大容量データ**: mzMLファイル群で数十GBになる場合があります
- **処理時間**: 数時間から一日程度かかる可能性があります
- **メモリ使用量**: 16GB以上のRAMを推奨
- **依存関係**: OpenMS 3.0以降が必要

## 次のステップ

このスクリプトの出力（OpenSWATH結果）は、次のステップでタンパク質定量マトリクスに変換されます：
- Step 4: `step_04_build_protein_matrix.py`で定量マトリクス構築
- Step 5: `step_05_overview_visualization.py`で結果可視化