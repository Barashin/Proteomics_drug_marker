# step_04_dia_analysis_sage.py の解説

## このコードの役割

**Sage-proteomics**を使用してDIA（Data-Independent Acquisition）質量分析データからタンパク質の同定と定量を行うスクリプトです。これがプロテオミクス解析の**最も重要なステップ**であり、mzMLファイルから「どのタンパク質が」「どれくらい」存在するかを決定します。

## 全体フロー

```
mzMLファイル群 + FASTAデータベース → Sage解析エンジン → タンパク質同定・定量結果
      ↑                                   ↑                        ↓
（32ファイル）                    （設定ファイル）         PSM・ペプチド・プロテインテーブル
```

| | ファイル | 内容 |
|-|---------|------|
| **入力** | `data/raw/raw_mzML/*.mzML` | 質量分析データ（32ファイル、16患者×2条件） |
| **入力** | `data/raw/human_proteome.fasta` | ヒトタンパク質配列データベース（約20,000種類） |
| **入力** | `scripts/sage_config.json` | Sage解析パラメータ設定ファイル |
| **出力** | `results/sage_output/results.sage.tsv` | PSM（Peptide-Spectrum Match）テーブル |
| **出力** | `results/sage_output/lfq.tsv` | ペプチドレベルの定量テーブル |
| **出力** | `results/sage_output/results.json` | 実行時の詳細ログ |

**所要時間の目安**: 30分〜1時間（CPUコア数とRAM容量による）

## セクション別解説

### セクション1: Sageコマンドの可用性確認

```python
sage_bin = shutil.which("sage")
if sage_bin is None:
    raise SystemExit("sage コマンドが見つかりません...")
```

**目的**: Sageがインストールされているか、PATHが通っているかを確認

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `shutil.which()` | コマンドがPATHに存在するか確認（which コマンドのPython版） |
| `raise SystemExit()` | エラーメッセージとともにプログラムを終了 |

**なぜこのチェックが必要**:
Sageは外部のCLI（Command Line Interface）ツールです。Pythonスクリプト内から`subprocess`で呼び出すため、事前にインストール状況を確認します。

**つまずきやすいポイント**:
- **仮想環境の問題**: `sage` が `crc-proteomics` 環境にインストールされているが、実行時に別の環境になっている
- **PATHの問題**: インストールはされているが、システムPATHに追加されていない

### セクション2: 入力ファイルの検証

```python
mzml_files = glob.glob(os.path.join(MZML_DIR, "*.mzML"))
if len(mzml_files) == 0:
    raise FileNotFoundError(f"mzMLファイルが見つかりません: {MZML_DIR}")
```

**目的**: 必要な入力ファイルがすべて存在することを確認

**引数の意味**:
- `glob.glob()`: ワイルドカード（*）を使ったファイル検索
- `*.mzML`: 拡張子が .mzML のファイルをすべて取得

**期待されるファイル構成**:
```
data/raw/raw_mzML/
├── CRC_DDA_1_T.mzML      # 患者1 腫瘍（Tumor）
├── CRC_DDA_1_N.mzML      # 患者1 正常（Normal）
├── CRC_DDA_2_T.mzML      # 患者2 腫瘍
├── CRC_DDA_2_N.mzML      # 患者2 正常
└── ...（合計32ファイル）
```

### セクション3: Sage設定ファイルの読み込み

```python
CONFIG_PATH = os.path.join(SCRIPT_DIR, "sage_config.json")
```

**目的**: Sage解析のパラメータを外部ファイルから読み込み

**sage_config.jsonの主要設定**:

| パラメータ | 値 | 生物学的意味 |
|-----------|---|------------|
| `enzyme` | `trypsin` | トリプシン酵素消化（K・R残基で切断） |
| `missed_cleavages` | `1` | 切断されずに残った結合を1個まで許容 |
| `min_peptide_len` | `7` | 最小ペプチド長（7アミノ酸） |
| `max_peptide_len` | `45` | 最大ペプチド長（45アミノ酸） |
| `precursor_tol` | `10 ppm` | プリカーサイオン質量精度 |
| `fragment_tol` | `10 ppm` | フラグメントイオン質量精度 |

**なぜJSONファイルに分離するか**:
パラメータ調整が頻繁にあるため、コードを変更せずに設定だけ変更できるよう設計されています。

### セクション4: Sage実行（メイン処理）

```python
cmd = [
    sage_bin, CONFIG_PATH,
    "--fasta", FASTA_PATH,
    "--output-directory", OUT_DIR,
] + mzml_files
```

**目的**: Sageコマンドを構築してサブプロセスで実行

**コマンド構造の解説**:
```bash
sage [設定ファイル] --fasta [データベース] --output-directory [出力先] [mzMLファイル1] [mzMLファイル2] ...
```

これは以下のようなコマンドライン実行と等価です：
```bash
sage sage_config.json \
  --fasta human_proteome.fasta \
  --output-directory results/sage_output/ \
  CRC_DDA_1_T.mzML CRC_DDA_1_N.mzML ...
```

### セクション5: 進行状況のモニタリング

```python
process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
for line in process.stdout:
    print(line.strip())
```

**目的**: Sageの実行進捗をリアルタイムで表示

**なぜリアルタイム表示が必要**:
- DIA解析は時間がかかる（30分-1時間）
- 途中でエラーが発生した場合の早期発見
- 解析の進捗確認によるユーザーの安心感

## 【深掘り】 DIA-MS（Data-Independent Acquisition）とは？

- **ひとことで**: すべての質量範囲を網羅的に測定する質量分析手法
- **従来手法（DDA）との違い**:
  - DDA（Data-Dependent）: 強いシグナルから順番に測定（見逃しあり）
  - DIA: 質量範囲を分割して全域を測定（見逃しなし）
- **利点**: 再現性が高い、低発現タンパク質も検出可能
- **課題**: データ解析が複雑（スペクトルが混合している）
- **参考**: Toyota et al. 論文の Methods セクション "DIA-MS Analysis"

## 【深掘り】 Sage-proteomicsとは？

- **ひとことで**: 次世代のオープンソースプロテオミクス検索エンジン
- **開発言語**: Rust（メモリ安全性と高速性を両立）
- **ライセンス**: MIT（商用利用制限なし）
- **特徴**:
  - Library-free DIA解析（事前のスペクトルライブラリが不要）
  - マルチスレッド対応で高速
  -論文標準のDIA-NNと同等の性能
- **なぜ選択したか**: DIA-NNの商用利用制限を回避するため

## 【深掘り】 PSM（Peptide-Spectrum Match）とは？

- **ひとことで**: 測定されたスペクトルとペプチド配列の対応関係
- **詳細**: 質量分析で得られた個別のスペクトルが、どのペプチド（タンパク質断片）由来かを同定した結果
- **FDR（False Discovery Rate）**: 間違った同定の割合（1%未満に設定）
- **品質評価**: qvalue、posterior error probability等の統計指標で評価
- **次のステップ**: PSMからペプチド定量、さらにプロテイン定量へ集約

## よくある質問

**Q: 実行中に "Segmentation fault" エラーが出ます**
A:
1. **メモリ不足**: 32ファイル同時処理には8GB以上のRAMが推奨
2. **ファイル破損**: mzMLファイルが不完全な可能性。step_03で検証済みか確認
3. **設定エラー**: sage_config.jsonの数値設定にエラーがある可能性

**Q: 解析時間を短縮したい**
A:
1. **CPUコア数**: 多コアCPUで並列処理が高速化されます
2. **ファイル数削減**: テスト用に2-4ファイルでの実行を試す
3. **tolerance設定**: 質量精度を少し緩めると高速化（精度とのトレードオフ）

**Q: 結果の品質をどう評価すればいい？**
A:
1. **同定数**: 数千〜数万のPSMが得られることを確認
2. **FDR**: qvalue < 0.01の同定が大部分を占めることを確認
3. **分布**: ペプチド長・質量分布が妥当な範囲にあることを確認

**Q: 他の生物種で使いたい場合は？**
A:
1. **FASTAファイル**: UniProtから該当種のプロテオームをダウンロード
2. **設定調整**: トリプシン以外の酵素を使用している場合は enzyme設定を変更
3. **質量精度**: 測定装置に応じてtoleranceを調整

**Q: DIA-NNとの結果の違いは？**
A:
- **同定数**: 通常はSageの方がやや少なめ（より保守的）
- **再現性**: Sageの方が安定している場合が多い
- **商用利用**: Sageは完全自由、DIA-NNは制限あり