# step_03_inspect_mzml.py の解説

## このコードの役割

mzMLファイルの内容を詳細に検査し、DIA解析に適したデータかどうかを確認するスクリプトです。各ファイルのスペクトル数、m/z範囲、保持時間、データ形式などのメタ情報を収集して表形式で出力します。

## 全体フロー

```
mzMLファイル群 → pymzmlで解析 → メタ情報収集 → CSV形式で出力
```

| | ファイル | 内容 |
|-|---------|------|
| 入力 | `data/raw/raw_mzML/*.mzML` | DIA-MSデータ（mzML形式） |
| 出力 | `results/tables/mzml_inventory.csv` | 各ファイルのメタ情報一覧 |

**所要時間の目安**: 約5-10分（32ファイル、ファイルサイズによる）

## セクション別解説

### セクション1: ヘルパー関数（Centroid判定）

```python
def _is_centroided(element) -> bool:
    for child in element.iter():
        acc = child.attrib.get("accession", "")
        if acc == "MS:1000127":  # centroid spectrum
            return True
        if acc == "MS:1000128":  # profile spectrum
            return False
    return False
```

**目的**: スペクトルデータがcentroid化（ピーク検出済み）かprofile（生データ）かを判定する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `element.iter()` | XML要素を再帰的に走査 |
| `attrib.get("accession")` | CV（Controlled Vocabulary）パラメータのアクセッション番号を取得 |

**引数の意味**:
- `MS:1000127`: 「centroid spectrum」を示すCV用語
- `MS:1000128`: 「profile spectrum」を示すCV用語

**つまずきやすいポイント**:
- mzMLはXML形式なので、名前空間の扱いに注意が必要
- 一部のファイルでCV パラメータが欠落している場合がある

### セクション2: mzMLファイル検査

```python
def inspect_mzml(path: str) -> dict:
    reader = pymzml.run.Reader(path)
    for spec in reader:
        total_spectra += 1
        ms_level = spec.ms_level
        rt = spec.scan_time_in_minutes()

        if ms_level == 1:
            ms1_count += 1
            # m/z範囲の計算
            mz_arr = spec.peaks("raw")[:, 0]

        elif ms_level == 2:
            ms2_count += 1
            # Isolation window情報の取得
```

**目的**: 各mzMLファイルを走査して、MS1/MS2スペクトル数、m/z範囲、保持時間範囲などを集計する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `pymzml.run.Reader()` | mzMLファイルを効率的に読み込むイテレータ |
| `spec.ms_level` | スペクトルのMSレベル（1=MS1、2=MS2） |
| `spec.scan_time_in_minutes()` | スペクトル取得時の保持時間（分） |
| `spec.peaks("raw")` | m/z値と強度のペア配列 |

**引数の意味**:
- `peaks("raw")`: ピークデータを生の数値配列として取得
- `ms_level`: 1=プリカーサー（親イオン）、2=フラグメント（子イオン）

**つまずきやすいポイント**:
- ファイルサイズが大きい場合、メモリ使用量に注意
- 一部のスペクトルでピークデータが空の場合がある

### セクション3: データ集計と出力

```python
# 全ファイルを処理
inventory = []
mzml_files = glob.glob(os.path.join(MZML_DIR, "*.mzML"))
for mzml_file in mzml_files:
    result = inspect_mzml(mzml_file)
    inventory.append(result)

# DataFrame作成とCSV出力
df = pd.DataFrame(inventory)
df.to_csv(OUTPUT_CSV, index=False)
```

**目的**: 各ファイルの検査結果を統合し、分析しやすい表形式で保存する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `glob.glob()` | パターンに一致するファイル一覧を取得 |
| `pd.DataFrame()` | 辞書のリストからデータフレームを作成 |
| `to_csv()` | CSVファイルとして保存 |

## 【深掘り】 Centroid vs Profile とは？

- **ひとことで**: 質量分析データの格納形式の違い
- **Profile**: 装置から得られる生の連続データ（アナログ信号をデジタル化）
- **Centroid**: ピーク検出処理を行い、各ピークを1点で代表した形式
- **どんなとき使う**: DIA解析ではcentroid化されたMS2データが必須
- **落とし穴**: Profile データは容量が大きく、解析ツールによっては対応していない

## 【深掘り】 MS1 vs MS2 とは？

- **ひとことで**: 質量分析の段階の違い
- **MS1**: プリカーサーイオン（親イオン）の質量を測定
- **MS2**: MS1で選択したイオンをフラグメント化して子イオンを測定
- **DIA-MS**: MS1で全範囲、MS2で連続的なウィンドウ（例：20Da幅）でフラグメント化
- **参考**: 論文 Methods セクションの「DIA-MS data acquisition」

## よくある質問

**Q: なぜmzMLファイルの検査が必要なの？**
A: RAW→mzML変換時のパラメータによって、centroid化の有無や質量範囲が変わります。DIA解析ツールは特定の形式を前提とするため、事前確認が重要です。

**Q: pymzmlとは何ですか？**
A: mzMLファイル（質量分析の標準フォーマット）を効率的に読み込むPythonライブラリです。MITライセンスで無料利用可能です。

**Q: 処理に時間がかかる場合は？**
A: mzMLファイルのサイズが大きい場合は時間がかかります。並列処理やより高速なライブラリ（pyOpenMS等）の利用も検討できます。

**Q: 自分のデータで試したい場合は？**
A: `MZML_DIR`のパスを変更し、自分のmzMLファイルを配置してください。ファイル拡張子が`.mzML`になっていることを確認してください。

## 出力されるCSVの内容

| 列名 | 内容 | 例 |
|------|------|-----|
| file | ファイル名 | `CRC001_Normal.mzML` |
| size_MB | ファイルサイズ（MB） | 245.7 |
| total_spectra | 総スペクトル数 | 15432 |
| ms1_count | MS1スペクトル数 | 2341 |
| ms2_count | MS2スペクトル数 | 13091 |
| ms1_mz_range | MS1のm/z範囲 | 400.0-1600.0 |
| ms2_mz_range | MS2のm/z範囲 | 200.0-1800.0 |
| rt_range_min | 保持時間範囲（分） | 85.2 |
| centroided_ms1 | MS1がcentroid化済みか | True |
| centroided_ms2 | MS2がcentroid化済みか | True |

## 重要な注意点

- **大容量ファイル**: 数百MB～数GBのファイルが多数あると時間がかかります
- **メモリ使用量**: 大きなmzMLファイルでは一時的にメモリを多く使用
- **依存関係**: pymzml>=2.5が必要
- **エラー処理**: 一部のファイルが破損している場合は該当ファイルをスキップ

## 次のステップ

このスクリプトの出力（`mzml_inventory.csv`）は以下で利用されます：
- ブログ記事での数値引用（ファイル数、スペクトル数など）
- DIA解析パラメータの妥当性確認
- データ品質レポートの作成