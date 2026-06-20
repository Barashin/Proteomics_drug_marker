# step_05_openms.py の解説

## このコードの役割

OpenMS + AlphaPeptDeep を組み合わせた深層学習ベースのDIA-MSデータ解析パイプラインです。sage（前のステップ）の代替として、AI予測を活用した高精度なタンパク質同定・定量を実現し、論文に近い検出性能を目指します。

## 全体フロー

```
FASTA → AlphaPeptDeep予測ライブラリ → OpenSWATH → PyProphet統計処理 → タンパク質マトリクス
```

| | ファイル | 内容 |
|-|---------|------|
| 入力 | `data/raw/raw_mzML/*.mzML` | DIA-MSの生データファイル |
| 入力 | `data/raw/human_proteome.fasta` | ヒトプロテオーム配列データベース |
| 出力 | `results/openms_output/library/predicted_library.tsv` | AlphaPeptDeep予測スペクトルライブラリ |
| 出力 | `results/openms_output/pyprophet/pyprophet_export.tsv` | PyProphet統計処理結果 |
| 出力 | `results/protein_matrix_from_openms.csv` | タンパク質×サンプルマトリクス |

**所要時間の目安**: 2-4時間（32サンプル、GPU環境でAlphaPeptDeep利用時）

## セクション別解説

### セクション1: ツール可用性チェック

```python
REQUIRED_TOOLS = ["peptdeep", "TargetedFileConverter", "OpenSwathWorkflow", "pyprophet"]

def ensure_tool(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(
            f"{name} が見つかりません。`micromamba install -c bioconda openms` および "
            f"`pip install peptdeep pyprophet` を実行してください。"
        )
```

**目的**: 必須ツールがインストール済みかを事前チェックし、不備があれば明確なインストール指示とともに終了する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `shutil.which()` | 実行パス上にコマンドが存在するかチェック |
| `RuntimeError` | インストール不備時に適切なエラーメッセージで終了 |

**引数の意味**:
- `peptdeep`: AlphaPeptDeepライブラリ（深層学習ベーススペクトル予測）
- `TargetedFileConverter`: OpenMSのライブラリ形式変換ツール
- `OpenSwathWorkflow`: OpenMSのDIA解析メインワークフロー
- `pyprophet`: 統計的FDR制御とスコアリングツール

**つまずきやすいポイント**:
- AlphaPeptDeepはpipインストール、OpenMSはcondaインストールが必要
- 環境によってはPyTorchの依存関係で競合することがある

### セクション2: AlphaPeptDeep予測ライブラリ生成

```python
# デフォルト設定を生成してから必要箇所のみ上書き
default_yaml = os.path.join(OPENMS_OUT, "peptdeep_default.yaml")
if subprocess.call(["peptdeep", "export-settings", default_yaml]) != 0:
    sys.exit(1)

# タスクと実行デバイス
settings["task_workflow"] = ["library"]
settings["torch_device"]["device_type"] = "cpu"

# ライブラリ予測設定
lib["fasta"]["protease"] = "trypsin"
lib["fasta"]["max_miss_cleave"] = 1
lib["decoy"] = "pseudo_reverse"
```

**目的**: FASTAファイルから深層学習モデルを使って理論スペクトルライブラリを生成する。保持時間・フラグメント強度・CCS値をTransformerで予測。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `peptdeep export-settings` | デフォルト設定YAMLファイルを生成 |
| `yaml.safe_load/safe_dump` | YAML設定ファイルの読み書き |
| `settings["task_workflow"] = ["library"]` | ライブラリ生成モードに設定 |

**引数の意味**:
- `protease: "trypsin"`: トリプシン酵素でペプチドを切断（K、Rの後ろ）
- `max_miss_cleave: 1`: 切断漏れを最大1箇所まで許容
- `decoy: "pseudo_reverse"`: 偽陽性制御用のデコイペプチドを生成
- `min_peptide_len: 7/max: 45`: MS検出可能な長さ範囲
- `min_precursor_charge: 2/max: 4`: 一般的なペプチド電荷状態

**つまずきやすいポイント**:
- GPU環境ではdevice_typeをGPUに設定すると大幅高速化
- メモリ不足時はペプチド長やプリカーサー範囲を制限する必要

### セクション3: OpenSWATH列名変換とライブラリ形式統一

```python
rename_map = {
    "RT": "NormalizedRetentionTime",
    "FragmentMz": "ProductMz",
    "RelativeIntensity": "LibraryIntensity",
    "StrippedPeptide": "PeptideSequence",
    "ModifiedPeptide": "ModifiedPeptideSequence",
    "ProteinID": "ProteinName",
}

# AlphaPeptDeep の `_PEPTIDE_` 形式から OpenSwath が読める形へ
mp = parts[mp_idx].strip("_")
for k, v in mod_map.items():
    mp = mp.replace(k, v)
```

**目的**: AlphaPeptDeep出力をOpenSWATH互換のTSV形式に変換し、修飾アミノ酸記法を統一する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `str.replace()` | 修飾記法の変換（[Carbamidomethyl] → (UniMod:4)） |
| `str.strip("_")` | AlphaPeptDeepの区切り文字を除去 |
| `TargetedFileConverter` | TSV → PQP（OpenMS標準ライブラリ）形式変換 |

**引数の意味**:
- `TransitionGroupId/TransitionId`: OpenSWATHが要求する一意識別子
- `(UniMod:4)`: カルバミドメチル化システインの標準記法
- `(UniMod:35)`: メチオニン酸化の標準記法
- `.pqp`: SQLiteベースのOpenMSライブラリ形式

**つまずきやすいポイント**:
- 修飾記法の不一致でライブラリ読み込みが失敗することが多い
- ProteinNameの複数エントリ（`Q9UL16;A0A087X182`）は先頭のみ採用

### セクション4: OpenSWATH DIA検索

```python
rc = run(
    [
        "OpenSwathWorkflow",
        "-in", mzml_path,
        "-tr", PQP_LIB,
        "-out_osw", osw_out,
        "-mz_extraction_window", str(MS2_PPM),
        "-mz_extraction_window_unit", "ppm",
        "-Scoring:stop_report_after_feature", "5",
        "-threads", str(args.threads),
    ],
    f"Step2b OpenSWATH {i+1}/{len(mzml_files)} {sample}",
)
```

**目的**: 各mzMLファイルを予測ライブラリと照合し、ペプチドピークを同定・定量する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `OpenSwathWorkflow` | OpenMSのDIA解析エンジン |
| `-mz_extraction_window` | 質量数抽出ウィンドウ（±10ppm） |
| `-Scoring:stop_report_after_feature` | 上位N個のピーク候補のみ報告 |
| `.osw` | OpenSWATH出力（SQLite形式） |

**引数の意味**:
- `-mz_extraction_window 10`: 理論値±10ppmの範囲でピーク抽出
- `-mz_extraction_window_ms1 10`: MS1（プリカーサー）の抽出ウィンドウ
- `-min_upper_edge_dist 1`: ピーク境界検出の最小距離
- `-min_peak_width 10`: 最小ピーク幅（秒）

**つまずきやすいポイント**:
- ppm値はデータ品質に応じて調整が必要（低分解能なら50ppm等）
- メモリ使用量が大きく、32サンプル同時処理では16GB以上推奨

### セクション5: PyProphet統計処理とFDR制御

```python
# OSWファイルをマージしてからまとめて統計処理
rc = run(["pyprophet", "merge", "--out", MERGED_OSW] + osw_files, "Step3 merge")

for stage in ["score", "peptide", "protein"]:
    cmd = ["pyprophet", stage, "--in", MERGED_OSW]
    if stage == "score":
        cmd += ["--level", "ms2", "--ss_initial_fdr", "0.15", "--ss_iteration_fdr", "0.05"]
    else:
        cmd += ["--context", "global"]
```

**目的**: 複数サンプルの結果を統合し、機械学習ベースのスコアリングでFDR<1%を維持しながら同定品質を向上。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `pyprophet merge` | 複数OSWファイルを統合してグローバル統計処理を可能に |
| `pyprophet score` | 機械学習ベースのピーク品質スコアリング |
| `pyprophet peptide/protein` | ペプチド・タンパク質レベルのFDR制御 |

**引数の意味**:
- `--ss_initial_fdr 0.15`: 初期FDR閾値（15%）
- `--ss_iteration_fdr 0.05`: 反復スコアリングでの目標FDR（5%）
- `--context global`: 全サンプル横断的な統計制御
- `--max_global_peptide_qvalue 0.01`: 最終ペプチドFDR<1%
- `--max_global_protein_qvalue 0.01`: 最終タンパク質FDR<1%

### セクション6: タンパク質マトリクス構築

```python
df["sample"] = df["filename"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
df["acc"] = df["ProteinName"].str.extract(r"sp\|(\w+)\|", expand=False)
df["gene"] = df["acc"].map(id_to_gene).fillna(df["acc"])

matrix = df.pivot_table(index="gene", columns="sample", values="Intensity", aggfunc="sum")
matrix = matrix.replace(0, np.nan).dropna(how="all")
```

**目的**: PyProphet出力から遺伝子シンボルごとに強度を集約し、統計解析で使用可能なタンパク質×サンプル行列を構築する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `str.extract(r"sp\|(\w+)\|")` | UniProtヘッダからアクセッション番号を抽出 |
| `pd.pivot_table()` | ロング形式からワイド形式（マトリクス）に変換 |
| `aggfunc="sum"` | 同じ遺伝子の複数ペプチド強度を合計 |
| `replace(0, np.nan)` | 0値を欠損値に変換（検出なし＝NA） |

**引数の意味**:
- `index="gene"`: 行方向に遺伝子シンボル
- `columns="sample"`: 列方向にサンプル名
- `values="Intensity"`: セルに強度値
- `dropna(how="all")`: 全サンプルで未検出のタンパク質は除外

## 【深掘り】 AlphaPeptDeepとは？

- **ひとことで**: 深層学習でペプチドのMS/MSスペクトルと保持時間を予測するツール
- **定義**: Transformerアーキテクチャを用いて、ペプチド配列から理論スペクトルライブラリを生成。実測データなしで高精度予測が可能
- **どんなとき使う**: 実測スペクトルライブラリが不足する場合や、新規サンプルでの網羅的プロテオーム解析
- **DIA-NNとの違い**: DIA-NNの深層学習予測部分を代替する商用利用可能ツール
- **参考**: Yang et al. Nature Communications 2023「AlphaPeptDeep」

## 【深掘り】 PyProphetとは？

- **ひとことで**: 機械学習を用いたDIA-MSデータの統計的スコアリングとFDR制御ツール
- **定義**: OpenSWATHの出力に対して、ランダムフォレスト等でピーク品質を学習し、偽陽性率を統計的に制御
- **どんなとき使う**: DIA-MSで高品質な同定結果を得たい場合（単純な閾値切りより高精度）
- **落とし穴**: サンプル数が少ないと機械学習の効果が限定的
- **参考**: Teleman et al. Nature Methods 2015「PyProphet」

## 【深掘り】 Match Between Runs (MBR)とは？

- **ひとことで**: サンプル間で保持時間を揃えて、片方でしか同定されなかったペプチドを補完する手法
- **定義**: 高信頼度で同定されたペプチドの保持時間を基準に、他サンプルでの同定を拡張する統計的補完
- **どんなとき使う**: サンプル間でマトリクス効果や装置コンディションが異なる場合
- **落とし穴**: 保持時間変動が大きいと偽陽性が増加する可能性
- **参考**: DIA-MSで標準的に使われる定量精度向上テクニック

## よくある質問

**Q: なぜAlphaPeptDeepを使うの？DIA-NNじゃダメ？**
A: DIA-NNは学術利用のみ無料で、商用利用には有料ライセンスが必要です。AlphaPeptDeep + OpenMSは完全に商用利用可能で、類似の深層学習予測性能を提供します。

**Q: OpenSWATHとDIA-NNの性能差は？**
A: 検出タンパク質数ではDIA-NNがやや有利ですが、OpenMS + AlphaPeptDeepでも論文の70-80%程度の検出性能は達成可能です。統計的FDR制御の信頼性ではOpenMSが優れています。

**Q: GPU環境がない場合は？**
A: CPU環境でも動作しますが、AlphaPeptDeepの予測に数時間～1日かかります。`device_type: "cpu"`設定で確実に動作します。

**Q: sage（前のステップ）との使い分けは？**
A: sageは軽量で高速、OpenMSは高精度で包括的です。リソースに余裕があり、より多くのタンパク質を検出したい場合はOpenMSを推奨します。

**Q: 自分のデータで試したい場合は？**
A: `--subset N`オプションで小規模テストから始めてください。FASTAファイルも対象生物種に変更が必要です。プロテアーゼがトリプシン以外の場合は設定変更も必要です。

## 出力されるマトリクスの形式

```
Protein,CRC001-N,CRC001-T,CRC002-N,CRC002-T,...
TP53,45234.2,67891.1,32456.7,58923.4,...
EGFR,23456.8,12789.5,45678.2,34567.9,...
MYC,12345.6,23456.7,17890.1,29834.5,...
```

- **行**: タンパク質（Gene Symbol）
- **列**: サンプル名（.mzML拡張子を除去）
- **値**: ペプチド強度の合計値（対数変換前の生データ）
- **欠損**: 検出されなかったタンパク質はNaN

## 重要な注意点

- **処理時間**: 32サンプルで2-4時間（CPU環境では8-12時間）
- **メモリ使用量**: 16GB以上推奨（大規模データセット時）
- **ディスク使用量**: 一時ファイルで20-50GB必要
- **依存関係**: OpenMS 3.0+、peptdeep 1.4+、pyprophet 2.2+が必要

## 次のステップ

このスクリプトの出力（`protein_matrix_from_openms.csv`）は以下のステップで利用されます：
- Step 5B: `step_05_overview_visualization.py`での可視化比較
- Step 6: `step_06_differential_expression.py`での差分発現解析
- OpenMS結果とsage結果の性能比較分析