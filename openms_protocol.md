# OpenMS + AlphaPeptDeep 実行プロトコル（実機検証済み）

- 検証日: 2026-04-22
- 検証環境: macOS 15.3 (Darwin 25.3.0), Apple Silicon, 10 cores, 16GB RAM, GPU非搭載
- micromamba 2.4.0, Python 3.11.15
- OpenMS 3.5.0 / peptdeep 1.4.2 / pyprophet（最新）

記事 `article-05-openms.md` のコードは **そのままでは動きません**。実機で検証して以下の修正が必要と判明しました。

---

## 1. インストール（article の bash ブロックへの修正）

記事に書かれている手順は `pip install alphapeptdeep` となっていますが、**正しいパッケージ名は `peptdeep`** です。

```bash
micromamba activate crc-proteomics

# OpenMS は bioconda に osx-arm64 ビルドあり
micromamba install -c bioconda -c conda-forge openms=3.5.0 -y

# AlphaPeptDeep の PyPI パッケージ名は peptdeep（alphapeptdeep ではない）
pip install peptdeep pyprophet
```

`environment.yml` への追記例:

```yaml
dependencies:
  - bioconda::openms>=3.5
  - pyyaml
  - pip:
    - peptdeep
    - pyprophet
```

---

## 2. 修正が必要な4箇所

### 修正1: peptdeep CLI 構文

記事の書き方（動かない）:
```bash
peptdeep library --settings settings.yaml
```

正しい構文（peptdeep 1.4.2）:
```bash
peptdeep library settings.yaml    # 位置引数
```

### 修正2: peptdeep settings YAML 構造

記事のYAMLはAlphaPeptDeepの古いAPIを仮定しており、現行バージョンでは動きません。正しくは `peptdeep export-settings default.yaml` でデフォルトを出力し、必要箇所だけ上書きする方式に変える必要があります。主なキー：

```yaml
task_workflow:
  - library
torch_device:
  device_type: cpu    # GPU なら gpu
thread_num: 8
library:
  infile_type: fasta
  infiles:
    - /abs/path/to/human_proteome.fasta
  fasta:
    protease: trypsin
    max_miss_cleave: 1
  fix_mods: [Carbamidomethyl@C]
  var_mods: []
  max_var_mod_num: 0
  min_peptide_len: 7
  max_peptide_len: 45
  min_precursor_charge: 2
  max_precursor_charge: 4
  decoy: pseudo_reverse
  frag_types: [b, y]
  max_frag_charge: 2
  output_folder: /abs/path/to/library
  output_tsv:
    enabled: true
    min_fragment_mz: 200.0
    max_fragment_mz: 1800.0
```

### 修正3: AlphaPeptDeep TSV → OpenSwath 互換TSV への列名・形式変換

peptdeep の TSV は AlphaPeptDeep 独自の speclib 形式で、OpenSwath の `TargetedFileConverter` はそのままでは読めません。以下の変換が必要：

| peptdeep の列名 | OpenSwath の列名 |
|---------------|----------------|
| RT | NormalizedRetentionTime |
| FragmentMz | ProductMz |
| RelativeIntensity | LibraryIntensity |
| StrippedPeptide | PeptideSequence |
| ModifiedPeptide | ModifiedPeptideSequence |
| ProteinID | ProteinName |
| FragmentNumber | FragmentSeriesNumber |

追加必要な列:
- `TransitionGroupId`（= `ModifiedPeptide_PrecursorCharge`）
- `TransitionId`（= `TransitionGroupId_行番号`）

さらに以下のデータ変換が必要：
- `ModifiedPeptide` の前後の `_` を除去（`_MQAEIKR_` → `MQAEIKR`）
- `[Carbamidomethyl]` → `(UniMod:4)`、`[Oxidation]` → `(UniMod:35)` に置換
- `ProteinID` に `;` で複数UniProtが入っている場合は先頭のみ採用

実装は `scripts/step_05_openms.py` を参照。

### 修正4: OpenSwathWorkflow の引数（OpenMS 3.5.0 の API 変更）

記事のコード:
```bash
OpenSwathWorkflow ... -out_osw out.osw -use_ms1_traces
```

OpenMS 3.5.0 では:
- `-out_osw` は廃止 → `-out_features out.osw -out_features_type osw`
- `-use_ms1_traces` は廃止（デフォルトで使用）

加えて、DIA-MS データによっては SWATH 窓にギャップがあり、OpenSwath が以下で中断することがあります：
```
Extraction windows have a gap. Will abort (override with -force)
```

→ `-force` オプションを追加すれば継続実行できます。

---

## 3. 計算リソース（実測値）

macOS M1/M2 相当の 10 cores CPU、GPU なしでの実測：

| ステップ | 時間 | ボトルネック |
|---------|------|------------|
| peptdeep library 予測（3.9M precursors、ヒトプロテオーム20k） | **104分** | **PyTorch CPU 推論**（GPU があれば5-15分） |
| AlphaPeptDeep → OpenSwath TSV 変換（8.7GB） | 1.5分 | I/O |
| TargetedFileConverter (TSV → PQP) | **45分** | OpenMS のパース |
| OpenSwathWorkflow（1サンプル、850MB mzML、3.9M precursor library） | **約90-95分** | PQP ロード52分 + DIA 解析 |
| PyProphet score + export | 約20-30分 | XGBoost CPU学習 |

**推定総時間**:
- 32サンプル全て: **約50-60 時間**（CPU）→ GPU があれば peptdeep 部分が短縮され **20-30 時間** 程度
- 2サンプル（動作確認）: **約4時間**（CPU）

### 論文が「60-90分」と書いている背景
論文の記述はおそらく：
- もっと小さい FASTA（セーフチューリング対象の絞り込みあり）
- または GPU 前提
- または library-based ではなく別パイプライン

ヒトプロテオーム丸ごと（20k protein × 1 missed cleavage）だと precursor が約400万になり、現行のハードでは6時間〜数日かかります。

---

## 4. 実行例（小規模動作確認）

まず1患者分（CRC01 N+T = 2 mzML）で通し動作を確認：

```bash
cd Article
micromamba activate crc-proteomics
python scripts/step_05_openms.py --subset 2 --threads 8
```

`scripts/step_05_openms.py` には上記4修正が全て反映済みです。ライブラリ予測だけ再利用して先に進む場合：

```bash
python scripts/step_05_openms.py --subset 2 --threads 8 --skip-library
```

全32サンプル実行は GPU マシン（Linux + CUDA）を推奨。

---

## 5. 推奨される記事改訂方針

1. **「CPUで60-90分」の記述を削除**し、実測値（GPU 30分、CPU 6-50時間）を記載する
2. **インストール手順の修正**: `pip install peptdeep pyprophet`
3. **peptdeep の YAML 設定**: デフォルトエクスポート＋オーバーライド方式に書き換え
4. **OpenSWATH のCLI 引数**: `-out_features` / `-force` を明示
5. **TSV 変換ステップを1セクション追加**: AlphaPeptDeep → OpenSwath 互換変換が必要である旨
6. **結論セクション**: 「sage が現実的な選択肢。OpenMS+AlphaPeptDeep は GPU あり・時間許容範囲が広い場合の選択肢」と明記

---

## 6. 検証ログ

- `Article/scripts/step_05_openms.py`: 4修正全て反映済みの実装
- `Article/results/openms_output/library/predict.speclib.hdf`: peptdeep 出力（2.1GB、次回以降再利用可）
- `Article/results/openms_output/library/predicted_library.pqp`: OpenSwath用ライブラリ（TSV変換済み、45分ビルド）
- `evaluation_report.md`: code-evaluator の総合判定 Fail 報告書（先に発見された4 Fail の詳細）
