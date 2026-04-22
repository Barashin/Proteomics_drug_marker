# Code Evaluation Report

- 対象: Proteomics_drug_marker — 記事 #5（OpenMS + AlphaPeptDeep）
- 評価日: 2026-04-22
- 対象ファイル数: 3（blog 1 本 / notebook 1 本 / environment.yml 1 本）
  - `Article/blog/article-05-openms.md`（Python 20 ブロック、bash 3 ブロック）
  - `Article/notebooks/step_05.ipynb`（code cell 17 個）
  - `Article/environment.yml`
- 参考: `Article/scripts/` は記事内で `scripts/step_06_openms.py` が参照されているが **不在**（関連スクリプトは一切存在しない）

## サマリー

| 指標 | Pass | Warn | Fail |
|------|:----:|:----:|:----:|
| 1. 構文 | 37 | 0 | 0 |
| 2. インポート整合 | 1 | 1 | 0 |
| 3. 依存宣言 | 0 | 0 | 1 |
| 4. ファイル参照 | 0 | 1 | 1 |
| 5. スタイル / Lint | 1 | 1 | 0 |
| 6. 再現性 | 0 | 1 | 0 |
| 7. 外部ツール可用性 | 0 | 0 | 1 |
| 8. スモーク実行 | 1 | 0 | 1 |

**総合判定**: **Fail**（Fail 4 件 / Warn 5 件）

Fail が 1 つでもあれば code-debug に渡さず、code-implement に差し戻して修正すること。

---

## Fail 一覧（要修正）

| # | ファイル:行 | 指標 | 内容 | 提案修正 |
|---|-----------|------|------|---------|
| 1 | `Article/environment.yml`:全体 | 3 依存宣言 | OpenMS (`openms` / `OpenSwathWorkflow` / `TargetedFileConverter`)、AlphaPeptDeep (`alphapeptdeep` / `peptdeep` CLI)、PyProphet が **一切宣言されていない**。記事 109-126 行の bash ブロックで install 手順は書かれているが、`environment.yml` からは `micromamba env create -f environment.yml` だけでは #5 を再現できない。 | `bioconda::openms` を dependencies に、`alphapeptdeep`・`pyprophet` を `pip:` セクションに追記（AlphaPeptDeep が PyTorch に依存する旨のコメントも付記）。 |
| 2 | `Article/blog/article-05-openms.md`:686 | 4 ファイル参照 | `python scripts/step_06_openms.py` と記載されているが、`Article/scripts/` に **当該ファイルは存在しない**（ディレクトリ内の step_* は 02/03/04/04/04/05/06/07/08 のみ、OpenMS 用スクリプトは皆無）。かつ番号自体 #5 記事なのに `step_06_` を参照しているのも不整合。 | 2 択：(a) `Article/scripts/step_05_openms.py` を新規実装して参照パスを一致させる、(b) 当該 bash ブロックを削除し Notebook 実行のみに誘導する。まずは (b) を即時修正、次サイクルで (a) を検討。 |
| 3 | `notebooks/step_05.ipynb` 全 subprocess 呼び出し / `article-05-openms.md`:246,315,356,434,454,477,497,519 | 7 外部ツール可用性 | `subprocess.Popen(["peptdeep", ...])` / `"OpenSwathWorkflow"` / `"pyprophet"` / `"TargetedFileConverter"` の 4 ツールすべてに対し `shutil.which(...)` による事前可用性ガードが無い。未インストール環境では `FileNotFoundError: [Errno 2] No such file or directory: 'peptdeep'` 等で即クラッシュする。現に crc-proteomics 環境で `shutil.which` を叩くと 4 ツールすべて `None`（下記指標 8 参照）。 | Step 1/2/3 冒頭にヘルパー `def ensure_tool(name): if shutil.which(name) is None: raise RuntimeError(f"{name} not installed. See article §インストール")` を追加し、各ステップ先頭で `ensure_tool("peptdeep")` 等を呼ぶ。 |
| 4 | 環境全体 | 8 スモーク実行 | `crc-proteomics` 環境は存在するが、`peptdeep`, `OpenSwathWorkflow`, `pyprophet`, `TargetedFileConverter` がいずれも未インストール（`shutil.which` は全て `None`）。記事の install 手順（`micromamba install -c bioconda openms` と `pip install alphapeptdeep pyprophet`）は手動実行が必要だが、`environment.yml` 再構築で自動化されていない。加えて出力先 `Article/results/openms_output/` は未生成（当然、未実行なので）。 | Fail 1 の修正で `environment.yml` を更新し、ユーザーが `micromamba env update -f environment.yml` で一発再現できるようにする。その後 CPU での小規模スモーク（例: FASTA を 100 タンパク質にサブセット）を別 flag で走らせるコードを追加推奨。 |

---

## Warn 一覧（推奨修正）

| # | ファイル:行 | 指標 | 内容 | 提案 |
|---|-----------|------|------|------|
| 1 | `notebooks/step_05.ipynb` cell 3 / `article-05-openms.md`:138 | 2 インポート整合 | `import os, glob, re, subprocess, sys, time` を 1 行に詰め込んでいる。`glob` は mzml_files 構築（cell 4）のみ、`re` は Step 4 のみで利用のため未使用ではないが、PEP 8 違反（複数モジュールを同一行で import）。 | 1 モジュール 1 行に分割。 |
| 2 | `article-05-openms.md`:160 / notebook cell 4 | 4 ファイル参照 | `MZML_DIR = "../data/raw/raw_mzML"` の相対パスは Notebook の working dir が `Article/notebooks/` である前提。`scripts/step_06_openms.py` から呼ぶ場合（Fail 2 の修正後）は `Article/scripts/` 基準になり壊れる。`Path(__file__).resolve().parent` 基準のパス解決が無い。 | `ROOT = Path(__file__).resolve().parent.parent` を定義し、そこから `data/raw/raw_mzML` を組み立てる（notebook 版のみ相対パスのまま OK とし、script 版で切り替え）。 |
| 3 | `notebooks/step_05.ipynb` 全セル | 5 スタイル/Lint | `for line in proc.stdout:` ループは `proc.stdout` を最後まで読むが、`proc.stderr` は別スレッドで吸わないと **パイプが埋まってハングする可能性**（PIPE バッファ満杯デッドロック）。現状は `stderr=subprocess.STDOUT` で統合しているため実害は低いが、`bufsize=1`（line-buffered）の明示は無い。 | `subprocess.Popen(..., bufsize=1)` を追加、または `subprocess.run` + `capture_output=True` に置換。 |
| 4 | `notebooks/step_05.ipynb` / `article-05-openms.md` | 6 再現性 | AlphaPeptDeep の深層学習推論は内部で PyTorch を使う。`torch.manual_seed` 等のシード固定が無いため、予測ライブラリが実行ごとに微細に変動する可能性がある（推論のみなら決定的だが、念の為）。`np.random.seed` も無い。 | Step 1 実行前に `import torch; torch.manual_seed(42); np.random.seed(42)` を追加。peptdeep 自体に seed オプションがあれば YAML へ追加。 |
| 5 | `article-05-openms.md`:109-126 | 5 スタイル/Lint | インストール bash ブロックが `pip install alphapeptdeep` としているが、パッケージ名は `alphapeptdeep`（CLI は `peptdeep`）。記事では明示なく「CLI は peptdeep」と小さい括弧で書かれているのみで、読者が混乱する恐れ。 | 「`pip install alphapeptdeep` → `peptdeep --help` で確認」のように、パッケージ名と CLI 名の対応を明示。 |

---

## Pass 項目の要点

- **指標 1 構文**: 記事 Python 20 ブロック、bash 3 ブロック、notebook 17 コードセル、いずれも `ast.parse` / `py_compile` 相当で **構文エラー 0**。f-string 内の YAML 生成（`settings_content`）も正しくパース可能。
- **指標 2 インポート整合（Pass 側）**: 記事と notebook の import 文は `os, glob, re, subprocess, sys, time, pandas, numpy` で一致。PyPI / conda に存在するパッケージのみで、import 解決不能なものは無い。
- **指標 5 スタイル（Pass 側）**: 未定義名・到達不能コード・無限ループは検出されず。
- **指標 8 スモーク（Pass 側）**: `micromamba env list` で `crc-proteomics` は確認済み。`python -c "import pandas, numpy"` も pandas 3.0.2 / numpy 2.4.3 で成功。
- **記事 vs notebook のコード一致**: 関数名（`parse_fasta_gene_map`, `run_openswath`）、定数名（`Q_THRESHOLD=0.01`, `MISSED_CLEAVAGES=1`, `PEPTIDE_MIN_LEN=7` 等）、パラメータ値（`--ss_initial_fdr 0.15`、`-min_upper_edge_dist 1`、`max_fragment_charge: 2` 等）が**完全一致**。記事の最終セル群が notebook の末尾 2 セルに集約されている差分はあるが、コード文字列そのものは同一。メモリ基準の「ブログとnotebookのコード一致」は Pass。

---

## 参考：tool 可用性実測

```
$ /Users/kshinba/micromamba/envs/crc-proteomics/bin/python -c "import shutil; ..."
peptdeep: None
OpenSwathWorkflow: None
pyprophet: None
TargetedFileConverter: None
```

環境には `crc-proteomics` が存在し、pandas/numpy の import は成功するが、記事 #5 で使う 4 つの外部 CLI はいずれも **未インストール**。記事本文の install 手順（bash ブロック 0）を人間が別途実行する必要があり、`environment.yml` だけでは再現できない。

---

## 次アクション

- [ ] **Fail 1**: `environment.yml` に `openms`（bioconda）・`alphapeptdeep`・`pyprophet`（pip）を追記
- [ ] **Fail 2**: `scripts/step_05_openms.py` を新規作成するか、記事 686 行の bash ブロックを削除
- [ ] **Fail 3**: notebook と記事に `shutil.which` ガード（`ensure_tool`）を追加
- [ ] **Fail 4**: 上記修正後に `micromamba env update -f environment.yml` で再構築し、`which peptdeep OpenSwathWorkflow pyprophet TargetedFileConverter` が全て解決することを確認
- [ ] Warn 項目（import 分割、パス解決、バッファ、seed 固定、パッケージ名/CLI 名対応）は Fail 修正後に一括対応
- [ ] 再評価（code-evaluator 再実行）で Fail 0 を確認してから code-debug へ
