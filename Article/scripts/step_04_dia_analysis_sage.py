#!/usr/bin/env python3
"""
Step 4: sage によるDIA解析本体
================================
mzMLファイルと human_proteome.fasta を入力として、sage-proteomics で
library-free DIA解析を実行します。sage は Rust 製の超高速検索エンジンで、
MIT ライセンスのため商用利用に制約がありません。

【論文との対応】
  論文 (Toyota et al. 2025) では DIA-NN v1.8.1 を使用していますが、
  DIA-NN は商用利用に制約があるため、本書では sage で代替しています。
  検索パラメータ（Trypsin、missed cleavage 1、peptide 7-45 AA、precursor m/z 495-865、
  fragment m/z 200-1800、MS1/MS2 tolerance 10 ppm、FDR < 1%）は論文とほぼ同じに設定。

【処理の流れ】
  1. sage_config.json を読み込み
  2. data/raw/raw_mzML/*.mzML を列挙
  3. sage CLI を呼び出してDIA解析を実行
  4. 出力を results/sage_output/ に保存
     - results.sage.tsv : PSM テーブル（FDR付き）
     - lfq.tsv          : ペプチドLFQ定量テーブル
     - results.json     : 実行時の全パラメータ

入力: data/raw/raw_mzML/*.mzML, data/raw/human_proteome.fasta
出力: results/sage_output/{results.sage.tsv, lfq.tsv, results.json}
"""

# ============================================================
# ライブラリの読み込み
# ============================================================
import os
import sys
import glob
import shutil
import subprocess
import time

# ============================================================
# パス設定
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, "..")
CONFIG_PATH = os.path.join(SCRIPT_DIR, "sage_config.json")
MZML_DIR = os.path.join(PROJECT_DIR, "data", "raw", "raw_mzML")
FASTA_PATH = os.path.join(PROJECT_DIR, "data", "raw", "human_proteome.fasta")
OUT_DIR = os.path.join(PROJECT_DIR, "results", "sage_output")
LOG_DIR = os.path.join(PROJECT_DIR, "results", "logs")


# ============================================================
# メイン処理
# ============================================================
def main():
    # sage コマンドが PATH にあるか確認
    sage_bin = shutil.which("sage")
    if sage_bin is None:
        raise SystemExit(
            "sage コマンドが見つかりません。micromamba run -n crc-proteomics python ... "
            "のように crc-proteomics 環境内で実行してください。"
        )
    print(f"sage binary: {sage_bin}")

    # 入力ファイルの存在確認
    if not os.path.exists(FASTA_PATH):
        raise SystemExit(f"FASTA not found: {FASTA_PATH}")
    mzml_files = sorted(glob.glob(os.path.join(MZML_DIR, "*.mzML")))
    if not mzml_files:
        raise SystemExit(f"No mzML files found in {MZML_DIR}")
    print(f"FASTA: {FASTA_PATH}")
    print(f"mzML files: {len(mzml_files)}")
    for f in mzml_files:
        print(f"  - {os.path.basename(f)}")

    # 出力ディレクトリを作成
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)

    # sage コマンドを組み立て
    cmd = [
        sage_bin,
        "--fasta", FASTA_PATH,
        "--output_directory", OUT_DIR,
        "--disable-telemetry-i-dont-want-to-improve-sage",
        CONFIG_PATH,
    ] + mzml_files

    print("\n===== Running sage =====")
    print(" ".join(cmd[:6] + ["...", f"({len(mzml_files)} mzML files)"]))

    t0 = time.time()
    # stdout/stderr を受け取りつつリアルタイム表示
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    try:
        for line in proc.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
    finally:
        rc = proc.wait()
    elapsed = time.time() - t0

    if rc != 0:
        raise SystemExit(f"sage exited with code {rc}")
    print(f"\n===== sage finished in {elapsed/60:.1f} min =====")
    print(f"Output directory: {OUT_DIR}")


if __name__ == "__main__":
    main()
