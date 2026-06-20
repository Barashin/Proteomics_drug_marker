#!/usr/bin/env python3
"""
OpenMS pipeline - Batched execution to prevent WSL crashes
"""
import os
import sys
import glob
import time
import subprocess
from pathlib import Path

# 設定
BATCH_SIZE = 4  # 同時処理サンプル数を4に制限
DATA_DIR = "data/raw/raw_mzML"
PQP_LIB = "results/openms_output/library/predicted_library.pqp"
OSW_DIR = "results/openms_output/openswath"

def run_openswath_batch(mzml_files, batch_num, total_batches):
    """単一バッチのOpenSWATH処理"""
    print(f"\n=== Batch {batch_num}/{total_batches} ===")
    print(f"Processing {len(mzml_files)} samples")

    for i, mzml_path in enumerate(mzml_files):
        sample = os.path.splitext(os.path.basename(mzml_path))[0]
        osw_out = os.path.join(OSW_DIR, f"{sample}.osw")

        # 既存ファイルスキップ
        if os.path.exists(osw_out) and os.path.getsize(osw_out) > 1024:
            print(f"[{i+1}/{len(mzml_files)}] SKIP {sample} (existing)")
            continue

        print(f"[{i+1}/{len(mzml_files)}] Processing {sample}...")
        cmd = [
            "OpenSwathWorkflow",
            "-in", mzml_path,
            "-tr", PQP_LIB,
            "-out_osw", osw_out,
            "-min_upper_edge_dist", "1",
            "-mz_extraction_window", "10.0",
            "-mz_extraction_window_unit", "ppm",
            "-mz_extraction_window_ms1", "10.0",
            "-mz_extraction_window_ms1_unit", "ppm",
            "-Scoring:stop_report_after_feature", "5",
            "-Scoring:TransitionGroupPicker:min_peak_width", "10",
            "-threads", "2",  # スレッド数を制限
            "-force",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"ERROR in {sample}: {result.stderr}")
            return False

        # メモリクリア待機
        time.sleep(2)

    print(f"Batch {batch_num} completed successfully")
    return True

def main():
    # mzMLファイル一覧
    mzml_files = sorted(glob.glob(f"{DATA_DIR}/*.mzML"))
    print(f"Found {len(mzml_files)} mzML files")

    # バッチ分割
    batches = [mzml_files[i:i+BATCH_SIZE] for i in range(0, len(mzml_files), BATCH_SIZE)]
    print(f"Divided into {len(batches)} batches of max {BATCH_SIZE} samples")

    # バッチ実行
    for i, batch in enumerate(batches, 1):
        success = run_openswath_batch(batch, i, len(batches))
        if not success:
            print(f"Batch {i} failed. Stopping.")
            sys.exit(1)

        # バッチ間休憩
        if i < len(batches):
            print(f"Waiting 10 seconds before next batch...")
            time.sleep(10)

    print("\nAll batches completed successfully!")

if __name__ == "__main__":
    main()