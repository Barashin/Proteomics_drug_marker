#!/usr/bin/env python3
"""
OpenSWATH実行（Step 2）- 32サンプルのDIA解析
"""

import os, glob, subprocess, sys, time

def run_openswath(mzml_path, pqp_lib, output_dir, ms1_ppm, ms2_ppm):
    """OpenSWATHを1つのmzMLファイルに対して実行"""
    sample_name = os.path.splitext(os.path.basename(mzml_path))[0]
    osw_out = os.path.join(output_dir, f"{sample_name}.osw")

    if os.path.exists(osw_out):
        print(f"スキップ: {sample_name} (既に存在)")
        return osw_out

    cmd = [
        "OpenSwathWorkflow",
        "-in", mzml_path,
        "-tr", pqp_lib,
        "-out_osw", osw_out,
        "-min_upper_edge_dist", "1",
        "-mz_extraction_window", str(ms2_ppm),
        "-mz_extraction_window_unit", "ppm",
        "-mz_extraction_window_ms1", str(ms1_ppm),
        "-mz_extraction_window_ms1_unit", "ppm",
        "-use_ms1_traces",
        "-Scoring:stop_report_after_feature", "5",
        "-Scoring:TransitionGroupPicker:min_peak_width", "10",
        "-threads", "4"
    ]

    print(f"実行: {sample_name}")
    start_time = time.time()

    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in proc.stdout:
            if "Progress" in line:  # 進行状況のみ表示
                sys.stdout.write(line)
                sys.stdout.flush()
        proc.wait()

        if proc.returncode == 0:
            elapsed = (time.time() - start_time) / 60
            print(f"完了: {sample_name} ({elapsed:.1f}分)")
            return osw_out
        else:
            print(f"エラー: {sample_name} (終了コード: {proc.returncode})")
            return None
    except Exception as e:
        print(f"例外: {sample_name} - {e}")
        return None

def main():
    print("=== OpenSWATH DIA解析開始 ===")

    # パス設定
    MZML_DIR = "data/raw/raw_mzML"
    PQP_LIB = "results/openms_output/library/predicted_library.pqp"
    OSWR_DIR = "results/openms_output/openswath"

    # Toyota parameters
    MS1_PPM = 10.0
    MS2_PPM = 10.0

    os.makedirs(OSWR_DIR, exist_ok=True)

    # mzMLファイル一覧
    mzml_files = sorted(glob.glob(os.path.join(MZML_DIR, "*.mzML")))
    print(f"対象ファイル: {len(mzml_files)}個")

    # PQPライブラリ確認
    if not os.path.exists(PQP_LIB):
        print(f"エラー: PQPファイルが見つかりません: {PQP_LIB}")
        return

    pqp_size = os.path.getsize(PQP_LIB) / 1024**2
    print(f"PQPライブラリ: {pqp_size:.1f} MB")

    # 全サンプルでOpenSWATH実行
    osw_files = []
    total_start = time.time()

    for i, mzml_path in enumerate(mzml_files):
        print(f"\n[{i+1}/{len(mzml_files)}] {os.path.basename(mzml_path)}")
        osw_out = run_openswath(mzml_path, PQP_LIB, OSWR_DIR, MS1_PPM, MS2_PPM)

        if osw_out:
            osw_files.append(osw_out)
        else:
            print(f"失敗: {os.path.basename(mzml_path)}")
            # 継続する（他のファイルも処理）

    total_elapsed = (time.time() - total_start) / 60
    print(f"\n=== OpenSWATH完了 ===")
    print(f"成功: {len(osw_files)}/{len(mzml_files)}ファイル")
    print(f"総時間: {total_elapsed:.1f}分")
    print(f"OSWファイル出力: {OSWR_DIR}")

if __name__ == "__main__":
    main()