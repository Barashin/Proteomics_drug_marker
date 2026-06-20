#!/usr/bin/env python3
"""
全32サンプル OpenSWATH解析実行
"""

import os, glob, subprocess, sys, time

def run_openswath(mzml_path, pqp_lib, output_dir, ms1_ppm, ms2_ppm):
    """OpenSWATH実行"""
    sample_name = os.path.splitext(os.path.basename(mzml_path))[0]
    osw_out = os.path.join(output_dir, f"{sample_name}.osw")

    if os.path.exists(osw_out):
        osw_size = os.path.getsize(osw_out) / 1024**2
        print(f"⏭️  スキップ: {sample_name} ({osw_size:.1f} MB)")
        return osw_out

    cmd = [
        "OpenSwathWorkflow",
        "-in", mzml_path,
        "-tr", pqp_lib,
        "-out_osw", osw_out,
        "-mz_extraction_window", str(ms2_ppm),
        "-mz_extraction_window_ms1", str(ms1_ppm),
        "-threads", "4"
    ]

    start_time = time.time()

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)

        if proc.returncode == 0:
            elapsed = (time.time() - start_time) / 60
            osw_size = os.path.getsize(osw_out) / 1024**2
            print(f"✅ {sample_name}: {elapsed:.1f}分 ({osw_size:.1f} MB)")
            return osw_out
        else:
            print(f"❌ {sample_name}: エラー{proc.returncode}")
            return None
    except Exception as e:
        print(f"❌ {sample_name}: {e}")
        return None

def main():
    print("=== 全32サンプル OpenSWATH解析 ===")
    print("Toyota et al. 2025論文再現 - GPU環境実行")

    # パス設定
    MZML_DIR = "data/raw/raw_mzML"
    PQP_LIB = "results/openms_output/library/predicted_library.pqp"
    OSWR_DIR = "results/openms_output/openswath"

    # Toyota parameters
    MS1_PPM = 10.0
    MS2_PPM = 10.0

    os.makedirs(OSWR_DIR, exist_ok=True)

    # 全32サンプル
    mzml_files = sorted(glob.glob(os.path.join(MZML_DIR, "*.mzML")))
    print(f"📁 対象ファイル: {len(mzml_files)}個")

    pqp_size = os.path.getsize(PQP_LIB) / 1024**2
    print(f"📚 PQPライブラリ: {pqp_size:.1f} MB")
    print(f"🎯 MS精度: {MS1_PPM}/{MS2_PPM} ppm")
    print()

    # 全サンプル実行
    osw_files = []
    total_start = time.time()

    for i, mzml_path in enumerate(mzml_files):
        sample_name = os.path.basename(mzml_path)
        print(f"[{i+1:2d}/32] {sample_name}", end=" ")
        sys.stdout.flush()

        osw_out = run_openswath(mzml_path, PQP_LIB, OSWR_DIR, MS1_PPM, MS2_PPM)
        if osw_out:
            osw_files.append(osw_out)

    total_elapsed = (time.time() - total_start) / 60
    print(f"\n=== OpenSWATH完了 ===")
    print(f"✅ 成功: {len(osw_files)}/32 ファイル")
    print(f"⏱️ 総時間: {total_elapsed:.1f}分")
    print(f"📊 出力: {OSWR_DIR}")

    if len(osw_files) == 32:
        print("🎉 全サンプル解析完了！PyProphetステップに進めます。")
    else:
        print(f"⚠️ {32 - len(osw_files)}個のファイルが失敗しました。")

    return osw_files

if __name__ == "__main__":
    osw_files = main()