#!/usr/bin/env python3
"""
修正版 OpenSWATH実行 - 正しいパラメータで32サンプル解析
"""

import os, glob, subprocess, sys, time

def run_openswath(mzml_path, pqp_lib, output_dir, ms1_ppm, ms2_ppm):
    """修正版OpenSWATH実行"""
    sample_name = os.path.splitext(os.path.basename(mzml_path))[0]
    osw_out = os.path.join(output_dir, f"{sample_name}.osw")

    if os.path.exists(osw_out):
        print(f"スキップ: {sample_name} (既に存在)")
        return osw_out

    # OpenMS 3.1.0対応の正しいパラメータ
    cmd = [
        "OpenSwathWorkflow",
        "-in", mzml_path,
        "-tr", pqp_lib,
        "-out_osw", osw_out,
        "-mz_extraction_window", str(ms2_ppm),
        "-mz_extraction_window_ms1", str(ms1_ppm),
        "-threads", "2"
    ]

    print(f"実行: {sample_name}")
    start_time = time.time()

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)

        if proc.returncode == 0:
            elapsed = (time.time() - start_time) / 60
            print(f"✅ 完了: {sample_name} ({elapsed:.1f}分)")
            return osw_out
        else:
            print(f"❌ エラー: {sample_name} (終了コード: {proc.returncode})")
            print(f"エラー詳細: {proc.stderr[:200]}...")
            return None
    except Exception as e:
        print(f"❌ 例外: {sample_name} - {e}")
        return None

def main():
    print("=== 修正版 OpenSWATH DIA解析 ===")

    # パス設定
    MZML_DIR = "data/raw/raw_mzML"
    PQP_LIB = "results/openms_output/library/predicted_library.pqp"
    OSWR_DIR = "results/openms_output/openswath"

    # Toyota parameters (ppm単位)
    MS1_PPM = 10.0
    MS2_PPM = 10.0

    os.makedirs(OSWR_DIR, exist_ok=True)

    # mzMLファイル一覧（最初の3個でテスト）
    mzml_files = sorted(glob.glob(os.path.join(MZML_DIR, "*.mzML")))[:3]
    print(f"テスト対象: {len(mzml_files)}個")

    # PQPライブラリ確認
    if not os.path.exists(PQP_LIB):
        print(f"❌ エラー: PQPファイルが見つかりません: {PQP_LIB}")
        return

    pqp_size = os.path.getsize(PQP_LIB) / 1024**2
    print(f"📚 PQPライブラリ: {pqp_size:.1f} MB")
    print(f"🎯 MS精度: {MS1_PPM}/{MS2_PPM} ppm")

    # 最初の3サンプルでテスト実行
    osw_files = []
    total_start = time.time()

    for i, mzml_path in enumerate(mzml_files):
        print(f"\n[{i+1}/{len(mzml_files)}] {os.path.basename(mzml_path)}")
        osw_out = run_openswath(mzml_path, PQP_LIB, OSWR_DIR, MS1_PPM, MS2_PPM)

        if osw_out:
            osw_files.append(osw_out)

    total_elapsed = (time.time() - total_start) / 60
    print(f"\n=== テスト完了 ===")
    print(f"✅ 成功: {len(osw_files)}/{len(mzml_files)}ファイル")
    print(f"⏱️ 総時間: {total_elapsed:.1f}分")

    if len(osw_files) > 0:
        print("🎉 OpenSWATH動作確認完了！")
        print("全32サンプルの実行準備ができました。")
    else:
        print("❌ OpenSWATH設定に問題があります。")

if __name__ == "__main__":
    main()