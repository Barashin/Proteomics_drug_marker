#!/usr/bin/env python3
"""
OpenSWATH緩いパラメータでテスト実行
10 ppmが厳しすぎる可能性があるため、より緩い条件でテスト
"""

import os, subprocess, time

def run_openswath_test(mzml_path, pqp_lib, output_dir, ms1_ppm, ms2_ppm):
    """OpenSWATHテスト実行（緩いパラメータ）"""
    sample_name = os.path.splitext(os.path.basename(mzml_path))[0]
    osw_out = os.path.join(output_dir, f"{sample_name}_relaxed.osw")

    # より緩いパラメータでテスト
    cmd = [
        "OpenSwathWorkflow",
        "-in", mzml_path,
        "-tr", pqp_lib,
        "-out_osw", osw_out,
        "-mz_extraction_window", str(ms2_ppm),
        "-mz_extraction_window_unit", "ppm",
        "-mz_extraction_window_ms1", str(ms1_ppm),
        "-mz_extraction_window_ms1_unit", "ppm",
        "-threads", "2"
    ]

    print(f"=== OpenSWATHテスト: {sample_name} ===")
    print(f"パラメータ: MS1={ms1_ppm}ppm, MS2={ms2_ppm}ppm")
    print(f"コマンド: {' '.join(cmd[:8])}...")

    start_time = time.time()

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        elapsed = (time.time() - start_time) / 60

        if proc.returncode == 0:
            print(f"✅ 成功: {elapsed:.1f}分")

            # OSWファイルサイズ確認
            if os.path.exists(osw_out):
                osw_size = os.path.getsize(osw_out) / 1024**2
                print(f"OSWファイル: {osw_size:.1f} MB")
                return osw_out
            else:
                print("❌ OSWファイル未生成")
                return None
        else:
            print(f"❌ エラー: 終了コード {proc.returncode}")
            print(f"エラー詳細: {proc.stderr[:500]}")
            return None

    except Exception as e:
        print(f"❌ 例外: {e}")
        return None

def main():
    print("=== OpenSWATH緩いパラメータテスト ===")

    # パス設定
    MZML_DIR = "data/raw/raw_mzML"
    PQP_LIB = "results/openms_output/library/predicted_library.pqp"
    TEST_DIR = "results/openms_output/test_relaxed"

    os.makedirs(TEST_DIR, exist_ok=True)

    # テスト用：最初の1サンプルのみ
    mzml_files = [os.path.join(MZML_DIR, "CRC01-N.mzML")]

    # パラメータテスト
    test_params = [
        (50, 50),   # 50 ppm (緩い)
        (100, 100), # 100 ppm (非常に緩い)
        (20, 20),   # 20 ppm (中間)
    ]

    for ms1_ppm, ms2_ppm in test_params:
        print(f"\n{'='*60}")
        print(f"テスト: MS1={ms1_ppm}ppm, MS2={ms2_ppm}ppm")

        for mzml_path in mzml_files:
            osw_out = run_openswath_test(mzml_path, PQP_LIB, TEST_DIR, ms1_ppm, ms2_ppm)

            if osw_out:
                # 簡単なFEATURE確認
                import sqlite3
                conn = sqlite3.connect(osw_out)
                cursor = conn.cursor()
                cursor.execute("SELECT COUNT(*) FROM FEATURE;")
                feature_count = cursor.fetchone()[0]
                print(f"FEATURE検出数: {feature_count}")
                conn.close()

                if feature_count > 0:
                    print("🎉 FEATUREが検出されました！このパラメータを使用できます")
                    return ms1_ppm, ms2_ppm

        print()

    print("❌ いずれのパラメータでもFEATUREが検出されませんでした")
    return None

if __name__ == "__main__":
    result = main()