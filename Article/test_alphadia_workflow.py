#!/usr/bin/env python3
"""
AlphaDIA検証スクリプト
OpenMS + AlphaPeptDeepの代替として、Apache 2.0ライセンスのAlphaDIAを試す
"""

import os
import subprocess
import time
import glob

def setup_alphadia():
    """AlphaDIA環境セットアップ"""
    print("=== AlphaDIA環境セットアップ ===")

    commands = [
        # AlphaDIA インストール
        ["pip", "install", "alphadia"],

        # バージョン確認
        ["alphadia", "--version"]
    ]

    for cmd in commands:
        print(f"実行: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"エラー: {result.stderr}")
            return False
        print(f"成功: {result.stdout.strip()}")

    return True

def run_alphadia_workflow():
    """AlphaDIA DIA解析ワークフロー"""

    # パス設定（実際のファイル構造に合わせて修正）
    FASTA_PATH = "Article/data/raw/human_proteome.fasta"
    MZML_DIR = "data_raw/raw_mzML"
    OUTPUT_DIR = "results/alphadia_output"

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # mzMLファイル一覧（テスト用に2ファイルから開始）
    mzml_files = sorted(glob.glob(os.path.join(MZML_DIR, "*.mzML")))[:2]

    if not mzml_files:
        print(f"mzMLファイルが見つかりません: {MZML_DIR}")
        return False

    print(f"=== AlphaDIA DIA解析開始 ===")
    print(f"入力: {len(mzml_files)} mzMLファイル")
    print(f"FASTA: {FASTA_PATH}")
    print(f"出力: {OUTPUT_DIR}")

    # AlphaDIA設定
    config = {
        "fasta": os.path.abspath(FASTA_PATH),
        "output": os.path.abspath(OUTPUT_DIR),
        "precursor_mz": [400, 1200],  # Toyota論文のm/z範囲
        "precursor_charge": [2, 4],
        "missed_cleavages": 1,
        "variable_mods": {
            "M": 15.9949,  # Oxidation
        },
        "fixed_mods": {
            "C": 57.021464,  # Carbamidomethyl
        }
    }

    # AlphaDIA実行コマンド構築
    cmd = [
        "alphadia",
        "--fasta", config["fasta"],
        "--output", config["output"],
        "--library", "predict",  # 予測ライブラリ使用
        "--quant", "lfq",        # Label-free quantification
        "--threads", "4"
    ] + mzml_files

    print(f"実行コマンド: {' '.join(cmd)}")

    # 実行時間測定
    start_time = time.time()

    try:
        # AlphaDIA実行（リアルタイム出力）
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        for line in proc.stdout:
            print(line.rstrip())

        proc.wait()

        elapsed = (time.time() - start_time) / 60

        if proc.returncode == 0:
            print(f"\n✅ AlphaDIA実行完了: {elapsed:.1f}分")
            analyze_alphadia_results(OUTPUT_DIR)
            return True
        else:
            print(f"\n❌ AlphaDIA実行エラー: 終了コード {proc.returncode}")
            return False

    except Exception as e:
        print(f"❌ 実行エラー: {e}")
        return False

def analyze_alphadia_results(output_dir):
    """AlphaDIA結果の分析"""

    print("\n=== AlphaDIA結果分析 ===")

    # 主要な出力ファイルを確認
    expected_files = [
        "*.tsv",  # 定量結果
        "*.parquet", # 生データ
        "*.log"   # ログファイル
    ]

    for pattern in expected_files:
        files = glob.glob(os.path.join(output_dir, pattern))
        if files:
            for file in files:
                size_mb = os.path.getsize(file) / 1024**2
                print(f"出力ファイル: {os.path.basename(file)} ({size_mb:.1f} MB)")
        else:
            print(f"ファイルなし: {pattern}")

def compare_with_existing():
    """既存のsage結果との比較準備"""

    print("\n=== 比較分析準備 ===")

    sage_result = "results/protein_matrix_from_sage.csv"
    alphadia_result = "results/alphadia_output/protein_groups.tsv"  # 予想パス

    if os.path.exists(sage_result):
        print(f"Sage結果ファイル存在: {sage_result}")

        if os.path.exists(alphadia_result):
            print(f"AlphaDIA結果ファイル存在: {alphadia_result}")
            print("→ 次のステップで検出タンパク質数・定量値の比較が可能")
        else:
            print(f"AlphaDIA結果待ち: {alphadia_result}")
    else:
        print("まずsage解析の完了が必要")

def main():
    """メイン実行関数"""

    print("🧬 AlphaDIA DIA解析ワークフロー 🧬")
    print("Apache 2.0ライセンス・商用利用可能・深層学習ベース")
    print()

    # Step 1: セットアップ
    if not setup_alphadia():
        print("❌ セットアップ失敗")
        return

    # Step 2: DIA解析実行
    if not run_alphadia_workflow():
        print("❌ 解析失敗")
        return

    # Step 3: 比較準備
    compare_with_existing()

    print("\n🎉 AlphaDIA検証完了！")
    print("次のステップ: sage結果との比較・性能評価")

if __name__ == "__main__":
    main()