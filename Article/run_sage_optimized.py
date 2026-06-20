#!/usr/bin/env python3
"""
Sage最適化版実行スクリプト
Toyota et al. 2025の検出数に近づけるための最適化パラメータでsageを実行
"""

import os
import subprocess
import time
import glob
import pandas as pd

def run_sage_optimized():
    """最適化パラメータでsageを実行"""

    # パス設定
    CONFIG_FILE = "scripts/sage_config_optimized.json"
    MZML_DIR = "../data_raw/raw_mzML"
    RESULTS_DIR = "results/sage_optimized"

    # 出力ディレクトリ作成
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # mzMLファイル一覧取得
    mzml_files = glob.glob(os.path.join(MZML_DIR, "*.mzML"))
    mzml_files.sort()

    print(f"=== Sage最適化版実行 ===")
    print(f"設定ファイル: {CONFIG_FILE}")
    print(f"入力ファイル: {len(mzml_files)} mzML")
    print(f"出力先: {RESULTS_DIR}")

    # sageコマンド構築
    cmd = ["sage", CONFIG_FILE] + mzml_files

    print(f"\\nコマンド: {' '.join(cmd[:3])} ... (+ {len(mzml_files)} mzML files)")
    print("実行開始...")

    # sage実行
    start_time = time.time()

    try:
        # リアルタイム出力でsage実行
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                 universal_newlines=True, bufsize=1)

        # 出力を1行ずつリアルタイム表示
        for line in process.stdout:
            print(line.rstrip())

        process.wait()

        # 実行時間計算
        elapsed_time = time.time() - start_time

        if process.returncode == 0:
            print(f"\\n✅ Sage実行完了: {elapsed_time/60:.1f}分")
            analyze_results()
        else:
            print(f"\\n❌ Sage実行エラー: 終了コード {process.returncode}")
            return False

    except Exception as e:
        print(f"❌ 実行中にエラー: {e}")
        return False

    return True

def analyze_results():
    """結果の分析と比較"""

    RESULTS_DIR = "results/sage_optimized"
    ORIGINAL_DIR = "results/sage_output"

    print("\\n=== 結果分析 ===")

    try:
        # 最適化版の結果確認
        optimized_files = glob.glob(os.path.join(RESULTS_DIR, "*.tsv"))

        if optimized_files:
            # タンパク質数の比較
            quant_file = None
            for f in optimized_files:
                if "quant" in f and "proteins" in f:
                    quant_file = f
                    break

            if quant_file and os.path.exists(quant_file):
                df_opt = pd.read_csv(quant_file, sep='\\t')
                protein_count_opt = len(df_opt)
                print(f"最適化版検出タンパク質数: {protein_count_opt}")

                # 元のsage結果と比較
                original_quant = os.path.join(ORIGINAL_DIR, "quant.proteins.tsv")
                if os.path.exists(original_quant):
                    df_orig = pd.read_csv(original_quant, sep='\\t')
                    protein_count_orig = len(df_orig)
                    improvement = ((protein_count_opt - protein_count_orig) / protein_count_orig) * 100

                    print(f"元のsage検出数: {protein_count_orig}")
                    print(f"改善率: {improvement:+.1f}% ({protein_count_opt - protein_count_orig:+d})")

                    # Toyota論文との比較
                    toyota_proteins = 10329
                    coverage = (protein_count_opt / toyota_proteins) * 100
                    print(f"Toyota et al.カバー率: {coverage:.1f}% ({protein_count_opt}/{toyota_proteins})")

                # プロテインマトリクス作成
                create_protein_matrix(df_opt)

            else:
                print("❌ タンパク質定量ファイルが見つかりません")
        else:
            print("❌ 結果ファイルが見つかりません")

    except Exception as e:
        print(f"❌ 結果分析エラー: {e}")

def create_protein_matrix(df):
    """タンパク質×サンプルマトリクスを作成"""

    print("\\nタンパク質マトリクス作成中...")

    try:
        # サンプル列を特定（CRC##-Nまたは-T パターン）
        sample_cols = [col for col in df.columns if 'CRC' in col and ('-N' in col or '-T' in col)]

        if sample_cols:
            # プロテインID + サンプル強度値のマトリクス作成
            protein_matrix = df[['protein_group'] + sample_cols].copy()
            protein_matrix.set_index('protein_group', inplace=True)

            # マトリクス保存
            output_file = "results/protein_matrix_from_sage_optimized.csv"
            protein_matrix.to_csv(output_file)

            print(f"✅ プロテインマトリクス保存: {output_file}")
            print(f"   サイズ: {protein_matrix.shape[0]} proteins × {protein_matrix.shape[1]} samples")

            return output_file
        else:
            print("❌ サンプル列が見つかりません")

    except Exception as e:
        print(f"❌ マトリクス作成エラー: {e}")

    return None

if __name__ == "__main__":
    print("Sage最適化版によるDIA解析")
    print("目標: Toyota et al. 2025の検出数に近づける")
    print("=" * 50)

    success = run_sage_optimized()

    if success:
        print("\\n🎉 最適化版sage実行完了!")
        print("次のステップ: 結果を article-06-preprocess.md で前処理")
    else:
        print("\\n💥 実行失敗。設定を確認してください。")