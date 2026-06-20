#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
完全版 OpenMS + AlphaPeptDeep パイプライン
GPU環境（RTX 4070）で実行 - Toyota et al. 2025論文再現
"""

import os, glob, re, subprocess, sys, time
import pandas as pd
import numpy as np

def wait_for_library(library_path, check_interval=30, max_wait=3600):
    """AlphaPeptDeepライブラリ生成の完了を待機"""
    print(f"ライブラリ生成を待機中: {library_path}")
    start_time = time.time()

    while time.time() - start_time < max_wait:
        if os.path.exists(library_path) and os.path.getsize(library_path) > 1000000:  # 1MB以上
            print("ライブラリ生成完了!")
            return True
        print(f"待機中... ({(time.time() - start_time)/60:.1f}分経過)")
        time.sleep(check_interval)

    print("タイムアウト: ライブラリ生成が完了しません")
    return False

def run_command(cmd, description):
    """コマンドを実行して結果を表示"""
    print(f"\n=== {description} ===")
    print(f"実行: {' '.join(cmd)}")

    start_time = time.time()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    for line in proc.stdout:
        sys.stdout.write(line)
        sys.stdout.flush()

    proc.wait()
    elapsed = (time.time() - start_time) / 60

    print(f"完了: {elapsed:.1f}分 (終了コード: {proc.returncode})")
    return proc.returncode == 0

def run_openswath(mzml_path, pqp_lib, output_dir, ms1_ppm, ms2_ppm):
    """OpenSWATHを1つのmzMLファイルに対して実行"""
    sample_name = os.path.splitext(os.path.basename(mzml_path))[0]
    osw_out = os.path.join(output_dir, f"{sample_name}.osw")

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

    success = run_command(cmd, f"OpenSWATH: {sample_name}")
    return osw_out if success else None

def parse_fasta_gene_map(fasta_path):
    """FASTAからGene Symbolマッピングを作成"""
    id_to_gene = {}
    gn_re = re.compile(r"\bGN=(\S+)")

    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue

            parts = line[1:].split("|", 2)
            if len(parts) >= 3:
                acc = parts[1]
                entry_name = parts[2].split()[0]
            else:
                acc = line[1:].split()[0]
                entry_name = acc

            m = gn_re.search(line)
            id_to_gene[acc] = m.group(1) if m else entry_name.split("_")[0]

    return id_to_gene

def main():
    print("=== 完全版 OpenMS + AlphaPeptDeep パイプライン ===")
    print("GPU環境（RTX 4070）+ Toyota et al. 2025パラメータ")

    # --- パス設定 ---
    MZML_DIR = "data/raw/raw_mzML"
    FASTA_PATH = "data/raw/human_proteome.fasta"
    RESULTS_DIR = "results"

    OPENMS_OUT = os.path.join(RESULTS_DIR, "openms_output")
    LIBRARY_DIR = os.path.join(OPENMS_OUT, "library")
    OSWR_DIR = os.path.join(OPENMS_OUT, "openswath")
    PYPROPHET_DIR = os.path.join(OPENMS_OUT, "pyprophet")

    for d in [OPENMS_OUT, LIBRARY_DIR, OSWR_DIR, PYPROPHET_DIR]:
        os.makedirs(d, exist_ok=True)

    # mzMLファイル一覧
    mzml_files = sorted(glob.glob(os.path.join(MZML_DIR, "*.mzML")))
    print(f"mzMLファイル: {len(mzml_files)}個")

    # --- Toyota et al. 2025パラメータ ---
    Q_THRESHOLD = 0.01
    MS1_PPM = 10.0
    MS2_PPM = 10.0

    # === Step 1: AlphaPeptDeep完了待ち ===
    print("\n=== Step 1: AlphaPeptDeep ライブラリ確認 ===")
    PREDICTED_LIB = os.path.join(LIBRARY_DIR, "predicted_library_toyota_params.tsv")

    if not wait_for_library(PREDICTED_LIB):
        print("エラー: AlphaPeptDeepライブラリが生成されませんでした")
        return

    lib_size = os.path.getsize(PREDICTED_LIB) / 1024**2
    print(f"ライブラリサイズ: {lib_size:.1f} MB")

    # === Step 2: OpenSWATH ===
    print("\n=== Step 2: OpenSWATH DIA検索 ===")

    # TSV → PQP変換
    PQP_LIB = os.path.join(LIBRARY_DIR, "predicted_library.pqp")
    if not os.path.exists(PQP_LIB):
        cmd_convert = [
            "TargetedFileConverter",
            "-in", PREDICTED_LIB,
            "-out", PQP_LIB
        ]
        if not run_command(cmd_convert, "TSV→PQP変換"):
            print("エラー: PQP変換に失敗")
            return

    # 全mzMLファイルに対してOpenSWATH実行
    print(f"OpenSWATH実行: {len(mzml_files)}ファイル")
    osw_files = []

    for i, mzml_path in enumerate(mzml_files):
        sample_name = os.path.basename(mzml_path)
        print(f"\n[{i+1}/{len(mzml_files)}] {sample_name}")

        osw_out = run_openswath(mzml_path, PQP_LIB, OSWR_DIR, MS1_PPM, MS2_PPM)
        if osw_out:
            osw_files.append(osw_out)
        else:
            print(f"エラー: {sample_name} の処理に失敗")
            return

    print(f"OpenSWATH完了: {len(osw_files)}ファイル")

    # === Step 3: PyProphet ===
    print("\n=== Step 3: PyProphet FDR推定 ===")

    MERGED_OSW = os.path.join(PYPROPHET_DIR, "merged.osw")

    # マージ
    cmd_merge = ["pyprophet", "merge", "--out", MERGED_OSW] + osw_files
    if not run_command(cmd_merge, "PyProphet merge"):
        return

    # MS2スコアリング
    cmd_score = [
        "pyprophet", "score",
        "--in", MERGED_OSW,
        "--level", "ms2",
        "--ss_initial_fdr", "0.15",
        "--ss_iteration_fdr", "0.05"
    ]
    if not run_command(cmd_score, "PyProphet MS2スコアリング"):
        return

    # ペプチドレベルFDR
    cmd_peptide = [
        "pyprophet", "peptide",
        "--in", MERGED_OSW,
        "--context", "global"
    ]
    if not run_command(cmd_peptide, "PyProphet ペプチドFDR"):
        return

    # タンパク質レベルFDR
    cmd_protein = [
        "pyprophet", "protein",
        "--in", MERGED_OSW,
        "--context", "global"
    ]
    if not run_command(cmd_protein, "PyProphet タンパク質FDR"):
        return

    # エクスポート
    EXPORT_TSV = os.path.join(PYPROPHET_DIR, "pyprophet_export.tsv")
    cmd_export = [
        "pyprophet", "export",
        "--in", MERGED_OSW,
        "--out", EXPORT_TSV,
        "--max_global_peptide_qvalue", str(Q_THRESHOLD),
        "--max_global_protein_qvalue", str(Q_THRESHOLD)
    ]
    if not run_command(cmd_export, "PyProphet エクスポート"):
        return

    # === Step 4: タンパク質マトリクス構築 ===
    print("\n=== Step 4: タンパク質マトリクス構築 ===")

    # データ読み込み
    df = pd.read_csv(EXPORT_TSV, sep="\t")
    print(f"PyProphetデータ: {len(df)} rows, {len(df.columns)} columns")

    # サンプル名抽出
    df["sample"] = df["filename"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
    print(f"サンプル数: {df['sample'].nunique()}")

    # Gene Symbol変換
    id_to_gene = parse_fasta_gene_map(FASTA_PATH)
    df["acc"] = df["ProteinName"].str.extract(r"sp\|(\w+)\|", expand=False)
    df["acc"] = df["acc"].fillna(df["ProteinName"])
    df["gene"] = df["acc"].map(id_to_gene).fillna(df["acc"])
    print(f"Gene Symbol: {df['gene'].nunique()} ユニーク遺伝子")

    # マトリクス構築
    matrix = df.pivot_table(
        index="gene",
        columns="sample",
        values="Intensity",
        aggfunc="sum"
    )
    matrix = matrix.replace(0, np.nan)
    matrix = matrix.dropna(how="all")
    matrix.index.name = "Protein"

    # 保存
    out_csv = os.path.join(RESULTS_DIR, "protein_matrix_from_openms.csv")
    matrix.to_csv(out_csv)

    print(f"\n=== パイプライン完了 ===")
    print(f"出力: {out_csv}")
    print(f"タンパク質マトリクス: {matrix.shape[0]} proteins × {matrix.shape[1]} samples")
    print("Toyota et al. 2025 論文再現用データ完成!")

if __name__ == "__main__":
    main()