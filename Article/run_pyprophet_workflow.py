#!/usr/bin/env python3
"""
PyProphet FDR推定ワークフロー（Step 3）
"""

import os, glob, subprocess, sys, time
import pandas as pd
import numpy as np

def run_command(cmd, description):
    """コマンド実行"""
    print(f"\n=== {description} ===")
    print(f"実行: {' '.join(cmd[:5])}...")

    start_time = time.time()

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        elapsed = (time.time() - start_time) / 60

        if proc.returncode == 0:
            print(f"✅ 完了: {elapsed:.1f}分")
            return True
        else:
            print(f"❌ エラー: 終了コード {proc.returncode}")
            print(f"詳細: {proc.stderr[:300]}...")
            return False
    except Exception as e:
        print(f"❌ 例外: {e}")
        return False

def parse_fasta_gene_map(fasta_path):
    """FASTAからGene Symbolマッピング作成"""
    import re
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
    print("=== PyProphet FDR推定ワークフロー ===")
    print("Toyota et al. 2025論文再現")

    # パス設定
    FASTA_PATH = "data/raw/human_proteome.fasta"
    RESULTS_DIR = "results"
    OSWR_DIR = "results/openms_output/openswath"
    PYPROPHET_DIR = "results/openms_output/pyprophet"

    os.makedirs(PYPROPHET_DIR, exist_ok=True)

    Q_THRESHOLD = 0.01  # 1% FDR

    # OSWファイル確認
    osw_files = sorted(glob.glob(os.path.join(OSWR_DIR, "*.osw")))
    print(f"📁 OSWファイル: {len(osw_files)}個")

    if len(osw_files) == 0:
        print("❌ OSWファイルが見つかりません。OpenSWATHを先に実行してください。")
        return

    # 最初の数個のファイルサイズを確認
    for i in range(min(3, len(osw_files))):
        size_mb = os.path.getsize(osw_files[i]) / 1024**2
        print(f"  {os.path.basename(osw_files[i])}: {size_mb:.1f} MB")

    # === Step 3.1: PyProphet merge ===
    MERGED_OSW = os.path.join(PYPROPHET_DIR, "merged.osw")

    # PyProphet 3.0.14の新しい構文を使用
    cmd_merge = [
        "pyprophet", "merge", "osw",
        "--out", MERGED_OSW,
        "--template", osw_files[0]  # 最初のファイルをテンプレートとして使用
    ] + osw_files
    if not run_command(cmd_merge, "PyProphet merge"):
        return

    if os.path.exists(MERGED_OSW):
        merged_size = os.path.getsize(MERGED_OSW) / 1024**2
        print(f"📊 統合ファイル: {merged_size:.1f} MB")

    # === Step 3.2: MS2スコアリング ===
    cmd_score = [
        "pyprophet", "score",
        "--in", MERGED_OSW,
        "--level", "ms2",
        "--ss_initial_fdr", "0.15",
        "--ss_iteration_fdr", "0.05"
    ]
    if not run_command(cmd_score, "PyProphet MS2スコアリング"):
        return

    # === Step 3.3: ペプチドレベルFDR ===
    cmd_peptide = [
        "pyprophet", "peptide",
        "--in", MERGED_OSW,
        "--context", "global"
    ]
    if not run_command(cmd_peptide, "PyProphet ペプチドFDR"):
        return

    # === Step 3.4: タンパク質レベルFDR ===
    cmd_protein = [
        "pyprophet", "protein",
        "--in", MERGED_OSW,
        "--context", "global"
    ]
    if not run_command(cmd_protein, "PyProphet タンパク質FDR"):
        return

    # === Step 3.5: エクスポート ===
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

    if os.path.exists(EXPORT_TSV):
        export_size = os.path.getsize(EXPORT_TSV) / 1024**2
        print(f"📊 エクスポート: {export_size:.1f} MB")

        # データ確認
        df = pd.read_csv(EXPORT_TSV, sep="\t", nrows=5)
        print(f"📋 データ: {len(df)} rows x {len(df.columns)} columns (サンプル)")
        print(f"📋 列例: {list(df.columns[:5])}")

    # === Step 4: タンパク質マトリクス構築 ===
    print("\n=== Step 4: タンパク質マトリクス構築 ===")

    df = pd.read_csv(EXPORT_TSV, sep="\t")
    print(f"📊 PyProphetデータ: {len(df)} rows x {len(df.columns)} columns")

    # サンプル名抽出
    df["sample"] = df["filename"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
    print(f"📊 サンプル数: {df['sample'].nunique()}")

    # Gene Symbol変換
    print("🧬 Gene Symbol変換中...")
    id_to_gene = parse_fasta_gene_map(FASTA_PATH)
    df["acc"] = df["ProteinName"].str.extract(r"sp\|(\w+)\|", expand=False)
    df["acc"] = df["acc"].fillna(df["ProteinName"])
    df["gene"] = df["acc"].map(id_to_gene).fillna(df["acc"])
    print(f"🧬 Gene Symbol: {df['gene'].nunique()} ユニーク遺伝子")

    # マトリクス構築
    print("🔧 タンパク質マトリクス構築中...")
    matrix = df.pivot_table(
        index="gene",
        columns="sample",
        values="Intensity",
        aggfunc="sum"
    )
    matrix = matrix.replace(0, np.nan)
    matrix = matrix.dropna(how="all")
    matrix.index.name = "Protein"

    print(f"📊 マトリクス: {matrix.shape[0]} proteins × {matrix.shape[1]} samples")

    # 保存
    out_csv = os.path.join(RESULTS_DIR, "protein_matrix_from_openms.csv")
    matrix.to_csv(out_csv)

    print(f"\n🎉 === Toyota et al. 2025論文再現完了！ ===")
    print(f"📄 出力: {out_csv}")
    print(f"📊 最終結果: {matrix.shape[0]} タンパク質 × {matrix.shape[1]} サンプル")
    print("🔬 OpenMS + AlphaPeptDeepパイプライン成功！")

if __name__ == "__main__":
    main()