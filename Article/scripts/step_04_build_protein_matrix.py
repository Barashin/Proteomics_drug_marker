#!/usr/bin/env python3
"""
Step 4.5: sage出力からタンパク質定量マトリクスを構築
====================================================
sage の lfq.tsv はペプチドレベルの定量テーブルなので、
これをタンパク質レベルに集約して、ダウンストリーム解析（差分発現、PCA等）で
使える形式に変換します。

【処理の流れ】
  1. sage の lfq.tsv を読み込み
  2. q_value < 0.01 でフィルタ（ペプチドFDR <1%）
  3. 各ペプチドを代表タンパク質（先頭のUniProt ID）に割り当て
  4. FASTA ヘッダから UniProt ID → Gene Symbol のマッピングを構築
  5. ペプチド強度を Gene Symbol ごとに合計してタンパク質強度に集約
  6. タンパク質 × サンプルのマトリクスを出力（CSV）

【設計上の判断】
  - タンパク質推論は razor protein 方式（ペプチドが複数タンパク質にヒットする場合、
    先頭のUniProt IDを代表とする）で単純化。論文のDIA-NNも同等の推論を行う。
  - サンプル名は lfq.tsv のカラムヘッダ（例: "CRC04-N.mzML"）から .mzML を除去。
  - 集約は単純な強度合計（sum）。MaxLFQ 的な手法は使わない（sage の段階で
    `combine_charge_states: true` により charge 間はすでに統合済み）。

入力: results/sage_output/lfq.tsv
      data/raw/human_proteome.fasta
出力: results/protein_matrix_from_sage.csv   (Gene Symbol × サンプル、強度値)
"""

# ============================================================
# ライブラリの読み込み
# ============================================================
import os
import re
import pandas as pd

# ============================================================
# パス設定
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, "..")
LFQ_TSV = os.path.join(PROJECT_DIR, "results", "sage_output", "lfq.tsv")
FASTA = os.path.join(PROJECT_DIR, "data", "raw", "human_proteome.fasta")
OUT_CSV = os.path.join(PROJECT_DIR, "results", "protein_matrix_from_sage.csv")

# ============================================================
# パラメータ
# ============================================================
# ペプチドFDRの閾値。論文は <1% なので 0.01。
Q_VALUE_THRESHOLD = 0.01


# ============================================================
# FASTA パース: UniProt ID -> Gene Symbol のマップ
# ============================================================
def build_id_to_gene_map(fasta_path: str) -> dict:
    """FASTAヘッダから UniProt accession -> Gene Symbol の辞書を作る。

    UniProt のヘッダ例:
        >sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=4

    ここから accession = "P04637"、gene symbol = "TP53" を抽出する。
    GN= が存在しないエントリは "Entry Name" の先頭語（例: P53_HUMAN -> P53）を使う。
    """
    id_to_gene = {}
    gn_re = re.compile(r"\bGN=([^\s]+)")
    with open(fasta_path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            # ヘッダ行: ">sp|ACC|NAME DESC..."
            header = line[1:].rstrip()
            first_space = header.find(" ")
            id_part = header[:first_space] if first_space > 0 else header
            desc_part = header[first_space + 1 :] if first_space > 0 else ""

            # id_part = "sp|ACC|NAME" or "tr|ACC|NAME"
            parts = id_part.split("|")
            if len(parts) >= 3:
                acc = parts[1]
                entry_name = parts[2]  # 例: "P53_HUMAN"
            else:
                acc = id_part
                entry_name = id_part

            # GN= フィールドがあればそれを使う
            m = gn_re.search(desc_part)
            if m:
                gene = m.group(1)
            else:
                # GN= がない場合は Entry Name の先頭語
                gene = entry_name.split("_")[0]

            id_to_gene[acc] = gene
    return id_to_gene


# ============================================================
# sage proteins カラムから代表UniProt accession を取り出す
# ============================================================
def extract_primary_acc(proteins_field: str) -> str:
    """sage の proteins カラムは "sp|ACC1|NAME;tr|ACC2|NAME;..." の形式。

    最初の要素の accession（2番目のフィールド）を返す。
    """
    if not isinstance(proteins_field, str) or not proteins_field:
        return ""
    first = proteins_field.split(";")[0]
    parts = first.split("|")
    if len(parts) >= 2:
        return parts[1]
    return first


# ============================================================
# メイン処理
# ============================================================
def main():
    # ----- 1. lfq.tsv を読み込み -----
    if not os.path.exists(LFQ_TSV):
        raise SystemExit(f"lfq.tsv not found: {LFQ_TSV} (先にstep_04_dia_analysis_sage.pyを実行してください)")
    print(f"読み込み: {LFQ_TSV}")
    df = pd.read_csv(LFQ_TSV, sep="\t")
    print(f"  {len(df)} rows, {len(df.columns)} columns")
    print(f"  columns: {df.columns.tolist()}")

    # ----- 2. q_value でフィルタ -----
    before = len(df)
    df = df[df["q_value"] < Q_VALUE_THRESHOLD].copy()
    print(f"  q_value < {Q_VALUE_THRESHOLD}: {before} -> {len(df)} peptides")

    # ----- 3. サンプル列（.mzML で終わる列）を特定 -----
    meta_cols = ["peptide", "charge", "proteins", "q_value", "score", "spectral_angle"]
    sample_cols = [c for c in df.columns if c not in meta_cols]
    print(f"  サンプル列: {len(sample_cols)}")
    for c in sample_cols:
        print(f"    - {c}")

    # .mzML サフィックスを除去してサンプル名を整形
    clean_name = lambda c: re.sub(r"\.mzML$", "", c)
    df = df.rename(columns={c: clean_name(c) for c in sample_cols})
    sample_cols = [clean_name(c) for c in sample_cols]

    # ----- 4. 代表 UniProt accession を取り出し -----
    df["primary_acc"] = df["proteins"].apply(extract_primary_acc)

    # ----- 5. FASTA から Gene Symbol マップを構築 -----
    print(f"\nFASTAから Gene Symbol を抽出: {FASTA}")
    id_to_gene = build_id_to_gene_map(FASTA)
    print(f"  {len(id_to_gene)} エントリ")

    # ペプチドごとに gene を割り当て（マップに無ければ accession 自体を使う）
    df["gene"] = df["primary_acc"].map(id_to_gene).fillna(df["primary_acc"])

    # ----- 6. Gene ごとに強度を合計 -----
    # ペプチド強度は log スケールではなく生の強度なので、合計が意味を持つ
    # ※ sage の lfq.tsv の値は peak integration の合計（Sum integration）
    matrix = df.groupby("gene")[sample_cols].sum(min_count=1)
    print(f"\nタンパク質マトリクス: {matrix.shape[0]} proteins × {matrix.shape[1]} samples")

    # 0 の値は NaN に置換（検出されなかった = 欠損値として扱う）
    matrix = matrix.replace(0, pd.NA)

    # 有効値のあるタンパク質のみ残す
    before = len(matrix)
    matrix = matrix.dropna(how="all")
    print(f"  全サンプルで0/欠損のタンパク質を除去: {before} -> {len(matrix)}")

    # ----- 7. 保存 -----
    matrix.index.name = "Protein"
    matrix.to_csv(OUT_CSV)
    print(f"\n保存: {OUT_CSV}")
    print(f"  ヘッダ例: {list(matrix.columns)[:5]}")
    print(f"  先頭行例:")
    print(matrix.head())


if __name__ == "__main__":
    main()
