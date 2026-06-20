#!/usr/bin/env python3
import sqlite3, pandas as pd, numpy as np, os, glob, re

def extract_precursor_data(osw_file, sample_name):
    try:
        conn = sqlite3.connect(osw_file)
        query = """
        SELECT p.PROTEIN_ACCESSION, COUNT(*) as FEATURE_COUNT
        FROM FEATURE f
        JOIN PRECURSOR prec ON f.PRECURSOR_ID = prec.ID
        JOIN PRECURSOR_PEPTIDE_MAPPING ppm ON prec.ID = ppm.PRECURSOR_ID
        JOIN PEPTIDE pep ON ppm.PEPTIDE_ID = pep.ID
        JOIN PEPTIDE_PROTEIN_MAPPING ppmap ON pep.ID = ppmap.PEPTIDE_ID
        JOIN PROTEIN p ON ppmap.PROTEIN_ID = p.ID
        GROUP BY p.PROTEIN_ACCESSION
        """
        df = pd.read_sql_query(query, conn)
        conn.close()
        if len(df) > 0:
            df['sample'] = sample_name
            df['INTENSITY'] = df['FEATURE_COUNT']
            return df
        return None
    except Exception as e:
        print(f"  Error: {e}")
        return None

def parse_fasta_gene_map(fasta_path):
    id_to_gene = {}
    gn_re = re.compile(r"\bGN=(\S+)")
    if not os.path.exists(fasta_path):
        return {}
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

# Main execution
osw_files = sorted([f for f in glob.glob("results/openms_output/openswath/*.osw") 
                   if not f.endswith("_fixed.osw")])
print(f"Processing {len(osw_files)} OSW files...")

all_data = []
for i, osw_file in enumerate(osw_files):
    sample_name = os.path.splitext(os.path.basename(osw_file))[0]
    print(f"[{i+1:2d}/{len(osw_files)}] {sample_name}")
    df = extract_precursor_data(osw_file, sample_name)
    if df is not None:
        all_data.append(df)
        print(f"  Proteins: {len(df)}")

if all_data:
    combined_df = pd.concat(all_data, ignore_index=True)
    print(f"\nTotal proteins: {combined_df['PROTEIN_ACCESSION'].nunique()}")
    print(f"Total samples: {combined_df['sample'].nunique()}")
    
    # Gene mapping
    id_to_gene = parse_fasta_gene_map("data/raw/human_proteome.fasta")
    combined_df["acc"] = combined_df["PROTEIN_ACCESSION"].str.extract(r"sp\|(\w+)\|", expand=False)
    combined_df["acc"] = combined_df["acc"].fillna(combined_df["PROTEIN_ACCESSION"])
    combined_df["gene"] = combined_df["acc"].map(id_to_gene).fillna(combined_df["acc"])
    
    # Create matrix
    matrix = combined_df.groupby(['gene', 'sample'])['INTENSITY'].sum().unstack(fill_value=0)
    matrix = matrix.replace(0, np.nan).dropna(how='all')
    matrix.index.name = "Protein"
    
    # Save
    matrix.to_csv("results/protein_matrix_openms_final.csv")
    print(f"\n🎉 SUCCESS: {matrix.shape[0]} proteins × {matrix.shape[1]} samples")
    print(f"📊 vs Sage (2,110): {matrix.shape[0]/2110:.1f}x improvement")
    print(f"🎯 vs Toyota (10,329): {matrix.shape[0]/10329:.1f}x")
else:
    print("❌ No data extracted")
