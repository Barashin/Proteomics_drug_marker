"""OpenMS + AlphaPeptDeep DIA解析 → タンパク質マトリクス構築。

notebook/step_05.ipynb と同一のパイプライン。`--subset N` で先頭Nサンプルに絞る。
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

REQUIRED_TOOLS = ["peptdeep", "TargetedFileConverter", "OpenSwathWorkflow", "pyprophet"]


def ensure_tool(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(
            f"{name} が見つかりません。`micromamba install -c bioconda openms` および "
            f"`pip install peptdeep pyprophet` を実行してください。"
        )


def run(cmd: list[str], label: str) -> int:
    print(f"\n[{label}] $ {' '.join(str(c) for c in cmd)}")
    t0 = time.time()
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
    )
    for line in proc.stdout:
        sys.stdout.write(line)
    proc.wait()
    print(f"[{label}] done in {(time.time()-t0)/60:.1f} min (exit={proc.returncode})")
    return proc.returncode


def parse_fasta_gene_map(fasta_path: str) -> dict[str, str]:
    id_to_gene: dict[str, str] = {}
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--subset", type=int, default=0, help="先頭Nサンプルのみ処理 (0=全件)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--skip-library", action="store_true", help="既存ライブラリを再利用")
    args = ap.parse_args()

    for t in REQUIRED_TOOLS:
        ensure_tool(t)

    np.random.seed(42)

    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    MZML_DIR = os.path.join(ROOT, "data/raw/raw_mzML")
    FASTA_PATH = os.path.join(ROOT, "data/raw/human_proteome.fasta")
    RESULTS_DIR = os.path.join(ROOT, "results")
    OPENMS_OUT = os.path.join(RESULTS_DIR, "openms_output")
    LIBRARY_DIR = os.path.join(OPENMS_OUT, "library")
    OSW_DIR = os.path.join(OPENMS_OUT, "openswath")
    PYPROPHET_DIR = os.path.join(OPENMS_OUT, "pyprophet")
    for d in [OPENMS_OUT, LIBRARY_DIR, OSW_DIR, PYPROPHET_DIR]:
        os.makedirs(d, exist_ok=True)

    mzml_files = sorted(glob.glob(os.path.join(MZML_DIR, "*.mzML")))
    if args.subset > 0:
        mzml_files = mzml_files[: args.subset]
    print(f"処理対象 mzML: {len(mzml_files)} files")
    for m in mzml_files:
        print(f"  - {os.path.basename(m)}")

    Q_THRESHOLD = 0.01
    MS1_PPM = 10.0
    MS2_PPM = 10.0

    # ---- Step 1: AlphaPeptDeep library prediction ----
    PREDICTED_LIB = os.path.join(LIBRARY_DIR, "predicted_library.tsv")
    if args.skip_library and os.path.exists(PREDICTED_LIB):
        print(f"\n[Step 1] SKIP: 既存ライブラリを使用 {PREDICTED_LIB}")
    else:
        import yaml
        PEPTDEEP_SETTINGS = os.path.join(OPENMS_OUT, "peptdeep_settings.yaml")
        # デフォルト設定を生成してから必要箇所のみ上書き
        default_yaml = os.path.join(OPENMS_OUT, "peptdeep_default.yaml")
        if subprocess.call(["peptdeep", "export-settings", default_yaml]) != 0:
            sys.exit(1)
        with open(default_yaml) as f:
            settings = yaml.safe_load(f)

        # タスクと実行デバイス
        settings["task_workflow"] = ["library"]
        settings["torch_device"]["device_type"] = "cpu"
        settings["thread_num"] = args.threads

        # ライブラリ予測設定
        lib = settings["library"]
        lib["infile_type"] = "fasta"
        lib["infiles"] = [os.path.abspath(FASTA_PATH)]
        lib["fasta"]["protease"] = "trypsin"
        lib["fasta"]["max_miss_cleave"] = 1
        lib["fix_mods"] = ["Carbamidomethyl@C"]
        lib["var_mods"] = []
        lib["min_var_mod_num"] = 0
        lib["max_var_mod_num"] = 0
        lib["min_peptide_len"] = 7
        lib["max_peptide_len"] = 45
        lib["min_precursor_charge"] = 2
        lib["max_precursor_charge"] = 4
        lib["min_precursor_mz"] = 200.0
        lib["max_precursor_mz"] = 2000.0
        lib["decoy"] = "pseudo_reverse"
        lib["frag_types"] = ["b", "y"]
        lib["max_frag_charge"] = 2
        lib["output_folder"] = os.path.abspath(LIBRARY_DIR)
        lib["output_tsv"]["enabled"] = True
        lib["output_tsv"]["min_fragment_mz"] = 200.0
        lib["output_tsv"]["max_fragment_mz"] = 1800.0

        with open(PEPTDEEP_SETTINGS, "w") as f:
            yaml.safe_dump(settings, f)
        rc = run(["peptdeep", "library", PEPTDEEP_SETTINGS], "Step1 peptdeep")
        if rc != 0:
            sys.exit(rc)
        # peptdeep は output_folder/predict.speclib.tsv を生成するので PREDICTED_LIB にコピー
        generated = os.path.join(LIBRARY_DIR, "predict.speclib.tsv")
        if os.path.exists(generated) and not os.path.exists(PREDICTED_LIB):
            shutil.copyfile(generated, PREDICTED_LIB)

    # ---- Step 2a: AlphaPeptDeep TSV → OpenSwath互換TSV → PQP ----
    OSW_TSV = os.path.join(LIBRARY_DIR, "predicted_library_oswath.tsv")
    if not os.path.exists(OSW_TSV) or os.path.getsize(OSW_TSV) < 1024:
        print(f"\n[Step2a-1] AlphaPeptDeep→OpenSwath 列名変換: {OSW_TSV}")
        t0 = time.time()
        rename_map = {
            "RT": "NormalizedRetentionTime",
            "FragmentMz": "ProductMz",
            "RelativeIntensity": "LibraryIntensity",
            "StrippedPeptide": "PeptideSequence",
            "ModifiedPeptide": "ModifiedPeptideSequence",
            "ProteinID": "ProteinName",
            "FragmentNumber": "FragmentSeriesNumber",
        }
        with open(PREDICTED_LIB) as fin, open(OSW_TSV, "w") as fout:
            header_line = fin.readline().rstrip("\n").split("\t")
            new_header = [rename_map.get(h, h) for h in header_line]
            # TransitionGroupId / TransitionId 列を追加
            fout.write("\t".join(new_header + ["TransitionGroupId", "TransitionId"]) + "\n")
            mp_idx = new_header.index("ModifiedPeptideSequence")
            pc_idx = new_header.index("PrecursorCharge")
            pn_idx = new_header.index("ProteinName")
            # AlphaPeptDeep [Mod] → UniMod ID マッピング
            mod_map = {
                "[Carbamidomethyl]": "(UniMod:4)",
                "[Oxidation]": "(UniMod:35)",
                "[Acetyl]": "(UniMod:1)",
                "[Phospho]": "(UniMod:21)",
            }
            for i, line in enumerate(fin):
                parts = line.rstrip("\n").split("\t")
                # AlphaPeptDeep の `_PEPTIDE_` 形式から OpenSwath が読める形へ
                mp = parts[mp_idx].strip("_")
                for k, v in mod_map.items():
                    mp = mp.replace(k, v)
                parts[mp_idx] = mp
                # ProteinName が `Q9UL16;A0A087X182` のように複数ある場合は先頭のみ採用
                if ";" in parts[pn_idx]:
                    parts[pn_idx] = parts[pn_idx].split(";")[0]
                tg = f"{parts[mp_idx]}_{parts[pc_idx]}"
                ti = f"{tg}_{i}"
                fout.write("\t".join(parts + [tg, ti]) + "\n")
        print(f"[Step2a-1] done in {(time.time()-t0)/60:.1f} min")

    PQP_LIB = os.path.join(LIBRARY_DIR, "predicted_library.pqp")
    if not os.path.exists(PQP_LIB) or os.path.getsize(PQP_LIB) < 1024:
        rc = run(
            ["TargetedFileConverter", "-in", OSW_TSV, "-out", PQP_LIB],
            "Step2a TargetedFileConverter",
        )
        if rc != 0:
            sys.exit(rc)

    # ---- Step 2b: OpenSWATH ----
    osw_files = []
    for i, mzml_path in enumerate(mzml_files):
        sample = os.path.splitext(os.path.basename(mzml_path))[0]
        osw_out = os.path.join(OSW_DIR, f"{sample}.osw")
        if os.path.exists(osw_out) and os.path.getsize(osw_out) > 1024:
            print(f"[Step2b {i+1}/{len(mzml_files)}] SKIP {sample} (既存)")
            osw_files.append(osw_out)
            continue
        rc = run(
            [
                "OpenSwathWorkflow",
                "-in", mzml_path,
                "-tr", PQP_LIB,
                "-out_osw", osw_out,
                "-min_upper_edge_dist", "1",
                "-mz_extraction_window", str(MS2_PPM),
                "-mz_extraction_window_unit", "ppm",
                "-mz_extraction_window_ms1", str(MS1_PPM),
                "-mz_extraction_window_ms1_unit", "ppm",
                "-Scoring:stop_report_after_feature", "5",
                "-Scoring:TransitionGroupPicker:min_peak_width", "10",
                "-threads", str(args.threads),
                "-force",
            ],
            f"Step2b OpenSWATH {i+1}/{len(mzml_files)} {sample}",
        )
        if rc != 0:
            sys.exit(rc)
        osw_files.append(osw_out)

    # ---- Step 3: PyProphet ----
    MERGED_OSW = os.path.join(PYPROPHET_DIR, "merged.osw")
    # merged.oswが既にあれば削除して新しいサブセットで再作成
    if os.path.exists(MERGED_OSW):
        os.remove(MERGED_OSW)
    rc = run(["pyprophet", "merge", "--out", MERGED_OSW] + osw_files, "Step3 merge")
    if rc != 0:
        sys.exit(rc)

    for stage in ["score", "peptide", "protein"]:
        cmd = ["pyprophet", stage, "--in", MERGED_OSW]
        if stage == "score":
            cmd += ["--level", "ms2", "--ss_initial_fdr", "0.15", "--ss_iteration_fdr", "0.05"]
        else:
            cmd += ["--context", "global"]
        rc = run(cmd, f"Step3 {stage}")
        if rc != 0:
            sys.exit(rc)

    EXPORT_TSV = os.path.join(PYPROPHET_DIR, "pyprophet_export.tsv")
    rc = run(
        [
            "pyprophet", "export",
            "--in", MERGED_OSW,
            "--out", EXPORT_TSV,
            "--max_global_peptide_qvalue", str(Q_THRESHOLD),
            "--max_global_protein_qvalue", str(Q_THRESHOLD),
        ],
        "Step3 export",
    )
    if rc != 0:
        sys.exit(rc)

    # ---- Step 4: matrix ----
    id_to_gene = parse_fasta_gene_map(FASTA_PATH)
    df = pd.read_csv(EXPORT_TSV, sep="\t")
    print(f"\nPyProphet export: {len(df)} rows, columns: {list(df.columns[:10])}...")

    df["sample"] = df["filename"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
    df["acc"] = df["ProteinName"].str.extract(r"sp\|(\w+)\|", expand=False)
    df["acc"] = df["acc"].fillna(df["ProteinName"])
    df["gene"] = df["acc"].map(id_to_gene).fillna(df["acc"])

    matrix = df.pivot_table(index="gene", columns="sample", values="Intensity", aggfunc="sum")
    matrix = matrix.replace(0, np.nan).dropna(how="all")
    matrix.index.name = "Protein"

    tag = f"_subset{args.subset}" if args.subset > 0 else ""
    out_csv = os.path.join(RESULTS_DIR, f"protein_matrix_from_openms{tag}.csv")
    matrix.to_csv(out_csv)
    print(f"\n=== 完了 ===")
    print(f"検出タンパク質数: {matrix.shape[0]}")
    print(f"サンプル数: {matrix.shape[1]}")
    print(f"保存: {out_csv}")


if __name__ == "__main__":
    main()
