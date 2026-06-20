#!/usr/bin/env python3
"""
DIA Pipeline Results Comparison Framework

Compares protein detection and quantification performance between:
1. Sage baseline (library-free, 2,110 proteins)
2. OpenMS + AlphaPeptDeep (deep-learning library prediction)

Metrics computed:
- Protein detection counts and overlap
- Quantification precision (CV analysis per condition group)
- Missing value rates
- Dynamic range comparison
- Statistical significance testing (paired Wilcoxon)
- Automated KPI report with figures

Usage:
    python compare_deep_learning_results.py                # full report
    python compare_deep_learning_results.py --sage-only    # baseline-only QC
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RESULTS_DIR = Path(__file__).resolve().parent / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
SAGE_MATRIX = RESULTS_DIR / "protein_matrix_from_sage.csv"
OPENMS_MATRIX = RESULTS_DIR / "protein_matrix_from_openms.csv"
CLINICAL_INFO = RESULTS_DIR / "clinical_info.csv"
SAMPLE_INFO = RESULTS_DIR / "sample_info.csv"
REPORT_OUT = RESULTS_DIR / "comparison_report.json"

# Toyota et al. 2025 reference target
TARGET_PROTEINS = 10_000

# Matplotlib defaults
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "font.size": 10,
})


# ===================================================================
# Data loading helpers
# ===================================================================

def load_matrix(path: Path, label: str) -> pd.DataFrame | None:
    """Load a protein x sample matrix CSV. Returns None if missing."""
    if not path.exists():
        print(f"[WARN] {label} matrix not found: {path}")
        return None
    df = pd.read_csv(path, index_col=0)
    print(f"[INFO] {label}: {df.shape[0]:,} proteins x {df.shape[1]} samples")
    return df


def load_sample_groups() -> dict[str, list[str]]:
    """Return {'Normal': [...], 'Tumor': [...]} from sample_info.csv."""
    if not SAMPLE_INFO.exists():
        return {}
    si = pd.read_csv(SAMPLE_INFO)
    groups: dict[str, list[str]] = {}
    for cond in si["Condition"].unique():
        groups[cond] = si.loc[si["Condition"] == cond, "Sample"].tolist()
    return groups


def load_clinical_stages() -> dict[str, str]:
    """Return {patient_pair_prefix: stage} from clinical_info.csv."""
    if not CLINICAL_INFO.exists():
        return {}
    ci = pd.read_csv(CLINICAL_INFO)
    return dict(zip(ci["Sample_N"].str.replace("-N", ""), ci["Stage"]))


# ===================================================================
# Metric functions
# ===================================================================

def compute_cv(matrix: pd.DataFrame, sample_groups: dict[str, list[str]]) -> pd.DataFrame:
    """Compute coefficient of variation (%) per protein within each condition group.

    Returns DataFrame with columns [Protein, Condition, CV_pct].
    """
    rows = []
    for cond, samples in sample_groups.items():
        cols = [s for s in samples if s in matrix.columns]
        if len(cols) < 2:
            continue
        sub = matrix[cols]
        mean = sub.mean(axis=1)
        std = sub.std(axis=1)
        cv = (std / mean) * 100
        for prot, val in cv.items():
            if np.isfinite(val):
                rows.append({"Protein": prot, "Condition": cond, "CV_pct": val})
    return pd.DataFrame(rows)


def compute_missing_rate(matrix: pd.DataFrame) -> dict:
    """Compute per-protein and global missing-value statistics."""
    total_cells = matrix.shape[0] * matrix.shape[1]
    missing_cells = matrix.isna().sum().sum()
    per_protein_missing = matrix.isna().sum(axis=1) / matrix.shape[1] * 100
    return {
        "global_missing_pct": missing_cells / total_cells * 100 if total_cells else 0,
        "per_protein_missing": per_protein_missing,
        "proteins_complete": int((per_protein_missing == 0).sum()),
        "proteins_gt50pct_missing": int((per_protein_missing > 50).sum()),
    }


def compute_dynamic_range(matrix: pd.DataFrame) -> dict:
    """Compute dynamic range stats (log2 scale assumed for preprocessed, raw otherwise)."""
    vals = matrix.values.flatten()
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if len(vals) == 0:
        return {"min": np.nan, "max": np.nan, "range_log2": np.nan}
    log2_vals = np.log2(vals) if vals.max() > 100 else vals  # heuristic: raw vs log2
    return {
        "min_log2": float(np.min(log2_vals)),
        "max_log2": float(np.max(log2_vals)),
        "range_log2": float(np.max(log2_vals) - np.min(log2_vals)),
        "median_log2": float(np.median(log2_vals)),
    }


def compute_overlap(sets: dict[str, set]) -> dict:
    """Compute pairwise overlap statistics for named protein sets."""
    names = list(sets.keys())
    results = {}
    for i, n1 in enumerate(names):
        for n2 in names[i + 1:]:
            s1, s2 = sets[n1], sets[n2]
            inter = s1 & s2
            union = s1 | s2
            results[f"{n1}_vs_{n2}"] = {
                "overlap": len(inter),
                "unique_to_first": len(s1 - s2),
                "unique_to_second": len(s2 - s1),
                "jaccard": len(inter) / len(union) if union else 0,
                "overlap_pct_of_first": len(inter) / len(s1) * 100 if s1 else 0,
                "overlap_pct_of_second": len(inter) / len(s2) * 100 if s2 else 0,
            }
    return results


def paired_intensity_test(sage_matrix: pd.DataFrame, other_matrix: pd.DataFrame) -> dict:
    """Wilcoxon signed-rank test on median protein intensities for shared proteins."""
    shared = sorted(set(sage_matrix.index) & set(other_matrix.index))
    if len(shared) < 10:
        return {"n_shared": len(shared), "statistic": np.nan, "p_value": np.nan}
    sage_medians = sage_matrix.loc[shared].median(axis=1)
    other_medians = other_matrix.loc[shared].median(axis=1)
    stat_result = stats.wilcoxon(sage_medians, other_medians, nan_policy="omit")
    return {
        "n_shared": len(shared),
        "statistic": float(stat_result.statistic),
        "p_value": float(stat_result.pvalue),
    }


# ===================================================================
# Visualization functions
# ===================================================================

def plot_detection_barplot(counts: dict[str, int], outpath: Path):
    """Horizontal bar chart of protein detection counts per pipeline."""
    fig, ax = plt.subplots(figsize=(7, max(3, len(counts) * 0.8)))
    methods = list(counts.keys())
    values = list(counts.values())
    colors = sns.color_palette("Set2", len(methods))
    bars = ax.barh(methods, values, color=colors, edgecolor="grey", linewidth=0.5)
    ax.axvline(TARGET_PROTEINS, color="red", linestyle="--", linewidth=1, label=f"Target ({TARGET_PROTEINS:,})")
    for bar, val in zip(bars, values):
        ax.text(bar.get_width() + max(values) * 0.01, bar.get_y() + bar.get_height() / 2,
                f"{val:,}", va="center", fontsize=9)
    ax.set_xlabel("Proteins detected")
    ax.set_title("Protein Detection: Pipeline Comparison")
    ax.legend(loc="lower right")
    fig.savefig(outpath)
    plt.close(fig)
    print(f"[FIG] Saved: {outpath}")


def plot_cv_distribution(cv_data: dict[str, pd.DataFrame], outpath: Path):
    """Overlay CV distributions for each pipeline."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, cond in zip(axes, ["Normal", "Tumor"]):
        for method, cv_df in cv_data.items():
            subset = cv_df.loc[cv_df["Condition"] == cond, "CV_pct"]
            if subset.empty:
                continue
            ax.hist(subset, bins=50, alpha=0.5, label=f"{method} (med={subset.median():.1f}%)",
                    range=(0, min(200, subset.quantile(0.99))))
        ax.set_xlabel("CV (%)")
        ax.set_ylabel("Protein count")
        ax.set_title(f"Quantification CV - {cond}")
        ax.legend(fontsize=8)
    fig.suptitle("Quantification Precision (CV Distribution)", y=1.02)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)
    print(f"[FIG] Saved: {outpath}")


def plot_missing_value_comparison(missing_data: dict[str, dict], outpath: Path):
    """Bar chart of missing value statistics per pipeline."""
    methods = list(missing_data.keys())
    global_miss = [missing_data[m]["global_missing_pct"] for m in methods]
    complete = [missing_data[m]["proteins_complete"] for m in methods]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    colors = sns.color_palette("Set2", len(methods))

    ax1.bar(methods, global_miss, color=colors, edgecolor="grey")
    ax1.set_ylabel("Missing values (%)")
    ax1.set_title("Global Missing Rate")
    for i, v in enumerate(global_miss):
        ax1.text(i, v + 0.3, f"{v:.1f}%", ha="center", fontsize=9)

    ax2.bar(methods, complete, color=colors, edgecolor="grey")
    ax2.set_ylabel("Number of proteins")
    ax2.set_title("Proteins with Complete Data (all 32 samples)")
    for i, v in enumerate(complete):
        ax2.text(i, v + 10, f"{v:,}", ha="center", fontsize=9)

    fig.suptitle("Missing Value Analysis", y=1.02)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)
    print(f"[FIG] Saved: {outpath}")


def plot_venn_two(set_a: set, set_b: set, label_a: str, label_b: str, outpath: Path):
    """Simple 2-set Venn-style overlap diagram using matplotlib (no venn library needed)."""
    only_a = len(set_a - set_b)
    only_b = len(set_b - set_a)
    both = len(set_a & set_b)

    fig, ax = plt.subplots(figsize=(7, 5))
    # Draw two overlapping circles
    circle_a = plt.Circle((-0.3, 0), 1, alpha=0.3, color="steelblue", label=label_a)
    circle_b = plt.Circle((0.3, 0), 1, alpha=0.3, color="coral", label=label_b)
    ax.add_patch(circle_a)
    ax.add_patch(circle_b)
    ax.set_xlim(-1.8, 1.8)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect("equal")
    ax.axis("off")

    ax.text(-0.85, 0, f"{only_a:,}\nonly", ha="center", va="center", fontsize=14, fontweight="bold")
    ax.text(0, 0, f"{both:,}\nshared", ha="center", va="center", fontsize=14, fontweight="bold")
    ax.text(0.85, 0, f"{only_b:,}\nonly", ha="center", va="center", fontsize=14, fontweight="bold")
    ax.text(-0.85, -1.2, label_a, ha="center", fontsize=11, color="steelblue")
    ax.text(0.85, -1.2, label_b, ha="center", fontsize=11, color="coral")
    ax.set_title("Protein Overlap", fontsize=13)

    fig.savefig(outpath)
    plt.close(fig)
    print(f"[FIG] Saved: {outpath}")


def plot_dynamic_range(matrices: dict[str, pd.DataFrame], outpath: Path):
    """Box/violin plot of log2-intensity distributions per pipeline."""
    plot_data = []
    for method, mat in matrices.items():
        vals = mat.values.flatten()
        vals = vals[np.isfinite(vals) & (vals > 0)]
        log2_vals = np.log2(vals) if vals.max() > 100 else vals
        for v in np.random.choice(log2_vals, size=min(50_000, len(log2_vals)), replace=False):
            plot_data.append({"Method": method, "log2_intensity": v})

    df_plot = pd.DataFrame(plot_data)
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.violinplot(data=df_plot, x="Method", y="log2_intensity", hue="Method",
                   ax=ax, palette="Set2", inner="quartile", legend=False)
    ax.set_ylabel("log2(intensity)")
    ax.set_title("Dynamic Range Comparison")
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)
    print(f"[FIG] Saved: {outpath}")


def plot_shared_protein_correlation(sage_matrix: pd.DataFrame, other_matrix: pd.DataFrame,
                                    other_label: str, outpath: Path):
    """Scatter plot of median intensities for shared proteins."""
    shared = sorted(set(sage_matrix.index) & set(other_matrix.index))
    if len(shared) < 10:
        print(f"[WARN] Too few shared proteins ({len(shared)}) for correlation plot")
        return
    sage_med = sage_matrix.loc[shared].median(axis=1)
    other_med = other_matrix.loc[shared].median(axis=1)

    # Use log2 if raw
    if sage_med.max() > 100:
        sage_med = np.log2(sage_med.replace(0, np.nan))
        other_med = np.log2(other_med.replace(0, np.nan))

    mask = np.isfinite(sage_med) & np.isfinite(other_med)
    r, p = stats.pearsonr(sage_med[mask], other_med[mask])

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(sage_med[mask], other_med[mask], alpha=0.3, s=8, color="steelblue")
    lims = [min(sage_med[mask].min(), other_med[mask].min()),
            max(sage_med[mask].max(), other_med[mask].max())]
    ax.plot(lims, lims, "k--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("Sage median log2(intensity)")
    ax.set_ylabel(f"{other_label} median log2(intensity)")
    ax.set_title(f"Shared proteins (n={mask.sum():,}), r={r:.3f}")
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)
    print(f"[FIG] Saved: {outpath}")


# ===================================================================
# KPI report builder
# ===================================================================

def build_kpi_report(
    detection_counts: dict[str, int],
    cv_summaries: dict[str, dict],
    missing_summaries: dict[str, dict],
    dynamic_range: dict[str, dict],
    overlap_stats: dict,
    stat_tests: dict,
) -> dict:
    """Assemble a structured KPI report dictionary."""
    report = {
        "generated_at": datetime.now().isoformat(),
        "target_proteins": TARGET_PROTEINS,
        "detection": {},
        "quantification_cv": {},
        "missing_values": {},
        "dynamic_range": {},
        "overlap": overlap_stats,
        "statistical_tests": stat_tests,
    }

    for method, count in detection_counts.items():
        report["detection"][method] = {
            "proteins_detected": count,
            "fold_vs_target": round(count / TARGET_PROTEINS, 3),
        }

    sage_count = detection_counts.get("Sage", 0)
    for method, count in detection_counts.items():
        if method != "Sage" and sage_count > 0:
            report["detection"][method]["fold_vs_sage"] = round(count / sage_count, 3)

    report["quantification_cv"] = cv_summaries
    report["missing_values"] = {
        m: {k: v for k, v in d.items() if k != "per_protein_missing"}
        for m, d in missing_summaries.items()
    }
    report["dynamic_range"] = dynamic_range

    return report


def print_kpi_summary(report: dict):
    """Print a human-readable KPI summary to stdout."""
    print("\n" + "=" * 65)
    print("  DIA PIPELINE COMPARISON -- KPI SUMMARY")
    print("=" * 65)

    print(f"\nGenerated: {report['generated_at']}")
    print(f"Target (Toyota et al. 2025): {report['target_proteins']:,} proteins")

    print("\n--- PROTEIN DETECTION ---")
    for method, d in report["detection"].items():
        line = f"  {method}: {d['proteins_detected']:,} proteins"
        line += f"  ({d['fold_vs_target']:.1%} of target)"
        if "fold_vs_sage" in d:
            line += f"  [{d['fold_vs_sage']:.2f}x vs Sage]"
        print(line)

    print("\n--- QUANTIFICATION CV (median %) ---")
    for method, cv_info in report["quantification_cv"].items():
        parts = [f"{cond}: {v:.1f}%" for cond, v in cv_info.items()]
        print(f"  {method}: {', '.join(parts)}")

    print("\n--- MISSING VALUES ---")
    for method, m in report["missing_values"].items():
        print(f"  {method}: {m['global_missing_pct']:.1f}% missing, "
              f"{m['proteins_complete']:,} complete, "
              f"{m['proteins_gt50pct_missing']:,} >50% missing")

    print("\n--- DYNAMIC RANGE (log2) ---")
    for method, dr in report["dynamic_range"].items():
        print(f"  {method}: [{dr['min_log2']:.1f}, {dr['max_log2']:.1f}] "
              f"range={dr['range_log2']:.1f}, median={dr['median_log2']:.1f}")

    if report["overlap"]:
        print("\n--- PROTEIN OVERLAP ---")
        for pair, o in report["overlap"].items():
            print(f"  {pair}: {o['overlap']:,} shared "
                  f"(Jaccard={o['jaccard']:.3f})")

    if report["statistical_tests"]:
        print("\n--- STATISTICAL TESTS ---")
        for pair, t in report["statistical_tests"].items():
            sig = "***" if t["p_value"] < 0.001 else "**" if t["p_value"] < 0.01 else "*" if t["p_value"] < 0.05 else "ns"
            print(f"  {pair}: Wilcoxon p={t['p_value']:.2e} ({sig}), n={t['n_shared']:,}")

    print("\n" + "=" * 65)


# ===================================================================
# Main entry point
# ===================================================================

def main():
    parser = argparse.ArgumentParser(description="DIA pipeline comparison framework")
    parser.add_argument("--sage-only", action="store_true",
                        help="Run baseline QC on Sage results only")
    args = parser.parse_args()

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    sage = load_matrix(SAGE_MATRIX, "Sage")
    if sage is None:
        print("[ERROR] Sage baseline matrix is required. Exiting.")
        sys.exit(1)

    openms = None if args.sage_only else load_matrix(OPENMS_MATRIX, "OpenMS+AlphaPeptDeep")

    sample_groups = load_sample_groups()

    # ---------------------------------------------------------------
    # Detection counts
    # ---------------------------------------------------------------
    detection_counts: dict[str, int] = {"Sage": sage.shape[0]}
    if openms is not None:
        detection_counts["OpenMS+AlphaPeptDeep"] = openms.shape[0]

    plot_detection_barplot(detection_counts, FIGURES_DIR / "cmp_detection_barplot.png")

    # ---------------------------------------------------------------
    # CV analysis
    # ---------------------------------------------------------------
    cv_data: dict[str, pd.DataFrame] = {}
    cv_summaries: dict[str, dict] = {}

    sage_cv = compute_cv(sage, sample_groups)
    cv_data["Sage"] = sage_cv
    cv_summaries["Sage"] = {
        cond: sage_cv.loc[sage_cv["Condition"] == cond, "CV_pct"].median()
        for cond in sample_groups
    }

    if openms is not None:
        openms_cv = compute_cv(openms, sample_groups)
        cv_data["OpenMS+AlphaPeptDeep"] = openms_cv
        cv_summaries["OpenMS+AlphaPeptDeep"] = {
            cond: openms_cv.loc[openms_cv["Condition"] == cond, "CV_pct"].median()
            for cond in sample_groups
        }

    if cv_data:
        plot_cv_distribution(cv_data, FIGURES_DIR / "cmp_cv_distribution.png")

    # ---------------------------------------------------------------
    # Missing value analysis
    # ---------------------------------------------------------------
    missing_summaries: dict[str, dict] = {"Sage": compute_missing_rate(sage)}
    if openms is not None:
        missing_summaries["OpenMS+AlphaPeptDeep"] = compute_missing_rate(openms)

    plot_missing_value_comparison(missing_summaries, FIGURES_DIR / "cmp_missing_values.png")

    # ---------------------------------------------------------------
    # Dynamic range
    # ---------------------------------------------------------------
    matrices_for_dr: dict[str, pd.DataFrame] = {"Sage": sage}
    if openms is not None:
        matrices_for_dr["OpenMS+AlphaPeptDeep"] = openms

    dr_stats: dict[str, dict] = {m: compute_dynamic_range(mat) for m, mat in matrices_for_dr.items()}
    plot_dynamic_range(matrices_for_dr, FIGURES_DIR / "cmp_dynamic_range.png")

    # ---------------------------------------------------------------
    # Protein overlap & correlation
    # ---------------------------------------------------------------
    overlap_stats: dict = {}
    stat_tests: dict = {}

    if openms is not None:
        protein_sets = {"Sage": set(sage.index), "OpenMS": set(openms.index)}
        overlap_stats = compute_overlap(protein_sets)
        plot_venn_two(set(sage.index), set(openms.index), "Sage", "OpenMS+AlphaPeptDeep",
                      FIGURES_DIR / "cmp_venn_overlap.png")
        plot_shared_protein_correlation(sage, openms, "OpenMS+AlphaPeptDeep",
                                        FIGURES_DIR / "cmp_correlation_shared.png")
        stat_tests["Sage_vs_OpenMS"] = paired_intensity_test(sage, openms)

    # ---------------------------------------------------------------
    # Build and save KPI report
    # ---------------------------------------------------------------
    report = build_kpi_report(
        detection_counts=detection_counts,
        cv_summaries=cv_summaries,
        missing_summaries=missing_summaries,
        dynamic_range=dr_stats,
        overlap_stats=overlap_stats,
        stat_tests=stat_tests,
    )
    print_kpi_summary(report)

    # Serialize (convert numpy types for JSON)
    def _convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        raise TypeError(f"Non-serializable: {type(obj)}")

    with open(REPORT_OUT, "w") as f:
        json.dump(report, f, indent=2, default=_convert)
    print(f"\n[INFO] Full report saved to: {REPORT_OUT}")
    print(f"[INFO] Figures saved to: {FIGURES_DIR}/cmp_*.png")


if __name__ == "__main__":
    main()
