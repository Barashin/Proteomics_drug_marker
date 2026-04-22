#!/usr/bin/env python3
"""
全図を論文(Toyota et al. 2025)スタイルで再生成するスクリプト。
データソース: mzML → sage → preprocess 由来の preprocessed_data.csv
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.patches import Ellipse, Patch
from matplotlib.gridspec import GridSpec
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests

np.random.seed(42)

# ============================================================
# Paths
# ============================================================
RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
TABLE_DIR = os.path.join(RESULTS_DIR, "tables")
BLOG_IMG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "blog", "images")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

# ============================================================
# Load data (mzML → sage → preprocess pipeline)
# ============================================================
print("=== Loading data (mzML-derived via sage) ===")
df = pd.read_csv(os.path.join(RESULTS_DIR, "preprocessed_data.csv"), index_col=0)
sample_info = pd.read_csv(os.path.join(RESULTS_DIR, "sample_info.csv"))
clinical = pd.read_csv(os.path.join(RESULTS_DIR, "clinical_info.csv"))
print(f"  {df.shape[0]} proteins x {df.shape[1]} samples")

# Build stage map
stage_map = {}
for _, row in clinical.iterrows():
    stage_map[row["Sample_N"]] = "Normal"
    stage_map[row["Sample_T"]] = row["Stage"]
sample_info["Stage"] = sample_info["Sample"].map(stage_map)

conditions = sample_info.set_index("Sample")["Condition"]
stages_series = sample_info.set_index("Sample")["Stage"]
normal_samples = sample_info[sample_info["Condition"] == "Normal"]["Sample"].tolist()
tumor_samples = sample_info[sample_info["Condition"] == "Tumor"]["Sample"].tolist()

# ============================================================
# Paper-matching color scheme
# ============================================================
# Histology colors (Fig 1b, 2a-e)
HIST_NONTUMOR = "#3498DB"  # Blue
HIST_TUMOR = "#E74C3C"     # Red

# Clinical stage colors (paper legend)
STAGE_COLORS = {
    "Normal": "#3498DB",
    "I":      "#66BB6A",
    "II":     "#FFD54F",
    "III":    "#FFA726",
    "IV":     "#AB47BC",
}

# Protein expression pattern (Fig 1b, 2c-e)
COLOR_UP_PROTEIN = "#E8A0BF"    # Pink = Up-regulated in tumor tissues
COLOR_DOWN_PROTEIN = "#A0C4E8"  # Light blue = Down-regulated in tumor tissues


def save_fig(fig_or_g, filename, dpi=150):
    """Save to both results/figures and blog/images."""
    for d in [FIG_DIR, BLOG_IMG_DIR]:
        path = os.path.join(d, filename)
        fig_or_g.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"    Saved: {filename}")


def confidence_ellipse(x, y, ax, n_std=2.0, **kwargs):
    """Draw a confidence ellipse on the given axes."""
    if len(x) < 2:
        return
    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, **kwargs)
    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_x, mean_y = np.mean(x), np.mean(y)
    transf = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)
    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


# ============================================================
# Figure 1a: Correlation matrix heatmap
# ============================================================
def plot_fig1a():
    print("\n[Fig 1a] Correlation matrix heatmap")
    corr = df.corr(method="pearson")

    # Row/col color = Histology
    histology_colors = conditions.map({
        "Normal": HIST_NONTUMOR, "Tumor": HIST_TUMOR
    }).reindex(corr.index)

    g = sns.clustermap(
        corr, method="average", metric="correlation",
        cmap="Reds", vmin=0.7, vmax=1.0,
        figsize=(10, 10),
        row_colors=histology_colors, col_colors=histology_colors,
        linewidths=0, xticklabels=True, yticklabels=True,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=6, rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=6)

    # Add Pearson r annotation on colorbar
    cbar = g.ax_cbar
    cbar.set_title("Pearson r", fontsize=8)

    legend_elements = [
        Patch(facecolor=HIST_NONTUMOR, label="Non-tumor tissue"),
        Patch(facecolor=HIST_TUMOR, label="Tumor tissue"),
    ]
    g.ax_heatmap.legend(handles=legend_elements, loc="upper left",
                        bbox_to_anchor=(1.05, 1.0), frameon=False, fontsize=8)

    save_fig(g, "fig1a_correlation.png")
    plt.close()
    corr.to_csv(os.path.join(TABLE_DIR, "correlation_matrix.csv"))


# ============================================================
# Figure 1b: Hierarchical clustering with Histology + Stage bars
# ============================================================
def plot_fig1b():
    print("\n[Fig 1b] Hierarchical clustering")

    # Dual row colors: Histology + Clinical stage
    hist_colors = conditions.map({
        "Normal": HIST_NONTUMOR, "Tumor": HIST_TUMOR
    }).reindex(df.columns)
    stage_colors = stages_series.map(STAGE_COLORS).reindex(df.columns)

    row_colors_df = pd.DataFrame({
        "Histology": hist_colors,
        "Clinical stage": stage_colors,
    })

    # Protein direction for column colors
    normal_mean = df[normal_samples].mean(axis=1)
    tumor_mean = df[tumor_samples].mean(axis=1)
    fc = tumor_mean - normal_mean
    protein_colors = fc.apply(lambda x: COLOR_UP_PROTEIN if x > 0 else COLOR_DOWN_PROTEIN)
    col_colors_df = pd.DataFrame({
        "Protein expression\npattern": protein_colors,
    })

    g = sns.clustermap(
        df.T, method="ward", metric="euclidean",
        cmap="RdBu_r", center=0, vmin=-3, vmax=3,
        figsize=(14, 8),
        row_colors=row_colors_df, col_colors=col_colors_df,
        col_cluster=True, row_cluster=True,
        xticklabels=False, yticklabels=True,
        z_score=1,
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=7)

    # Combined legend
    legend_elements = [
        Patch(facecolor=HIST_NONTUMOR, label="Non-tumor tissue"),
        Patch(facecolor=HIST_TUMOR, label="Tumor tissue"),
    ]
    for stage_name in ["I", "II", "III", "IV"]:
        legend_elements.append(Patch(facecolor=STAGE_COLORS[stage_name],
                                     label=f"Stage {stage_name}"))
    legend_elements.append(Patch(facecolor=COLOR_UP_PROTEIN,
                                 label="Up-regulated in tumor tissues"))
    legend_elements.append(Patch(facecolor=COLOR_DOWN_PROTEIN,
                                 label="Down-regulated in tumor tissues"))

    g.ax_heatmap.legend(handles=legend_elements, loc="lower right",
                        bbox_to_anchor=(1.4, -0.25), frameon=False, fontsize=7, ncol=2)

    save_fig(g, "fig1b_clustering.png")
    plt.close()


# ============================================================
# Figure 1c: PCA with Histology + Clinical stage + dashed ellipses
# ============================================================
def plot_fig1c():
    print("\n[Fig 1c] PCA")

    pca = PCA(n_components=2)
    scores = pca.fit_transform(df.T)

    fig, ax = plt.subplots(figsize=(8, 6))

    normal_idx = [i for i, s in enumerate(df.columns) if s in normal_samples]
    tumor_idx = [i for i, s in enumerate(df.columns) if s in tumor_samples]

    # Non-tumor: all green (paper style)
    ax.scatter(scores[normal_idx, 0], scores[normal_idx, 1],
               c="#2E8B57", marker="o", s=100, alpha=0.85,
               edgecolors="white", linewidth=0.5, zorder=3)

    # Tumor: colored by clinical stage (paper style)
    plotted_stages = set()
    for i in tumor_idx:
        sample = df.columns[i]
        stage = stage_map.get(sample, "III")
        color = STAGE_COLORS.get(stage, "#E74C3C")
        ax.scatter(scores[i, 0], scores[i, 1],
                   c=color, marker="o", s=100, alpha=0.85,
                   edgecolors="white", linewidth=0.5, zorder=3)
        plotted_stages.add(stage)

    # Dashed confidence ellipses (paper style)
    confidence_ellipse(scores[normal_idx, 0], scores[normal_idx, 1], ax, n_std=2.0,
                       facecolor="none", edgecolor="#2E8B57", linestyle="--", linewidth=1.5)
    confidence_ellipse(scores[tumor_idx, 0], scores[tumor_idx, 1], ax, n_std=2.0,
                       facecolor="none", edgecolor="#E74C3C", linestyle="--", linewidth=1.5)

    ax.set_xlabel(f"Component 1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"Component 2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")

    # Histology legend
    hist_handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#2E8B57",
                    markersize=8, label="Non-tumor tissue"),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#E74C3C",
                    markersize=8, label="Tumor tissue"),
    ]
    legend1 = ax.legend(handles=hist_handles, loc="upper left", frameon=False,
                        fontsize=8, title="Histology", title_fontsize=9)
    ax.add_artist(legend1)

    # Clinical stage legend
    stage_handles = []
    for stage_name in ["I", "II", "III", "IV"]:
        stage_handles.append(
            plt.Line2D([0], [0], marker='o', color='w',
                       markerfacecolor=STAGE_COLORS[stage_name],
                       markersize=6, label=f"Stage {stage_name}"))
    ax.legend(handles=stage_handles, loc="lower left", frameon=False,
              fontsize=7, title="Clinical stage", title_fontsize=8)
    ax.add_artist(legend1)

    ax.grid(True, color="gray", alpha=0.3, linestyle="-", linewidth=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    save_fig(fig, "fig1c_pca.png")
    plt.close()

    var_df = pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(len(pca.explained_variance_ratio_))],
        "Variance_Ratio": pca.explained_variance_ratio_,
        "Cumulative": np.cumsum(pca.explained_variance_ratio_),
    })
    var_df.to_csv(os.path.join(TABLE_DIR, "pca_variance.csv"), index=False)


# ============================================================
# Welch's t-test (differential expression)
# ============================================================
def run_welch_ttest():
    print("\n[DE analysis] Welch's t-test")
    results = []
    for protein in df.index:
        nv = df.loc[protein, normal_samples].dropna()
        tv = df.loc[protein, tumor_samples].dropna()
        if len(nv) < 2 or len(tv) < 2:
            continue
        t_stat, p_val = stats.ttest_ind(tv, nv, equal_var=False)
        log2fc = tv.mean() - nv.mean()
        results.append({
            "Protein": protein,
            "Mean_Normal": nv.mean(), "Mean_Tumor": tv.mean(),
            "Log2FC": log2fc, "T_statistic": t_stat,
            "P_value": p_val,
            "Neg_log10_P": -np.log10(max(p_val, 1e-300)),
        })
    result_df = pd.DataFrame(results)
    result_df["Significant"] = "NS"
    result_df.loc[(result_df["P_value"] < 0.05) & (result_df["Log2FC"] > 1.0), "Significant"] = "Up"
    result_df.loc[(result_df["P_value"] < 0.05) & (result_df["Log2FC"] < -1.0), "Significant"] = "Down"

    n_up = (result_df["Significant"] == "Up").sum()
    n_down = (result_df["Significant"] == "Down").sum()
    print(f"  Significant: Up={n_up}, Down={n_down}, NS={len(result_df)-n_up-n_down}")

    result_df.to_csv(os.path.join(TABLE_DIR, "differential_proteins.csv"), index=False)
    return result_df


# ============================================================
# Figure Bonus: Volcano plot
# ============================================================
def plot_volcano(result_df):
    print("\n[Volcano] Bonus figure")
    fig, ax = plt.subplots(figsize=(8, 6))

    for sig, color, alpha in [("NS", "#CCCCCC", 0.3), ("Up", "#E74C3C", 0.6), ("Down", "#3498DB", 0.6)]:
        mask = result_df["Significant"] == sig
        ax.scatter(result_df.loc[mask, "Log2FC"], result_df.loc[mask, "Neg_log10_P"],
                   c=color, s=10, alpha=alpha, label=f"{sig} ({mask.sum()})")

    ax.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.5)
    ax.axvline(1.0, color="gray", linestyle="--", linewidth=0.5)
    ax.axvline(-1.0, color="gray", linestyle="--", linewidth=0.5)

    ax.set_xlabel("Log2 Fold Change (Tumor / Normal)")
    ax.set_ylabel("-Log10(P-value)")
    ax.set_title("Differential Protein Abundance: Tumor vs Non-tumor")
    ax.legend(frameon=False, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    save_fig(fig, "fig_bonus_volcano.png")
    plt.close()


# ============================================================
# Figure 2a: All DE proteins heatmap
# ============================================================
def plot_fig2a(result_df):
    print("\n[Fig 2a] All DE proteins heatmap")

    sig_proteins = result_df[result_df["Significant"] != "NS"]["Protein"].values
    sig_data = df.loc[df.index.isin(sig_proteins)]
    if len(sig_data) == 0:
        print("  Skipped: no significant proteins")
        return

    # Row colors: Histology + Clinical stage
    hist_colors = conditions.map({
        "Normal": HIST_NONTUMOR, "Tumor": HIST_TUMOR
    }).reindex(df.columns)
    stage_colors = stages_series.map(STAGE_COLORS).reindex(df.columns)
    row_colors_df = pd.DataFrame({
        "Histology": hist_colors,
        "Clinical stage": stage_colors,
    })

    # Column colors: Protein direction
    protein_fc = result_df.set_index("Protein")["Log2FC"]
    col_color_vals = sig_data.index.map(
        lambda p: COLOR_UP_PROTEIN if protein_fc.get(p, 0) > 0 else COLOR_DOWN_PROTEIN
    )
    col_colors = pd.Series(col_color_vals, index=sig_data.index, name="Expression")

    g = sns.clustermap(
        sig_data.T, method="ward", metric="euclidean",
        cmap="RdBu_r", center=0, vmin=-3, vmax=3,
        figsize=(12, 8),
        row_colors=row_colors_df, col_colors=col_colors,
        col_cluster=True, row_cluster=True,
        xticklabels=False, yticklabels=True,
        z_score=1,
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=7)

    save_fig(g, "fig2a_heatmap_all.png")
    plt.close()


# ============================================================
# Figure 2b: All DE proteins PCA
# ============================================================
def plot_fig2b(result_df):
    print("\n[Fig 2b] All DE proteins PCA")

    sig_proteins = result_df[result_df["Significant"] != "NS"]["Protein"].values
    sig_data = df.loc[df.index.isin(sig_proteins)]

    pca = PCA(n_components=2)
    scores = pca.fit_transform(sig_data.T)

    fig, ax = plt.subplots(figsize=(8, 6))

    normal_idx = [i for i, s in enumerate(df.columns) if s in normal_samples]
    tumor_idx = [i for i, s in enumerate(df.columns) if s in tumor_samples]

    # Paper style: Non-tumor=blue, Tumor=red
    ax.scatter(scores[normal_idx, 0], scores[normal_idx, 1],
               c=HIST_NONTUMOR, marker="o", s=80, alpha=0.8, label="Non-tumor tissue",
               edgecolors="white", linewidth=0.5, zorder=3)
    ax.scatter(scores[tumor_idx, 0], scores[tumor_idx, 1],
               c=HIST_TUMOR, marker="o", s=80, alpha=0.8, label="Tumor tissue",
               edgecolors="white", linewidth=0.5, zorder=3)

    # Dashed ellipses
    confidence_ellipse(scores[normal_idx, 0], scores[normal_idx, 1], ax, n_std=2.0,
                       facecolor="none", edgecolor=HIST_NONTUMOR, linestyle="--", linewidth=1.5)
    confidence_ellipse(scores[tumor_idx, 0], scores[tumor_idx, 1], ax, n_std=2.0,
                       facecolor="none", edgecolor=HIST_TUMOR, linestyle="--", linewidth=1.5)

    ax.set_xlabel(f"Component 1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"Component 2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title("Histology")
    ax.legend(frameon=False)
    ax.grid(True, color="gray", alpha=0.3, linestyle="-", linewidth=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    save_fig(fig, "fig2b_pca_all.png")
    plt.close()


# ============================================================
# Figure 2c-e: Top N clustering & Figure 2f-h: Top N PCA
# ============================================================
def plot_fig2_topN(result_df):
    sig_df = result_df[result_df["Significant"] != "NS"].copy()
    sig_up = sig_df[sig_df["Significant"] == "Up"].sort_values("Log2FC", ascending=False)
    sig_down = sig_df[sig_df["Significant"] == "Down"].sort_values("Log2FC", ascending=True)

    hist_colors = conditions.map({
        "Normal": HIST_NONTUMOR, "Tumor": HIST_TUMOR
    }).reindex(df.columns)
    stage_colors = stages_series.map(STAGE_COLORS).reindex(df.columns)
    row_colors_df = pd.DataFrame({
        "Histology": hist_colors,
        "Clinical stage": stage_colors,
    })

    normal_idx = [i for i, s in enumerate(df.columns) if s in normal_samples]
    tumor_idx = [i for i, s in enumerate(df.columns) if s in tumor_samples]

    for n, fig_c, fig_f in [(50, "fig2c_clustering_top50.png", "fig2f_pca_top50.png"),
                             (100, "fig2d_clustering_top100.png", "fig2g_pca_top100.png"),
                             (200, "fig2e_clustering_top200.png", "fig2h_pca_top200.png")]:
        print(f"\n[Fig 2] Top {n} clustering + PCA")

        # Select top N/2 up + top N/2 down (or available)
        n_half = n // 2
        top_up = sig_up.head(min(n_half, len(sig_up)))
        top_down = sig_down.head(min(n_half, len(sig_down)))
        top_combined = pd.concat([top_up, top_down])
        top_proteins = top_combined["Protein"].values
        top_data = df.loc[df.index.isin(top_proteins)]

        if len(top_data) == 0:
            print(f"  Skipped: insufficient data for top {n}")
            continue

        # Column colors: protein direction
        protein_fc = result_df.set_index("Protein")["Log2FC"]
        col_vals = top_data.index.map(
            lambda p: COLOR_UP_PROTEIN if protein_fc.get(p, 0) > 0 else COLOR_DOWN_PROTEIN
        )
        col_colors = pd.Series(col_vals, index=top_data.index, name="Expression")

        # Clustering (Fig 2c-e)
        g = sns.clustermap(
            top_data.T, method="ward", metric="euclidean",
            cmap="RdBu_r", center=0, vmin=-3, vmax=3,
            figsize=(10, 8),
            row_colors=row_colors_df, col_colors=col_colors,
            col_cluster=True, row_cluster=True,
            xticklabels=False, yticklabels=True,
            z_score=1,
        )
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=7)
        save_fig(g, fig_c)
        plt.close()

        # PCA (Fig 2f-h)
        pca = PCA(n_components=2)
        scores = pca.fit_transform(top_data.T)

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(scores[normal_idx, 0], scores[normal_idx, 1],
                   c=HIST_NONTUMOR, marker="o", s=80, alpha=0.8, label="Non-tumor tissue",
                   edgecolors="white", linewidth=0.5, zorder=3)
        ax.scatter(scores[tumor_idx, 0], scores[tumor_idx, 1],
                   c=HIST_TUMOR, marker="o", s=80, alpha=0.8, label="Tumor tissue",
                   edgecolors="white", linewidth=0.5, zorder=3)

        confidence_ellipse(scores[normal_idx, 0], scores[normal_idx, 1], ax, n_std=2.0,
                           facecolor="none", edgecolor=HIST_NONTUMOR, linestyle="--", linewidth=1.5)
        confidence_ellipse(scores[tumor_idx, 0], scores[tumor_idx, 1], ax, n_std=2.0,
                           facecolor="none", edgecolor=HIST_TUMOR, linestyle="--", linewidth=1.5)

        ax.set_xlabel(f"Component 1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        ax.set_ylabel(f"Component 2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        ax.legend(frameon=False)
        ax.grid(True, color="gray", alpha=0.3, linestyle="-", linewidth=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        save_fig(fig, fig_f)
        plt.close()


# ============================================================
# COSMIC coverage figure
# ============================================================
def plot_cosmic():
    print("\n[COSMIC] Coverage figure")

    COSMIC_CRC_GENES = [
        "APC", "TP53", "KRAS", "BRAF", "PIK3CA", "SMAD4", "FBXW7",
        "CTNNB1", "AMER1", "ARID1A", "MSH6", "MSH2", "MLH1", "PMS2",
        "RNF43", "ZNRF3", "AXIN2", "DCC", "SMAD2",
        "ERBB2", "ERBB3", "MET", "IGF1R", "CDH1", "POLE", "POLD1",
        "NRAS", "SOX9", "TCF7L2", "ACVR2A", "TGFBR2", "ATM", "BRCA2",
        "PTCH1", "STK11", "MAP2K1", "CREBBP", "EP300",
        "KMT2A", "KMT2D", "KMT2C", "ARID2", "SMARCA4",
        "NOTCH1", "NOTCH2", "NOTCH3", "JAK1", "JAK2",
        "FGFR1", "FGFR2", "FGFR3", "PDGFRA", "KIT",
        "RET", "ALK", "ROS1", "NTRK1", "NTRK2", "NTRK3",
        "IDH1", "IDH2", "EZH2", "NF1", "NF2",
    ]

    COSMIC_ALL_GENES = list(set(COSMIC_CRC_GENES + [
        "ABL1", "AKT1", "ALK", "AR", "BRCA1", "BRCA2", "CDK4", "CDKN2A",
        "EGFR", "EWSR1", "FANCA", "FLT3", "GATA3", "GNA11", "GNAQ",
        "HIF1A", "HRAS", "KDR", "MAP2K2", "MAP3K1", "MDM2", "MDM4",
        "MTOR", "MYC", "MYCN", "NPM1", "PAX5", "PBRM1", "PDGFRB",
        "PHF6", "PTEN", "RAC1", "RAF1", "RB1", "RUNX1", "SETD2",
        "SF3B1", "SOCS1", "SPOP", "STAG2", "STAT3", "SUZ12", "TET2",
        "TSHR", "U2AF1", "VHL", "WT1", "XPO1",
        "AKAP9", "ALDH2", "B2M", "CALR", "CARS1", "CHD4", "CMPK1",
        "COL1A1", "CTNNA1", "DDX3X", "DICER1", "EIF4A2", "FUBP1",
        "GNAS", "HSP90AA1", "HSP90AB1", "LARP4B", "LRP1B", "MAP2K4",
        "MAPK1", "MAX", "MEN1", "MUTYH", "MYD88", "NCOR1", "PPP2R1A",
        "PTPN11", "RAD21", "RECQL4", "RHOA", "SDHB", "SDHD",
        "SMARCB1", "SRC", "SRSF2", "TAF15", "TCF3", "TOP1",
        "TP63", "TRIM33", "TSC1", "TSC2", "WHSC1",
        "ABI1", "ACKR3", "ACSL3", "ACSL6", "AFDN",
        "AFF1", "AFF3", "AFF4", "AJUBA", "ATIC",
        "ATP1A1", "ATP2B3", "AXIN1", "BAP1", "BCL2",
        "BCL6", "BCL9", "BCOR", "BCORL1", "BCR",
        "BRD3", "BRD4", "BTK", "CAMTA1", "CARD11",
        "CASP8", "CBFA2T3", "CBFB", "CBL", "CCND1",
        "CCND2", "CCND3", "CCNE1", "CD274", "CD79A",
        "CD79B", "CDC73", "CDH11", "CDK12", "CDK6",
        "CDKN1A", "CDKN1B", "CDKN2B", "CDKN2C", "CEBPA",
    ]))

    identified = set(df.index)
    overlap_all = identified & set(COSMIC_ALL_GENES)
    overlap_crc = identified & set(COSMIC_CRC_GENES)

    coverage_all = len(overlap_all) / len(set(COSMIC_ALL_GENES)) * 100
    coverage_crc = len(overlap_crc) / len(set(COSMIC_CRC_GENES)) * 100

    print(f"  All cancer: {len(overlap_all)}/{len(set(COSMIC_ALL_GENES))} ({coverage_all:.1f}%)")
    print(f"  CRC: {len(overlap_crc)}/{len(set(COSMIC_CRC_GENES))} ({coverage_crc:.1f}%)")

    # Save overlap table
    overlap_data = []
    for gene in sorted(overlap_all):
        gene_type = "CRC_Specific" if gene in set(COSMIC_CRC_GENES) else "All_Cancer"
        overlap_data.append({"Gene": gene, "Type": gene_type})
    pd.DataFrame(overlap_data).to_csv(os.path.join(TABLE_DIR, "cosmic_overlap.csv"), index=False)

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # All cancer
    bars1 = axes[0].bar(
        ["COSMIC\n(All Cancer)", "Identified\nin This Study"],
        [len(set(COSMIC_ALL_GENES)), len(overlap_all)],
        color=["#FFA726", "#66BB6A"], width=0.6, alpha=0.8
    )
    axes[0].set_title(f"Cancer-Associated Proteins\n(Coverage: {coverage_all:.1f}%)")
    axes[0].set_ylabel("Number of Proteins")
    axes[0].bar_label(bars1, fontweight="bold", fontsize=11)

    # CRC
    bars2 = axes[1].bar(
        ["COSMIC\n(CRC)", "Identified\nin This Study"],
        [len(set(COSMIC_CRC_GENES)), len(overlap_crc)],
        color=["#FFA726", "#66BB6A"], width=0.6, alpha=0.8
    )
    axes[1].set_title(f"CRC-Associated Proteins\n(Coverage: {coverage_crc:.1f}%)")
    axes[1].set_ylabel("Number of Proteins")
    axes[1].bar_label(bars2, fontweight="bold", fontsize=11)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    save_fig(fig, "fig_cosmic_coverage.png")
    plt.close()


# ============================================================
# Figure 3: Stage analysis (paper style: RdBu_r + profile)
# ============================================================
def plot_fig3():
    print("\n[Fig 3] Stage-specific analysis")

    stage_order = ["Normal", "I", "II", "III", "IV"]
    available_stages = [s for s in stage_order if s in sample_info["Stage"].unique()]

    # Compute stage medians
    medians = {}
    for stage in available_stages:
        samples = sample_info[sample_info["Stage"] == stage]["Sample"].tolist()
        valid = [s for s in samples if s in df.columns]
        if valid:
            medians[stage] = df[valid].median(axis=1)
    median_df = pd.DataFrame(medians)

    # One-way ANOVA
    results = []
    for protein in df.index:
        groups = []
        for stage in available_stages:
            samples = sample_info[sample_info["Stage"] == stage]["Sample"].tolist()
            vals = df.loc[protein, [s for s in samples if s in df.columns]].dropna()
            if len(vals) >= 2:
                groups.append(vals.values)
        if len(groups) < 2:
            continue
        f_stat, p_val = stats.f_oneway(*groups)
        results.append({"Protein": protein, "F_statistic": f_stat, "P_value": p_val})

    anova_df = pd.DataFrame(results)
    if len(anova_df) == 0:
        print("  Skipped: no ANOVA results")
        return

    reject, fdr, _, _ = multipletests(anova_df["P_value"], method="fdr_bh")
    anova_df["FDR"] = fdr
    anova_df["Significant"] = fdr < 0.01

    n_sig = anova_df["Significant"].sum()
    print(f"  ANOVA significant (FDR<0.01): {n_sig}")
    anova_df.to_csv(os.path.join(TABLE_DIR, "anova_results.csv"), index=False)

    # Get significant proteins
    sig_proteins = anova_df[anova_df["Significant"]]["Protein"].tolist()
    if not sig_proteins:
        sig_proteins = anova_df.nsmallest(min(200, len(anova_df)), "FDR")["Protein"].tolist()

    sig_data = median_df.loc[median_df.index.isin(sig_proteins)]

    # Z-score normalization
    z_data = sig_data.apply(lambda x: (x - x.mean()) / x.std(), axis=1).dropna()
    if len(z_data) == 0:
        print("  Skipped: no valid Z-score data")
        return

    # Hierarchical clustering into 30 clusters
    n_clust = min(30, len(z_data))
    Z = linkage(z_data.values, method="ward")
    clusters = fcluster(Z, t=n_clust, criterion="maxclust")

    cluster_df = pd.DataFrame({"Protein": z_data.index, "Cluster": clusters})
    cluster_df.to_csv(os.path.join(TABLE_DIR, "cluster_assignments.csv"), index=False)

    # Find interesting clusters (increasing and decreasing patterns)
    cluster_profiles = {}
    for c in range(1, n_clust + 1):
        mask = clusters == c
        if mask.sum() == 0:
            continue
        mean_profile = z_data.iloc[mask.nonzero()[0]].mean(axis=0)
        cluster_profiles[c] = {
            "mean": mean_profile,
            "n_proteins": mask.sum(),
            "trend": mean_profile.iloc[-1] - mean_profile.iloc[0],
        }

    # Select 2 increasing + 2 decreasing clusters
    sorted_by_trend = sorted(cluster_profiles.items(), key=lambda x: x[1]["trend"])
    decreasing = [(k, v) for k, v in sorted_by_trend[:5] if v["n_proteins"] >= 3][:2]
    increasing = [(k, v) for k, v in sorted_by_trend[-5:][::-1] if v["n_proteins"] >= 3][:2]

    # If not enough, just pick the first available
    if len(increasing) < 2:
        increasing = [(k, v) for k, v in sorted(cluster_profiles.items(),
                       key=lambda x: -x[1]["trend"]) if (k, v) not in decreasing][:2]
    if len(decreasing) < 2:
        decreasing = [(k, v) for k, v in sorted(cluster_profiles.items(),
                       key=lambda x: x[1]["trend"]) if (k, v) not in increasing][:2]

    selected = increasing + decreasing
    if len(selected) < 4:
        selected = list(cluster_profiles.items())[:4]

    # Stage color bar for heatmap
    stage_color_list = [STAGE_COLORS.get(s, "#999999") for s in available_stages]

    # Paper style Figure 3: 2x2 grid, each cell = heatmap + profile
    fig = plt.figure(figsize=(14, 12))
    outer_gs = GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)

    labels = ["(a)", "(b)", "(c)", "(d)"]
    titles_prefix = ["Increasing", "Increasing", "Decreasing", "Decreasing"]

    for panel_idx, (clust_id, info) in enumerate(selected[:4]):
        inner_gs = outer_gs[panel_idx].subgridspec(2, 1, height_ratios=[3, 2], hspace=0.15)

        mask = clusters == clust_id
        cluster_z = z_data.iloc[mask.nonzero()[0]]
        n_prot = len(cluster_z)

        # --- Heatmap ---
        ax_heat = fig.add_subplot(inner_gs[0])

        # Reorder columns to match stage order
        plot_data = cluster_z[available_stages].values

        im = ax_heat.imshow(plot_data, aspect="auto", cmap="RdBu_r",
                            vmin=-2.5, vmax=2.5, interpolation="nearest")

        # Stage color bar on top
        for j, color in enumerate(stage_color_list):
            ax_heat.add_patch(plt.Rectangle((j - 0.5, -1.8), 1, 1.2,
                              facecolor=color, edgecolor="none",
                              clip_on=False, transform=ax_heat.transData))

        ax_heat.set_xticks(range(len(available_stages)))
        ax_heat.set_xticklabels([])
        ax_heat.set_ylabel("Protein", fontsize=8)
        ax_heat.set_title(f"Cluster {clust_id}: {n_prot} proteins {labels[panel_idx]}",
                          fontsize=10, fontweight="bold")

        if n_prot <= 30:
            ax_heat.set_yticks(range(n_prot))
            ax_heat.set_yticklabels(cluster_z.index, fontsize=4)
        else:
            ax_heat.set_yticks([])

        # Colorbar
        cbar = plt.colorbar(im, ax=ax_heat, fraction=0.03, pad=0.02)
        cbar.set_label("Z-score", fontsize=7)
        cbar.ax.tick_params(labelsize=6)

        # --- Profile plot ---
        ax_prof = fig.add_subplot(inner_gs[1])

        # Individual protein lines (gray)
        for _, row in cluster_z.iterrows():
            ax_prof.plot(range(len(available_stages)), row.values,
                         color="gray", alpha=0.15, linewidth=0.5)

        # Mean profile (thick colored)
        mean_vals = cluster_z.mean(axis=0).values
        for j in range(len(available_stages)):
            ax_prof.scatter(j, mean_vals[j], c=stage_color_list[j], s=50, zorder=5,
                            edgecolors="white", linewidth=0.5)
        ax_prof.plot(range(len(available_stages)), mean_vals,
                     color="#333333", linewidth=2.5, zorder=4)

        ax_prof.set_xticks(range(len(available_stages)))
        ax_prof.set_xticklabels(["Non-tumor"] + [f"Stage {s}" for s in available_stages[1:]],
                                fontsize=7, rotation=20)
        ax_prof.set_ylabel("Z-score", fontsize=8)
        ax_prof.axhline(0, color="gray", linestyle="--", linewidth=0.5)
        ax_prof.spines["top"].set_visible(False)
        ax_prof.spines["right"].set_visible(False)

    # Bottom legend
    legend_elements = []
    for stage, color in STAGE_COLORS.items():
        label = "Non-tumor tissue" if stage == "Normal" else f"Tumor tissue: Stage {stage}"
        legend_elements.append(Patch(facecolor=color, label=label))

    fig.legend(handles=legend_elements, loc="lower center", ncol=5,
               frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.02))

    save_fig(fig, "fig3_stage_heatmap.png")
    plt.close()


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("Regenerating all figures (paper style)")
    print("Data source: mzML -> sage -> preprocess pipeline")
    print("=" * 60)

    plot_fig1a()
    plot_fig1b()
    plot_fig1c()

    result_df = run_welch_ttest()
    plot_volcano(result_df)

    plot_fig2a(result_df)
    plot_fig2b(result_df)
    plot_fig2_topN(result_df)

    plot_cosmic()
    plot_fig3()

    print("\n" + "=" * 60)
    print("All figures regenerated successfully!")
    print("=" * 60)
