"""
Analysis Module — Generates publication-quality visualizations for normalized omics data.
Features:
  1. Sample Correlation Heatmap (Pearson / Spearman)
  2. Principal Component Analysis (PCA) — Raw vs Normalized side-by-side
  3. Distribution Comparisons — Box plots and Histograms (Raw vs Normalized)
"""

import io
import base64
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for server-side rendering
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import spearmanr
from sklearn.decomposition import PCA


def _fig_to_base64(fig) -> str:
    """Convert a matplotlib figure to a base64-encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def generate_correlation_heatmap(df: pd.DataFrame, sample_cols: list, method: str = "spearman") -> dict:
    """
    Generate a correlation heatmap across sample columns.
    Returns base64 image and the correlation matrix as a list-of-lists.
    """
    subset = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    if method == "spearman":
        corr_matrix, _ = spearmanr(subset.values)
        if len(sample_cols) == 1:
            corr_matrix = np.array([[1.0]])
        elif corr_matrix.ndim == 0:
            corr_matrix = np.array([[float(corr_matrix)]])
    else:
        corr_matrix = subset.corr(method="pearson").values

    fig, ax = plt.subplots(figsize=(max(6, len(sample_cols) * 1.2), max(5, len(sample_cols) * 1.0)))
    fig.patch.set_facecolor("white")

    cmap = LinearSegmentedColormap.from_list("custom", ["#4169E1", "#FFFFFF", "#DC143C"])
    im = ax.imshow(corr_matrix, cmap=cmap, vmin=corr_matrix.min() - 0.01, vmax=1.0, aspect="auto")

    ax.set_xticks(range(len(sample_cols)))
    ax.set_yticks(range(len(sample_cols)))
    ax.set_xticklabels(sample_cols, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(sample_cols, fontsize=9)

    ax.set_title(f"Sample Correlation Heatmap ({method.capitalize()})", fontsize=14, fontweight="bold", pad=15)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.ax.tick_params(labelsize=9)

    fig.tight_layout()

    return {
        "image": _fig_to_base64(fig),
        "matrix": corr_matrix.tolist(),
        "columns": sample_cols,
    }


def generate_pca_plots(df: pd.DataFrame, raw_cols: list, norm_cols: list, gene_id_col: str) -> str:
    """
    Generate side-by-side PCA scatter plots: Raw (log2) vs Normalized.
    Returns base64 image string.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor("white")

    datasets = [
        ("Raw (log2)", raw_cols, "#4169E1", axes[0]),
        ("Normalized", norm_cols, "#228B22", axes[1]),
    ]

    for label, cols, color, ax in datasets:
        subset = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0).values

        # Apply log2(x+1) transform for raw counts to compress dynamic range
        if label.startswith("Raw"):
            subset = np.log2(subset + 1)

        # Transpose: we want samples as observations (rows), genes as features (columns)
        data_T = subset.T

        n_components = min(2, data_T.shape[0], data_T.shape[1])
        if n_components < 2 or data_T.shape[0] < 2:
            ax.text(0.5, 0.5, "Not enough samples\nfor PCA", ha="center", va="center", transform=ax.transAxes, fontsize=12)
            ax.set_title(f"PCA - {label}", fontsize=12, fontweight="bold")
            continue

        pca = PCA(n_components=2)
        transformed = pca.fit_transform(data_T)
        ev = pca.explained_variance_ratio_

        ax.scatter(transformed[:, 0], transformed[:, 1], c=color, s=60, edgecolors="black", linewidth=0.5, zorder=5)
        ax.set_xlabel("PC1", fontsize=11)
        ax.set_ylabel("PC2", fontsize=11)
        ax.set_title(f"PCA - {label}\nEV: {ev[0]:.2f}, {ev[1]:.2f}", fontsize=12, fontweight="bold")
        ax.grid(True, alpha=0.3)

    fig.tight_layout(pad=3)
    return _fig_to_base64(fig)


def generate_distribution_plots(df: pd.DataFrame, raw_cols: list, norm_cols: list) -> str:
    """
    Generate distribution comparison plots: box plots + histograms for raw vs normalized data.
    Returns base64 image string.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.patch.set_facecolor("white")

    raw_data = df[raw_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    norm_data = df[norm_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    # Apply log2(x+1) transform
    raw_log2 = np.log2(raw_data + 1)
    norm_log2 = np.log2(norm_data + 1)

    # Row 1: Box plots
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

    # Raw box plot
    bp1 = axes[0, 0].boxplot(
        [raw_log2[c].dropna().values for c in raw_cols],
        labels=[c[:12] for c in raw_cols],
        patch_artist=True,
    )
    for i, patch in enumerate(bp1["boxes"]):
        patch.set_facecolor(colors[i % len(colors)])
    axes[0, 0].set_title("Raw Counts (log2) Distribution", fontsize=11, fontweight="bold")
    axes[0, 0].tick_params(axis="x", rotation=45)

    # Normalized box plot
    norm_labels = [c[:12] for c in norm_cols]
    bp2 = axes[0, 1].boxplot(
        [norm_log2[c].dropna().values for c in norm_cols],
        labels=norm_labels,
        patch_artist=True,
    )
    for i, patch in enumerate(bp2["boxes"]):
        patch.set_facecolor(colors[i % len(colors)])
    axes[0, 1].set_title("Normalized (log2) Distribution", fontsize=11, fontweight="bold")
    axes[0, 1].tick_params(axis="x", rotation=45)

    # Row 2: Histograms
    all_raw = raw_log2.values.flatten()
    all_raw = all_raw[np.isfinite(all_raw)]
    axes[1, 0].hist(all_raw, bins=40, color="#1f77b4", alpha=0.7, edgecolor="white", linewidth=0.3)
    # KDE-like line using histogram smoothing
    if len(all_raw) > 5:
        from scipy.ndimage import gaussian_filter1d
        counts, bin_edges = np.histogram(all_raw, bins=40)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        smooth = gaussian_filter1d(counts.astype(float), sigma=1.5)
        axes[1, 0].plot(bin_centers, smooth, color="#4169E1", linewidth=2)
    axes[1, 0].set_title("Raw Counts (log2) Histogram", fontsize=11, fontweight="bold")
    axes[1, 0].set_xlabel("Value")
    axes[1, 0].set_ylabel("Count")

    all_norm = norm_log2.values.flatten()
    all_norm = all_norm[np.isfinite(all_norm)]
    axes[1, 1].hist(all_norm, bins=40, color="#1f77b4", alpha=0.7, edgecolor="white", linewidth=0.3)
    if len(all_norm) > 5:
        from scipy.ndimage import gaussian_filter1d
        counts, bin_edges = np.histogram(all_norm, bins=40)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        smooth = gaussian_filter1d(counts.astype(float), sigma=1.5)
        axes[1, 1].plot(bin_centers, smooth, color="#4169E1", linewidth=2)
    axes[1, 1].set_title("Normalized (log2) Histogram", fontsize=11, fontweight="bold")
    axes[1, 1].set_xlabel("Value")
    axes[1, 1].set_ylabel("Count")

    fig.tight_layout(pad=3)
    return _fig_to_base64(fig)


def run_full_analysis(df: pd.DataFrame, gene_id_col: str, method: str = "spearman") -> dict:
    """
    Orchestrate all three analysis pipelines and return results as a single dict.
    Automatically detects raw sample columns vs TPM/RPKM normalized columns.
    """
    all_cols = df.columns.tolist()

    # Identify column types
    tpm_cols = [c for c in all_cols if c.startswith("TPM_")]
    rpkm_cols = [c for c in all_cols if c.startswith("RPKM_")]
    
    # Normalized columns: prefer TPM, fall back to RPKM
    norm_cols = tpm_cols if tpm_cols else rpkm_cols

    # Raw sample columns: numeric columns that are NOT gene_id, gene_length, or normalized
    exclude = {gene_id_col, "gene_length_bp"} | set(tpm_cols) | set(rpkm_cols)
    # Also exclude any columns with "length" in the name
    exclude |= {c for c in all_cols if "length" in c.lower()}
    
    raw_cols = [
        c for c in all_cols
        if c not in exclude and pd.api.types.is_numeric_dtype(df[c])
    ]

    results = {}

    # 1. Correlation Heatmap (use raw sample cols)
    if len(raw_cols) >= 2:
        results["correlation"] = generate_correlation_heatmap(df, raw_cols, method)
    else:
        results["correlation"] = None

    # 2. PCA
    if len(raw_cols) >= 2 and len(norm_cols) >= 2:
        results["pca"] = generate_pca_plots(df, raw_cols, norm_cols, gene_id_col)
    else:
        results["pca"] = None

    # 3. Distribution
    if raw_cols and norm_cols:
        results["distribution"] = generate_distribution_plots(df, raw_cols, norm_cols)
    else:
        results["distribution"] = None

    return results
