"""
OmicsNorm — One-Click RNA-seq Normalization & Analysis Tool
============================================================
Supports: CSV / TSV / .gz compressed matrices
Gene identifiers: gene symbols OR Ensembl IDs
Normalizations: TPM · RPKM · FPKM · CPM
Analytics: Correlation Heatmap · PCA · Distribution · Outlier Detection
"""

import io, gzip, re, warnings
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings("ignore")

# ─── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="OmicsNorm",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─── Custom CSS ────────────────────────────────────────────────────────────────
st.markdown("""
<style>
  /* ── Base ── */
  [data-testid="stAppViewContainer"] { background: #0d1117; }
  [data-testid="stSidebar"] { background: #161b22; border-right: 1px solid #21262d; }
  .block-container { padding-top: 1.2rem; max-width: 1100px; }

  /* ── Typography ── */
  h1 { color: #e6edf3 !important; font-weight: 700 !important; letter-spacing: -0.02em; }
  h2, h3 { color: #c9d1d9 !important; font-weight: 600 !important; }

  /* ── Metrics ── */
  div[data-testid="stMetric"] {
    background: #161b22; border: 1px solid #21262d;
    border-radius: 8px; padding: 14px 18px;
  }
  div[data-testid="stMetricValue"] { color: #e6edf3 !important; font-size: 1.4rem !important; }
  div[data-testid="stMetricLabel"] { color: #8b949e !important; font-size: 0.78rem !important; text-transform: uppercase; letter-spacing: 0.04em; }

  /* ── Alerts ── */
  .stAlert { border-radius: 6px; }

  /* ── Buttons ── */
  div.stButton > button {
    border-radius: 6px; font-weight: 600;
    background: #21262d; border: 1px solid #30363d; color: #c9d1d9;
    transition: background 0.15s ease;
  }
  div.stButton > button:hover { background: #30363d; color: #e6edf3; }
  [data-testid="stDownloadButton"] button {
    background: #238636; border: 1px solid #2ea043; color: #fff; border-radius: 6px;
  }
  [data-testid="stDownloadButton"] button:hover { background: #2ea043; }

  /* ── Badges ── */
  .badge-raw  { background: #0d419d; color: #79c0ff; padding: 3px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: 600; }
  .badge-norm { background: #0f3d2e; color: #3fb950; padding: 3px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: 600; }

  /* ── Section headers ── */
  .step-header {
    border-left: 3px solid #30363d; padding-left: 12px;
    margin: 2rem 0 1rem 0; color: #e6edf3;
    font-size: 1.05rem; font-weight: 600;
  }

  /* ── Tables ── */
  [data-testid="stDataFrame"] { border-radius: 6px; overflow-x: auto; }

  /* ── Tabs ── */
  button[data-baseweb="tab"] { font-weight: 600 !important; font-size: 0.85rem !important; }
</style>
""", unsafe_allow_html=True)


# ═══════════════════════════════════════════════════════════════════════════════
#  UTILITIES
# ═══════════════════════════════════════════════════════════════════════════════

@st.cache_data(show_spinner=False)
def load_gene_lengths() -> dict:
    """
    Load exonic gene lengths from bundled files.
    Merges two sources into one dict so lookups work for both
    Ensembl IDs (ENSG...) and gene symbols (TSPAN6, TP53, etc.).
    """
    import os
    merged = {}

    # Source 1: Ensembl ID → length  (from group_2 backend)
    ensg_path = os.path.join(os.path.dirname(__file__), "gene_lengths_exonic.csv")
    if os.path.exists(ensg_path):
        df = pd.read_csv(ensg_path)
        ensg_map = dict(zip(
            df.iloc[:, 0].astype(str).str.split(".").str[0],
            df.iloc[:, 1],
        ))
        merged.update(ensg_map)
        print(f"[LENGTHS] Loaded {len(ensg_map)} Ensembl ID → length mappings")

    # Source 2: Gene symbol → length  (pre-built mapping)
    sym_path = os.path.join(os.path.dirname(__file__), "gene_symbol_lengths.csv")
    if os.path.exists(sym_path):
        df2 = pd.read_csv(sym_path)
        sym_map = dict(zip(
            df2.iloc[:, 0].astype(str),
            df2.iloc[:, 1],
        ))
        merged.update(sym_map)
        print(f"[LENGTHS] Loaded {len(sym_map)} gene symbol → length mappings")

    print(f"[LENGTHS] Total merged: {len(merged)} entries")
    return merged


def read_uploaded_file(uploaded) -> pd.DataFrame:
    """Read CSV / TSV / .gz compressed file into a DataFrame."""
    name = uploaded.name.lower()
    raw = uploaded.read()

    if name.endswith(".gz"):
        raw = gzip.decompress(raw)

    # Auto-detect delimiter
    sample = raw[:4096].decode("utf-8", errors="replace")
    sep = "\t" if sample.count("\t") > sample.count(",") else ","

    df = pd.read_csv(io.BytesIO(raw), sep=sep, engine="python")
    # Drop unnamed index columns
    df = df.loc[:, ~df.columns.str.contains("^Unnamed")]
    return df


def detect_gene_id_col(df: pd.DataFrame) -> str | None:
    """Return the most likely gene identifier column."""
    for col in df.columns:
        vals = df[col].astype(str)
        # Ensembl pattern
        if vals.str.match(r"^ENSG\d+").mean() > 0.5:
            return col
        # Gene symbol pattern (mostly letters, some digits)
        if vals.str.match(r"^[A-Za-z][A-Za-z0-9\-_\.]+$").mean() > 0.7:
            return col
    return df.columns[0]


def detect_data_type(df: pd.DataFrame, gene_col: str) -> dict:
    """Heuristically decide if data is raw counts or already normalised."""
    num = df.drop(columns=[gene_col], errors="ignore").select_dtypes(include=np.number)
    if num.empty:
        return {"is_raw": True, "reason": "No numeric columns"}

    flat = num.values.flatten()
    flat = flat[np.isfinite(flat) & (flat >= 0)]

    frac_integers = np.sum(flat == np.round(flat)) / len(flat)
    median_val     = np.median(flat)
    max_val        = np.max(flat)
    variance       = np.var(flat)

    if frac_integers > 0.85 and max_val > 100:
        reason = f"≥85% integer values, median={median_val:.1f} → likely raw counts"
        return {"is_raw": True, "reason": reason}
    elif median_val < 50 and max_val < 1000 and frac_integers < 0.5:
        reason = f"Small decimal values, median={median_val:.3f} → likely already normalised"
        return {"is_raw": False, "reason": reason}
    else:
        reason = f"Mixed signal (median={median_val:.2f}, int_frac={frac_integers:.0%}) — treating as raw"
        return {"is_raw": True, "reason": reason}


def preprocess_dataframe(df: pd.DataFrame, gene_col: str) -> tuple[pd.DataFrame, dict]:
    """
    Clean and validate input data before normalization.

    Steps:
      1. Gene ID column is excluded from all numeric operations.
      2. All sample columns coerced to numeric (errors → NaN).
      3. NaN replaced with 0.
      4. Negative values clipped to 0.
      5. All sample columns cast to float64.

    Validation:
      - Warns if all values across the entire matrix are zero.
      - Errors if >50% of cells across sample columns are NaN before filling.

    Returns:
      (cleaned_df, report_dict) where report_dict contains:
        nan_pct, neg_count, all_zero, n_sample_cols
    """
    work = df.copy()

    # Step 1 — isolate sample columns (gene col is strictly excluded)
    sample_cols = [c for c in work.columns if c != gene_col]

    if not sample_cols:
        return work, {"nan_pct": 0.0, "neg_count": 0, "all_zero": True, "n_sample_cols": 0}

    # Step 2 — coerce to numeric (non-parseable → NaN)
    raw_numeric = work[sample_cols].apply(pd.to_numeric, errors="coerce")

    # Measure NaN rate BEFORE filling (for validation report)
    total_cells = raw_numeric.size
    nan_count   = int(raw_numeric.isna().sum().sum())
    nan_pct     = nan_count / total_cells if total_cells > 0 else 0.0

    # Step 3 — replace NaN with 0
    raw_numeric = raw_numeric.fillna(0)

    # Step 4 — clip negatives to 0
    neg_count = int((raw_numeric < 0).sum().sum())
    raw_numeric = raw_numeric.clip(lower=0)

    # Step 5 — enforce float64
    raw_numeric = raw_numeric.astype(np.float64)

    # Write cleaned sample columns back; gene col stays untouched
    work[sample_cols] = raw_numeric

    all_zero = bool(raw_numeric.values.sum() == 0)

    print(
        f"[PREPROCESS] gene_col={gene_col!r}  sample_cols={len(sample_cols)}  "
        f"nan_pct={nan_pct:.2%}  neg_clipped={neg_count}  all_zero={all_zero}"
    )

    return work, {
        "nan_pct": nan_pct,
        "neg_count": neg_count,
        "all_zero": all_zero,
        "n_sample_cols": len(sample_cols),
    }

# ─── Normalisation math ───────────────────────────────────────────────────────
# All functions guarantee:
#   • Input coerced to float64, NaN → 0
#   • No division by zero (epsilon guard)
#   • Output is float64, same shape as input, contains no NaN / None
#   • Debug min/max printed to stdout

_EPS = 1e-9  # tiny constant to prevent division by zero


def _safe_array(arr) -> np.ndarray:
    """Coerce any array-like to a clean float64 ndarray with NaN replaced by 0."""
    a = pd.to_numeric(pd.Series(arr), errors="coerce").fillna(0).values.astype(np.float64)
    return a


def _debug(label: str, arr: np.ndarray) -> None:
    print(f"[DEBUG] {label:<20s}  min={arr.min():.6g}  max={arr.max():.6g}  "
          f"nan={np.isnan(arr).sum()}  shape={arr.shape}")


def _assert_clean(label: str, arr: np.ndarray, expected_len: int) -> None:
    assert arr.shape == (expected_len,), \
        f"{label}: shape mismatch — expected ({expected_len},), got {arr.shape}"
    assert not np.any(np.isnan(arr)), \
        f"{label}: output contains NaN values"
    assert arr.dtype == np.float64, \
        f"{label}: expected float64, got {arr.dtype}"


def cpm(counts) -> np.ndarray:
    """
    CPM = counts / total_counts_per_sample * 1e6
    Epsilon guard prevents division by zero when library size is 0.
    """
    c = _safe_array(counts)
    total = c.sum()
    result = c / (total + _EPS) * 1e6
    result = np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)
    _debug("CPM", result)
    _assert_clean("CPM", result, len(c))
    return result


def rpkm(counts, lengths_bp) -> np.ndarray:
    """
    RPKM = counts / (length_kb) / (total_millions)
    Epsilon guards on both length_kb and total to avoid /0.
    """
    c = _safe_array(counts)
    L = _safe_array(lengths_bp)
    total = c.sum()
    length_kb = L / 1e3 + _EPS           # avoid /0 on zero-length genes
    total_M   = (total / 1e6) + _EPS     # avoid /0 on empty library
    result = (c / length_kb) / total_M
    result = np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)
    _debug("RPKM", result)
    _assert_clean("RPKM", result, len(c))
    return result


def tpm(counts, lengths_bp) -> np.ndarray:
    """
    TPM:
      1. RPK  = counts / gene_length_kb
      2. sf   = sum(RPK per sample) / 1e6
      3. TPM  = RPK / sf
    Epsilon guard on sf prevents /0 when all RPKs are 0.
    """
    c = _safe_array(counts)
    L = _safe_array(lengths_bp)
    length_kb = L / 1e3 + _EPS
    rpk = c / length_kb
    sf  = rpk.sum() / 1e6 + _EPS
    result = rpk / sf
    result = np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)
    _debug("TPM", result)
    _assert_clean("TPM", result, len(c))
    return result


def fpkm(counts, lengths_bp) -> np.ndarray:
    """FPKM = RPKM for standard single-fragment counting."""
    result = rpkm(counts, lengths_bp)
    _debug("FPKM", result)
    return result


def rpkm_to_tpm(rpkm_vals) -> np.ndarray:
    """Exact mathematical conversion: RPKM → TPM."""
    r = _safe_array(rpkm_vals)
    s = r.sum() + _EPS
    result = r / s * 1e6
    result = np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)
    _debug("RPKM→TPM", result)
    return result


def tpm_to_rpkm_approx(tpm_vals, median_lib: float = 2e7) -> np.ndarray:
    """Approximate TPM → RPKM using assumed total read depth."""
    t = _safe_array(tpm_vals)
    result = (t / 1e6) * (median_lib / 1e6)
    result = np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)
    _debug("TPM→RPKM(approx)", result)
    return result


def normalize_dataframe(
    df: pd.DataFrame,
    gene_col: str,
    gene_length_db: dict,
    do_tpm: bool,
    do_rpkm: bool,
    do_fpkm: bool,
    do_cpm: bool,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """
    Run selected normalisation methods on raw count columns.
    Returns (result_df, raw_cols, normalised_prefix_cols).
    """
    work = df.copy()
    # Strip version suffixes from gene IDs (ENSG000001.3 → ENSG000001)
    work[gene_col] = work[gene_col].astype(str).str.split(".").str[0]

    # Identify numeric sample columns; coerce each to float64 and fill NaN→0
    num_cols = [c for c in work.columns
                if c != gene_col and pd.api.types.is_numeric_dtype(work[c])]
    for col in num_cols:
        work[col] = pd.to_numeric(work[col], errors="coerce").fillna(0).astype(np.float64)

    # Map gene lengths
    lengths_series = work[gene_col].map(gene_length_db)
    has_length = lengths_series.notna().any()

    if (do_tpm or do_rpkm or do_fpkm) and not has_length:
        st.warning("⚠️ No gene lengths found for your gene IDs. CPM only (length-free) will be calculated for TPM/RPKM/FPKM.")

    work["__gene_length__"] = lengths_series.fillna(np.nan)
    valid_mask = (
        work["__gene_length__"].notna()
        if has_length
        else pd.Series([True] * len(work), index=work.index)
    )

    norm_prefixes = []
    for col in num_cols:
        # Always use _safe_array — extra defence against any upstream oddities
        counts  = _safe_array(work.loc[valid_mask, col].values)
        lengths = _safe_array(work.loc[valid_mask, "__gene_length__"].values) if has_length else None

        if do_cpm:
            work.loc[valid_mask, f"CPM_{col}"] = cpm(counts)
            if "CPM" not in norm_prefixes:
                norm_prefixes.append("CPM")
        if do_tpm and lengths is not None:
            work.loc[valid_mask, f"TPM_{col}"] = tpm(counts, lengths)
            if "TPM" not in norm_prefixes:
                norm_prefixes.append("TPM")
        if do_rpkm and lengths is not None:
            work.loc[valid_mask, f"RPKM_{col}"] = rpkm(counts, lengths)
            if "RPKM" not in norm_prefixes:
                norm_prefixes.append("RPKM")
        if do_fpkm and lengths is not None:
            work.loc[valid_mask, f"FPKM_{col}"] = fpkm(counts, lengths)
            if "FPKM" not in norm_prefixes:
                norm_prefixes.append("FPKM")

    work.drop(columns=["__gene_length__"], inplace=True)

    # Guarantee: no None/NaN anywhere in the returned DataFrame before display
    work.fillna(0, inplace=True)

    # Final assertion — must always hold before returning
    assert work is not None, "normalize_dataframe: returned DataFrame is None"
    assert work.shape[0] > 0, "normalize_dataframe: returned DataFrame has 0 rows"

    print(f"[NORMALIZE] output shape={work.shape}  prefixes={norm_prefixes}")
    return work, num_cols, norm_prefixes


# ─── Plotting Helpers ─────────────────────────────────────────────────────────

DARK_BG = "#0d1117"
DARK_AX = "#161b22"
GRID_CLR = "#21262d"
TEXT_CLR = "#c9d1d9"


def _style_ax(ax, title=""):
    ax.set_facecolor(DARK_AX)
    ax.tick_params(colors=TEXT_CLR, labelsize=8)
    for spine in ax.spines.values():
        spine.set_edgecolor(GRID_CLR)
    ax.xaxis.label.set_color(TEXT_CLR)
    ax.yaxis.label.set_color(TEXT_CLR)
    if title:
        ax.set_title(title, color=TEXT_CLR, fontsize=10, fontweight="bold", pad=8)
    ax.grid(True, color=GRID_CLR, linewidth=0.5, alpha=0.6)


def _dark_fig(w, h):
    fig, ax = plt.subplots(figsize=(w, h))
    fig.patch.set_facecolor(DARK_BG)
    ax.set_facecolor(DARK_AX)
    return fig, ax


def plot_correlation_heatmap(df_num: pd.DataFrame, method: str) -> plt.Figure:
    if method == "spearman":
        corr, _ = spearmanr(df_num.values)
        if df_num.shape[1] == 1:
            corr = np.array([[1.0]])
        elif np.ndim(corr) == 0:
            corr = np.array([[float(corr)]])
    else:
        corr = df_num.corr(method="pearson").values

    n = len(df_num.columns)

    # Scale figure size: minimum 10×8, grow with sample count
    w = max(10, min(n * 0.7, 18))
    h = max(8, min(n * 0.6, 15))
    fig, ax = plt.subplots(figsize=(w, h))
    fig.patch.set_facecolor(DARK_BG)
    ax.set_facecolor(DARK_AX)

    # Coolwarm colormap on dark background
    im = ax.imshow(corr, cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")

    # Readable labels — truncate long names, scale font with matrix size
    max_label_len = 18 if n <= 15 else 12
    label_fontsize = max(6, min(9, 120 // n))
    labels = [c[:max_label_len] for c in df_num.columns]

    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=label_fontsize, color=TEXT_CLR)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=label_fontsize, color=TEXT_CLR)

    # Annotate cells only when readable (≤15 samples)
    if n <= 15:
        annot_fontsize = max(5.5, min(8, 100 // n))
        for i in range(n):
            for j in range(n):
                val = corr[i, j]
                text_color = "#e6edf3" if abs(val) > 0.6 else "#21262d"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=annot_fontsize, color=text_color, fontweight="medium")

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.03)
    cbar.ax.tick_params(colors=TEXT_CLR, labelsize=8)
    cbar.outline.set_edgecolor(GRID_CLR)

    # Title
    ax.set_title(
        f"Sample Correlation — {method.capitalize()}",
        color="#e6edf3", fontsize=13, fontweight="bold", pad=14,
    )

    # Clean spines
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.tight_layout(rect=[0, 0, 0.95, 1])
    return fig, corr


def plot_pca(df_raw: pd.DataFrame, df_norm: pd.DataFrame | None) -> plt.Figure:
    ncols = 2 if df_norm is not None else 1
    fig, axes = plt.subplots(1, ncols, figsize=(6 * ncols, 5))
    fig.patch.set_facecolor(DARK_BG)
    if ncols == 1:
        axes = [axes]

    datasets = [("Raw (log₂+1)", df_raw, "#58a6ff")]
    if df_norm is not None:
        datasets.append(("Normalised (log₂+1)", df_norm, "#3fb950"))

    for (label, df, color), ax in zip(datasets, axes):
        ax.set_facecolor(DARK_AX)
        data = np.log2(df.values.astype(float) + 1).T  # samples × genes
        n_comp = min(2, *data.shape)
        if n_comp < 2 or data.shape[0] < 2:
            ax.text(0.5, 0.5, "Need ≥2 samples", ha="center", va="center",
                    transform=ax.transAxes, color=TEXT_CLR)
            ax.set_title(f"PCA — {label}", color=TEXT_CLR, fontsize=10, fontweight="bold")
            continue

        pca_model = PCA(n_components=2)
        coords = pca_model.fit_transform(StandardScaler().fit_transform(data))
        ev = pca_model.explained_variance_ratio_

        ax.scatter(coords[:, 0], coords[:, 1], c=color, s=80,
                   edgecolors="#30363d", linewidth=0.8, zorder=5)

        for i, name in enumerate(df.columns):
            ax.annotate(name[:10], (coords[i, 0], coords[i, 1]),
                        fontsize=6.5, color=TEXT_CLR,
                        xytext=(4, 4), textcoords="offset points")

        ax.set_xlabel(f"PC1 ({ev[0]:.1%})", color=TEXT_CLR, fontsize=9)
        ax.set_ylabel(f"PC2 ({ev[1]:.1%})", color=TEXT_CLR, fontsize=9)
        ax.set_title(f"PCA — {label}", color=TEXT_CLR, fontsize=10, fontweight="bold")
        ax.tick_params(colors=TEXT_CLR, labelsize=7)
        for spine in ax.spines.values():
            spine.set_edgecolor(GRID_CLR)
        ax.grid(True, color=GRID_CLR, linewidth=0.5, alpha=0.5)

    fig.tight_layout(pad=2)
    return fig


def plot_distribution(df_raw: pd.DataFrame, df_norm: pd.DataFrame | None):
    ncols = 2 if df_norm is not None else 1
    fig, axes = plt.subplots(2, ncols, figsize=(6 * ncols, 8))
    fig.patch.set_facecolor(DARK_BG)
    if ncols == 1:
        axes = axes.reshape(2, 1)

    palette = ["#58a6ff", "#3fb950", "#f78166", "#d2a8ff",
               "#ffa657", "#79c0ff", "#56d364", "#ff7b72"]

    def _draw_panel(col_idx, df, label):
        # Boxplot
        ax_box = axes[0, col_idx]
        ax_box.set_facecolor(DARK_AX)
        log_data = [np.log2(df[c].values.astype(float) + 1) for c in df.columns]
        bp = ax_box.boxplot(log_data, patch_artist=True,
                             medianprops={"color": "white", "linewidth": 1.5},
                             whiskerprops={"color": GRID_CLR},
                             capprops={"color": GRID_CLR},
                             flierprops={"marker": "o", "markersize": 2,
                                         "markerfacecolor": "#8b949e", "alpha": 0.3})
        for i, patch in enumerate(bp["boxes"]):
            patch.set_facecolor(palette[i % len(palette)])
            patch.set_alpha(0.85)

        ax_box.set_xticklabels([c[:10] for c in df.columns], rotation=45,
                                ha="right", fontsize=7, color=TEXT_CLR)
        ax_box.set_ylabel("log₂(x+1)", color=TEXT_CLR, fontsize=9)
        ax_box.set_title(f"Boxplot — {label}", color=TEXT_CLR, fontsize=10, fontweight="bold")
        for spine in ax_box.spines.values():
            spine.set_edgecolor(GRID_CLR)
        ax_box.tick_params(axis="y", colors=TEXT_CLR)
        ax_box.grid(True, color=GRID_CLR, linewidth=0.5, alpha=0.5, axis="y")

        # Histogram
        ax_hist = axes[1, col_idx]
        ax_hist.set_facecolor(DARK_AX)
        all_vals = np.log2(df.values.astype(float) + 1).flatten()
        all_vals = all_vals[np.isfinite(all_vals)]
        ax_hist.hist(all_vals, bins=60, color="#58a6ff", alpha=0.75, edgecolor="none")
        ax_hist.set_xlabel("log₂(x+1)", color=TEXT_CLR, fontsize=9)
        ax_hist.set_ylabel("Gene count", color=TEXT_CLR, fontsize=9)
        ax_hist.set_title(f"Histogram — {label}", color=TEXT_CLR, fontsize=10, fontweight="bold")
        ax_hist.tick_params(colors=TEXT_CLR, labelsize=7)
        for spine in ax_hist.spines.values():
            spine.set_edgecolor(GRID_CLR)
        ax_hist.grid(True, color=GRID_CLR, linewidth=0.5, alpha=0.5, axis="y")

    _draw_panel(0, df_raw, "Raw Counts")
    if df_norm is not None:
        _draw_panel(1, df_norm, "Normalised")

    fig.tight_layout(pad=2.5)
    return fig


def detect_outliers_pca(df_num: pd.DataFrame, z_threshold: float = 2.0) -> pd.DataFrame:
    data = np.log2(df_num.values.astype(float) + 1).T
    if data.shape[0] < 3:
        return pd.DataFrame()

    pca_model = PCA(n_components=min(2, *data.shape))
    coords = pca_model.fit_transform(StandardScaler().fit_transform(data))

    centroid = coords.mean(axis=0)
    dists = np.linalg.norm(coords - centroid, axis=1)
    z = (dists - dists.mean()) / (dists.std() + 1e-9)

    result = pd.DataFrame({
        "Sample": df_num.columns,
        "PCA Distance": dists.round(4),
        "Z-Score": z.round(3),
        "Outlier": z > z_threshold,
    }).sort_values("Z-Score", ascending=False)
    return result


# ═══════════════════════════════════════════════════════════════════════════════
#  APP  LAYOUT
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    # ── Header ─────────────────────────────────────────────────────────────────
    st.markdown("# OmicsNorm")
    st.caption("One-click RNA-seq normalization & analysis  ·  TPM · RPKM · FPKM · CPM")
    st.markdown("")  # spacer

    gene_length_db = load_gene_lengths()

    # ── Sidebar ────────────────────────────────────────────────────────────────
    with st.sidebar:
        st.markdown("### Upload Data")
        uploaded = st.file_uploader(
            "Gene expression matrix",
            type=["csv", "tsv", "txt", "gz"],
            help="Rows = genes, Columns = samples. Supports .gz compressed files.",
        )
        st.divider()
        st.markdown("### Normalization")
        do_tpm  = st.checkbox("TPM", value=True)
        do_rpkm = st.checkbox("RPKM", value=True)
        do_fpkm = st.checkbox("FPKM", value=False)
        do_cpm  = st.checkbox("CPM", value=True)

        st.divider()
        st.markdown("### Analysis")
        corr_method = st.selectbox("Correlation method", ["spearman", "pearson"])
        outlier_z = st.slider("Outlier Z-score threshold", 1.5, 4.0, 2.0, 0.1)

        st.divider()
        st.caption("OmicsNorm v1.0")

    # ── Guard ──────────────────────────────────────────────────────────────────
    if uploaded is None:
        st.info("Upload a gene expression matrix using the sidebar to get started.")
        st.markdown("""
        **Accepted formats**
        - Rows = genes (gene symbols or Ensembl IDs)
        - Columns = samples
        - Values = raw counts or normalized floats
        - Compressed `.gz` files supported
        """)
        return

    # ── Load ───────────────────────────────────────────────────────────────────
    with st.spinner("Reading and parsing file..."):
        try:
            df_raw = read_uploaded_file(uploaded)
        except Exception as e:
            st.error(f"Failed to read file: {e}")
            return

        gene_col = detect_gene_id_col(df_raw)
        dtype_info = detect_data_type(df_raw, gene_col)

    # ── Preprocess ─────────────────────────────────────────────────────────────
    with st.spinner("Preprocessing data..."):
        df_raw, prep_report = preprocess_dataframe(df_raw, gene_col)

    # Validation gates — show issues but continue where possible
    if prep_report["nan_pct"] > 0.50:
        st.error(
            f"⛔ {prep_report['nan_pct']:.1%} of cells were non-numeric and replaced with 0. "
            "Check that your file has the correct delimiter and that sample columns contain numbers."
        )
        return
    elif prep_report["nan_pct"] > 0.10:
        st.warning(
            f"⚠️ {prep_report['nan_pct']:.1%} of cells were non-numeric and replaced with 0."
        )

    if prep_report["neg_count"] > 0:
        st.warning(
            f"⚠️ {prep_report['neg_count']} negative value(s) were clipped to 0 before normalization."
        )

    if prep_report["all_zero"]:
        st.warning(
            "⚠️ All sample values are zero after preprocessing. "
            "Normalization results will also be zero — please check your input file."
        )
    sample_cols = [c for c in df_raw.columns
                   if c != gene_col and pd.api.types.is_numeric_dtype(df_raw[c])]
    n_genes    = len(df_raw)
    n_samples  = len(sample_cols)

    # ── File summary ──────────────────────────────────────────────────────────
    badge_cls = "badge-raw" if dtype_info["is_raw"] else "badge-norm"
    badge_lbl = "Raw Counts" if dtype_info["is_raw"] else "Normalised"
    st.markdown(
        f"`{uploaded.name}` &nbsp; "
        f"<span class='{badge_cls}'>{badge_lbl}</span>",
        unsafe_allow_html=True,
    )
    st.markdown("")  # spacer

    c1, c2, c3 = st.columns(3)
    c1.metric("Genes", f"{n_genes:,}")
    c2.metric("Samples", n_samples)
    c3.metric("Gene ID Column", gene_col)

    st.markdown("")  # spacer
    if dtype_info["is_raw"]:
        st.success(f"{dtype_info['reason']} — normalization recommended.")
    else:
        st.info(f"{dtype_info['reason']} — you may still re-normalize.")

    with st.expander("Override gene ID column"):
        gene_col = st.selectbox("Gene identifier column", df_raw.columns.tolist(),
                                index=df_raw.columns.tolist().index(gene_col))

    with st.expander("Data preview"):
        st.dataframe(df_raw.head(5), use_container_width=True)

    st.markdown("")  # spacer
    st.divider()

    # ═══════════════════ NORMALIZATION ════════════════════════════════════════
    st.markdown('<p class="step-header">Step 1 — Normalization</p>', unsafe_allow_html=True)

    if not any([do_tpm, do_rpkm, do_fpkm, do_cpm]):
        st.warning("Select at least one normalization method in the sidebar.")
        return

    length_avail   = bool(gene_length_db)
    needs_length   = do_tpm or do_rpkm or do_fpkm
    length_methods = [m for m, flag in [("TPM", do_tpm), ("RPKM", do_rpkm), ("FPKM", do_fpkm)] if flag]

    if needs_length and not length_avail:
        st.warning(
            f"Gene lengths not found for {length_methods}. "
            "Only CPM (length-free) will be computed unless you upload a gene length file."
        )

    with st.spinner("Processing normalization..."):
        df_out, raw_cols, norm_prefixes = normalize_dataframe(
            df_raw, gene_col, gene_length_db,
            do_tpm, do_rpkm, do_fpkm, do_cpm,
        )

    if not norm_prefixes:
        st.error("No normalizations could be computed. Check gene ID column and method selection.")
        return

    # Validate df_out before any rendering
    assert df_out is not None, "df_out is None after normalization"
    assert df_out.shape[0] > 0, "df_out has 0 rows after normalization"
    df_out.fillna(0, inplace=True)

    st.success(f"Computed **{', '.join(norm_prefixes)}** across {len(raw_cols)} sample columns")

    # Show normalised preview
    all_norm_cols = [c for c in df_out.columns if any(c.startswith(p+"_") for p in norm_prefixes)]
    preview_df = df_out[[gene_col] + all_norm_cols].head(5)
    preview_df = preview_df.fillna(0)  # extra safety for display
    with st.expander("Normalized data preview"):
        st.dataframe(preview_df, use_container_width=True)

    # Download — fill any residual NaN before encoding
    download_df = df_out.fillna(0)
    csv_bytes = download_df.to_csv(index=False).encode()
    st.download_button(
        "Download Normalized CSV",
        data=csv_bytes,
        file_name="normalized_output.csv",
        mime="text/csv",
    )

    # ═══════════════════ LOG TRANSFORMATION ════════════════════════════════════
    st.markdown("")  # spacer
    st.divider()
    st.markdown('<p class="step-header">Step 2 — Log Transformation</p>', unsafe_allow_html=True)

    log_method = st.radio(
        "Transform method",
        ["log2(x + 1)  — recommended, handles zeros safely",
         "log2(x)  — strict, zeros become -∞ → replaced with 0"],
        index=0,
        horizontal=True,
    )

    use_pseudocount = log_method.startswith("log2(x + 1)")

    # Build log-transformed DataFrame from the normalized columns
    log_df = df_out[[gene_col]].copy()

    for col in all_norm_cols:
        vals = pd.to_numeric(df_out[col], errors="coerce").fillna(0).values.astype(np.float64)

        if use_pseudocount:
            # log2(x + 1): always safe, 0 → 0
            log_vals = np.log2(vals + 1)
        else:
            # log2(x): mask zeros and negatives to avoid -inf/NaN
            safe = vals.copy()
            safe[safe <= 0] = np.nan
            log_vals = np.log2(safe)
            log_vals = np.nan_to_num(log_vals, nan=0.0, posinf=0.0, neginf=0.0)

        # Column naming: TPM_R-20 → log2_TPM_R-20
        log_df[f"log2_{col}"] = log_vals.astype(np.float64)

    log_df.fillna(0, inplace=True)
    label = "log₂(x+1)" if use_pseudocount else "log₂(x)"
    st.success(f"Applied **{label}** to {len(all_norm_cols)} normalized columns")

    with st.expander("Log-transformed preview"):
        st.dataframe(log_df.head(5).fillna(0), use_container_width=True)

    log_csv = log_df.fillna(0).to_csv(index=False).encode()
    st.download_button(
        f"Download {label} Transformed CSV",
        data=log_csv,
        file_name=f"log2_transformed_output.csv",
        mime="text/csv",
    )

    st.markdown("")  # spacer
    st.divider()

    # ═══════════════════ ANALYTICS ═══════════════════════════════════════════
    st.markdown('<p class="step-header">Step 3 — Analysis</p>', unsafe_allow_html=True)

    raw_num  = df_raw.set_index(gene_col)[raw_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    # First normalised prefix for comparison
    first_prefix = norm_prefixes[0]
    norm_sample_cols = [f"{first_prefix}_{c}" for c in raw_cols if f"{first_prefix}_{c}" in df_out.columns]
    norm_num = (
        df_out.set_index(gene_col)[norm_sample_cols]
        .apply(pd.to_numeric, errors="coerce").fillna(0)
        if norm_sample_cols else None
    )

    tab_corr, tab_pca, tab_dist, tab_outlier = st.tabs(
        ["Correlation", "PCA", "Distribution", "Outliers"]
    )

    with tab_corr:
        st.caption(f"Sample-to-sample {corr_method} correlation on raw counts")
        if len(raw_cols) < 2:
            st.info("Need at least 2 sample columns for correlation.")
        else:
            with st.spinner("Running correlation analysis..."):
                fig_corr, corr_matrix = plot_correlation_heatmap(raw_num, corr_method)
            st.pyplot(fig_corr, use_container_width=False)

            corr_df = pd.DataFrame(corr_matrix, index=raw_cols, columns=raw_cols)
            st.download_button(
                f"Download Correlation Matrix",
                data=corr_df.to_csv().encode(),
                file_name=f"correlation_matrix_{corr_method}.csv",
                mime="text/csv",
            )

    with tab_pca:
        if len(raw_cols) < 2:
            st.info("Need at least 2 samples for PCA.")
        else:
            with st.spinner("Computing principal components..."):
                fig_pca = plot_pca(raw_num, norm_num)
            st.pyplot(fig_pca, use_container_width=True)

    with tab_dist:
        with st.spinner("Generating distribution plots..."):
            fig_dist = plot_distribution(raw_num, norm_num)
        st.pyplot(fig_dist, use_container_width=True)

    with tab_outlier:
        st.caption(f"Threshold: Z-score > {outlier_z} (Euclidean distance from PCA centroid)")
        if len(raw_cols) < 3:
            st.info("Need 3 or more samples for outlier detection.")
        else:
            with st.spinner("Running outlier detection..."):
                outlier_df = detect_outliers_pca(raw_num, z_threshold=outlier_z)

            flagged = outlier_df[outlier_df["Outlier"]]
            if flagged.empty:
                st.success("No outlier samples detected at this threshold.")
            else:
                st.warning(f"{len(flagged)} outlier sample(s) found:")
                st.dataframe(flagged, use_container_width=True, hide_index=True)

            with st.expander("All sample Z-scores"):
                st.dataframe(outlier_df, use_container_width=True, hide_index=True)


if __name__ == "__main__":
    main()
