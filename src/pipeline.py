import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge

# Optional interactive graph
try:
    from pyvis.network import Network
    PYVIS_AVAILABLE = True
except ImportError:
    PYVIS_AVAILABLE = False


# ============================================================
# CONFIG
# ============================================================
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE_PATH = os.path.join(PROJECT_ROOT, "data")
OUT_DIR = os.path.join(PROJECT_ROOT, "results")
NORMAL_FILE = "BC_TCGA-Normal.txt"
TUMOR_FILE = "BC_TCGA-Tumor.txt"
COMBINED_FILE = "combined_expression.csv"
DEG_FILE = "DEG_filtered.csv"
TF_FILE = "Homo_sapiens_TF.txt"
TF_LIST_FILE = "tf_list.csv"
FILTERED_TF_EXPR_FILE = "filtered_TF_expression.csv"
CORR_NETWORK_FILE = "TF_network_correlation.csv"


# ============================================================
# UTILS
# ============================================================
def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def file_path(base_path, filename):
    return os.path.join(base_path, filename)


def safe_read_csv(path, **kwargs):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    return pd.read_csv(path, **kwargs)


def standardize_gene_index(df):
    df.index = df.index.astype(str).str.strip().str.upper()
    df = df.groupby(df.index).mean()
    return df


# ============================================================
# STEP 1: LOAD / CREATE COMBINED EXPRESSION
# ============================================================
def load_or_create_combined_expression(base_path=BASE_PATH):
    print("\n[1] Loading or creating combined expression matrix...")

    combined_path = file_path(base_path, COMBINED_FILE)
    normal_path = file_path(base_path, NORMAL_FILE)
    tumor_path = file_path(base_path, TUMOR_FILE)

    if os.path.exists(combined_path):
        print(f"    Found existing {COMBINED_FILE}")
        combined = pd.read_csv(combined_path, index_col=0)
        combined = standardize_gene_index(combined)
        print(f"    Combined shape: {combined.shape}")
        return combined

    # Otherwise create from tumor + normal
    print(f"    {COMBINED_FILE} not found. Creating from tumor + normal files...")

    normal = safe_read_csv(normal_path, sep="\t", index_col=0)
    tumor = safe_read_csv(tumor_path, sep="\t", index_col=0)

    normal = standardize_gene_index(normal)
    tumor = standardize_gene_index(tumor)

    common_genes = tumor.index.intersection(normal.index)
    tumor = tumor.loc[common_genes]
    normal = normal.loc[common_genes]

    tumor.columns = [f"Tumor_{c}" for c in tumor.columns]
    normal.columns = [f"Normal_{c}" for c in normal.columns]

    combined = pd.concat([tumor, normal], axis=1)
    combined.to_csv(combined_path)

    print(f"    Combined created: {combined.shape}")
    print(f"    Saved: {combined_path}")

    return combined


# ============================================================
# STEP 2: LOAD TF LIST
# ============================================================
def load_tf_list(base_path=BASE_PATH):
    print("\n[2] Loading TF list...")

    tf_list_path = file_path(base_path, TF_LIST_FILE)
    tf_file_path = file_path(base_path, TF_FILE)

    # Prefer tf_list.csv if exists (matches notebook)
    if os.path.exists(tf_list_path):
        try:
            tf_list = pd.read_csv(tf_list_path, header=None)[0].astype(str).str.strip().str.upper().tolist()
            print(f"    Loaded TFs from {TF_LIST_FILE}: {len(tf_list)}")
            return list(dict.fromkeys(tf_list))
        except Exception:
            print(f"    Warning: Could not parse {TF_LIST_FILE}, falling back to {TF_FILE}")

    # Fallback to Homo_sapiens_TF.txt
    if not os.path.exists(tf_file_path):
        raise FileNotFoundError(f"Neither {TF_LIST_FILE} nor {TF_FILE} found in {base_path}")

    # Try full table format first
    try:
        tf_df = pd.read_csv(tf_file_path, sep="\t")
        if "Symbol" in tf_df.columns:
            tf_list = tf_df["Symbol"].dropna().astype(str).str.strip().str.upper().unique().tolist()
            print(f"    Loaded TFs from {TF_FILE} (Symbol column): {len(tf_list)}")

            # Save tf_list.csv for notebook compatibility
            pd.Series(tf_list).to_csv(tf_list_path, index=False, header=False)
            print(f"    Saved derived {TF_LIST_FILE}")

            return tf_list
    except Exception:
        pass

    # Fallback single-column format
    try:
        tf_list = pd.read_csv(tf_file_path, header=None)[0].dropna().astype(str).str.strip().str.upper().tolist()
        tf_list = list(dict.fromkeys(tf_list))
        print(f"    Loaded TFs from {TF_FILE} (single-column): {len(tf_list)}")

        # Save tf_list.csv
        pd.Series(tf_list).to_csv(tf_list_path, index=False, header=False)
        print(f"    Saved derived {TF_LIST_FILE}")

        return tf_list
    except Exception as e:
        raise ValueError(f"Could not parse TF file: {e}")


# ============================================================
# STEP 3: SAVE FILTERED TF EXPRESSION (NOTEBOOK STYLE)
# ============================================================
def create_filtered_tf_expression(combined, tf_list, base_path=BASE_PATH):
    print("\n[3] Creating filtered TF expression matrix...")

    tf_in_data = list(set(tf_list).intersection(combined.index))
    tf_data = combined.loc[tf_in_data].copy()

    out_path = file_path(base_path, FILTERED_TF_EXPR_FILE)
    tf_data.to_csv(out_path)

    print(f"    TFs present in data: {len(tf_in_data)}")
    print(f"    Saved: {out_path}")

    return tf_data, tf_in_data


# ============================================================
# STEP 4: LOAD / COMPUTE DEGs
# ============================================================
def load_or_compute_deg(combined, base_path=BASE_PATH, variance_thresh=0.1, logfc_thresh=0.5, adj_p_thresh=0.05):
    print("\n[4] Loading or computing DEGs...")

    deg_path = file_path(base_path, DEG_FILE)

    if os.path.exists(deg_path):
        print(f"    Found existing {DEG_FILE}")
        deg = pd.read_csv(deg_path)
        deg["gene"] = deg["gene"].astype(str).str.strip().str.upper()
        print(f"    Loaded DEGs: {len(deg)}")
        return deg

    print(f"    {DEG_FILE} not found. Computing DEGs from combined expression...")

    tumor_cols = [c for c in combined.columns if str(c).startswith("Tumor_")]
    normal_cols = [c for c in combined.columns if str(c).startswith("Normal_")]

    if len(tumor_cols) == 0 or len(normal_cols) == 0:
        raise ValueError("Combined expression does not contain Tumor_ / Normal_ columns required for DEG computation.")

    tumor = combined[tumor_cols].fillna(0)
    normal = combined[normal_cols].fillna(0)

    # Variance filter (notebook style)
    variance = combined.var(axis=1)
    filtered_genes = variance[variance > variance_thresh].index

    combined_f = combined.loc[filtered_genes]
    tumor_f = tumor.loc[filtered_genes]
    normal_f = normal.loc[filtered_genes]

    print(f"    Genes after variance filtering: {len(filtered_genes)}")

    results = []

    for gene in combined_f.index:
        tumor_vals = tumor_f.loc[gene]
        normal_vals = normal_f.loc[gene]

        if np.isnan(tumor_vals).any() or np.isnan(normal_vals).any():
            continue

        try:
            t_stat, p_val = ttest_ind(tumor_vals, normal_vals, equal_var=False)
            log_fc = tumor_vals.mean() - normal_vals.mean()  # notebook style
            results.append([gene, log_fc, p_val])
        except Exception:
            continue

    deg = pd.DataFrame(results, columns=["gene", "log2FC", "pval"])

    if deg.empty:
        raise ValueError("DEG computation produced no results.")

    deg["adj_pval"] = multipletests(deg["pval"], method="fdr_bh")[1]

    deg_filtered = deg[
        (deg["log2FC"].abs() > logfc_thresh) &
        (deg["adj_pval"] < adj_p_thresh)
    ].copy()

    deg_filtered.to_csv(deg_path, index=False)

    print(f"    DEGs found: {len(deg_filtered)}")
    print(f"    Saved: {deg_path}")

    return deg_filtered


# ============================================================
# STEP 5: PREPARE TF + TARGET SETS (NOTEBOOK STYLE)
# ============================================================
def prepare_tf_target_sets(combined, deg_filtered, tf_list):
    print("\n[5] Preparing TF and target sets...")

    deg_genes = deg_filtered["gene"].astype(str).str.strip().str.upper().tolist()
    deg_genes = [g for g in deg_genes if g in combined.index]

    # NOTEBOOK STYLE: TFs = TF list ∩ DEG
    tf_genes = list(set(tf_list).intersection(deg_genes))
    tf_genes = [g for g in tf_genes if g in combined.index]

    print(f"    TFs (TF ∩ DEG): {len(tf_genes)}")
    print(f"    DEGs in expression: {len(deg_genes)}")

    if len(tf_genes) == 0:
        raise ValueError("No TF genes found after TF ∩ DEG intersection.")
    if len(deg_genes) == 0:
        raise ValueError("No DEG genes found in expression matrix.")

    return tf_genes, deg_genes


# ============================================================
# STEP 6: CORRELATION BASELINE (NOTEBOOK STYLE)
# ============================================================
def infer_correlation_network(combined, tf_genes, deg_genes, threshold=0.7, base_path=BASE_PATH):
    print("\n[6] Building correlation baseline network...")

    X = combined.loc[tf_genes]
    Y = combined.loc[deg_genes]

    print(f"    X shape (TFs x samples): {X.shape}")
    print(f"    Y shape (DEGs x samples): {Y.shape}")

    X_t = X.T
    Y_t = Y.T

    corr_matrix = np.corrcoef(X_t.values.T, Y_t.values.T)

    n_tf = X.shape[0]
    corr_tf_gene = corr_matrix[:n_tf, n_tf:]

    edges = []

    for i, tf in enumerate(tf_genes):
        for j, gene in enumerate(deg_genes):
            if tf == gene:
                continue

            corr = corr_tf_gene[i, j]

            if not np.isnan(corr) and abs(corr) > threshold:
                edges.append([tf, gene, corr])

    edges_df = pd.DataFrame(edges, columns=["TF", "Target", "Weight"])

    out_path = file_path(base_path, CORR_NETWORK_FILE)
    edges_df.to_csv(out_path, index=False)

    print(f"    Correlation edges: {len(edges_df)}")
    print(f"    Saved: {out_path}")

    return edges_df, corr_tf_gene


# ============================================================
# STEP 7: RIDGE NETWORK (NOTEBOOK STYLE)
# ============================================================
def infer_ridge_network(
    combined,
    tf_genes,
    deg_genes,
    output_dir=OUT_DIR,
    ridge_alpha=1.0,
    coef_thresh=1e-3,
    use_gpu=True
):
    print("\n[7] Building Ridge TF-gene network...")

    ensure_dir(output_dir)

    # Notebook-aligned matrices
    X = combined.loc[tf_genes].T.copy()   # samples × TFs
    Y = combined.loc[deg_genes].T.copy()  # samples × genes

    print(f"    Raw X shape: {X.shape}")
    print(f"    Raw Y shape: {Y.shape}")

    # NaN report
    print(f"    NaNs in X: {np.isnan(X.values).sum()}")
    print(f"    NaNs in Y: {np.isnan(Y.values).sum()}")

    # Fill NaNs column-wise mean
    X = X.fillna(X.mean())
    Y = Y.fillna(Y.mean())

    # Final cleanup
    X = X.fillna(0)
    Y = Y.fillna(0)

    # Remove constant TF columns
    X = X.loc[:, X.std() > 0]

    print(f"    Cleaned X shape: {X.shape}")
    print(f"    Cleaned Y shape: {Y.shape}")

    if X.shape[1] == 0:
        raise ValueError("No non-constant TF columns remain after filtering.")

    # Normalize X only
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    Y_values = Y.values

    print(f"    Final X shape: {X_scaled.shape}")
    print(f"    Final Y shape: {Y_values.shape}")

    coef = None
    model_used = "sklearn Ridge (CPU)"

    # Try GPU Ridge first
    if use_gpu:
        try:
            import cupy as cp
            from cuml.linear_model import Ridge as cuRidge

            print("    Using cuML Ridge (GPU)...")

            X_gpu = cp.asarray(X_scaled, dtype=cp.float32)
            Y_gpu = cp.asarray(Y_values, dtype=cp.float32)

            model = cuRidge(alpha=ridge_alpha)
            model.fit(X_gpu, Y_gpu)

            coef_gpu = model.coef_   # (n_genes, n_tfs)
            coef = cp.asnumpy(coef_gpu)

            model_used = "cuML Ridge (GPU)"
        except Exception as e:
            print(f"    GPU Ridge unavailable, falling back to CPU Ridge. Reason: {e}")

    # CPU fallback
    if coef is None:
        print("    Using sklearn Ridge (CPU)...")

        model = Ridge(alpha=ridge_alpha)
        model.fit(X_scaled, Y_values)

        coef = model.coef_   # (n_genes, n_tfs)

    print(f"    Model used: {model_used}")
    print(f"    Coefficient matrix shape: {coef.shape}")

    # Build edges
    edges = []
    tf_names = list(X.columns)
    gene_names = list(Y.columns)

    for i, gene in enumerate(gene_names):
        for j, tf in enumerate(tf_names):
            weight = float(coef[i, j])

            if abs(weight) > coef_thresh:
                if tf != gene:   # avoid self-loop
                    edges.append([tf, gene, weight])

        if (i + 1) % 500 == 0 or (i + 1) == len(gene_names):
            print(f"    Processed {i + 1}/{len(gene_names)} genes")

    network_df = pd.DataFrame(edges, columns=["TF", "Target", "Weight"])

    out_file = os.path.join(output_dir, "TF_network_ridge.csv")
    network_df.to_csv(out_file, index=False)

    print(f"    Total edges: {len(network_df)}")
    print(f"    Saved: {out_file}")

    return network_df


# ============================================================
# STEP 8: CLEAN NETWORK (NOTEBOOK STYLE)
# ============================================================
def clean_network(network_df, k=3, threshold=0.05):
    print("\n[8] Cleaning network (top-k + threshold)...")

    if network_df.empty:
        print("    Warning: network_df is empty.")
        return network_df.copy()

    network_df = network_df.copy()
    network_df["abs_weight"] = network_df["Weight"].abs()

    # Vectorized notebook-style clean
    network_sorted = network_df.sort_values(["Target", "abs_weight"], ascending=[True, False])
    network_clean = network_sorted.groupby("Target").head(k)
    network_clean = network_clean[network_clean["abs_weight"] > threshold]
    network_clean = network_clean.drop(columns=["abs_weight"])

    print(f"    Before cleaning: {len(network_df)} edges")
    print(f"    After cleaning: {len(network_clean)} edges")

    return network_clean


# ============================================================
# STEP 9: SAVE NETWORK FILES
# ============================================================
def save_network_outputs(network_df, network_clean, out_dir=OUT_DIR):
    print("\n[9] Saving network outputs...")

    ensure_dir(out_dir)

    raw_path = os.path.join(out_dir, "tf_gene_network_raw.csv")
    clean_path = os.path.join(out_dir, "tf_gene_network_clean.csv")

    network_df.to_csv(raw_path, index=False)
    network_clean.to_csv(clean_path, index=False)

    print(f"    Saved raw network: {raw_path}")
    print(f"    Saved clean network: {clean_path}")


# ============================================================
# STEP 10: STATIC VISUALS
# ============================================================
def plot_correlation_distribution(corr_tf_gene, out_dir=OUT_DIR):
    print("\n[10A] Plotting correlation distribution...")

    ensure_dir(out_dir)

    plt.figure(figsize=(8, 5))
    plt.hist(corr_tf_gene.flatten(), bins=100)
    plt.title("Correlation Distribution")
    plt.xlabel("Correlation")
    plt.ylabel("Count")
    out_path = os.path.join(out_dir, "correlation_distribution.png")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"    Saved: {out_path}")


def plot_top_tfs_bar(edges_df, out_dir=OUT_DIR, title="Top TFs (Correlation)"):
    print("\n[10B] Plotting top TF bar chart...")

    ensure_dir(out_dir)

    if edges_df.empty:
        print("    No edges to plot.")
        return

    top_tfs = edges_df["TF"].value_counts().head(10)

    plt.figure(figsize=(8, 5))
    top_tfs.plot(kind="bar")
    plt.title(title)
    plt.xlabel("TF")
    plt.ylabel("Edge Count")
    out_path = os.path.join(out_dir, "top_tfs_bar.png")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"    Saved: {out_path}")


def plot_small_network(edges_df, out_dir=OUT_DIR, max_edges=200, filename="small_network.png"):
    print("\n[10C] Plotting small network...")

    ensure_dir(out_dir)

    if edges_df.empty:
        print("    No edges to plot.")
        return

    sample = edges_df.head(max_edges)
    G = nx.from_pandas_edgelist(sample, "TF", "Target", create_using=nx.DiGraph())

    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(G, pos, node_size=20, with_labels=False, arrows=True)
    plt.title("Small TF-Gene Network")
    out_path = os.path.join(out_dir, filename)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"    Saved: {out_path}")


# ============================================================
# STEP 11: CLEANED NETWORK VISUALIZATION
# ============================================================
def build_focus_graph(network_clean, top_tf_count=10):
    if network_clean.empty:
        return nx.DiGraph(), [], [], [], [], []

    top_tfs = network_clean["TF"].value_counts().head(top_tf_count).index.tolist()
    network_focus = network_clean[network_clean["TF"].isin(top_tfs)].copy()

    G = nx.DiGraph()
    for _, row in network_focus.iterrows():
        G.add_edge(row["TF"], row["Target"], weight=row["Weight"])

    tf_set = set(network_focus["TF"])
    degree_dict = dict(G.degree())

    node_colors = []
    node_sizes = []

    for node in G.nodes():
        if node in tf_set:
            node_colors.append("blue")
        else:
            node_colors.append("orange")
        node_sizes.append(degree_dict[node] * 80 + 100)

    edge_colors = []
    edge_widths = []

    for (_, _, d) in G.edges(data=True):
        if d["weight"] > 0:
            edge_colors.append("green")
        else:
            edge_colors.append("red")
        edge_widths.append(abs(d["weight"]) * 5)

    return G, node_colors, node_sizes, edge_colors, edge_widths, top_tfs


def plot_clean_network(network_clean, out_dir=OUT_DIR):
    print("\n[11] Plotting cleaned biological network...")

    ensure_dir(out_dir)

    if network_clean.empty:
        print("    No cleaned edges to plot.")
        return

    G, node_colors, node_sizes, edge_colors, edge_widths, top_tfs = build_focus_graph(network_clean)

    if G.number_of_edges() == 0:
        print("    Focus graph is empty.")
        return

    try:
        pos = nx.kamada_kawai_layout(G)
    except Exception:
        pos = nx.spring_layout(G, seed=42)

    labels = {node: node for node in top_tfs if node in G.nodes()}

    plt.figure(figsize=(14, 12))

    nx.draw(
        G, pos,
        with_labels=False,
        node_color=node_colors,
        node_size=node_sizes,
        edge_color=edge_colors,
        width=edge_widths,
        alpha=0.8,
        arrows=True
    )

    nx.draw_networkx_labels(G, pos, labels, font_size=10)

    plt.title(
        "Clean TF-Gene Regulatory Network\n"
        "Blue=TF, Orange=Gene | Green=Activation, Red=Repression | Size=Importance"
    )

    out_path = os.path.join(out_dir, "clean_tf_gene_network.png")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"    Saved: {out_path}")


# ============================================================
# STEP 12: INTERACTIVE HTML NETWORK
# ============================================================
def create_interactive_network(network_clean, out_dir=OUT_DIR, max_edges=200):
    print("\n[12] Creating interactive HTML network...")

    ensure_dir(out_dir)

    if not PYVIS_AVAILABLE:
        print("    PyVis not installed. Skipping interactive HTML.")
        print("    Install with: pip install pyvis")
        return

    if network_clean.empty:
        print("    No cleaned edges to visualize.")
        return

    # Use strongest edges only for responsiveness
    network_small = network_clean.reindex(
        network_clean["Weight"].abs().sort_values(ascending=False).index
    ).head(max_edges)

    net = Network(height="750px", width="100%", directed=True, bgcolor="#ffffff", font_color="black")

    tf_set = set(network_small["TF"])

    for _, row in network_small.iterrows():
        tf = row["TF"]
        gene = row["Target"]
        weight = row["Weight"]

        if tf not in net.node_ids:
            net.add_node(tf, label=tf, color="blue", title=f"TF: {tf}")

        if gene not in net.node_ids:
            gene_color = "blue" if gene in tf_set else "orange"
            net.add_node(gene, label=gene, color=gene_color, title=f"Gene: {gene}")

        edge_color = "green" if weight > 0 else "red"
        edge_width = max(1, min(10, abs(weight) * 5))

        net.add_edge(tf, gene, value=edge_width, color=edge_color, title=f"Weight: {weight:.4f}")

    net.repulsion(node_distance=150, central_gravity=0.2, spring_length=150, spring_strength=0.05)

    out_path = os.path.join(out_dir, "interactive_tf_gene_network.html")
    net.save_graph(out_path)

    print(f"    Saved: {out_path}")


# ============================================================
# STEP 13: NETWORK SUMMARY
# ============================================================
def summarize_network(network_df, network_clean, out_dir=OUT_DIR):
    print("\n[13] Summarizing network...")

    ensure_dir(out_dir)

    summary_lines = []

    # Raw network summary
    summary_lines.append("===== RAW NETWORK SUMMARY =====")
    summary_lines.append(f"Total raw edges: {len(network_df)}")

    if not network_df.empty:
        activation = (network_df["Weight"] > 0).sum()
        repression = (network_df["Weight"] < 0).sum()

        summary_lines.append(f"Activation edges: {activation}")
        summary_lines.append(f"Repression edges: {repression}")

        top_tfs = network_df["TF"].value_counts().head(10)
        summary_lines.append("\nTop TFs by edge count:")
        for tf, cnt in top_tfs.items():
            summary_lines.append(f"  {tf}: {cnt}")

        strongest = network_df.reindex(
            network_df["Weight"].abs().sort_values(ascending=False).index
        ).head(10)

        summary_lines.append("\nStrongest TF-Gene interactions:")
        for _, row in strongest.iterrows():
            summary_lines.append(f"  {row['TF']} -> {row['Target']} ({row['Weight']:.4f})")

    # Clean network summary
    summary_lines.append("\n===== CLEAN NETWORK SUMMARY =====")
    summary_lines.append(f"Total clean edges: {len(network_clean)}")

    if not network_clean.empty:
        G = nx.DiGraph()
        for _, row in network_clean.iterrows():
            G.add_edge(row["TF"], row["Target"], weight=row["Weight"])

        out_degree = dict(G.out_degree())
        top_regs = sorted(out_degree.items(), key=lambda x: x[1], reverse=True)[:10]

        summary_lines.append("Top Master Regulators:")
        for tf, deg in top_regs:
            summary_lines.append(f"  {tf}: regulates {deg} genes")

    summary_text = "\n".join(summary_lines)

    print(summary_text)

    out_path = os.path.join(out_dir, "network_summary.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(summary_text)

    print(f"\n    Saved summary: {out_path}")


# ============================================================
# MAIN PIPELINE
# ============================================================
def run_pipeline(
    base_path=BASE_PATH,
    out_dir=OUT_DIR,
    corr_threshold=0.7,
    ridge_coef_thresh=1e-3,
    clean_k=3,
    clean_threshold=0.05
):
    print("=" * 70)
    print("        TF-GENE REGULATORY NETWORK PIPELINE (NOTEBOOK-ALIGNED)")
    print("=" * 70)

    ensure_dir(out_dir)

    # 1. Combined expression
    combined = load_or_create_combined_expression(base_path)

    # 2. TF list
    tf_list = load_tf_list(base_path)

    # 3. Filtered TF expression
    tf_data, tf_in_data = create_filtered_tf_expression(combined, tf_list, base_path)

    # 4. DEGs
    deg_filtered = load_or_compute_deg(combined, base_path)

    # 5. TF + target sets
    tf_genes, deg_genes = prepare_tf_target_sets(
        combined,
        deg_filtered,
        tf_list
    )

    print(f"    FINAL deg_genes count passed to RIDGE: {len(deg_genes)}")

    # 6. Correlation baseline
    corr_edges_df, corr_tf_gene = infer_correlation_network(
        combined,
        tf_genes,
        deg_genes,
        threshold=corr_threshold,
        base_path=base_path
    )

    # 7. RIDGE network
    network_df = infer_ridge_network(
        combined,
        tf_genes,
        deg_genes,
        output_dir=out_dir,
        ridge_alpha=1.0,
        coef_thresh=ridge_coef_thresh,
        use_gpu=True
    )

    # Save notebook-style raw file too
    raw_network_path = file_path(base_path, "tf_gene_network.csv")
    network_df.to_csv(raw_network_path, index=False)
    print(f"\n    Saved notebook-style raw RIDGE network: {raw_network_path}")

    # 8. Clean network
    network_clean = clean_network(
        network_df,
        k=clean_k,
        threshold=clean_threshold
    )

    clean_network_path = file_path(base_path, "tf_gene_network_clean.csv")
    network_clean.to_csv(clean_network_path, index=False)
    print(f"    Saved notebook-style clean network: {clean_network_path}")

    # 9. Save outputs
    save_network_outputs(network_df, network_clean, out_dir)

    # 10. Plots
    plot_correlation_distribution(corr_tf_gene, out_dir)
    plot_top_tfs_bar(corr_edges_df, out_dir, title="Top TFs (Correlation Baseline)")
    plot_small_network(corr_edges_df, out_dir, max_edges=200, filename="correlation_small_network.png")

    # 11. Cleaned network plot
    plot_clean_network(network_clean, out_dir)

    # 12. Interactive HTML
    create_interactive_network(network_clean, out_dir, max_edges=200)

    # 13. Summary
    summarize_network(network_df, network_clean, out_dir)

    print("\n" + "=" * 70)
    print("Pipeline completed successfully ✅")
    print("=" * 70)

    return {
        "combined": combined,
        "deg_filtered": deg_filtered,
        "corr_network": corr_edges_df,
        "ridge_network_raw": network_df,
        "ridge_network_clean": network_clean
    }


# ============================================================
# ENTRY POINT
# ============================================================
if __name__ == "__main__":
    run_pipeline() 
