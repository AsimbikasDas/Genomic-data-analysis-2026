import streamlit as st
import pandas as pd
import polars as pl
import numpy as np
import plotly.express as px
from data_engine import run_deseq2_analysis
from pathway_enrichment import get_enriched_pathways
from llm_integration import generate_biological_summary

# PAGE CONFIG
st.set_page_config(page_title="GenePath AI Dashboard", layout="wide")
st.title("🧬 GenePath AI: Production Dashboard")

# --- CACHED FUNCTIONS ---

@st.cache_data(show_spinner=False)
def load_counts_with_polars(file):
    """Uses the Rust-based Polars engine to instantly load massive CSVs."""
    pl_df = pl.read_csv(file)
    pd_df = pl_df.to_pandas()
    pd_df.set_index(pd_df.columns[0], inplace=True)
    return pd_df

@st.cache_data
def cached_deseq2(counts, metadata, design_factors, p_val, lfc, treated, control):
    return run_deseq2_analysis(
        counts, metadata, 
        design_factors=design_factors, 
        treated_name=treated, 
        control_name=control, 
        p_value_threshold=p_val, 
        lfc_threshold=lfc
    )

@st.cache_data
def convert_ensembl_to_symbol(gene_list):
    """Automatically detects Ensembl IDs and converts them to Gene Symbols."""
    if not gene_list: 
        return []
    if str(gene_list[0]).upper().startswith("ENS"):
        import mygene
        mg = mygene.MyGeneInfo()
        clean_ids = [str(g).split('.')[0] for g in gene_list]
        results = mg.querymany(clean_ids, scopes='ensembl.gene', fields='symbol', species='human', verbose=False)
        symbols = [res['symbol'] for res in results if 'symbol' in res]
        return list(set(symbols))
    return gene_list

@st.cache_data
def add_gene_symbols_to_df(df):
    """Takes a DESeq2 dataframe and adds a human-readable Gene Symbol column for downloads."""
    if df.empty or not str(df.index[0]).upper().startswith("ENS"):
        return df 
        
    import mygene
    mg = mygene.MyGeneInfo()
    clean_ids = [str(g).split('.')[0] for g in df.index]
    results = mg.querymany(clean_ids, scopes='ensembl.gene', fields='symbol', species='human', as_dataframe=True, verbose=False)
    
    df_out = df.copy()
    if 'symbol' in results.columns:
        results = results[~results.index.duplicated(keep='first')] 
        df_out['Gene_Symbol'] = results['symbol'].reindex(clean_ids).values
    else:
        df_out['Gene_Symbol'] = "Unknown"
        
    cols = ['Gene_Symbol'] + [c for c in df_out.columns if c != 'Gene_Symbol']
    return df_out[cols]

def auto_format_datasets(counts_df, metadata_df):
    """An invisible butler that automatically fixes common bioinformatics CSV errors."""
    
    # 1. Clean Metadata Whitespace (e.g., turns " treated " into "treated")
    metadata_df.columns = metadata_df.columns.str.strip().str.lower() # forces lowercase 'condition'
    if 'condition' in metadata_df.columns:
        metadata_df['condition'] = metadata_df['condition'].astype(str).str.strip()
    if 'batch' in metadata_df.columns:
        metadata_df['batch'] = metadata_df['batch'].astype(str).str.strip()

    # 2. THE TRANSPOSE CHECK (Crucial!)
    # If the metadata sample names match the counts COLUMNS, the user uploaded it genes=rows. 
    # We silently auto-transpose it to samples=rows so PyDESeq2 doesn't crash!
    if set(metadata_df.index).intersection(set(counts_df.columns)):
        counts_df = counts_df.T 

    # 3. Align the Datasets
    # Drops any samples that are in the counts but missing from the metadata (or vice versa)
    common_samples = list(set(counts_df.index).intersection(set(metadata_df.index)))
    if len(common_samples) == 0:
        st.error("🚨 Critical Error: The sample names in your metadata DO NOT match the sample names in your counts matrix!")
        st.stop()
        
    counts_df = counts_df.loc[common_samples]
    metadata_df = metadata_df.loc[common_samples]

    # 4. Clean the Math
    counts_df = counts_df.fillna(0) # Convert any empty cells to 0
    
    # Force everything to integers (PyDESeq2 crashes if it sees decimals)
    # We use a try-except block just in case a user left a string (like a gene name) inside the matrix
    try:
        counts_df = counts_df.astype(int)
    except ValueError:
        st.error("🚨 Critical Error: Your counts matrix contains text/letters inside the data cells. It must be numbers only!")
        st.stop()

    return counts_df, metadata_df

# --- SIDEBAR (WRAPPED IN A FORM TO PREVENT RELOADS) ---
with st.sidebar.form(key="settings_form"):
    st.header("⚙️ Analysis Parameters")
    
    # 🧪 Context is now dataset-agnostic
    tissue_type = st.text_input("Experiment Context (Tissue/Disease):", "Airway smooth muscle cells treated with dexamethasone")
    
    p_val_thresh = st.slider("P-value Threshold", 0.001, 0.10, 0.05, 0.005)
    lfc_thresh = st.slider("Log2 Fold Change", 0.0, 3.0, 0.58, 0.1)

    st.markdown("---")
    st.header("🧬 Pathway & AI Settings")
    db_choice = st.selectbox("Pathway Database", ["KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2021", "MSigDB_Hallmark_2020"])
    top_n_pathways = st.slider("Pathways to Analyze", 5, 20, 10)
    detail_level = st.selectbox("AI Report Detail", ["Concise Summary", "Detailed Mechanistic Report", "Clinical/Therapeutic Focus"])

    st.markdown("---")
    st.header("📊 UI Display Settings")
    preview_rows = st.slider("Samples to preview in tables:", 1, 20, 5)
    top_genes_display = st.slider("Top genes for Heatmap:", 10, 100, 30)
    
    # This button locks in the settings without running the heavy math!
    st.form_submit_button("💾 Save Settings")

# --- UI: DATA FORMAT GUIDELINES ---
with st.expander("ℹ️ Data Format Guidelines (Read Before Uploading)"):
    st.markdown("""
    **1. Gene Counts (CSV):** Must have Sample IDs as rows and Gene Symbols/Ensembl as columns. 
    **2. Metadata (CSV):** Must have Sample IDs as rows. Must include a column named `condition`. 
    *Optional:* You can include a `batch` column to correct for technical noise (Batch Effects).
    """)

# --- DATA UPLOAD ---
col1, col2 = st.columns(2)
with col1:
    counts_file = st.file_uploader("Upload Gene Counts (CSV)", type="csv")
with col2:
    metadata_file = st.file_uploader("Upload Metadata (CSV)", type="csv")

if counts_file and metadata_file:
    with st.spinner("Loading massive dataset instantly with Polars..."):
        raw_counts = load_counts_with_polars(counts_file)
        raw_metadata = pd.read_csv(metadata_file, index_col=0) 
        
        #  PASS IT THROUGH THE AUTO-FORMATTER!
        counts_df, metadata_df = auto_format_datasets(raw_counts, raw_metadata)
    
    if "condition" not in metadata_df.columns:
        st.error("🚨 Your metadata MUST contain a column named exactly 'condition'. Please fix your CSV and re-upload.")
        st.stop()

    # BATCH EFFECT DETECTION
    design_factors = ["condition"]
    if "batch" in metadata_df.columns:
        design_factors = ["batch", "condition"]
        st.info("✅ Batch column detected! The math engine will automatically correct for batch effects.")

    # DYNAMIC TREATMENT/CONTROL SELECTION
    st.markdown("### 🧪 Define Your Experiment")
    unique_conditions = metadata_df['condition'].unique().tolist()
    
    t_col1, t_col2 = st.columns(2)
    with t_col1:
        control_group = st.selectbox("Select CONTROL Group (Baseline):", unique_conditions, index=0)
    with t_col2:
        treated_group = st.selectbox("Select TREATED Group (Experimental):", unique_conditions, index=len(unique_conditions)-1)

    
# STEP 1: DYNAMIC RAW DATA INSPECTION
    with st.expander("🔍 Step 1: Raw Data Inspection", expanded=True):
        st.write(f"**Full Dataset Loaded:** Detected {counts_df.shape[0]} samples and {counts_df.shape[1]} genes.")
        
        # Highly visible message explaining the sample view
        st.info(f"💡 **Data Preview:** Displaying the first {preview_rows} samples and first 20 genes to keep your browser lightning fast. Your full dataset is safely loaded in the background engine!")
        
        # Display the sliced dataframe
        st.dataframe(counts_df.iloc[:preview_rows, :20])

    # --- EXECUTION BUTTON ---
    # The pipeline ONLY runs when this button is clicked
    if st.button("🚀 Run Full Pipeline", use_container_width=True):
        
        # STEP 2: MATH & QC VISUALIZATIONS
        st.subheader("🧮 Step 2: DESeq2 Math, QC, & Filtering")
        
        with st.status("Initializing Bioinformatics Pipeline...", expanded=True) as status:
            st.write("⚙️ Filtering low-count noise genes...")
            st.write(f"🧮 Running PyDESeq2 contrasting '{treated_group}' vs '{control_group}'...")
            
            up_genes, down_genes, sig_genes, full_results, norm_counts = cached_deseq2(
                counts_df, metadata_df, design_factors, p_val_thresh, lfc_thresh, treated_group, control_group
            )
            
            st.write("📊 Normalizing read counts for visualizations...")
            status.update(label="Math Pipeline Complete!", state="complete", expanded=False)

        st.success(f"Filtered to {len(up_genes)} Up-regulated and {len(down_genes)} Down-regulated genes.")
        
        # --- PRO VISUALIZATIONS ---
        v_col1, v_col2 = st.columns(2)
        with v_col1:
            st.write("🎯 **PCA Plot (Quality Control)**")
            from sklearn.decomposition import PCA
            pca = PCA(n_components=2)
            pca_results = pca.fit_transform(norm_counts)
            pca_df = pd.DataFrame(pca_results, columns=['PC1', 'PC2'], index=norm_counts.index).join(metadata_df)
            fig_pca = px.scatter(pca_df, x='PC1', y='PC2', color='condition', hover_name=pca_df.index)
            st.plotly_chart(fig_pca, use_container_width=True)
        
        with v_col2:
            st.write("🌋 **Volcano Plot**")
            full_results['-log10(padj)'] = -np.log10(full_results['padj'].fillna(1))
            conditions = [
                (full_results['padj'] < p_val_thresh) & (full_results['log2FoldChange'] > lfc_thresh),
                (full_results['padj'] < p_val_thresh) & (full_results['log2FoldChange'] < -lfc_thresh)
            ]
            full_results['Status'] = np.select(conditions, ['Up-regulated', 'Down-regulated'], default='Not Significant')
            fig_volc = px.scatter(full_results, x='log2FoldChange', y='-log10(padj)', color='Status', hover_name=full_results.index, color_discrete_map={'Up-regulated': 'red', 'Down-regulated': 'blue', 'Not Significant': 'gray'})
            st.plotly_chart(fig_volc, use_container_width=True)
        
        st.write(f"🔥 **Expression Heatmap (Top {top_genes_display} Significant Genes)**")
        if not sig_genes.empty:
            top_sig_genes = sig_genes.sort_values('padj').head(top_genes_display).index
            heatmap_data = norm_counts[top_sig_genes].T 
            fig_heat = px.imshow(heatmap_data, color_continuous_scale='RdBu_r', aspect="auto")
            st.plotly_chart(fig_heat, use_container_width=True)
        else:
            st.warning("No significant genes found to plot on heatmap.")

        # --- TRANSLATE GENES FOR BIOLOGY/AI ---
        with st.status("Translating Ensembl IDs to Gene Symbols...", expanded=True) as t_status:
            translated_up_genes = convert_ensembl_to_symbol(up_genes)
            translated_down_genes = convert_ensembl_to_symbol(down_genes)
            t_status.update(label=f"Translated {len(translated_up_genes)} Up and {len(translated_down_genes)} Down genes!", state="complete", expanded=False)

        # STEP 3: PATHWAY ENRICHMENT
        st.subheader("🧬 Step 3: Pathway Enrichment")
        with st.status(f"Querying {db_choice} Database...", expanded=True) as p_status:
            pathways_df = get_enriched_pathways(translated_up_genes, gene_sets=db_choice, top_n=top_n_pathways)
            p_status.update(label="Pathway Enrichment Complete!", state="complete", expanded=False)
            
        if not isinstance(pathways_df, str) and not pathways_df.empty:
            fig_path = px.bar(pathways_df, x='Adjusted P-value', y='Term', orientation='h', color='Adjusted P-value', color_continuous_scale='Viridis_r')
            fig_path.update_layout(yaxis={'categoryorder':'total descending'})
            st.plotly_chart(fig_path, use_container_width=True)

# STEP 4: LLM INTERPRETATION
        st.subheader("🧠 Step 4: AI Biological Interpretation")
        
        # Combine the user's tissue context with the exact condition names!
        # Example output: "Airway smooth muscle cells (Comparing Asthma vs Healthy)"
        full_context = f"{tissue_type} (Comparing '{treated_group}' vs '{control_group}')"
        
        with st.status(f"Generating {detail_level}...", expanded=True) as ai_status:
            interpretation = generate_biological_summary(
                translated_up_genes, 
                translated_down_genes, 
                pathways_df, 
                tissue_type=full_context, # <-- Pass the combined context here!
                detail_level=detail_level
            )
            ai_status.update(label="AI Report Generated!", state="complete", expanded=False)
        st.markdown(interpretation)
        
        # STEP 5: DOWNLOADS
        st.markdown("---")
        st.subheader("💾 Download Final Reports")
        
        with st.spinner("Formatting files for download..."):
            download_sig_genes = add_gene_symbols_to_df(sig_genes)
            
        d_col1, d_col2, d_col3 = st.columns(3)
        d_col1.download_button("Download Significant Genes", download_sig_genes.to_csv(), "significant_genes.csv", "text/csv")
        if not isinstance(pathways_df, str) and not pathways_df.empty:
            d_col2.download_button("Download Pathways", pathways_df.to_csv(), "enriched_pathways.csv", "text/csv")
        d_col3.download_button("Download AI Report", interpretation, "ai_biological_report.txt", "text/plain")