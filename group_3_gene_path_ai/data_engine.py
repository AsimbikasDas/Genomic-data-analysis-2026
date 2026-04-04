import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet 
from pydeseq2.ds import DeseqStats 

def run_deseq2_analysis(counts_df, metadata_df, design_factors=["condition"], treated_name="B", control_name="A", p_value_threshold=0.05, lfc_threshold=0.58):
    """
    Runs PyDESeq2 with filtering, batch effect correction, and extracts normalized counts.
    """
    # 1. PRE-FILTERING: Drop genes with less than 10 total reads across all samples
    counts_df = counts_df.loc[:, counts_df.sum(axis=0) >= 10]
    
    # 2. RUN MATH & BATCH CORRECTION
    # PyDESeq2 handles batch effects if you pass multiple columns (e.g., ['batch', 'condition'])
    dds = DeseqDataSet(counts=counts_df, metadata=metadata_df, design_factors=design_factors, n_cpus=1) 
    dds.deseq2() 
    
    # The primary biological condition must be the one we calculate the contrast for
    primary_condition = design_factors[-1] if isinstance(design_factors, list) else design_factors
    
    stat_res = DeseqStats(dds, contrast=[primary_condition, treated_name, control_name], n_cpus=1) 
    stat_res.summary() 
    results_df = stat_res.results_df 
    
    # 3. EXTRACT NORMALIZED COUNTS (For PCA and Heatmap)
    # Get the mathematically normalized counts from PyDESeq2
    base_norm_counts = dds.layers['normed_counts']
    
    # Apply a log1p transformation (log(x + 1)) to scale the data perfectly for visual charts
    log_norm_counts = np.log1p(base_norm_counts)
    
    # Convert it back into a nice Pandas DataFrame
    normalized_counts = pd.DataFrame(log_norm_counts, index=dds.obs_names, columns=dds.var_names)
    
    # 4. FILTER SIGNIFICANT GENES
    significant_genes = results_df[(results_df.padj < p_value_threshold) & (abs(results_df.log2FoldChange) > lfc_threshold)] 
    
    up_regulated = significant_genes[significant_genes.log2FoldChange > 0].index.tolist()
    down_regulated = significant_genes[significant_genes.log2FoldChange < 0].index.tolist()
    
    return up_regulated, down_regulated, significant_genes, results_df, normalized_counts