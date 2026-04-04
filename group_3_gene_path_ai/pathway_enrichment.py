import gseapy as gp
import pandas as pd

# Updated to accept gene_sets and top_n as arguments
def get_enriched_pathways(gene_list, gene_sets='KEGG_2021_Human', top_n=10):
    """
    Takes a list of significant genes and finds enriched biological pathways dynamically.
    """
    if not gene_list:
        return "No significant genes provided for pathway enrichment."

    try:
        enr = gp.enrichr(gene_list=gene_list, 
                         gene_sets=gene_sets, # Uses the user's choice
                         organism='human', 
                         outdir=None) 
        
        results_df = enr.results
        significant_pathways = results_df[results_df['Adjusted P-value'] < 0.05]
        top_pathways = significant_pathways.head(top_n)[['Term', 'Adjusted P-value', 'Overlap']]
        
        return top_pathways
        
    except Exception as e:
        return f"Error during pathway enrichment: {e}"