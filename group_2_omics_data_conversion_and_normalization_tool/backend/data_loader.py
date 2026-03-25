import pandas as pd

def load_gene_lengths(file_path: str = "gene_lengths_exonic.csv") -> dict:
    """In-memory loaded singleton database mapping Ensembl Gene IDs to physical length in base-pairs."""
    try:
        df = pd.read_csv(file_path)
        return dict(zip(df["gene_id"], df["gene_length_bp"]))
    except Exception as e:
        print(f"Warning: Failed to load gene lengths database -> {e}")
        return {}

gene_length_db = load_gene_lengths()
