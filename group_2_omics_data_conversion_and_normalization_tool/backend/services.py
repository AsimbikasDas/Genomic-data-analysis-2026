import io
import pandas as pd
import numpy as np

from math_utils import compute_rpkm, compute_tpm
from data_loader import gene_length_db

def parse_csv_file(content: bytes) -> pd.DataFrame:
    """Load raw bytes into a Pandas DataFrame and sanitise standard unnamed anomalies."""
    df = pd.read_csv(io.BytesIO(content))
    return df.loc[:, ~df.columns.str.contains("^Unnamed")]

def get_preview_data(df: pd.DataFrame) -> dict:
    """Isolate metadata and sample rows to send back a light payload before intensive math."""
    gene_id_col = next((c for c in df.columns if df[c].astype(str).str.startswith("ENS").any()), None)
    return {
        "columns": df.columns.tolist(),
        "gene_id_col": gene_id_col,
        "preview": df.head(5).fillna("").to_dict(orient="records")
    }

def process_normalization(df: pd.DataFrame, gene_id_col: str, is_tpm: bool, is_rpkm: bool) -> list:
    """Optimized domain pipeline with memory-safe operations for large-scale omics files."""
    if gene_id_col not in df.columns:
        raise ValueError(f"Selected gene column '{gene_id_col}' is missing.")

    # 1. Clean Gene IDs efficiently
    df[gene_id_col] = df[gene_id_col].astype(str).str.split(".").str[0]
    
    # 2. Map lengths and handle fallback without full dataframe copies
    df["gene_length_bp"] = df[gene_id_col].map(gene_length_db)
    
    length_fallback_col = next((c for c in df.columns if "length" in c.lower() and c != "gene_length_bp"), None)
    if length_fallback_col:
        df["gene_length_bp"] = df["gene_length_bp"].fillna(df[length_fallback_col])

    # 3. Filter in-place (mostly) by dropping rows with missing length dependencies
    # We drop from the original 'df' to save memory rather than creating a 'df_valid' copy
    df.dropna(subset=["gene_length_bp"], inplace=True)
    
    # Cast lengths to float32 to halve memory usage compared to float64
    lengths = df["gene_length_bp"].values.astype(np.float32)
    
    # Identification of numeric columns for vectorization
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if "gene_length_bp" in numeric_cols:
        numeric_cols.remove("gene_length_bp")
    if length_fallback_col and length_fallback_col in numeric_cols:
        numeric_cols.remove(length_fallback_col)

    for col in numeric_cols:
        # Cast input counts to float32 for processing
        counts = df[col].fillna(0).values.astype(np.float32)
        
        if is_tpm:
            # Resulting columns also stored as float32
            df[f"TPM_{col}"] = compute_tpm(counts, lengths).astype(np.float32)
        if is_rpkm:
            df[f"RPKM_{col}"] = compute_rpkm(counts, lengths).astype(np.float32)

    # 4. Final preparation: fill missing values before returning
    return df.fillna("")


