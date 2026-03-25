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
    """Core domain pipeline managing data transformation across normalisation strategies."""
    if gene_id_col not in df.columns:
        raise ValueError(f"Selected gene column '{gene_id_col}' is missing from the file.")

    # Strip suffixes off Ensembl IDs
    df[gene_id_col] = df[gene_id_col].astype(str).str.split(".").str[0]
    
    # Try to map lengths from internal DB
    df["gene_length_bp"] = df[gene_id_col].map(gene_length_db)
    
    # [NEW] [Refinement] Fallback Mechanism
    # If lookup failed but the user's file already has a column named GeneLength or Length, use that!
    length_fallback_col = next((c for c in df.columns if "length" in c.lower()), None)
    if length_fallback_col:
        df["gene_length_bp"] = df["gene_length_bp"].fillna(df[length_fallback_col])

    df_valid = df.dropna(subset=["gene_length_bp"]).copy()
    lengths = df_valid["gene_length_bp"].values.astype(float)
    
    numeric_cols = df_valid.select_dtypes(include=np.number).columns.tolist()
    if "gene_length_bp" in numeric_cols:
        numeric_cols.remove("gene_length_bp")
    if length_fallback_col in numeric_cols:
        numeric_cols.remove(length_fallback_col)

    for col in numeric_cols:
        counts = df_valid[col].fillna(0).values.astype(float)
        if is_tpm:
            df_valid[f"TPM_{col}"] = compute_tpm(counts, lengths)
        if is_rpkm:
            df_valid[f"RPKM_{col}"] = compute_rpkm(counts, lengths)

    return df_valid.fillna("").to_dict(orient="records")
