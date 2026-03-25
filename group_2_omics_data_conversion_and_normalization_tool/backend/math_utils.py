import numpy as np

def compute_rpkm(counts: np.ndarray, lengths_bp: np.ndarray) -> np.ndarray:
    """Calculates Reads Per Kilobase Million (RPKM) utilizing fast NumPy vectorization."""
    total = counts.sum()
    return (counts / (lengths_bp / 1e3)) / (total / 1e6) if total else np.zeros_like(counts)

def compute_tpm(counts: np.ndarray, lengths_bp: np.ndarray) -> np.ndarray:
    """Calculates Transcripts Per Million (TPM) utilizing fast NumPy vectorization."""
    rpk = counts / (lengths_bp / 1e3)
    scale = rpk.sum() / 1e6
    return rpk / scale if scale else np.zeros_like(counts)
