import streamlit as st
import pandas as pd
import numpy as np

# ─────────────────────────────────────────────
# ⚙️ CONFIG
# ─────────────────────────────────────────────
st.set_page_config(page_title="OmicsForge", layout="wide")
st.title("🧬 OmicsForge — Ultra Fast RNA Normalisation")

# ─────────────────────────────────────────────
# LOAD DB
# ─────────────────────────────────────────────
@st.cache_data
def load_gene_lengths():
    df = pd.read_csv("gene_lengths_exonic.csv")
    return dict(zip(df["gene_id"], df["gene_length_bp"]))

gene_length_db = load_gene_lengths()

# ─────────────────────────────────────────────
# FUNCTIONS
# ─────────────────────────────────────────────
def compute_rpkm(counts, lengths_bp):
    total = counts.sum()
    return (counts / (lengths_bp / 1e3)) / (total / 1e6) if total else np.zeros_like(counts)

def compute_tpm(counts, lengths_bp):
    rpk = counts / (lengths_bp / 1e3)
    scale = rpk.sum() / 1e6
    return rpk / scale if scale else np.zeros_like(counts)

def rpkm_to_tpm(rpkm):
    total = rpkm.sum()
    return (rpkm / total) * 1e6 if total else np.zeros_like(rpkm)

def tpm_to_rpkm(tpm):
    total = tpm.sum()
    return (tpm / 1e6) * total if total else np.zeros_like(tpm)

# ─────────────────────────────────────────────
# UPLOAD
# ─────────────────────────────────────────────
uploaded = st.file_uploader("Upload CSV", type=["csv"])

if uploaded:
    df = pd.read_csv(uploaded)
    df = df.loc[:, ~df.columns.str.contains("^Unnamed")]

    st.dataframe(df.head())

    # Detect gene column
    gene_id_col = next((c for c in df.columns if df[c].astype(str).str.startswith("ENS").any()), None)

    if not gene_id_col:
        st.error("No gene column found")
        st.stop()

    df[gene_id_col] = df[gene_id_col].astype(str).str.split(".").str[0]

    numeric_cols = df.select_dtypes(include=np.number).columns.tolist()

    method = st.radio(
        "Method",
        ["TPM only", "RPKM only", "Both TPM & RPKM"],
        horizontal=True
    )

    do_tpm = "TPM" in method
    do_rpkm = "RPKM" in method

    # ─────────────────────────────────────────
    # RUN NORMALISATION
    # ─────────────────────────────────────────
    if st.button("🚀 Run Normalisation"):

        df["gene_length_bp"] = df[gene_id_col].map(gene_length_db)

        df_valid = df.dropna(subset=["gene_length_bp"]).copy()
        lengths = df_valid["gene_length_bp"].values.astype(float)

        for col in numeric_cols:
            counts = df_valid[col].fillna(0).values.astype(float)

            if do_tpm:
                df_valid[f"TPM_{col}"] = compute_tpm(counts, lengths)

            if do_rpkm:
                df_valid[f"RPKM_{col}"] = compute_rpkm(counts, lengths)

        # ✅ STORE RESULT IN SESSION
        st.session_state["df_valid"] = df_valid
        st.session_state["gene_id_col"] = gene_id_col

    # ─────────────────────────────────────────
    # SHOW RESULTS (PERSISTENT)
    # ─────────────────────────────────────────
    if "df_valid" in st.session_state:

        df_valid = st.session_state["df_valid"]
        gene_id_col = st.session_state["gene_id_col"]

        st.success("✅ Normalisation complete")

        st.subheader("📊 Results")
        st.dataframe(df_valid.head())

        st.download_button(
            "⬇ Download Results",
            df_valid.to_csv(index=False),
            "normalized.csv"
        )

        # ─────────────────────────────────────
        # 🔄 CONVERSION (NOW FIXED)
        # ─────────────────────────────────────
        st.markdown("---")
        st.subheader("🔄 Conversion")

        tpm_cols = [c for c in df_valid.columns if c.startswith("TPM_")]
        rpkm_cols = [c for c in df_valid.columns if c.startswith("RPKM_")]

        conversion_type = st.radio(
            "Conversion Type",
            ["RPKM → TPM", "TPM → RPKM"],
            horizontal=True,
            key="conv_type"
        )

        if conversion_type == "RPKM → TPM":

            selected = st.multiselect("Select RPKM columns", rpkm_cols, default=rpkm_cols)

            if st.button("Convert RPKM → TPM"):
                conv_df = df_valid[[gene_id_col]].copy()

                for col in selected:
                    conv_df[col.replace("RPKM", "TPM_conv")] = rpkm_to_tpm(
                        df_valid[col].values
                    )

                st.dataframe(conv_df.head())

        else:

            selected = st.multiselect("Select TPM columns", tpm_cols, default=tpm_cols)

            if st.button("Convert TPM → RPKM"):
                conv_df = df_valid[[gene_id_col]].copy()

                for col in selected:
                    conv_df[col.replace("TPM", "RPKM_conv")] = tpm_to_rpkm(
                        df_valid[col].values
                    )

                st.dataframe(conv_df.head())