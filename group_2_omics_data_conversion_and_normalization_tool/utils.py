import numpy as np
import pandas as pd

def counts_to_tpm(df_counts, gene_length):

    rpk = df_counts.div(gene_length/1000, axis=0)

    scaling = rpk.sum()

    tpm = rpk.div(scaling) * 1e6

    return tpm


def counts_to_rpkm(df_counts, gene_length):

    total_reads = df_counts.sum()

    rpkm = (1e9 * df_counts).div(gene_length, axis=0).div(total_reads)

    return rpkm


def counts_to_cpm(df_counts):

    total_reads = df_counts.sum()

    cpm = df_counts.div(total_reads) * 1e6

    return cpm


def log_transform(df):

    return np.log2(df + 1)