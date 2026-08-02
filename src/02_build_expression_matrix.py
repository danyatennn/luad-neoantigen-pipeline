"""
Gene-Expression Dataset (Lung Cancer, GSE81089)
"""

import pandas as pd
import numpy as np
import requests
import gzip
import os
import matplotlib.pyplot as plt
import config

FPKM_URL = config.GEO_FPKM_URL
RAW_FILE = config.GEO_FPKM
OUTPUT_FILE = config.TPM_MATRIX

TUMOR_SUFFIX = "T" # sample columns ending in this = tumor
NORMAL_SUFFIX = "N" # sample columns ending in this = matched normal
KEEP_ONLY_TUMOR = True  # deliverable matrix should be only tumor samples


def download_file(url, path):
    if os.path.exists(path):
        print(f"Already downloaded: {path}")
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print(f"Downloading {url} ...")
    r = requests.get(url, stream=True)
    r.raise_for_status()
    with open(path, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    print(f"Saved to {path}")


def main():
    # download FPKM matrix
    download_file(FPKM_URL, RAW_FILE)

    # load matrix
    with gzip.open(RAW_FILE, "rt") as f:
        expr = pd.read_csv(f, sep="\t", index_col=0)

    print(f"Loaded raw FPKM matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")
    print("First few columns:", list(expr.columns[:5]))
    print("Index name (gene identifier column):", expr.index.name)

    # Split tumor vs normal using the sample-name suffix convention
    tumor_cols, normal_cols, unmatched = [], [], []
    for col in expr.columns:
        base = col.split("_")[0]  # strip suffixes like _2122, _1
        if base.endswith(TUMOR_SUFFIX):
            tumor_cols.append(col)
        elif base.endswith(NORMAL_SUFFIX):
            normal_cols.append(col)
        else:
            unmatched.append(col)

    print(f"\nTumor samples: {len(tumor_cols)}")
    print(f"Normal samples: {len(normal_cols)}")
    if unmatched:
        print(f"Unmatched columns (check manually): {unmatched}")

    working = expr[tumor_cols] if KEEP_ONLY_TUMOR else expr[tumor_cols + normal_cols]

    # Convert FPKM --> TPM
    print("\nConverting FPKM to TPM...")
    tpm = working.div(working.sum(axis=0), axis=1) * 1e6

    # Handle duplicated gene symbols: collapse by MEDIAN across rows mapping to the same gene symbol
    tpm.index.name = "GeneName"
    tpm = tpm.groupby("GeneName").median()

    # QC printout 
    n_zero_genes = (tpm.sum(axis=1) == 0).sum()
    missing_pct = tpm.isna().mean().mean() * 100
    print(f"\nFinal TPM matrix: {tpm.shape[0]} genes x {tpm.shape[1]} samples")
    print(f"Genes with zero expression across all samples: {n_zero_genes}")
    print(f"Missing values: {missing_pct:.2f}%")

    # Write csv
    tpm.reset_index().to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"\nWrote {OUTPUT_FILE}")

    # Make histogram & calculate distribution of TPM vals
    tpm_values = tpm.values.flatten()

    median_tpm = pd.Series(tpm_values).median()
    variance_tpm = pd.Series(tpm_values).var()
    print(pd.Series(tpm_values).describe(percentiles=[0.5, 0.9, 0.95, 0.99]))
    print("Median TPM:", median_tpm)
    print("Variance of TPM:", variance_tpm)

    # histogram
    log_tpm_values = np.log2(tpm_values + 1)

    plt.figure(figsize=(8, 5))
    plt.hist(log_tpm_values, bins=100)
    plt.xlabel("log2(TPM + 1)")
    plt.ylabel("Frequency")
    plt.title("Distribution of GeneLevelTPM values (log2 scale)")
    plt.savefig(config.FIG_TPM_HISTOGRAM, dpi=150)
    plt.show()

if __name__ == "__main__":
    main()