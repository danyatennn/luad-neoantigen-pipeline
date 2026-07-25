"""
Part 6 - Gene-Expression Dataset (Lung Cancer, GSE81089)

Install requirements:
    pip install pandas requests --break-system-packages
"""

import pandas as pd
import numpy as np
import requests
import gzip
import os
import matplotlib.pyplot as plt

FPKM_URL = ("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81089"
            "&format=file&file=GSE81089%5FFPKM%5Fcufflinks%2Etsv%2Egz")
DESTDIR = "./geo_data"
RAW_FILE = os.path.join(DESTDIR, "GSE81089_FPKM_cufflinks.tsv.gz")
OUTPUT_FILE = "02_gene_by_sample_TPM.tsv"

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
    # 1. download FPKM matrix
    download_file(FPKM_URL, RAW_FILE)

    # 2. load matrix
    with gzip.open(RAW_FILE, "rt") as f:
        expr = pd.read_csv(f, sep="\t", index_col=0)

    print(f"Loaded raw FPKM matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")
    print("First few columns:", list(expr.columns[:5]))
    print("Index name (gene identifier column):", expr.index.name)

    # 3. Split tumor vs normal using the sample-name suffix convention
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

    # 4. Convert FPKM --> TPM (per-sample renormalization to sum to 1e6)
    print("\nConverting FPKM to TPM...")
    tpm = working.div(working.sum(axis=0), axis=1) * 1e6

    # 5. Handle duplicated gene symbols: collapse by MEDIAN across rows mapping to the same gene symbol
    tpm.index.name = "GeneName"
    tpm = tpm.groupby("GeneName").median()

    # 6. QC printout 
    n_zero_genes = (tpm.sum(axis=1) == 0).sum()
    missing_pct = tpm.isna().mean().mean() * 100
    print(f"\nFinal TPM matrix: {tpm.shape[0]} genes x {tpm.shape[1]} samples")
    print(f"Genes with zero expression across all samples: {n_zero_genes}")
    print(f"Missing values: {missing_pct:.2f}%")

    # 7. Write csv
    tpm.reset_index().to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"\nWrote {OUTPUT_FILE}")

    # 8. Make histogram & calculate distribution of TPM vals
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
    plt.savefig("tpm_distribution_log.png", dpi=150)
    plt.show()

if __name__ == "__main__":
    main()