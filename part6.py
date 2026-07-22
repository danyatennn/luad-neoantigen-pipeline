"""
Install requirements first:
    pip install GEOparse pandas --break-system-packages

Usage:
    python 02_process_expression.py
"""

import GEOparse
import pandas as pd
import numpy as np
import os

GEO_ACCESSION = "GSE81089" # dataset number
DESTDIR = "./geo_data" # where raw GEO files will be cached
OUTPUT_FILE = "02_gene_by_sample_TPM.tsv"

TUMOR_SAMPLES_ONLY = True

def main():
    os.makedirs(DESTDIR, exist_ok=True)

    # 1. Download series (metadata + supplementary files) from GEO
    print(f"Downloading {GEO_ACCESSION} from GEO ...")
    gse = GEOparse.get_GEO(geo=GEO_ACCESSION, destdir=DESTDIR)

    # 2. Inspect sample metadata to find the tumor/normal label. Print a sample's metadata to see what field marks tumor vs normal.
    example_gsm = list(gse.gsms.values())[0]
    print("\nExample sample metadata (use this to find the tumor/normal field):")
    for key, val in example_gsm.metadata.items():
        print(f"  {key}: {val}")

    # Adjust the keyword ("tumor" / "normal") below if the actual wording differs once you've looked at the printout above.
    sample_group = {}
    for gsm_name, gsm in gse.gsms.items():
        chars = " ".join(gsm.metadata.get("characteristics_ch1", [])).lower()
        title = gsm.metadata.get("title", [""])[0].lower()
        text = chars + " " + title
        if "normal" in text or "adjacent" in text:
            sample_group[gsm_name] = "normal"
        else:
            sample_group[gsm_name] = "tumor"

    n_tumor = sum(v == "tumor" for v in sample_group.values())
    n_normal = sum(v == "normal" for v in sample_group.values())
    print(f"\nDetected {n_tumor} tumor samples and {n_normal} normal samples.")

    # 3. Build the raw expression matrix (genes x samples) from each GSM's table. 
    if len(example_gsm.table) > 0:
        print("\nUsing per-sample GSM tables...")
        expr = gse.pivot_samples("VALUE")   # genes x samples, FPKM values
    else:
        print("\nNo per-sample tables found — check gse.gpls / gse.metadata")
        print("for a supplementary series matrix / FPKM file and load it")
        print("with pandas.read_csv() instead, e.g.:")
        print('  expr = pd.read_csv("<downloaded_supp_file>", sep="\\t", index_col=0)')
        raise SystemExit(
            "Manual step needed: locate the supplementary FPKM file. "
            "Check ./geo_data/ after this run, or the GEO page's "
            "'Series Matrix File(s)' / 'Supplementary file' links."
        )

    print(f"Raw expression matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")

    # 4. Restrict to tumor samples (optional, see flag above)
    if TUMOR_SAMPLES_ONLY:
        tumor_cols = [s for s in expr.columns if sample_group.get(s) == "tumor"]
        expr = expr[tumor_cols]
        print(f"Restricted to {len(tumor_cols)} tumor samples")

    # 5. Convert FPKM --> TPM 
    print("\nConverting FPKM to TPM...")
    tpm = expr.div(expr.sum(axis=0), axis=1) * 1e6

    # 6. Handle duplicated gene symbols: take the MEDIAN across probes/
    #    rows mapping to the same gene symbol (report this choice in
    #    your README, per Part 6 requirements)
    tpm["GeneName"] = tpm.index  # assumes index is already gene symbol;
    tpm = tpm.groupby("GeneName").median()

    # 7. Basic QC printout (some of what Part 8 will ask for later)
    missing_pct = tpm.isna().mean().mean() * 100
    print(f"\nFinal TPM matrix: {tpm.shape[0]} genes x {tpm.shape[1]} samples")
    print(f"Missing values: {missing_pct:.2f}%")

    # 8. Write tab-delimited output
    tpm.reset_index().rename(columns={"GeneName": "GeneName"}).to_csv(
        OUTPUT_FILE, sep="\t", index=False
    )
    print(f"\nWrote {OUTPUT_FILE}")


if __name__ == "__main__":
    main()