"""
Part 5 - Mutation-by-Sample Matrix

Install requirements:
    pip install pandas requests --break-system-packages
"""

import pandas as pd
import numpy as np
import requests
import gzip
import io
import json
import os

# change as you wish
PROJECT = "TCGA-LUAD" 
WORKFLOW = "MuTect2 Variant Aggregation and Masking"
DESTDIR = "./gdc_data"
RAW_MAF = os.path.join(DESTDIR, f"{PROJECT}_mutect2.maf.gz")
OUTPUT_FILE = "01_mutation_by_sample_matrix.tsv"

# paste MAF filepath here
LOCAL_MAF_PATH = r"C:\Users\lyssa\Documents\GitHub\project130_lung_cancer\cohortMAF.2026-07-20.maf.gz"

GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
GDC_DATA_ENDPOINT = "https://api.gdc.cancer.gov/data"

def find_maf_file_id():
    # Query the GDC API for the project's aggregated MuTect2 MAF file.
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [PROJECT]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Masked Somatic Mutation"]}},
            {"op": "in", "content": {"field": "data_format", "value": ["MAF"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": [WORKFLOW]}},
        ],
    }
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,file_size",
        "format": "JSON",
        "size": "10",
    }
    r = requests.get(GDC_FILES_ENDPOINT, params=params)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    if not hits:
        raise RuntimeError(
            f"No MAF file found for {PROJECT} / {WORKFLOW}. "
            "Check https://portal.gdc.cancer.gov/repository and search manually."
        )
    print(f"Found file: {hits[0]['file_name']} ({hits[0]['file_size']/1e6:.1f} MB)")
    return hits[0]["file_id"]


def download_maf(file_id, path):
    if os.path.exists(path):
        print(f"Already downloaded: {path}")
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print(f"Downloading file {file_id} ...")
    r = requests.get(f"{GDC_DATA_ENDPOINT}/{file_id}", stream=True)
    r.raise_for_status()
    with open(path, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    print(f"Saved to {path}")


def load_maf(path):
    """MAF files have '#'-prefixed comment/version lines before the header."""
    with gzip.open(path, "rt") as f:
        lines = [l for l in f if not l.startswith("##")]
    maf = pd.read_csv(io.StringIO("".join(lines)), sep="\t", low_memory=False)
    return maf


def main():
    # 1. Locate + download the aggregated project MAF, UNLESS a local
    #    file path was provided above, in which case use that instead.
    if LOCAL_MAF_PATH:
        maf_path = LOCAL_MAF_PATH
        print(f"Using local MAF file: {maf_path}")
        if not os.path.exists(maf_path):
            raise FileNotFoundError(
                f"Could not find {maf_path} — check the path (and that "
                "it's reachable from wherever you're running this script)."
            )
    else:
        file_id = find_maf_file_id()
        download_maf(file_id, RAW_MAF)
        maf_path = RAW_MAF

    # 2. Load
    maf = load_maf(maf_path)
    n_total = len(maf)
    print(f"\nLoaded MAF: {n_total} total variant records")

    # 3. Filter (documented per Part 4)
    filt = maf.copy()
    if "FILTER" in filt.columns:
        filt = filt[filt["FILTER"] == "PASS"]
    filt = filt[filt["Variant_Classification"] == "Missense_Mutation"]
    filt = filt[filt["Gene"].notna() | filt["Hugo_Symbol"].notna()]
    filt = filt[filt["HGVSp_Short"].notna()]

    n_filtered = len(filt)
    print(f"After filtering (PASS, missense, protein-coding): {n_filtered} records")
    print(f"Unique genes: {filt['Hugo_Symbol'].nunique()}")
    print(f"Unique tumor samples: {filt['Tumor_Sample_Barcode'].nunique()}")

    # 4. Build long-format table with the assignment's required columns.
    long_df = filt[[
        "Hugo_Symbol", "HGVSc", "HGVSp_Short", "Tumor_Sample_Barcode"
    ]].rename(columns={
        "Hugo_Symbol": "Gene_Name",
        "HGVSc": "Mutation",
        "HGVSp_Short": "AminoAcid_Change",
        "Tumor_Sample_Barcode": "Sample",
    }).drop_duplicates()

    # 5. Pivot to binary Gene x Sample matrix. 
    long_df["present"] = 1
    matrix = long_df.pivot_table(
        index=["Gene_Name", "Mutation", "AminoAcid_Change"],
        columns="Sample",
        values="present",
        fill_value=0,
    ).reset_index()

    print(f"\nFinal mutation-by-sample matrix: {matrix.shape[0]} mutations x "
          f"{matrix.shape[1]-3} samples")

    # 6. Write tab-delimited output
    matrix.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"Wrote {OUTPUT_FILE}")

    # 7. for part 8?
    print(f"Mutations before filtering: {n_total}")
    print(f"Mutations after filtering: {n_filtered}")
    print(f"Unique genes: {filt['Hugo_Symbol'].nunique()}")
    print(f"Tumor samples: {filt['Tumor_Sample_Barcode'].nunique()}")
    print("Reference genome assembly: GRCh38/hg38 (GDC harmonized)")


if __name__ == "__main__":
    main()