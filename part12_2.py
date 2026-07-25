"""
Part 12.2 - Peptide-MHC Binding Prediction: Class II, 15-mers (IEDB REST API / NetMHCIIpan)
Calls IEDB's public RESTful API, which runs the real NetMHCIIpan tool on
their servers -- no local install, no TensorFlow/Keras, no OS-specific
issues. Runs directly from any machine with internet access, including
plain Windows + Python.

Docs: https://tools.iedb.org/main/tools-api/

Install requirements:
    pip install pandas requests tqdm --break-system-packages
"""

import pandas as pd
import requests
import time
import io
from tqdm import tqdm

PEPTIDE_FILE = r"C:\Users\lyssa\Downloads\peptides_15mers.tsv" # filtered part 10 output (15-mers only)
OUTPUT_FILE = "06_classII_15mer_predictions.tsv"
NEOANTIGEN_CANDIDATES_FILE = "06_classII_neoantigen_candidates.tsv"

HLA_CLASS_II = [
    "HLA-DRB1*07:01",
    "HLA-DRB1*15:01",
]

IEDB_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
METHOD = "netmhciipan_ba"   # binding affinity mode, consistent with Part 12.1

BATCH_SIZE = 200            # peptides per request -- keeps requests reasonably sized
REQUEST_DELAY = 1.0         # seconds between requests -- be polite to the shared server

STRONG_BINDER_RANK = 2.0
WEAK_BINDER_RANK = 10.0
CANDIDATE_MUT_RANK_MAX = 10.0
CANDIDATE_WT_RANK_MIN = 10.0

def classify_binder(rank):
    if pd.isna(rank):
        return "NA"
    if rank <= STRONG_BINDER_RANK:
        return "Strong"
    elif rank <= WEAK_BINDER_RANK:
        return "Weak"
    else:
        return "Non-binder"


def build_fasta(peptides, id_prefix="pep"):
    lines = []
    for i, p in enumerate(peptides):
        lines.append(f">{id_prefix}{i}")
        lines.append(p)
    return "\n".join(lines)


def call_iedb(peptides, allele):
    fasta = build_fasta(peptides)
    resp = requests.post(IEDB_URL, data={
        "method": METHOD,
        "sequence_text": fasta,
        "allele": allele,
        "length": "asis",
    })
    resp.raise_for_status()
    return resp.text


def main():
    # 1. Load 15-mers
    df15 = pd.read_csv(PEPTIDE_FILE, sep="\t")
    print(f"Loaded {len(df15)} 15-mer rows")

    key_cols = ["GeneName", "TranscriptID", "ProteinChange", "ProteinPosition", "MutPos"]
    pair_counts = df15.groupby(key_cols)["Type"].apply(lambda x: set(x))
    bad_pairs = pair_counts[pair_counts.apply(lambda s: s != {"Mutant", "WildType"})]
    print(f"Mutation-position groups: {len(pair_counts)}, incomplete pairs: {len(bad_pairs)}")

    unique_peptides = df15["Peptide"].unique().tolist()
    print(f"Unique 15-mer peptides: {len(unique_peptides)}")

    # Sanity check- 3 peptides, 1 allele, look at the raw response ---
    print("\nSanity check: calling IEDB API with 3 peptides...")
    raw = call_iedb(unique_peptides[:3], HLA_CLASS_II[0])
    print(raw[:2000])
    print("\n^ Confirm this looks like a real results table (tab-delimited), "
          "not an error message, before continuing.\n")
    input("Press Enter to continue with the full run, or Ctrl+C to stop and fix something...")

    # 2. Full run: batch peptides, loop per allele
    n_batches = (len(unique_peptides) + BATCH_SIZE - 1) // BATCH_SIZE
    print(f"\nRunning {len(HLA_CLASS_II)} alleles x {n_batches} batches "
          f"({BATCH_SIZE} peptides/batch)...")

    all_results = []
    for allele in HLA_CLASS_II:
        for i in tqdm(range(0, len(unique_peptides), BATCH_SIZE),
                      desc=f"Allele {allele}"):
            batch = unique_peptides[i:i + BATCH_SIZE]
            raw = call_iedb(batch, allele)
            try:
                result = pd.read_csv(io.StringIO(raw), sep="\t")
                result["Allele"] = allele
                all_results.append(result)
            except Exception as e:
                print(f"WARNING: failed to parse batch at index {i} for {allele}: {e}")
                print(raw[:500])
            time.sleep(REQUEST_DELAY)

    preds = pd.concat(all_results, ignore_index=True)
    print(f"\nRaw prediction columns returned: {preds.columns.tolist()}")
    preds.to_csv("raw_iedb_response_debug.tsv", sep="\t", index=False)
    print("Saved raw response to raw_iedb_response_debug.tsv for inspection.")

    # 3. NOTE: adjust these column names once you've seen the sanity-check
    #    output above. Update the rename map below to match.
    rename_map = {
        "peptide": "Peptide",
        "ic50": "Affinity_nM",
        "percentile_rank": "Percentile_Rank",
        "rank": "Percentile_Rank",
    }
    preds = preds.rename(columns={k: v for k, v in rename_map.items() if k in preds.columns})

    keep_cols = [c for c in ["Peptide", "Allele", "Affinity_nM", "Percentile_Rank"] if c in preds.columns]
    preds = preds[keep_cols].drop_duplicates()

    # 4. Join back onto full peptide table
    grid = df15.merge(preds, on="Peptide", how="left")

    if "Percentile_Rank" in grid.columns:
        grid["Binder"] = grid["Percentile_Rank"].apply(classify_binder)
    grid["Software"] = "NetMHCIIpan (via IEDB API)"
    grid["Prediction_Mode"] = METHOD

    print(f"\nBinder classification counts:\n{grid.get('Binder', pd.Series(dtype=str)).value_counts()}")

    grid.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"Wrote {OUTPUT_FILE} ({len(grid)} rows)")

    # 5. Mutant vs wild-type comparison
    if "Percentile_Rank" in grid.columns:
        compare_cols = key_cols + ["Allele"]
        mut = grid[grid["Type"] == "Mutant"][compare_cols + ["Peptide", "Percentile_Rank", "Affinity_nM"]]
        wt = grid[grid["Type"] == "WildType"][compare_cols + ["Peptide", "Percentile_Rank", "Affinity_nM"]]
        merged = mut.merge(wt, on=compare_cols, suffixes=("_Mutant", "_WildType"))
        merged["Rank_Improvement"] = merged["Percentile_Rank_WildType"] - merged["Percentile_Rank_Mutant"]

        candidates = merged[
            (merged["Percentile_Rank_Mutant"] <= CANDIDATE_MUT_RANK_MAX)
            & (merged["Percentile_Rank_WildType"] > CANDIDATE_WT_RANK_MIN)
        ].sort_values("Rank_Improvement", ascending=False)

        print(f"\nCandidate class II neoantigens: {len(candidates)}")
        candidates.to_csv(NEOANTIGEN_CANDIDATES_FILE, sep="\t", index=False)
        print(f"Wrote {NEOANTIGEN_CANDIDATES_FILE}")


if __name__ == "__main__":
    main()