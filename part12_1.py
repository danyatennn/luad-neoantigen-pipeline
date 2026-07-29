"""
Part 12.1 - Peptide-MHC Binding Prediction: Class I, 9-mers (MHCflurry)
Loads the 9-mer peptides from Part 10, checks mutant/wild-type pairing
integrity, scores every peptide against the fixed Class I HLA panel
from Part 11 using MHCflurry, and computes the mutant-vs-wildtype
binding delta to flag candidate neoantigens.

Install requirements:
    pip install mhcflurry --break-system-packages
    mhcflurry-downloads fetch models_class1_presentation

"""

import pandas as pd
from mhcflurry import Class1AffinityPredictor
import mhcflurry
from tqdm import tqdm

PEPTIDE_FILE = r"C:\Users\lyssa\Downloads\mutant_peptides.tsv" # Part 10 output filepath
OUTPUT_FILE = "9mer_predictions.tsv"
NEOANTIGEN_CANDIDATES_FILE = "classI_neoantigen_candidates.tsv"

HLA_CLASS_I = [
    "HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01",
]

# Percentile-rank thresholds for binder classification
STRONG_BINDER_RANK = 0.5
WEAK_BINDER_RANK = 2.0

# A mutation/allele pair is flagged as a candidate neoantigen if the
# mutant percentile rank is <= this AND the wild-type rank is not
# already a strong binder itself (i.e. binding is newly gained/improved
# because of the mutation, not just generally strong for that region)
CANDIDATE_MUT_RANK_MAX = 2.0
CANDIDATE_WT_RANK_MIN = 2.0

def classify_binder(rank):
    if rank <= STRONG_BINDER_RANK:
        return "Strong"
    elif rank <= WEAK_BINDER_RANK:
        return "Weak"
    else:
        return "Non-binder"


def main():
    # 1. Load Part 10 output, keep only 9-mers
    df = pd.read_csv(PEPTIDE_FILE, sep="\t")
    df9 = df[df["Length"] == 9].copy()
    print(f"Loaded {len(df)} total peptide rows, {len(df9)} are 9-mers")

    # Remove peptides with 'non-standard' amino acids (U)
    STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

    def has_nonstandard_aa(peptide):
        return not set(peptide).issubset(STANDARD_AA)

    df9["Nonstandard_AA"] = df9["Peptide"].apply(has_nonstandard_aa)

    n_excluded = df9["Nonstandard_AA"].sum()
    excluded_genes = df9.loc[df9["Nonstandard_AA"], "GeneName"].unique()
    print(f"Excluding {n_excluded} peptides with nonstandard amino acids "
        f"(selenocysteine 'U'), affecting genes: {list(excluded_genes)}")

    df9 = df9[~df9["Nonstandard_AA"]].copy()

    # 2. Pairing check: every mutation/position should have exactly one Mutant and one WildType peptide
    key_cols = ["GeneName", "TranscriptID", "ProteinChange", "ProteinPosition", "MutPos"]
    pair_counts = df9.groupby(key_cols)["Type"].apply(lambda x: set(x))
    bad_pairs = pair_counts[pair_counts.apply(lambda s: s != {"Mutant", "WildType"})]
    print(f"\nMutation-position groups: {len(pair_counts)}")
    print(f"Incomplete pairs: {len(bad_pairs)}")
    if len(bad_pairs) > 0:
        print("First few incomplete pairs:")
        print(bad_pairs.head())

    # 3. Get the list of unique peptides to predict (avoid redundant calls)
    unique_peptides = df9["Peptide"].unique().tolist()
    print(f"\nUnique 9-mer peptides: {len(unique_peptides)}")
    print(f"Total predictions to run: {len(unique_peptides)} peptides x "
          f"{len(HLA_CLASS_I)} alleles = {len(unique_peptides) * len(HLA_CLASS_I)}")

    # 4. Load MHCflurry's affinity predictor. 
    affinity_predictor = Class1AffinityPredictor.load()

    per_allele_frames = []
    for allele in tqdm(HLA_CLASS_I, desc="Alleles"):
        pred_df = affinity_predictor.predict_to_dataframe(
            peptides=unique_peptides,
            allele=allele,
        )
        pred_df["Allele"] = allele
        per_allele_frames.append(pred_df)

    preds = pd.concat(per_allele_frames, ignore_index=True)
    preds = preds.rename(columns={
        "peptide": "Peptide",
        "prediction": "Affinity_nM",
        "prediction_percentile": "Percentile_Rank",
    })[["Peptide", "Allele", "Affinity_nM", "Percentile_Rank"]]

    # 5. Join predictions back onto the full peptide table (this
    #    reconstructs the peptide x allele grid, now with real scores)
    grid = df9.merge(preds, on="Peptide", how="left")

    grid["Binder"] = grid["Percentile_Rank"].apply(classify_binder)
    grid["Software"] = "MHCflurry"
    grid["Software_Version"] = mhcflurry.__version__
    grid["Prediction_Mode"] = "affinity"

    print(f"\nBinder classification counts:\n{grid['Binder'].value_counts()}")

    # 6. Save full results table
    grid.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"\nWrote {OUTPUT_FILE} ({len(grid)} rows)")

    # 7. Mutant vs wild-type comparison -> candidate neoantigens
    compare_cols = key_cols + ["Allele"]
    mut = grid[grid["Type"] == "Mutant"][compare_cols + ["Peptide", "Percentile_Rank", "Affinity_nM"]]
    wt = grid[grid["Type"] == "WildType"][compare_cols + ["Peptide", "Percentile_Rank", "Affinity_nM"]]

    merged = mut.merge(
        wt, on=compare_cols, suffixes=("_Mutant", "_WildType")
    )
    merged["Rank_Improvement"] = merged["Percentile_Rank_WildType"] - merged["Percentile_Rank_Mutant"]

    candidates = merged[
        (merged["Percentile_Rank_Mutant"] <= CANDIDATE_MUT_RANK_MAX)
        & (merged["Percentile_Rank_WildType"] > CANDIDATE_WT_RANK_MIN)
    ].sort_values("Rank_Improvement", ascending=False)

    print(f"\nCandidate neoantigens (mutant binds, wild-type does not): {len(candidates)}")
    candidates.to_csv(NEOANTIGEN_CANDIDATES_FILE, sep="\t", index=False)
    print(f"Wrote {NEOANTIGEN_CANDIDATES_FILE}")

    if len(candidates) > 0:
        print("\nTop 10 candidates by rank improvement:")
        print(candidates.head(10)[
            ["GeneName", "ProteinChange", "Allele", "Peptide_Mutant",
             "Percentile_Rank_Mutant", "Percentile_Rank_WildType"]
        ])


if __name__ == "__main__":
<<<<<<< Updated upstream
    main()
=======
    main()  
>>>>>>> Stashed changes
