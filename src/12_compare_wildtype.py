"""
Part 14 - Comparison with Wild-Type Peptide: 

"""
# Part 14

import pandas as pd

import config

class_1 = pd.read_csv(config.CLASS1_CANDIDATES, sep="\t")
class_2 = pd.read_csv(config.CLASS2_CANDIDATES, sep="\t")
mutation_by_sample_matrix = pd.read_csv(config.INTEGRATED_MATRIX, sep="\t")
immunogenicity_class_1 = pd.read_csv(config.IMMUNOGENICITY, sep="\t")


#equation 1: DeltaAffinity =WildTypeAffinity- MutantAffinity
    #positive DeltaAffinity = mutant has lower IC50 and binds more strongly;
#equation 2: AffinityFoldChange = WildTypeAffinity / MutantAffinity
    #AffinityFoldChange > 1 = mutant binds more strongly than wild type.


class_1["DeltaAffinity"] = (
    class_1["Affinity_nM_WildType"]
    - class_1["Affinity_nM_Mutant"]
)

class_1["AffinityFoldChange"] = (
    class_1["Affinity_nM_WildType"]
    / class_1["Affinity_nM_Mutant"]
)

# for class 2, since the classii_neoantigen_candidates.csv file table did not contain the 
# necessary columns to calculate delta affinity of affinity fold change,
# delta Rank EL was used instead
#deltaRank EL equation: WildType_Rank_EL - Mut_Rank_EL
    #A positive value means that the mutant has a lower, better presentation rank
#Rank fold change can also be calculated
#RankFoldChange = WildType_Rank_EL / Mut_Rank_EL

class_2["DeltaRank_EL"] = (
    class_2["WT_Rank_EL"]
    - class_2["Mut_Rank_EL"]
)

class_2["RankFoldChange_EL"] = (
    class_2["WT_Rank_EL"]
    / class_2["Mut_Rank_EL"]
)

print(class_1.head())
print("Class 1 columns:", class_1.columns.tolist())

print(class_2.head())
print("Class 2 columns:", class_2.columns.tolist())

#in order to deduce whether a neoantigen should be prioretised or not, a transparent feature-count score can be used.
#This will be based on 6 criteria:
#Priority neoantigen candidates will have:
    # Strong mutant peptide–HLA binding: High Mutant binding affinity
    # Better mutant binding than wild-type binding: We want a positive DeltaAffinity
    # Favourable immunogenicity score: we want %Rank to be ≤ 0.5
    # Expression of the mutated gene: GeneLevelTPM >= 15 (moderately expressed; the
    #   threshold is a project choice and is stated in the README)
    # Occurrence in one or more tumour samples: TumourCount >= 1
    # Low or absent presentation predicted for the wild-type peptide: WT binding %Rank > 2, meaning that WT is not a predicted binder.

#For Class 1 HLA

sample_columns = [
    column
    for column in mutation_by_sample_matrix.columns
    if column.startswith("TCGA-")
]

mutation_by_sample_matrix["TumourCount"] = (
    mutation_by_sample_matrix[sample_columns]
    .eq(1)
    .sum(axis=1)
)

# The immunogenicity table is deduplicated first: the same peptide/allele appears twice
# when one protein change is reachable from two nucleotide changes, and PRIME scores it
# identically both times. Without this the join adds 264 spurious rows.
IMMUNO_KEY = ["GeneName", "ProteinChange", "Peptide_Mutant", "Allele"]

first_merge = class_1.merge(
    immunogenicity_class_1.drop_duplicates(subset=IMMUNO_KEY),
    on=IMMUNO_KEY,
    how="left"
)

mutation_by_sample_matrix = mutation_by_sample_matrix.rename(
    columns={"Gene_Name": "GeneName", "AminoAcid_Change": "ProteinChange"}
)

# GeneLevelTPM is one value per gene --> dedupe on GeneName before joining
tpm_lookup = (
    mutation_by_sample_matrix[["GeneName", "GeneLevelTPM"]]
    .drop_duplicates(subset=["GeneName"])
)

tumour_count_lookup = (
    mutation_by_sample_matrix[["GeneName", "ProteinChange", "TumourCount"]]
    .drop_duplicates(subset=["GeneName", "ProteinChange"])
)

HLA_1 = first_merge.merge(tpm_lookup, on="GeneName", how="left")
HLA_1 = HLA_1.merge(tumour_count_lookup, on=["GeneName", "ProteinChange"], how="left")

print(HLA_1.head())
print("HLA_1 columns:", HLA_1.columns.tolist())

HLA_1["PrioritisationScore"] = (
   # 1 point: strong mutant binding rank
   (pd.to_numeric(
       HLA_1["Percentile_Rank_Mutant"], errors="coerce"
   ) <= 0.5).astype(int)
   +

   # 1 point: mutant affinity is at least 2-fold better than WT
   (pd.to_numeric(
       HLA_1["AffinityFoldChange"], errors="coerce"
   ) >= 2).astype(int)
   +

   # 1 point: favourable PRIME immunogenicity rank
   (pd.to_numeric(
       HLA_1["ImmunogenicityScore"], errors="coerce"
   ) <= 0.5).astype(int)
   +

   # 1 point: mutated gene is expressed
   (pd.to_numeric(
       HLA_1["GeneLevelTPM"], errors="coerce"
   ) >= 15).astype(int)
   +

   # 1 point: wild-type peptide has low predicted presentation
   (pd.to_numeric(
       HLA_1["Percentile_Rank_WildType"], errors="coerce"
   ) > 2)
   +
   (pd.to_numeric(
       HLA_1["TumourCount"], errors="coerce"
   ) >= 1)
)

print(HLA_1.head())
print("HLA_1 columns:", HLA_1.columns.tolist())
HLA_1.to_csv(
    config.WT_VS_MUTANT,
    sep="\t",
    index=False,
    na_rep="NA"
)


#for finding those with a specific score
# Keep only candidates with a prioritisation score of exactly 6
candidates_score_6 = HLA_1[
    HLA_1["PrioritisationScore"] == 6
].copy()

# Save as a tab-delimited file
candidates_score_6 = candidates_score_6.drop_duplicates().copy()

candidates_score_6.to_csv(
    config.CANDIDATES_SCORE_6,
    sep="\t",
    index=False,
    na_rep="NA"
)

#print(
#    f"Saved {len(candidates_score_6)} candidates to "
#    "'candidates_prioritisation_score_6.tsv'"
#)


#print(
#    selected_candidates
#    .sort_values("GeneName")
#    .to_string(index=False)
#)

#for finding those with a score greater than what is set.
minimum_score = 4

high_priority_candidates = HLA_1[
    HLA_1["PrioritisationScore"] >= minimum_score
].copy()

#print(
#    high_priority_candidates
#    .sort_values("GeneName")
#    .to_string(index=False)
#)

score_summary = (
    HLA_1
    .groupby("PrioritisationScore")
    .agg(
        CandidateCount=("PrioritisationScore", "size"),
        UniqueGeneCount=("GeneName", "nunique")
    )
    .reindex(range(1, 7), fill_value=0)
    .reset_index()
)

print(score_summary.to_string(index=False))



































