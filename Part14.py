#Part 14

import pandas as pd

class_1 = pd.read_csv('/Users/naigelu/Desktop/Imperial College London/Year 2/STJU Summer Internship/Project 130/Lung Cancer/Google Drive Downloads/Part 11 + 12.1/classI_neoantigen_candidates.tsv', sep="\t")
class_2 = pd.read_csv('/Users/naigelu/Desktop/Imperial College London/Year 2/STJU Summer Internship/Project 130/Lung Cancer/Google Drive Downloads/Part 11 + 12.1/classii_neoantigen_candidates.csv', sep=",")
mutation_by_sample_matrix = pd.read_csv('/Users/naigelu/Desktop/Imperial College London/Year 2/STJU Summer Internship/Project 130/Lung Cancer/Google Drive Downloads/integrated_mutation_expression_matrix.tsv', sep="\t")
immunogenicity_class_1 = pd.read_csv('/Users/naigelu/Desktop/Imperial College London/Year 2/STJU Summer Internship/Project 130/Lung Cancer/Part 13/part13_immunogenicity_results.tsv', sep="\t")


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

#for class 2, since the classii_neoantigen_candidates.csv file table did not contain the necessary columns to calculate delta affinity of affinity fold change, delta Rank EL was used instead
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
    # Expression of the mutated gene: GeneLevelTPM > 0
    # Occurrence in one or more tumour samples: Mutation count more than 1 tumour
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

first_merge = class_1.merge(
    immunogenicity_class_1,
    on=[
        "GeneName",
        "ProteinChange",
        "Peptide_Mutant",
        "Allele"
    ],
    how="left"
)

mutation_by_sample_matrix = mutation_by_sample_matrix.rename(
    columns={"Gene_Name": "GeneName"}
)

tpm_lookup = mutation_by_sample_matrix[
    ["GeneName", "GeneLevelTPM", "TumourCount"]
].copy()

HLA_1 = first_merge.merge(
    tpm_lookup,
    on="GeneName",
    how="left",
)

print(HLA_1.head())
print("HLA_1 columns:", HLA_1.columns.tolist())

HLA_1["PrioritisationScore"] = (
   # 1 point: strong mutant binding rank
   (pd.to_numeric(
       HLA_1["Percentile_Rank_Mutant"],
   ) <= 0.5).astype(int)
   +

   # 1 point: mutant affinity is at least 2-fold better than WT
   (pd.to_numeric(
       HLA_1["AffinityFoldChange"],
   ) >= 2).astype(int)
   +

   # 1 point: favourable PRIME immunogenicity rank
   (pd.to_numeric(
       HLA_1["ImmunogenicityScore"],
   ) <= 0.5).astype(int)
   +

   # 1 point: mutated gene is expressed
   (pd.to_numeric(
       HLA_1["GeneLevelTPM"],
   ) >= 15).astype(int)
   +

   # 1 point: wild-type peptide has low predicted presentation
   (pd.to_numeric(
       HLA_1["Percentile_Rank_WildType"],
   ) > 2)
   +
   (pd.to_numeric(
       HLA_1["TumourCount"],
   ) >= 1)
)

print(HLA_1.head())
print("HLA_1 columns:", HLA_1.columns.tolist())
HLA_1.to_csv(
    "/Users/naigelu/Desktop/Imperial College London/Year 2/STJU Summer Internship/Project 130/Lung Cancer/All_Combined_tables_HLA_1.tsv",
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
    "/Users/naigelu/Desktop/Imperial College London/Year 2/STJU Summer Internship/Project 130/Lung Cancer/Part 14/candidates_prioritisation_score_6.tsv",
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



































