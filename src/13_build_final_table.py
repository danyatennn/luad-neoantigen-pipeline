#Part 15

import pandas as pd

import config

HLA1_immunogenicity = pd.read_csv(config.IMMUNOGENICITY, sep='\t')
ninemers_presentation = pd.read_csv(config.CLASS1_CANDIDATES, sep='\t')
fifteenmers_presentation = pd.read_csv(config.CLASS2_CANDIDATES, sep=',')
genelevelTPM = pd.read_csv(config.INTEGRATED_MATRIX, sep='\t')

#in order to find the mutation frequency, as well as
sample_columns = [
    column
    for column in genelevelTPM.columns
    if column.startswith("TCGA-")
]

# Ensure the sample values are numeric
genelevelTPM[sample_columns] = genelevelTPM[sample_columns].apply(
    pd.to_numeric,
    errors="coerce"
)

# Count how many tumour samples contain the mutation
genelevelTPM["TumourMutationCount"] = (
    genelevelTPM[sample_columns]
    .eq(1)
    .sum(axis=1)
)

# Total number of tumour samples in the matrix
total_tumour_samples = len(sample_columns)

# Fraction of tumour samples containing the mutation
genelevelTPM["MutationFrequency"] = (
    genelevelTPM["TumourMutationCount"]
    / total_tumour_samples
)

class1_9mers = HLA1_immunogenicity.merge(
    ninemers_presentation,
    on=["GeneName", "ProteinChange", "Peptide_Mutant", "Allele"],
    how="left"
)

print("HLA1_immunogenicity columns:")
print(HLA1_immunogenicity.columns.tolist())

print("\n9-mers_presentation columns:")
print(ninemers_presentation.columns.tolist())

print("\n15-mers_presentation columns:")
print(fifteenmers_presentation.columns.tolist())

print("\ngenelevelTPM columns:")
print(genelevelTPM.columns.tolist())

print("class1_9mers columns:")
print(class1_9mers.columns.tolist())

#GeneName, Chromosome, Position, Ref, Alt, TranscriptID, ProteinChange, GeneLevelTPM, MutationFrequency, PeptideType, Peptide, PeptideLength, MutationPosition, HLAAllele, BindingAffinity, BindingRank, BindingScore, ImmunogenicityScore
#, Chromosome, , , , , , , , PeptideType, Peptide, PeptideLength, , , BindingAffinity, BindingRank, BindingScore,


final_columns = [
    "GeneName",
    "Chromosome",
    "Position",
    "TranscriptID",
    "ProteinChange",
    "GeneLevelTPM",
    "MutationFrequency",
    "PeptideType",
    "Peptide",
    "PeptideLength",
    "MutationPosition",
    "HLAAllele",
    "BindingAffinity",
    "BindingRank",
    "ImmunogenicityScore"
]

# Optional blank table
final_table = pd.DataFrame(columns=final_columns)

#what we want from HLA1_immunogenicity: GeneName, ProteinChange, ImmunogenicityScore -- this should first be combined with the ninemers_presenation
#what we want from ninemers_presentation: GeneName, ProteinPosition, Peptide_WildType (Ref), Peptide_Mutant (Alt), TranscriptID, MutPos, Allele, Affinity_nM_Mutant, Percentile_Rank_Mutant
#what we want from fifteenmers_presentation: GeneName, ProteinChange, Peptide, MutPos, Allele, Mut_Binding_Call, Mut_Rank_EL
#what we want from genelevelTPM: GeneLevelTPM, MutationFrequency

#for the mutant 9mers:
final_9mers = pd.DataFrame({
    "GeneName": class1_9mers["GeneName"],
    "Chromosome": "hg38/GRCh38",
    "Position": class1_9mers["ProteinPosition"],
    "TranscriptID": class1_9mers["TranscriptID"],
    "ProteinChange": class1_9mers["ProteinChange_x"],

    "PeptideType": "Mutant",
    "Peptide": class1_9mers["Peptide_Mutant_x"],
    "PeptideLength": 9,
    "MutationPosition": class1_9mers["MutPos"],
    "HLAAllele": class1_9mers["Allele_x"],

    "BindingAffinity": class1_9mers["Affinity_nM_Mutant"],
    "BindingRank": class1_9mers["Percentile_Rank_Mutant"],
    "ImmunogenicityScore": class1_9mers["ImmunogenicityScore"]
})

#for the 15mers:
final_15mers = pd.DataFrame({
    "GeneName": fifteenmers_presentation["GeneName"],
    "Chromosome": "hg38/GRCh38",
    "Position": fifteenmers_presentation["ProteinPosition"],
    "TranscriptID": pd.NA,
    "ProteinChange": fifteenmers_presentation["ProteinChange"],

    "PeptideType": "Mutant",
    "Peptide": fifteenmers_presentation["Peptide"],
    "PeptideLength": 15,
    "MutationPosition": fifteenmers_presentation["MutPos"],
    "HLAAllele": fifteenmers_presentation["Allele"],

    "BindingAffinity": pd.NA,
    "BindingRank": fifteenmers_presentation["Mut_Rank_EL"],

    # PRIME was not performed for class II peptides
    "ImmunogenicityScore": pd.NA
})

final_table = pd.concat(
    [final_9mers, final_15mers],
    ignore_index=True
)

##GeneLevelTPMs and MutationFrequency need to be added last since they have both the 9 and 15 mers mixed together based on the gene name.
##please check this when you come back, it still needs some work
#TODO: whoever wrote this is this still true or

genelevelTPM = genelevelTPM.rename(
    columns={"Gene_Name": "GeneName", "AminoAcid_Change": "ProteinChange"}
)

# GeneLevelTPM is one value per gene -> dedupe on GeneName 
tpm_lookup = (
    genelevelTPM[["GeneName", "GeneLevelTPM"]]
    .drop_duplicates(subset=["GeneName"])
)

freq_lookup = (
    genelevelTPM[["GeneName", "ProteinChange", "MutationFrequency"]]
    .drop_duplicates(subset=["GeneName", "ProteinChange"])
)

final_table = final_table.merge(tpm_lookup, on="GeneName", how="left")
final_table = final_table.merge(freq_lookup, on=["GeneName", "ProteinChange"], how="left")
final_table = final_table[final_columns]
print("Final table shape:", final_table.shape)
print(final_table.head(10).to_string(index=False))

print("final_table columns:")
print(final_table.columns.tolist())

final_table.to_csv(
    config.NEOANTIGEN_TABLE,
    sep="\t",
    index=False,
    na_rep="NA"
)
