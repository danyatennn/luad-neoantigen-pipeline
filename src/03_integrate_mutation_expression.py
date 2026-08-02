"""
Part 7 - Integration of mutations and gene expression.

Converted from task7.ipynb (the original notebook is kept in bin/).
Adds GeneLevelTPM (median TPM across tumour samples) to every mutation row and
writes data/03_integrated_mutation_expression.tsv.
"""

import pandas as pd

import config

MAF_PATH = config.MAF
MUT_MATRIX = config.MUTATION_MATRIX
TPM_PATH = config.TPM_MATRIX
OUT_PATH = config.INTEGRATED_MATRIX

mut = pd.read_csv(MUT_MATRIX, sep="\t")
print("mutation matrix shape:", mut.shape)

tpm = pd.read_csv(TPM_PATH, sep="\t", index_col="GeneName")
print("TPM table shape:", tpm.shape, "(genes, samples)")

# aggregation: median TPM per gene across all tumour samples
gene_level_tpm = tpm.median(axis=1).rename("GeneLevelTPM")

# build Hugo -> Ensembl mapping directly from the MAF (columns Hugo_Symbol and Gene)
mapping = pd.read_csv(MAF_PATH, sep="\t", usecols=["Hugo_Symbol", "Gene"], low_memory=False)
# one Ensembl ID per Hugo symbol (take the most frequent one if MAF lists several)
mapping = (
    mapping.dropna()
    .groupby("Hugo_Symbol")["Gene"]
    .agg(lambda s: s.mode().iat[0])
)
print("unique Hugo symbols in mapping:", len(mapping))

# attach Ensembl ID -> GeneLevelTPM to every mutation row
mut["Ensembl"] = mut["Gene_Name"].map(mapping)
mut["GeneLevelTPM"] = mut["Ensembl"].map(gene_level_tpm)

# put GeneLevelTPM right after AminoAcid_Change, drop helper Ensembl column
sample_cols = [c for c in mut.columns if c.startswith("TCGA-")]
result = mut[["Gene_Name", "Mutation", "AminoAcid_Change", "GeneLevelTPM", *sample_cols]]

print("integrated matrix shape:", result.shape)
print("rows with missing TPM:", result["GeneLevelTPM"].isna().sum())

# save integrated matrix as tab-delimited file
result.to_csv(OUT_PATH, sep="\t", index=False)
print("saved to:", OUT_PATH)

# sanity check: GeneLevelTPM for well-known LUAD drivers and the most frequent mutations
print(
    result.assign(freq=result[sample_cols].sum(axis=1))
          .nlargest(10, "freq")[["Gene_Name", "Mutation", "AminoAcid_Change", "GeneLevelTPM", "freq"]]
)
