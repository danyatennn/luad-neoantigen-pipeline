"""
Somatic mutation dataset and the binary mutation-by-sample matrix
Reads the TCGA-LUAD MAF and writes data/01_mutation_by_sample.tsv
"""

import pandas as pd
import config

MAF_PATH = config.MAF
OUT_PATH = config.MUTATION_MATRIX

cols = [
    "Hugo_Symbol", # gene name
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Variant_Classification",
    "Variant_Type",
    "HGVSc", # coding DNA mutation, f.e. c.35G>A
    "HGVSp", # protein change, e.g. p.Gly12Asp
    "Tumor_Sample_Barcode", # sample
    "GDC_FILTER", # for filtering PASS
]

maf = pd.read_csv(MAF_PATH, sep="\t", usecols=cols, low_memory=False)
print("number of rows in MAF:", len(maf))

# leave only missense SNV; PASS = empty GDC_FILTER (quality reads)
mask = (
    (maf["Variant_Classification"] == "Missense_Mutation")
    & (maf["Variant_Type"] == "SNP")
    & (maf["GDC_FILTER"].isna() | (maf["GDC_FILTER"] == ""))
)
muts = maf.loc[mask].copy()
print("number of rows after filtering:", len(muts))
print("number of aliquot barcodes:", muts["Tumor_Sample_Barcode"].nunique())

# A TCGA barcode identifies an aliquot (one extraction), not a tumour: the same
# tumour can be sequenced on several plates and appear under several barcodes.
# The first 16 characters identify the sample, so truncate before counting.
muts["Tumor_Sample_Barcode"] = muts["Tumor_Sample_Barcode"].str[:16]
print("number of unique samples:", muts["Tumor_Sample_Barcode"].nunique())

# each unique mutation = (gene, HGVSc, HGVSp). Different mutations in one gene — different rows.
muts["MutID"] = muts["Hugo_Symbol"] + "|" + muts["HGVSc"] + "|" + muts["HGVSp"]

# binary matrix: rows — mutations, columns — samples
matrix = (
    pd.crosstab(muts["MutID"], muts["Tumor_Sample_Barcode"])
    .clip(upper=1)  # if a mutation occurs multiple times in one sample, we still want to count it as 1
)

# restore separate columns Gene_Name / Mutation / AminoAcid_Change from MutID
meta = matrix.index.to_series().str.split("|", expand=True)
meta.columns = ["Gene_Name", "Mutation", "AminoAcid_Change"]

result = pd.concat([meta.reset_index(drop=True), matrix.reset_index(drop=True)], axis=1)
print("Size of matrix:", result.shape, "(rows, columns with meta)")

# save to TSV
result.to_csv(OUT_PATH, sep="\t", index=False)
print("Saved to:", OUT_PATH)
