import pandas as pd

import config

# VEP web output: tab-separated, one row per variant
vep = pd.read_csv(config.VEP_OUTPUT, sep="\t")

# Enforce a single, consistent transcript-selection rule: MANE Select only.
# Variants without a MANE Select transcript are dropped so the strategy stays uniform.
before = len(vep)
vep = vep[vep["MANE"] == "MANE_Select"]

# Keep missense only
vep = vep[vep["Consequence"].str.contains("missense_variant")]
print(f"Kept {len(vep)} of {before} variants (MANE Select + missense).")

annot = pd.DataFrame({
    "Reference_Assembly": "GRCh38",
    "Gene_Symbol": vep["SYMBOL"],
    "Transcript_ID": vep["Feature"],
    "Protein_ID": vep["ENSP"],
    "Genomic_Variant": vep["Location"] + " " + vep["REF_ALLELE"] + ">" + vep["Allele"],
    "Coding_DNA_Change": vep["HGVSc"],
    "Protein_Change": vep["HGVSp"],
    "Variant_Consequence": vep["Consequence"],
    "Protein_Position": vep["Protein_position"],
    "Transcript_Selection_Rule": "MANE Select",
})

annot.to_csv(config.VARIANT_ANNOTATION, sep="\t", index=False)
print(f"Saved {len(annot)} annotated missense variants to {config.VARIANT_ANNOTATION}")
print(annot.head())
