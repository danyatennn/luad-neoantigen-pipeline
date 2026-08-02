import pandas as pd

import config

# Read only the columns that define a variant and its class
cols = ["Chromosome", "Start_Position", "Reference_Allele",
        "Tumor_Seq_Allele2", "Variant_Classification"]
maf = pd.read_csv(config.MAF, sep="\t", comment="#",
                  usecols=cols, low_memory=False)

# Keep missense only
maf = maf[maf["Variant_Classification"] == "Missense_Mutation"]

# One variant may appear in many samples, VEP only needs the unique variant list
maf = maf.drop_duplicates(subset=["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"])

# Build a VCF. MAF coordinates are 1-based, same as VCF.
# Missense variants are all SNP/ONP/TNP (REF and ALT same length), so no anchor base is needed.
vcf = pd.DataFrame({
    "#CHROM": maf["Chromosome"].str.replace("chr", "", regex=False),
    "POS": maf["Start_Position"],
    "ID": ".",
    "REF": maf["Reference_Allele"],
    "ALT": maf["Tumor_Seq_Allele2"],
    "QUAL": ".",
    "FILTER": ".",
    "INFO": ".",
})

with open(config.VEP_INPUT, "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    vcf.to_csv(f, sep="\t", index=False)

print(f"Wrote {len(vcf)} unique missense variants to {config.VEP_INPUT}")
