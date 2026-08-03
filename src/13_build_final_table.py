"""
Part 15 - Advanced output format.

Builds the submitted neoantigen table: one row per peptide x HLA allele, carrying the
genomic variant, expression, mutation frequency, binding and immunogenicity scores.
"""

import pandas as pd

import config

HLA1_immunogenicity = pd.read_csv(config.IMMUNOGENICITY, sep='\t')
ninemers_presentation = pd.read_csv(config.CLASS1_CANDIDATES, sep='\t')
fifteenmers_presentation = pd.read_csv(config.CLASS2_CANDIDATES, sep='\t')
genelevelTPM = pd.read_csv(config.INTEGRATED_MATRIX, sep='\t')

# --- mutation frequency, keyed on the genomic variant ---
# It cannot be keyed on the protein change: the mutation matrix numbers residues on the
# transcript the MAF was annotated with, while the peptides use VEP's MANE Select
# transcript, so the same variant can be p.Asp567Tyr in one and p.Asp608Tyr in the other.
# The genomic coordinate is transcript-independent, so it is used instead. The filter
# below is the one applied in step 01, so the two stay consistent.
maf = pd.read_csv(
    config.MAF, sep="\t", low_memory=False,
    usecols=["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2",
             "Variant_Classification", "Variant_Type", "GDC_FILTER", "Tumor_Sample_Barcode"],
)
maf = maf[(maf["Variant_Classification"] == "Missense_Mutation")
          & (maf["Variant_Type"] == "SNP")
          & (maf["GDC_FILTER"].isna() | (maf["GDC_FILTER"] == ""))].copy()
maf["Chromosome"] = maf["Chromosome"].str.replace("chr", "", regex=False)
maf["Sample"] = maf["Tumor_Sample_Barcode"].str[:16]   # aliquot barcode -> sample
total_tumour_samples = maf["Sample"].nunique()

freq_lookup = (
    maf.groupby(["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"])["Sample"]
    .nunique()
    .rename("MutationFrequency")
    .reset_index()
    .rename(columns={"Start_Position": "Position",
                     "Reference_Allele": "Ref", "Tumor_Seq_Allele2": "Alt"})
)
freq_lookup["MutationFrequency"] /= total_tumour_samples
print(f"tumour samples: {total_tumour_samples}, distinct variants: {len(freq_lookup):,}")

# --- class I: attach immunogenicity to the candidate table ---
# The candidate table is the spine, so the result keeps one row per candidate x allele.
# The immunogenicity table is deduplicated first: the same peptide/allele appears twice
# when one protein change is reachable from two different nucleotide changes, and PRIME
# scores it identically both times.
KEY = ["GeneName", "ProteinChange", "Peptide_Mutant", "Allele"]
class1_9mers = ninemers_presentation.merge(
    HLA1_immunogenicity.drop_duplicates(subset=KEY), on=KEY, how="left"
)
print(f"class I candidates: {len(ninemers_presentation)} -> after merge: {len(class1_9mers)}")

# --- genomic variant and transcript, from the VEP annotation of step 06 ---
# The candidate tables carry protein-level coordinates only, so the genomic ones are
# looked up here rather than left blank.
annot = pd.read_csv(
    config.VARIANT_ANNOTATION, sep="\t",
    usecols=["Gene_Symbol", "Protein_Change", "Genomic_Variant", "Transcript_ID"],
).rename(columns={"Gene_Symbol": "GeneName", "Transcript_ID": "TranscriptID"})

# "ENSP00000493376.2:p.Ile239Met" -> "p.Ile239Met"
annot["ProteinChange"] = annot["Protein_Change"].str.split(":").str[-1]
# "1:69744-69744 C>G" -> chromosome, position, reference allele, alternate allele
annot[["Chromosome", "Position", "Ref", "Alt"]] = annot["Genomic_Variant"].str.extract(
    r"^([^:]+):(\d+)-\d+ (\S+)>(\S+)$"
)
# A few protein changes are reachable from two nucleotide changes; keep the first.
annot = annot.drop_duplicates(subset=["GeneName", "ProteinChange"])[
    ["GeneName", "ProteinChange", "Chromosome", "Position", "Ref", "Alt", "TranscriptID"]
]

# Tool name and version are already recorded by step 09; read them from there so the
# two cannot drift apart.
class1_tool = pd.read_csv(
    config.CLASS1_PREDICTIONS, sep="\t", nrows=1, usecols=["Software", "Software_Version"]
).iloc[0]


def class1_rows(peptide_type):
    """One row per class I candidate x allele, for the mutant or the wild-type peptide."""
    mutant = peptide_type == "Mutant"
    return pd.DataFrame({
        "GeneName": class1_9mers["GeneName"],
        "ProteinChange": class1_9mers["ProteinChange"],
        "PeptideType": peptide_type,
        "Peptide": class1_9mers["Peptide_Mutant" if mutant else "Peptide_WildType"],
        "PeptideLength": 9,
        "MutationPosition": class1_9mers["MutPos"],
        "HLAAllele": class1_9mers["Allele"],
        "BindingAffinity": class1_9mers["Affinity_nM_Mutant" if mutant else "Affinity_nM_WildType"],
        "BindingRank": class1_9mers["Percentile_Rank_Mutant" if mutant else "Percentile_Rank_WildType"],
        # MHCflurry ran in affinity mode: it reports an IC50 and a percentile rank,
        # there is no separate presentation score.
        "BindingScore": pd.NA,
        # PRIME was run on the mutant peptides only.
        "ImmunogenicityScore": class1_9mers["ImmunogenicityScore"] if mutant else pd.NA,
        "PredictionTool": class1_tool["Software"],
        "ToolVersion": class1_tool["Software_Version"],
    })


def class2_rows():
    """One row per class II candidate x allele. Mutant only - see the note below."""
    return pd.DataFrame({
        "GeneName": fifteenmers_presentation["GeneName"],
        "ProteinChange": fifteenmers_presentation["ProteinChange"],
        "PeptideType": "Mutant",
        "Peptide": fifteenmers_presentation["Peptide"],
        "PeptideLength": 15,
        "MutationPosition": fifteenmers_presentation["MutPos"],
        "HLAAllele": fifteenmers_presentation["Allele"],
        # NetMHCIIpan ran in eluted-ligand mode: it reports %Rank_EL, not an IC50.
        "BindingAffinity": pd.NA,
        "BindingRank": fifteenmers_presentation["Mut_Rank_EL"],
        "BindingScore": pd.NA,
        # PRIME is a class I method and was not run for class II peptides.
        "ImmunogenicityScore": pd.NA,
        "PredictionTool": "NetMHCIIpan",
        "ToolVersion": "4.2",
    })


# Wild-type rows are emitted for class I, where the wild-type peptide and its scores are
# both available. The class II candidate table keeps only the wild-type rank and not the
# wild-type peptide sequence, so no class II wild-type rows are produced.
final_table = pd.concat(
    [class1_rows("Mutant"), class1_rows("WildType"), class2_rows()],
    ignore_index=True,
)

final_table = final_table.merge(annot, on=["GeneName", "ProteinChange"], how="left")

# --- expression and mutation frequency ---
genelevelTPM = genelevelTPM.rename(
    columns={"Gene_Name": "GeneName", "AminoAcid_Change": "ProteinChange"}
)

# GeneLevelTPM is one value per gene.
tpm_lookup = genelevelTPM[["GeneName", "GeneLevelTPM"]].drop_duplicates(subset=["GeneName"])

final_table["Position"] = pd.to_numeric(final_table["Position"], errors="coerce")
final_table = final_table.merge(tpm_lookup, on="GeneName", how="left")
final_table = final_table.merge(freq_lookup, on=["Chromosome", "Position", "Ref", "Alt"], how="left")

final_columns = [
    "GeneName", "Chromosome", "Position", "Ref", "Alt", "TranscriptID", "ProteinChange",
    "GeneLevelTPM", "MutationFrequency", "PeptideType", "Peptide", "PeptideLength",
    "MutationPosition", "HLAAllele", "BindingAffinity", "BindingRank", "BindingScore",
    "ImmunogenicityScore", "PredictionTool", "ToolVersion",
]
final_table = final_table[final_columns]

print("Final table shape:", final_table.shape)
print(final_table["PeptideType"].value_counts().to_string())
print(final_table.head(5).to_string(index=False))

final_table.to_csv(config.NEOANTIGEN_TABLE, sep="\t", index=False, na_rep="NA")
print(f"\nWrote {config.NEOANTIGEN_TABLE}")

# Stronger literature-based filtering & ranking

EXPR_TPM_MIN = 3034 # expression fallback threshold (RNA read-depth arm not applicable here)
IC50_STRONG_MAX = 50.0 # nM, Teku & Vihinen 2018
IC50_WEAK_MAX = 500.0 # nM, Teku & Vihinen 2018
PRIME_RANK_MAX = 0.5 # PRIME %Rank, Schmidt et al. 2021
TOP_N  = 20

final_table["BindingAffinity"]  = pd.to_numeric(final_table["BindingAffinity"], errors="coerce")
final_table["BindingRank"] = pd.to_numeric(final_table["BindingRank"], errors="coerce")
final_table["ImmunogenicityScore"] = pd.to_numeric(final_table["ImmunogenicityScore"], errors="coerce")
final_table["GeneLevelTPM"] = pd.to_numeric(final_table["GeneLevelTPM"], errors="coerce")

is_class1 = final_table["PeptideLength"] == 9
is_class2 = final_table["PeptideLength"] == 15

expressed = final_table["GeneLevelTPM"] >= EXPR_TPM_MIN
strong_binder = final_table["BindingAffinity"] < IC50_WEAK_MAX # Class I only — NA for Class II fails safely
# PRIME wasn't run for Class II, so don't let a real NA there disqualify 15-mers
immunogenic_ok = (final_table["ImmunogenicityScore"] <= PRIME_RANK_MAX) | is_class2

hard_filtered = final_table[
    expressed
    & (strong_binder | is_class2) # IC50 filter is Class-I-only
    & immunogenic_ok
].copy()

print(f"{len(final_table)} rows -> {len(hard_filtered)} pass hard literature filters "
      f"(expression >= {EXPR_TPM_MIN} TPM, IC50 < {IC50_WEAK_MAX} nM for Class I, "
      f"PRIME rank <= {PRIME_RANK_MAX} for Class I)")

# Tiered ranking: high-affinity binders first, then by immunogenicity, then by raw affinity
hard_filtered["HighAffinityTier"] = (hard_filtered["BindingAffinity"] < IC50_STRONG_MAX).fillna(False).astype(int)

ranked = hard_filtered.sort_values(
    by=["HighAffinityTier", "ImmunogenicityScore", "BindingAffinity"],
    ascending=[False, True, True],
    na_position="last",
)

top_candidates = ranked.head(TOP_N)
top_candidates.to_csv(config.TOP_CANDIDATES, sep="\t", index=False, na_rep="NA")
print(f"Wrote top {len(top_candidates)} candidates to {config.TOP_CANDIDATES}")
