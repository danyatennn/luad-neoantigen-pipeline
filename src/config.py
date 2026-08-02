"""
Shared configuration
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

DATA = ROOT / "data"
RAW = DATA / "raw"
INTERIM = DATA / "interim"
FIGURES = ROOT / "figures"

INTERIM.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)


# raw inputs
MAF = RAW / "cohortMAF.2026-07-20.maf.gz"
GEO_FPKM = RAW / "GSE81089_FPKM_cufflinks.tsv.gz"
PROTEOME = RAW / "Homo_sapiens.GRCh38.pep.all.fa.gz"
VEP_OUTPUT = RAW / "vep_output_GRCh38.txt"  # downloaded from the VEP web service

GEO_FPKM_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81089"
    "&format=file&file=GSE81089%5FFPKM%5Fcufflinks%2Etsv%2Egz"
)
ENSEMBL_PEP_URL = (
    "https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/pep/"
    "Homo_sapiens.GRCh38.pep.all.fa.gz"
)


# deliverables
MUTATION_MATRIX = DATA / "01_mutation_by_sample.tsv"
TPM_MATRIX = DATA / "02_gene_by_sample_TPM.tsv"
INTEGRATED_MATRIX = DATA / "03_integrated_mutation_expression.tsv"
NEOANTIGEN_TABLE = DATA / "04_neoantigen_predictions.tsv"


# intermediate tables
VEP_INPUT = INTERIM / "vep_input.vcf"
VARIANT_ANNOTATION = INTERIM / "variant_protein_annotation.tsv"
PEPTIDES = INTERIM / "mutant_peptides.tsv"
PEPTIDES_15MER_FILTERED = INTERIM / "15mers_filtered.tsv"

CLASS1_PREDICTIONS = INTERIM / "9mer_predictions.tsv"
CLASS1_CANDIDATES = INTERIM / "classI_neoantigen_candidates.tsv"

CLASS2_RAW = INTERIM / "netmhciipan_raw.txt"
CLASS2_PREDICTIONS = INTERIM / "15mer_preds.tsv"
CLASS2_CANDIDATES = INTERIM / "classii_neoantigen_candidates.tsv"

MUTANT_9MERS_FASTA = INTERIM / "mutant_9mers.fasta"
PRIME_OUTPUT = INTERIM / "prime_mutant_9mers.txt"  # produced by running PRIME
IMMUNOGENICITY = INTERIM / "part13_immunogenicity_results.tsv"

WT_VS_MUTANT = INTERIM / "wt_vs_mutant_scored.tsv"
CANDIDATES_SCORE_6 = INTERIM / "candidates_prioritisation_score_6.tsv"


# figures
FIG_TPM_HISTOGRAM = FIGURES / "tpm_distribution_log.png"
FIG_MUTATION_VS_EXPRESSION = FIGURES / "mutation_vs_expression.png"
FIG_MUTATION_HEATMAP = FIGURES / "mutation_heatmap.png"
FIG_TOP_MUTATED_GENES = FIGURES / "top_mutated_genes.png"


# HLA panel
HLA_CLASS_I = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01"]
HLA_CLASS_II = ["DRB1_0101", "DRB1_0701"]

PUBLIC_URLS = {
    GEO_FPKM: GEO_FPKM_URL,
    PROTEOME: ENSEMBL_PEP_URL,
}

DRIVE_IDS = {
    MAF: "1tZ3FVMv__hWnGByNSAIgevIc_mGRQY3C",
    VEP_OUTPUT: "1lTw3Ldde-Vh97GNBDIWHSE4SH7WKbgFE",
    PEPTIDES_15MER_FILTERED: "1j2cZj-ZcgWoOwTWJdU7YMBs0cSP-xB3y",
    CLASS1_PREDICTIONS: "1FZdtoWRDdxyikoiYCZKanAYr8-d5bNT7",
    CLASS1_CANDIDATES: "1QrXeUIzLMMSI21Xe97k_NW_DjYMcYSCn",
    CLASS2_RAW: "1nxxaOkk4H3be2_uJrd4G26UVfI5Y6dMf",
    CLASS2_PREDICTIONS: "1MtO3uYwLXWLc4gVopkMOEuCgKjVFdbvx",
    CLASS2_CANDIDATES: "1Pmxoo8ffjsUwW4YY075VEMAS1f588_U6",
    PRIME_OUTPUT: "1YlJScPrmcBUgSzRUIfk3i9RDpL5P_1Hc",
    IMMUNOGENICITY: "1ceTAFoK6zfg3aUFtmhzFjLsQH09y2V63",
}
