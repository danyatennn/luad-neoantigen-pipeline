# Project 130: Integrating Cancer Mutations, Gene Expression, and Neoantigen Prediction

**Group members:** Naige Lu, Daniil Ten, Hae-gyung (Tessa) Han

**Selected cancer type:** Lung adenocarcinoma (LUAD)

---

## Data sources

### Somatic mutations — NCI Genomic Data Commons

| | |
|---|---|
| Project ID | **TCGA-LUAD** (Lung Adenocarcinoma) |
| Program | The Cancer Genome Atlas (TCGA) |
| dbGaP Study Accession | phs000178 |
| Portal | https://portal.gdc.cancer.gov/ |
| Access date | 20 July 2026, 11:17 |
| Data Category | Simple Nucleotide Variation |
| Data Type | Masked Somatic Mutation |
| Experimental Strategy | Whole Exome Sequencing |
| Sequencing centre | Broad Institute (BI) |
| Cohort size | 194,729 records · 557 patients · 575 samples · 616 aliquots |


### Gene expression — NCBI Gene Expression Omnibus

| | |
|---|---|
| Accession | **GSE81089** |
| Series | Next Generation Sequencing (RNAseq) of non-small cell lung cancer |
| Platform | **GPL16791** — Illumina HiSeq 2500 |
| Download date | 22 July 2026 |
| Samples | 199 tumour, 19 matched normal (tumour only retained) |
| Values as supplied | FPKM (cufflinks) — **converted to TPM by us** |
| Gene identifier type | Ensembl gene ID (ENSG) |


### Reference genome assembly

**GRCh38 / hg38 throughout.** All 194,729 MAF records are annotated to GRCh38, verified
programmatically; VEP was run against GRCh38.p14 (Ensembl release 116); the expression matrix
is indexed by Ensembl gene identifiers and carries no coordinates

---

## Software

**Environment:** Python 3.11.5 on Ubuntu 26.04 LTS.

**Python packages** (see `requirements.txt`): pandas, numpy, matplotlib, requests, tqdm,
gdown, mhcflurry 2.2.1.

**External tools:**

| Tool | Version | Used in | Notes |
|---|---|---|---|
| Ensembl Variant Effect Predictor | release 116 (GRCh38.p14) | Part 9 | Web service; maps genomic variants to transcript and protein consequences. Accessed 23 July 2026 |
| MHCflurry | 2.2.1 | Part 12.1 | Class I binding, affinity mode. Chosen over NetMHCpan because it runs locally and installs via pip |
| NetMHCIIpan | 4.2 | Part 12.2 | Class II presentation, eluted-ligand mode. Free for academic use; requires local installation |
| PRIME | 2.1 | Part 13 | Class I immunogenicity. Requires local installation |

---

## How to run the pipeline

Nothing has to be downloaded by hand. From a fresh clone:

```bash
pip install -r requirements.txt
python src/00_fetch_inputs.py     # pulls every input (public sources + Google Drive)
python src/01_build_mutation_matrix.py
python src/02_build_expression_matrix.py
python src/03_integrate_mutation_expression.py
python src/04_qc_figures.py
python src/05_prepare_vep_input.py
# 06 needs the VEP web output, which 00_fetch_inputs.py has already downloaded
python src/06_parse_vep_output.py
python src/07_generate_peptides.py
python src/08_filter_peptides.py
python src/09_predict_mhc_class1.py
python src/11_predict_immunogenicity.py
python src/12_compare_wildtype.py
python src/13_build_final_table.py
```

`00_fetch_inputs.py` skips files that are already present, so it is safe to re-run.

Two stages are not run by the commands above because they need licensed software installed
locally (NetMHCIIpan 4.2 and PRIME 2.1). Their outputs are downloaded from Google Drive by
step 00, so the rest of the pipeline still runs end to end. To reproduce them yourself:

```bash
python src/10_predict_mhc_class2.py --netmhciipan-bin /path/to/netMHCIIpan
```

All paths and the HLA panel live in `src/config.py`; no script hard-codes a path.

### Pipeline overview

| Step | Assignment part | Produces |
|---|---|---|
| `00_fetch_inputs` | — | Downloads all inputs |
| `01_build_mutation_matrix` | 4–5 | `data/01_mutation_by_sample.tsv` |
| `02_build_expression_matrix` | 6 | `data/02_gene_by_sample_TPM.tsv` |
| `03_integrate_mutation_expression` | 7 | `data/03_integrated_mutation_expression.tsv` |
| `04_qc_figures` | 8 | Four figures in `figures/` |
| `05_prepare_vep_input` | 9 | VCF for the VEP web service |
| *(manual)* | 9 | VEP output, downloaded by step 00 |
| `06_parse_vep_output` | 9 | Variant-to-protein annotation, MANE Select only |
| `07_generate_peptides` | 10 | Mutant and wild-type 9-mers and 15-mers |
| `08_filter_peptides` | 10 | 15-mers restricted to the ten most mutated genes |
| `09_predict_mhc_class1` | 12.1 | Class I binding predictions and candidates |
| `10_predict_mhc_class2` | 12.2 | Class II predictions *(needs NetMHCIIpan)* |
| `11_predict_immunogenicity` | 13 | PRIME immunogenicity *(needs PRIME output)* |
| `12_compare_wildtype` | 14 | Wild-type comparison and prioritisation score |
| `13_build_final_table` | 15 | `data/04_neoantigen_predictions.tsv` |

Directory layout: `src/` code · `data/` the four submitted tables · `data/raw/` downloaded
inputs · `data/interim/` intermediates, not submitted · `figures/` · `docs/` · `notebooks/`
exploratory analysis · `bin/` superseded drafts.

---

## Output columns

### `01_mutation_by_sample.tsv` — binary mutation-by-sample matrix

| Column | Meaning |
|---|---|
| `Gene_Name` | HGNC gene symbol |
| `Mutation` | Coding DNA change, HGVS (e.g. `c.35G>A`) |
| `AminoAcid_Change` | Protein change, HGVS 3-letter (e.g. `p.Gly12Asp`) |
| `TCGA-…` (575 columns) | 1 if the mutation is present in that tumour sample, 0 otherwise |

One row per distinct mutation, defined by the triple (gene, coding change, protein change);
several mutations in one gene are separate rows. TCGA barcodes identify aliquots rather than
tumours, so they are truncated to their first 16 characters, collapsing repeated sequencing of
the same specimen into a single column.

### `02_gene_by_sample_TPM.tsv` — gene-by-sample expression matrix

| Column | Meaning |
|---|---|
| `GeneName` | **Ensembl gene identifier** as supplied by GSE81089 (see limitations) |
| `L400T`, `L401T`, … (199 columns) | TPM for that gene in that tumour sample |

FPKM values were converted to TPM by renormalising each sample to a column sum of 10⁶.

### `03_integrated_mutation_expression.tsv` — integrated matrix

Columns of `01` plus `GeneLevelTPM`, inserted after `AminoAcid_Change`: the **median TPM
across the 199 expression tumour samples**, chosen for robustness to outliers. `NA` where the
gene is absent from the expression dataset.

### `04_neoantigen_predictions.tsv` — neoantigen table

| Column | Meaning |
|---|---|
| `GeneName` | HGNC gene symbol |
| `Chromosome`, `Position`, `Ref`, `Alt` | Genomic variant on GRCh38 |
| `TranscriptID` | Ensembl MANE Select transcript |
| `ProteinChange` | Protein change, HGVS 3-letter |
| `GeneLevelTPM` | Median TPM of the gene across tumour samples |
| `MutationFrequency` | Fraction of tumour samples carrying this variant |
| `PeptideType` | `Mutant` or `WildType` |
| `Peptide` | Peptide sequence |
| `PeptideLength` | 9 (class I) or 15 (class II) |
| `MutationPosition` | 1-based position of the substituted residue within the peptide |
| `HLAAllele` | Allele the scores below refer to |
| `BindingAffinity` | Predicted IC50 in nM (class I only; lower is stronger) |
| `BindingRank` | Percentile rank (lower is stronger) |
| `BindingScore` | `NA` — neither predictor was run in a mode reporting a separate score |
| `ImmunogenicityScore` | PRIME %Rank (lower is more immunogenic); `NA` for wild-type and class II |
| `PredictionTool`, `ToolVersion` | Predictor that produced the binding values in that row |

Missing tool outputs are written as `NA`, never as zero.

---

## Known limitations

**Predictions only.** Every result is computational. No peptide was tested experimentally, so
nothing here establishes antigen presentation or T-cell recognition.

**Two independent cohorts.** Mutations come from TCGA-LUAD and expression from GSE81089, so
no patient-level matching is possible. `GeneLevelTPM` is a cohort-level summary, and a high
value shows that a gene is transcriptionally active in LUAD generally — not that the mutant
allele is transcribed in the tumour carrying it.

**Fixed HLA panel.** Patient-specific typing (Option B) via the TCGA PanImmune Atlas is
controlled-access under dbGaP and was not obtainable within the assignment timeframe, so a
fixed panel of common alleles was used. Coverage is skewed towards Western European
populations, and results do not generalise to patients carrying other alleles.

**Class II immunogenicity was not predicted.** PRIME is a class I method. No equivalent class
II predictor was integrated, so CD4⁺ candidates are prioritised on binding alone. The class II
analysis is also restricted to the ten most frequently mutated genes, because NetMHCIIpan is
too slow for all 1.67 million mutant 15-mers.

**Missense SNVs only.** In-frame indels and frameshifts, which can generate long novel
peptide stretches, are excluded. This also means driver frequencies reported here fall below
published figures for genes inactivated mainly by truncating mutations — STK11 and RBM10 in
particular.

**Recurrence is confounded by gene length.** The most frequently mutated genes (TTN, MUC16,
CSMD3, LRP1B) are long and weakly expressed; mutation frequency correlates with protein length
across all genes (r = 0.59). Formal driver identification would require length- and
covariate-corrected methods such as MutSigCV or dNdScv, which were outside scope.

**Expression coverage.** 1,313 mutation rows (1.1 %, 195 genes) have no `GeneLevelTPM` because
those genes are absent from GSE81089. A further 409 variants carry different gene symbols in
VEP (Ensembl 116) and in the MAF owing to nomenclature updates, for example
*BLTP3B*/*UHRF1BP1L*; their coordinates and alleles agree exactly.
