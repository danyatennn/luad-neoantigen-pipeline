# Code review — Project 130 (Lung Cancer / LUAD)

Review of all scripts and notebooks against `Assignment_Project130.docx`.
No code was modified — the following are findings only, each one verified against the actual data in `data/`.

Tags: **[BLOCK]** — breaks the result, **[REPRO]** — reproducibility (15 points),
**[LOGIC]** — logic error, **[FORMAT]** — does not match the required submission format, **[NIT]** — minor.

---

## 1. Blockers

### 1.1 [BLOCK] `Part15.py:39-43` — merging on `GeneName` alone pairs peptides with the wrong affinities

```python
class1_9mers = HLA1_immunogenicity.merge(ninemers_presentation, on="GeneName", how="left")
```

Both tables have 51,750 rows, and each gene accounts for dozens of rows in each (`TTN`: 226 and 226).
Joining on gene alone produces a Cartesian product within every gene: **~623,000 rows instead of 51,750 (×12)**.

The final table is then assembled (`Part15.py:91-107`) as follows:

- `Peptide` ← `Peptide_Mutant_x`, `HLAAllele` ← `Allele_x` — from the **left** table;
- `BindingAffinity` ← `Affinity_nM_Mutant`, `BindingRank` ← `Percentile_Rank_Mutant`,
  `MutationPosition` ← `MutPos`, `TranscriptID` — from the **right** table.

So a single row ends up holding a peptide from one record and the affinity of a **different peptide**.
Verified on `TTN` (first 2×3 rows):

```
Peptide_Mutant_x   Peptide_Mutant_y   Affinity_nM_Mutant
KAQLARQQY          KAQLARQQY           867.18   <- matches
KAQLARQQY          KAGVGEHAY           745.38   <- peptide and affinity from different rows
KAQLARQQY          VATVVAAVY          1542.89
KAGVGEHAY          KAQLARQQY           867.18
...
4 of 6 rows mismatched
```

This is the final deliverable (`04_neoantigen_predictions.tsv`) — the artefact the "MHC-binding prediction"
and "Candidate prioritization" criteria are graded on. Right now it contains incorrect peptide↔score pairs.

**How to fix:** merge on a key that uniquely identifies a row:
`on=["GeneName", "ProteinChange", "Peptide_Mutant", "Allele"]`.
All four columns exist in both tables. The `_x`/`_y` suffixes will then disappear on their own —
their presence in the code is precisely the signal that the join is ambiguous.

---

### 1.2 [BLOCK] `Part14.py:88-96` — `tpm_lookup` is not deduplicated, rows multiply 24-fold

```python
tpm_lookup = mutation_by_sample_matrix[["GeneName", "GeneLevelTPM", "TumourCount"]].copy()
HLA_1 = first_merge.merge(tpm_lookup, on="GeneName", how="left")
```

`mutation_by_sample_matrix` is `integrated_mutation_expression_matrix.tsv`, where **one row = one mutation**
(115,180 rows across 16,207 genes, ~7 mutations per gene on average). Merging on `GeneName` without
aggregation gives:

```
classI_neoantigen_candidates.tsv       51,750 rows
Part 14 All_Combined_tables_HLA_1.tsv   1,268,351 rows  (×24.5, 295 MB file)
```

Beyond the bloat, this **makes `TumourCount` wrong**: each peptide is assigned the sample count of an
*arbitrary* mutation in the same gene rather than the one the peptide was derived from.
`GeneLevelTPM` is merely duplicated with an identical value (harmless), but `PrioritisationScore` is
computed against 24 different `TumourCount` values for one and the same candidate — the score stops
being well defined.

**How to fix:** split into two separate joins.
- `GeneLevelTPM` — by gene, but with `.drop_duplicates("GeneName")` first (there is only one value per gene anyway).
- `TumourCount` — by the **specific mutation**: `on=["GeneName", "ProteinChange"]`.
  The formats already match, verified: in the integrated matrix `AminoAcid_Change` = `p.Phe81Leu`,
  in the candidates `ProteinChange` = `p.Asp1341Tyr` — both 3-letter HGVS. Only a rename is needed.

### 1.3 [BLOCK] `Part15.py:138-152` — the same mistake a second time

`expression_lookup` is taken from the same per-mutation table and merged on `GeneName` without
deduplication. It stacks on top of the already-inflated join from 1.1 — the two blow-ups multiply.

---

### 1.4 [BLOCK] `Project 130 Lung Cancer.ipynb`, cell[12] — the `GDC_FILTER == "PASS"` filter yields an empty table

```python
filtered_maf = maf[... & (maf["GDC_FILTER"] == "PASS")].copy()
```

The `GDC_FILTER` values in this MAF (from the output of cell[7] in that same notebook):

```
NaN          191643
NonExonic      3075
gdc_pon          11
```

There is no `"PASS"` string at all — PASS is encoded as an empty value. `filtered_maf` comes out **empty**.

`task5.ipynb` cell[2] does this correctly (`GDC_FILTER.isna() | == ""`), so the pipeline is unaffected.
But the notebook is part of the submission, and if any QC numbers for the report were taken from here,
they are wrong. Either fix it or delete the cell.

---

## 2. Reproducibility

### 2.1 [REPRO] Absolute paths pointing at other people's machines — the pipeline runs nowhere

| File | Lines | Path |
|---|---|---|
| `Part13.py` | 5, 6, 66, 156 | `/Users/naigelu/Desktop/Imperial College London/...` |
| `Part14.py` | 5-8, 139, 156 | same |
| `Part15.py` | 5-8, 162 | same |
| `part12_1.py` | 19 | `C:\Users\lyssa\Downloads\mutant_peptides.tsv` |
| `opengzip.py` | 3 | `C:\Users\lyssa\Documents\GitHub\...` |

Assignment, Part V.3: *"All major steps must be reproducible from submitted code"*.
As it stands, 5 of 11 scripts will not start on a clean machine. The minimum is relative paths from the
repository root (`data/...`); better still, a shared module of path constants or `argparse`, as already
done in `run_netmhciipan.py`.

### 2.2 [REPRO] No script produces `15mers_filtered.tsv`

```
data/peptides_15mers.tsv    3,337,140 rows   (= 1,668,570 WT + 1,668,570 Mutant, full Part 10 output)
data/15mers_filtered.tsv        63,322 rows   (53× fewer)
```

`run_netmhciipan.py` takes `--tsv 15mers_filtered.tsv`, but **what reduced 3.3 M rows to 63 K, and by what
criterion**, is nowhere in the repository. Likewise there is no script extracting `peptides_15mers.tsv`
from `mutant_peptides.tsv`. This is a direct gap against Part V.3 and against the "Organised and readable
code" criterion. The filtering needs to be written up as a script and described in the README
(how many rows remained and why).

### 2.3 [REPRO] File names do not line up between steps — Part 15 will not run

```
Part13.py:156  writes    part13_HLA1_immunogenicity_results.tsv
Part14.py:8    reads     part13_immunogenicity_results.tsv        <- different name
Part15.py:5    reads     part13_HLA1_immunogenicity_results.tsv
data/          contains  part13_immunogenicity_results.tsv        <- only this one
```

There is no `part13_HLA1_immunogenicity_results.tsv` in `data/` — `Part15.py` will fail with
`FileNotFoundError`. One name is needed in all three places.

### 2.4 [REPRO] `.idea/` is committed, `data/` is not

`git ls-files` shows `.idea/workspace.xml` and `.idea/runConfigurations/Part15.xml` in the repository,
while the entire `data/` folder holding the deliverables is untracked. These should be swapped:
`.idea/` into `.gitignore`, and the results (or at least the small final tables) into the repository.

---

## 3. Logic and correctness

### 3.1 [LOGIC] `Part14.py:101-134` — 2 of the 6 PrioritisationScore criteria are always true

Criterion 5 is `Percentile_Rank_WildType > 2`. But the candidates in `classI_neoantigen_candidates.tsv`
were already selected in `part12_1.py:127-130` by exactly that condition (`CANDIDATE_WT_RANK_MIN = 2.0`).
Verified against the file:

```
Percentile_Rank_WildType > 2   -> 100.0 %  of rows   (constant, always +1)
TumourCount >= 1               -> effectively always true (a candidate exists => the mutation is in the cohort)
Percentile_Rank_Mutant <= 0.5  ->  14.0 %
AffinityFoldChange >= 2        ->  87.1 %
```

Net effect: the scale is presented as 6-point, but it actually discriminates candidates on only 2-4
features, and the minimum attainable score is 2, not 0. "Score = 6" does not mean what the comment claims.
Either compute the score before candidate filtering (on the full `9mer_predictions.tsv` table), or state
plainly in the report that two criteria hold by construction.

### 3.2 [LOGIC] `Part14.py:120-133` — comments disagree with the code

| Comment (lines 55-57) | Code |
|---|---|
| `GeneLevelTPM > 0` | `GeneLevelTPM >= 15` (lines 121-123) |
| `Mutation count more than 1 tumour` | `TumourCount >= 1` (lines 131-133) |

The 15 TPM threshold is nowhere justified and is not mentioned in the report. Either bring the code in
line with the comment or justify the threshold (either way, in the README/report, since it affects the
final selection).

### 3.3 [LOGIC] `run_netmhciipan.py:166` — class II candidates are sorted as strings

`parse_output()` (lines 67-93) puts **strings** from `split()` into the dict, so `pred_df` is entirely
dtype `object`. `classify()` casts to float only for the classification; the column itself stays a string.
Then:

```python
candidates = paired[paired["Neoepitope_Candidate"]].sort_values("Mut_Rank_EL")
```

— a lexicographic sort, where `"10.5" < "2.3"`. Verified on `classii_neoantigen_candidates.csv`:
as strings the column is monotonic, as numbers it is not (first violation at row 142).
The "top candidates" cannot be taken from this file.

**Fix:** `pred_df["Rank_EL"] = pd.to_numeric(pred_df["Rank_EL"], errors="coerce")` after parsing.

### 3.4 [LOGIC] `task5.ipynb` cell[7] — `savefig` after `show`, the heatmap is never saved

```python
plt.show()
plt.savefig("130heatmap.png")
```

`plt.show()` clears the current figure, so `savefig` writes a blank image.
And indeed `130heatmap.png` is not in the repository. The bar plot in cell[8] has no `savefig` at all.

Assignment Part 8 requires at least 2 figures — `tpm_distribution_log.png` and `mutation_vs_expression.png`
are actually saved, so the minimum is formally met, but the heatmap and the bar plot (two of the four
figures suggested in the assignment) never reach disk.

### 3.5 [LOGIC] Missing values are treated as "bad" rather than as NA

Assignment Part V.10 / Part 15: *"Represent missing tool outputs as NA, not zero"*.

- `part12_1.py:37-43` — `classify_binder(NaN)`: every comparison returns `False`, so `"Non-binder"` is
  returned. It should be `"NA"`.
- `Part14.py:103-133` — `pd.to_numeric(...)` without `errors="coerce"`, and `NaN <= 0.5` → `False`.
  So "immunogenicity could not be computed" and "immunogenicity is poor" contribute the same 0 to the score.

### 3.6 [LOGIC] `02_gene_by_sample_TPM.tsv`: the column is called `GeneName` but holds Ensembl IDs

```
GeneName          L400T      L401T ...
ENSG00000000003   75.169...  35.472...
```

`part6.py:75` renames the index to `GeneName`, although the source column in GSE81089 is `Ensembl_gene_id`.
The assignment (Part 6) requires a matrix keyed by `GeneName` (examples: ARID1A, TP53, KRAS) and separately
requires reporting the *"Method used to map identifiers to gene symbols"*.

Consequence: `part6.py:76` `tpm.groupby("GeneName").median()` is a **no-op**, because ENSG IDs are unique
(verified: 0 duplicated index entries). So the *"Method used to handle duplicated gene symbols"* requirement
is effectively unmet — there was nothing to collapse.

The mapping happens later, in `task7.ipynb` cell[3-4], via the `Gene` column from the MAF. It works, but
coverage is incomplete: **1,313 rows (1.1 %, 195 genes) end up without `GeneLevelTPM`**.
Either name the column honestly (`EnsemblGeneID`) or map to symbols directly in Part 6 — and in either
case document the method and the losses in the README.

### 3.7 [LOGIC] The matrix is built over aliquots (616), not samples (575)

`task5.ipynb` cell[3] builds the `crosstab` over the full `Tumor_Sample_Barcode`
(`TCGA-86-A4D0-01A-11D-A24D-08`) — that is an **aliquot**, not a sample. This gives 616 columns.

Meanwhile `Project 130 Lung Cancer.ipynb` cell[10-11] already computes
`Sample_ID = barcode[:16]` → **575 unique samples**. But `Sample_ID` is never used afterwards.

One sample can enter the matrix as several columns → the mutation-frequency denominator is inflated by
~7 %, and `MutationFrequency` in Part 15 is systematically understated. The assignment requirement
("At least 20 tumor samples") is met with an enormous margin either way, but the sample count in the
report must be stated correctly (575 samples / 616 aliquots, not "616 samples").

### 3.8 [LOGIC] Part 5 and Part 9 filter the MAF differently

| Step | Rule | Unique variants |
|---|---|---|
| `task5.ipynb:2` (Part 5) | missense **+ SNP + PASS** | 115,180 |
| `09a_prepare_vep_input.py:10` (Part 9) | missense only | 115,195 |

The difference is 15 variants: 11 flagged `gdc_pon`, 4 ONP/TNP. Trivial in volume, but it means some
peptides in Parts 10-15 originate from mutations that are **not** in the Part 5 matrix, and no
`TumourCount`/`MutationFrequency` will be found for them. Better to apply the same filter in both places
(or describe the discrepancy explicitly in the README).

---

## 4. Submission-format compliance

### 4.1 [FORMAT] Deliverable names do not match the required ones (Part III.A)

| Required | Present in `data/` |
|---|---|
| `01_mutation_by_sample.tsv` | ✅ present |
| `02_gene_by_sample_TPM.tsv` | ✅ present |
| `03_integrated_mutation_expression.tsv` | ❌ named `integrated_mutation_expression_matrix.tsv` |
| `04_neoantigen_predictions.tsv` | ❌ **missing** |

The Part 15 final table is not in `data/` at all — `Part15.py` writes it to
`/Users/naigelu/.../Part 15/final_neoantigen_table.tsv`. This is the headline artefact of the advanced
component; it needs to land in `data/04_neoantigen_predictions.tsv`.

### 4.2 [FORMAT] Some tables are CSV where TSV is required

Part V.1: *"All primary output tables must be tab-delimited (.tsv)"*.

- `run_netmhciipan.py:151-168` writes `predictions_full.csv` and `neoantigen_candidates.csv`;
- `data/` holds `classii_neoantigen_candidates.csv` and `15mer_preds.csv`;
- `Part14.py:6` and `Part15.py:7` read them as CSV.

Switching to `sep="\t"` and the matching extension is enough (with the reads updated in step).

### 4.3 [FORMAT] `Part15.py:91-129` — the final table does not match the Part 15 specification

Required columns vs what is actually assembled:

| Column | Status |
|---|---|
| `Chromosome` | ❌ filled with the string `"hg38/GRCh38"` (lines 93, 112) — that is the assembly, not a chromosome |
| `Position` | ❌ populated with the **protein** position `ProteinPosition`, not a genomic coordinate |
| `Ref` / `Alt` | ❌ columns absent entirely |
| `TranscriptID` | ⚠️ hard-coded to `pd.NA` for 15-mers (line 114), although it exists in `15mers_filtered.tsv` |
| `BindingScore` | ❌ absent |
| `PredictionTool` | ❌ absent (even though it is recorded in `9mer_predictions.tsv`) |
| `ToolVersion` | ❌ absent (MHCflurry 2.2.1 is recorded in `9mer_predictions.tsv`) |
| `PeptideType` | ⚠️ only `"Mutant"` — no WildType rows, although the assignment example shows both |

`Chromosome`/`Position`/`Ref`/`Alt` are recoverable from `GenomicVariant` — it is present in the source
tables as `1:69744-69744 C>G` and parses with a single regular expression.

### 4.4 [NIT] Column names in the integrated matrix differ from the Part 7 example

`task7.ipynb` cell[4] emits `Gene_Name / Mutation / AminoAcid_Change / GeneLevelTPM / …`, whereas the
Part 7 example in the assignment uses `GeneName / Mutation / AminoAcidChange / GeneLevelTPM / …`
(no underscores). The assignment is itself inconsistent (the Part 5 example does use underscores), so this
is not an error — but since the Part 15 final table requires `GeneName`, it is simpler to standardise on
the underscore-free variant; the joins in 1.2 would then need no `rename` either.

The column order is correct (metadata, then `GeneLevelTPM`, then samples), and the aggregation is
`median` across tumour samples (`task7.ipynb` cell[2]), exactly as the assignment requires.

### 4.5 [FORMAT] Scripts are not numbered in execution order (Part III.B)

The assignment asks for `01_download_mutations.py`, `02_build_mutation_matrix.py`, … Currently:

```
task5.ipynb  part6.py  task7.ipynb  09a_...py  09b_...py  10_peptides.py
part12_1.py  run_netmhciipan.py  Part13.py  Part14.py  Part15.py  plot.py  opengzip.py
```

Three different naming styles (`task5` / `part6` / `Part13` / `09a_`), a mix of `.py` and `.ipynb`, and
numbers that do not follow the run order. Standardise on a single `NN_description.py` scheme — this sits
directly under the "Reproducibility and presentation" criteria.

### 4.6 [FORMAT] The README is not filled in

Required (Part III.C), currently missing:

- dataset accession numbers (present only in the report draft: TCGA-LUAD + GSE81089);
- reference genome assembly (GRCh38 — in the code, not in the README);
- pipeline run instructions (section empty);
- an explanation of **all** output columns (section empty, with the leftover note
  `i was kind of unclear... - tessa` — the assignment wants the columns of the Part 15 final table,
  per "Explanation of all output columns");
- known limitations (section empty, although they are already written in the report draft — move them over);
- software versions: Python 3.11.5 / Ubuntu are given, but there are no versions for NetMHCIIpan (4.2),
  PRIME (2.1), MHCflurry (2.2.1 — recorded in `9mer_predictions.tsv`), VEP (Ensembl 116), nor access dates
  for the web tools (Part V.4 requires them explicitly; the VEP job date is visible on the screenshots in
  `scratchpad_pdf/` — 23.07.2026).

### 4.7 [FORMAT] Class II alleles disagree between the code and the report

- `run_netmhciipan.py:50` — default `DRB1_0101,DRB1_0701`;
- the report draft says "DRB1\*01:01, DRB1\*07:01" in one place and "HLA-DRB1\*07:01 and HLA-DRB1\*15:01"
  in another.

Part V.7: *"Do not report peptide-MHC scores without specifying the HLA allele"*. One consistent statement
is needed in the README and the report, matching what was actually run.

---

## 5. Minor items

- **`opengzip.py`** — dead code: someone else's Windows path, just prints the MAF to stdout. Delete it or
  rework it into `01_download_mutations.py`.
- **`exp.ipynb`** — a single cell with `maf.head()`; does not belong in the repository.
- **`part6.py:107`** — `plt.show()` in a script blocks execution in a headless environment; paired with
  `savefig` it is redundant.
- **Duplicated large files in the root and in `data/`.** `02_gene_by_sample_TPM.tsv`, `integrated_...tsv`,
  `mutant_peptides.tsv`, `variant_protein_annotation.tsv` all exist in two places. Notably the two copies
  of `02_gene_by_sample_TPM.tsv` **differ only in line endings** (CRLF in `data/`, LF in the root);
  the data is identical. Have the scripts write straight into `data/` and keep the root clean.
- **`QgqZf4GgO7KDUMli.txt` / `QgqZf4GgO7KDUMli.vep.txt`** — VEP job names that say nothing. Also, `.vep.txt`
  is the older output format **without** the `SYMBOL`/`MANE`/`ENSP` columns; the one actually used is
  `.txt` (`09b_vep_annotation.py:4`). The second file only causes confusion — drop it and rename the first
  (e.g. `vep_output_GRCh38.txt`).
- **`part12_1.py:114`** — the helper column `Nonstandard_AA` leaks into `9mer_predictions.tsv`.
  Better dropped from the final output.
- **`run_netmhciipan.py:69`** — the docstring says "netMHCIIpan 4.3" while everything else states 4.2.
- **Different logarithm bases across figures.** `part6.py:99` uses `np.log2` with the label
  "log2(TPM + 1)"; `plot.py:15` uses `np.log1p` (natural) with the label "log(TPM + 1)". One base is
  better for the report.
- **`Part14.py:200-234`** — ~35 blank lines at the end of the file.
- **`Part13.py:47-54`** — the "read 10 FASTA lines" block is entirely commented out inside, so the loop
  does nothing.
- **`09a_prepare_vep_input.py:28-30`** — the VCF is written without coordinate sorting. VEP web accepts
  this, but for reuse of the file (`bcftools`, tabix) it is worth sorting.

---

## 6. What is done well (leave alone)

- **`10_peptides.py`** — the cleanest file in the project. `windows()` (lines 33-36) correctly generates
  every window containing the mutation position without running past the protein boundary; the reference
  amino-acid check against the FASTA (lines 55-57) and the four `assert`s (lines 63-66) cover exactly the
  four requirements the assignment lists in Part 10. Skip counters are printed. The result checks out:
  1,009,767 9-mers and 1,668,570 15-mers, WT and Mutant strictly balanced.
- **`09b_vep_annotation.py`** — a single, explicitly documented transcript-selection strategy (MANE Select);
  all 10 required Part 9 fields are present; no duplicates on `(Genomic_Variant, Transcript_ID)`; no gaps
  in `Protein_ID`/`Protein_Change`.
- **`part6.py`** — FPKM→TPM is done correctly (per-column renormalisation to 1e6); verified: all 199 columns
  of the resulting matrix sum to exactly 1,000,000. The tumour/normal split by suffix works correctly
  (199 T + 19 N, no "other" columns). The assignment's rule against labelling FPKM as TPM without
  conversion is satisfied.
- **`run_netmhciipan.py`** — the only script with a proper CLI (`argparse`), required-column checks,
  peptide-length validation and a clear docstring. The other scripts should be modelled on it.
- **`part12_1.py`** — predictions are computed over unique peptides and then joined back (an order of
  magnitude saved); software/version/mode are recorded (exactly what Part 12.1 asks for); there is an
  explicit Mutant/WildType pairing check and non-standard amino acids are filtered out.
- Assemblies are nowhere mixed: MAF is GRCh38 (verified, 194,729/194,729 rows), VEP is GRCh38.p14, and
  expression uses Ensembl IDs without coordinates. Part V.5 is satisfied.

---

## 7. Priority of actions

1. Fix the joins: `Part15.py:39` and `Part14.py:92` (items 1.1-1.3) — without this the final table is wrong.
2. Generate and place `data/04_neoantigen_predictions.tsv` with the full Part 15 column set (items 4.1, 4.3).
3. Remove the absolute paths and reconcile file names across Parts 13/14/15 (items 2.1, 2.3).
4. Add the 15-mer filtering script (item 2.2).
5. Fill in the README (item 4.6) and renumber the scripts in execution order (item 4.5).
6. Everything else — per the list above.
