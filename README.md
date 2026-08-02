# Project 130: Integrating Cancer Mutations, Gene Expression, and Neoantigen Prediction
Group Members: Naige Lu, Daniil Ten, Hae-gyung (Tessa) Han

Selected cancer type: Lung Adenocarcinoma (LUAD)

Resources: 

GDC Project ID: TCGA-LUAD

Project name: Lung Adenocarcinoma

Program: The Cancer Genome Atlas (TCGA)

Access date: Monday, 20 July 2026 at 11:17

Filters: 
	Data category – simplified nucleotide variation
	Data Type – Masked Somatic Mutations
	Experimental Strategy – Whole Exome Sequencing


Download dates: July 22, 2026

Reference genome assembly:

Software dependencies:
- Python libraries: numpy, pandas, requests, matplotlib, gzip, os, re, mhcflurry, argparse, subprocess, tempfile, sys, pathlib
- Ensembl Variant Effect Predictor (Ensembl VEP) release 116: WRITE WHAT IT DOES
- NetMHCIIpan 4.2: This is an HLA class II predictor which can be obtained for free for academic purposes which we used for part 12.2. 
- PRedictor of IMmunogenic Epitopes (PRIME) Vers. 2.1: This is a class I immunogenicity predictor used for part 13.

Software versions: Python 3.11.5, Ubuntu 26.04 LTS

## How to run the pipeline:

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
python src/09_predict_mhc_class1.py
python src/11_predict_immunogenicity.py
python src/12_compare_wildtype.py
python src/13_build_final_table.py
```

`00_fetch_inputs.py` skips files that are already present, so it is safe to re-run.

Two steps are not run by the commands above because they need licensed software
installed locally (NetMHCIIpan 4.2 and PRIME 2.1). Their outputs are downloaded
from Google Drive by step 00, so the rest of the pipeline still runs end to end.
To reproduce them yourself:

```bash
python src/10_predict_mhc_class2.py --netmhciipan-bin /path/to/netMHCIIpan
```

All paths and the HLA panel live in `src/config.py`; no script hard-codes a path.

### Pipeline overview:

## Output columns
i was kind of unclear on if this is like the final output columns or the columns for each part of the pipeline... - tessa

## Known limitations
