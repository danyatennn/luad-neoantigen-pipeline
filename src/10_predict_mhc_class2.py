#!/usr/bin/env python3
"""
Part 12.2 - Peptide-MHC Binding Prediction: Class II, 15-mers (NetMHCIIpan)

Evaluate WildType/Mutant 15-mer peptide pairs (neoantigen candidate windows)
for HLA class II presentation using a local NetMHCIIpan 4.2 installation.

Expects a tab-delimited TSV with these columns (header row required):
    GenomicVariant  GeneName  TranscriptID  ProteinID  ProteinChange
    ProteinPosition  Length  MutPos  Peptide  Type
where Type is "WildType" or "Mutant". Each (GenomicVariant, ProteinPosition,
MutPos) triple identifies one sliding-window pair: the same 15-mer window
in wildtype vs mutant sequence.

Usage:
    python3 src/10_predict_mhc_class2.py \
        --netmhciipan-bin /path/to/netMHCIIpan-4.2/netMHCIIpan

    Input file and alleles default to the values in config.py; override with
    --tsv / --alleles if needed.

Output (paths come from config.py, all under data/interim/):
    netmhciipan_raw.txt              - raw netMHCIIpan stdout
    15mer_preds.csv                  - every input row x allele, with %Rank_EL + binding call
    classii_neoantigen_candidates.csv - Mutant peptides that bind (Strong/Weak) where the
                                  matched WildType peptide does not - i.e. peptides
                                  where the mutation appears to create new HLA-II
                                  presentation, sorted by strongest mutant rank

Requires pandas (pip install pandas --break-system-packages).
"""

import argparse
import subprocess
import tempfile
import os
import sys

import pandas as pd

import config


def parse_args():
    p = argparse.ArgumentParser(description="Run NetMHCIIpan 4.2 on WT/Mutant 15-mer peptide TSV")
    p.add_argument("--tsv", default=config.PEPTIDES_15MER_FILTERED,
                    help="Input TSV (GenomicVariant, GeneName, ..., Peptide, Type)")
    p.add_argument("--alleles", default=",".join(config.HLA_CLASS_II),
                    help="Comma-separated allele list (default: the Part 11 class II panel)")
    p.add_argument("--netmhciipan-bin", required=True, help="Path to the netMHCIIpan executable")
    p.add_argument("--strong-threshold", type=float, default=2.0, help="%%Rank_EL cutoff for 'Strong' binder")
    p.add_argument("--weak-threshold", type=float, default=10.0, help="%%Rank_EL cutoff for 'Weak' binder")
    return p.parse_args()


def run_netmhciipan(binary, pep_file, alleles, out_path):
    cmd = [binary, "-f", str(pep_file), "-inptype", "1", "-a", alleles]
    with open(out_path, "w") as out_fh:
        result = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"netMHCIIpan failed:\n{result.stderr}")


def parse_output(raw_path):
    """
    Parses netMHCIIpan 4.3 whitespace-delimited stdout into a list of dict rows.
    """
    fixed_cols = ["Pos", "MHC", "Peptide", "Of", "Core", "Core_Rel",
                  "Inverted", "Identity", "Score_EL", "Rank_EL", "Exp_Bind"]
    rows = []
    seen_header = False
    with open(raw_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or set(stripped) == {"-"}:
                continue
            if stripped.startswith("Pos") and "MHC" in stripped:
                seen_header = True
                continue
            if not seen_header:
                continue
            fields = stripped.split()
            if len(fields) < len(fixed_cols):
                continue
            row = dict(zip(fixed_cols, fields[:len(fixed_cols)]))
            extra = fields[len(fixed_cols):]
            row["BindLevel"] = " ".join(extra).replace("<=", "").strip() if extra else ""
            rows.append(row)
    return rows


def classify(rank, strong_cut, weak_cut):
    try:
        r = float(rank)
    except (TypeError, ValueError):
        return "NA"
    if r <= strong_cut:
        return "Strong"
    elif r <= weak_cut:
        return "Weak"
    return "Non-binder"


def main():
    args = parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    required_cols = {"GenomicVariant", "GeneName", "ProteinChange", "ProteinPosition", "MutPos", "Peptide", "Type"}
    missing = required_cols - set(df.columns)
    if missing:
        sys.exit(f"Input TSV is missing expected columns: {missing}")

    bad_len = df[df["Peptide"].astype(str).str.len() != 15]
    if len(bad_len):
        print(f"warning: dropping {len(bad_len)} rows where Peptide is not a 15-mer", file=sys.stderr)
        df = df[df["Peptide"].astype(str).str.len() == 15]

    unique_peptides = df["Peptide"].drop_duplicates().tolist()
    print(f"{len(df)} input rows -> {len(unique_peptides)} unique peptides to score "
          f"against alleles [{args.alleles}]")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".pep", delete=False) as tmp:
        tmp.write("\n".join(unique_peptides) + "\n")
        pep_file = tmp.name

    raw_out = config.CLASS2_RAW
    print("Running NetMHCIIpan")
    run_netmhciipan(args.netmhciipan_bin, pep_file, args.alleles, raw_out)
    os.unlink(pep_file)

    rows = parse_output(raw_out)
    if not rows:
        sys.exit(f"No parseable prediction rows found in {raw_out}. "
                  f"Open that file and check the header row or adjust parse_output().")

    pred_df = pd.DataFrame(rows)
    print(f"Parsed {len(pred_df)} prediction rows from {raw_out}")

    pred_df["Binding_Call"] = pred_df["Rank_EL"].apply(lambda r: classify(r, args.strong_threshold, args.weak_threshold))
    pred_df = pred_df.rename(columns={"MHC": "Allele"})
    pred_df = pred_df[["Peptide", "Allele", "Rank_EL", "Binding_Call"]]

    # merge predictions (one row per peptide x allele) onto every metadata row
    full = df.merge(pred_df, on="Peptide", how="left")
    full_csv = config.CLASS2_PREDICTIONS
    full.to_csv(full_csv, index=False)

    # pair WildType vs Mutant within the same mutation window, per allele
    key_cols = ["GenomicVariant", "GeneName", "ProteinChange", "ProteinPosition", "MutPos", "Allele"]
    wt = full[full["Type"] == "WildType"][key_cols + ["Rank_EL", "Binding_Call"]].rename(
        columns={"Rank_EL": "WT_Rank_EL", "Binding_Call": "WT_Binding_Call"})
    mut = full[full["Type"] == "Mutant"][key_cols + ["Peptide", "Rank_EL", "Binding_Call"]].rename(
        columns={"Rank_EL": "Mut_Rank_EL", "Binding_Call": "Mut_Binding_Call"})
    paired = mut.merge(wt, on=key_cols, how="left")

    def is_candidate(row):
        return row["Mut_Binding_Call"] in ("Strong", "Weak") and row["WT_Binding_Call"] not in ("Strong", "Weak")

    paired["Neoepitope_Candidate"] = paired.apply(is_candidate, axis=1)
    candidates = paired[paired["Neoepitope_Candidate"]].sort_values("Mut_Rank_EL")
    cand_csv = config.CLASS2_CANDIDATES
    candidates.to_csv(cand_csv, index=False)

    print(f"Raw output:            {raw_out}")
    print(f"Full predictions:      {full_csv}")
    print(f"Neoantigen candidates: {cand_csv}  ({len(candidates)} mutant windows where matched WT does not bind)")


if __name__ == "__main__":
    main()
