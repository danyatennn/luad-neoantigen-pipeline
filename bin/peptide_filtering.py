import pandas as pd

INPUT_FILE = "peptides_15mers.tsv" # full 15mer list
OUTPUT_FILE = "peptides_15mers_filtered.tsv"

# TO EDIT: top mutated genes from Part 8
TOP_GENES = ["KRAS", "EGFR", "TP53", "BRAF"]

df = pd.read_csv(INPUT_FILE, sep="\t")
print(f"Starting rows: {len(df)}")
print(f"Starting unique peptides: {df['Peptide'].nunique()}")

filtered = df[df["GeneName"].isin(TOP_GENES)]

print(f"\nAfter filtering to {TOP_GENES}:")
print(f"Rows: {len(filtered)}")
print(f"Unique peptides: {filtered['Peptide'].nunique()}")

filtered.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"\nWrote {OUTPUT_FILE}")