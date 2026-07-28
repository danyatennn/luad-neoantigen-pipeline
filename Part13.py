#Part 13 (Immunogenicity prediction)

import pandas as pd

file_1 = pd.read_csv('/Users/naigelu/Desktop/Imperial College London/Year 2/STJU Summer Internship/Project 130/Lung Cancer/Google Drive Downloads/Part 11 + 12.1/9mer_predictions.tsv', sep="\t")
file_2 = pd.read_csv('/Users/naigelu/Desktop/Imperial College London/Year 2/STJU Summer Internship/Project 130/Lung Cancer/Google Drive Downloads/Part 11 + 12.1/classI_neoantigen_candidates.tsv', sep="\t")

mutant_9mers = file_2[
    [
        "Peptide_Mutant",
        "Peptide_WildType",
        "Allele"
    ]
].copy()

print(
    mutant_9mers["Allele"]
    .dropna()
    .unique()
)

mutant_9mers["Peptide_Mutant"] = (
    mutant_9mers["Peptide_Mutant"]
    .astype(str)
    .str.strip()
    .str.upper()
)

mutant_9mers = mutant_9mers.drop_duplicates(
    subset=["Peptide_Mutant"]
).reset_index(drop=True)

#print(mutant_9mers[["Peptide_Mutant"]].head())
#print("Number of unique mutant 9-mers:", len(mutant_9mers))


with open("mutant_9mers.fasta", "w") as fasta_file:
    for index, peptide in enumerate(
        mutant_9mers["Peptide_Mutant"],
        start=1
    ):
        fasta_file.write(f">mutant_9mer_{index}\n")
        fasta_file.write(f"{peptide}\n")

print("FASTA file created: mutant_9mers.fasta")

with open("mutant_9mers.fasta", "r") as fasta_file:
    for _ in range(10):
        line = fasta_file.readline()

        if not line:
            break

        #print(line.strip())

#for PRIME, the lower the % rank, the better the immunogenicity. For the report:

#Immunogenicity score = PRIME %Rank
#Score direction = Lower is better
#Raw PRIME score: higher is better
#PRIME %Rank: lower is better

#PRediciton of IMmunogenic Epitopes (PRIME2.1) was used.

prime = pd.read_csv(
    '/Users/naigelu/Desktop/Imperial College London/Year 2/STJU Summer Internship/Project 130/Lung Cancer/project130_lung_cancer/mutant_9mers_res (1).txt',
    sep=r"\s+",
    comment="#"
)

file_2["PRIMEAllele"] = (
    file_2["Allele"]
    .astype("string")
    .str.replace("HLA-", "", regex=False)
    .str.replace("*", "", regex=False)
    .str.replace(":", "", regex=False)
)

print(
    file_2[
        ["Allele", "PRIMEAllele"]
    ].drop_duplicates()
)

combined = file_2.merge(
    prime,
    left_on="Peptide_Mutant",
    right_on="Peptide",
    how="left",
    validate="many_to_one"
)



def get_prime_rank(row):
    allele = row["PRIMEAllele"]

    if pd.isna(allele):
        return pd.NA

    column_name = f"%Rank_{allele}"

    if column_name not in row.index:
        return pd.NA

    return row[column_name]


combined["ImmunogenicityScore"] = combined.apply(
    get_prime_rank,
    axis=1
)

def get_prime_raw_score(row):
    allele = row["PRIMEAllele"]

    if pd.isna(allele):
        return pd.NA

    column_name = f"Score_{allele}"

    if column_name not in row.index:
        return pd.NA

    return row[column_name]


combined["PRIMERawScore"] = combined.apply(
    get_prime_raw_score,
    axis=1
)

combined["PredictedImmunogenic"] = (
    pd.to_numeric(
        combined["ImmunogenicityScore"],
        errors="coerce"
    ) <= 0.5
)


print(combined.head())
print(combined.columns.tolist())

part13_columns = [
    "GeneName",
    "ProteinChange",
    "Peptide_Mutant",
    "Allele",
    "ImmunogenicityScore",
    "PredictedImmunogenic",
]

part13 = combined[part13_columns].copy()

print(part13.head())
print(part13.columns.tolist())


part13.to_csv(
    "/Users/naigelu/Desktop/Imperial College London/Year 2/STJU Summer Internship/Project 130/Lung Cancer/Part 13/part13_immunogenicity_results.tsv",
    sep="\t",
    index=False,
    na_rep="NA"
)



