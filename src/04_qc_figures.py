"""
Basic data quality control figures
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import config


# Mutation heatmap and most frequently mutated genes
result = pd.read_csv(config.MUTATION_MATRIX, sep="\t")
sample_columns = result.columns[3:]

print("Top 10 most frequent mutations:")
print(
    result.assign(freq=result[sample_columns].sum(axis=1))
          .nlargest(10, "freq")[["Gene_Name", "Mutation", "AminoAcid_Change", "freq"]]
)

top_mutations = (
    result
    .assign(Mutation_Frequency=result[sample_columns].sum(axis=1))
    .sort_values("Mutation_Frequency", ascending=False)
    .head(30)
)

plt.imshow(top_mutations[sample_columns], aspect="auto", interpolation="nearest")
plt.yticks(
    range(len(top_mutations)),
    top_mutations["Gene_Name"] + " " + top_mutations["AminoAcid_Change"]
)
plt.xlabel("Tumour samples")
plt.ylabel("Mutation")
plt.title("Thirty Most Frequent LUAD Mutations")
plt.tight_layout()
# savefig must come before show(): show() clears the current figure, so saving
# afterwards writes an empty canvas.
plt.savefig(config.FIG_MUTATION_HEATMAP, dpi=150)
plt.show()


gene_by_sample = (
    result
    .groupby("Gene_Name")[sample_columns]
    .max()
)

# Count the number of tumour samples mutated for each gene
gene_frequency = (
    gene_by_sample
    .sum(axis=1)
    .sort_values(ascending=False)
)

# Select the ten most frequently mutated genes
top_10_genes = gene_frequency.head(10)

print(top_10_genes)

# Create the bar plot
plt.figure(figsize=(9, 6))
top_10_genes.sort_values().plot(kind="barh")

plt.xlabel("Number of tumour samples")
plt.ylabel("Gene")
plt.title("Ten Most Frequently Mutated Genes in TCGA-LUAD")
plt.tight_layout()
plt.savefig(config.FIG_TOP_MUTATED_GENES, dpi=150)
plt.show()


# Mutation frequency vs gene expression
df = pd.read_csv(config.INTEGRATED_MATRIX, sep="\t")
sample_cols = df.columns[4:]
gene_any_mut = df.groupby("Gene_Name")[sample_cols].max()
mutated_samples = gene_any_mut.sum(axis=1)
mutation_frequency = mutated_samples / len(sample_cols)
expression = df.groupby("Gene_Name")["GeneLevelTPM"].first()
summary = pd.DataFrame({
    "mutation_frequency": mutation_frequency,
    "expression_tpm": expression,
})
summary["log_expression"] = np.log1p(summary["expression_tpm"])
fig, ax = plt.subplots(figsize=(9, 6))
ax.scatter(summary["mutation_frequency"], summary["log_expression"],
           s=8, alpha=0.3, color="steelblue")
top = summary.nlargest(15, "mutation_frequency")
ax.scatter(top["mutation_frequency"], top["log_expression"],
           s=30, color="crimson", label="top 15 most mutated")
for gene, row in top.iterrows():
    ax.annotate(gene, (row["mutation_frequency"], row["log_expression"]),
                fontsize=8, xytext=(3, 3), textcoords="offset points")
ax.set_xlabel("Mutation frequency (fraction of samples with a mutation)")
ax.set_ylabel("Gene expression, log(TPM + 1)")
ax.set_title("Mutation frequency vs. gene expression")
ax.legend()
fig.tight_layout()
fig.savefig(config.FIG_MUTATION_VS_EXPRESSION, dpi=150)
print("Saved", config.FIG_MUTATION_VS_EXPRESSION)
print("Top mutated genes:")
print(top)
