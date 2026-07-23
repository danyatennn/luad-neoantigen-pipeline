import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("integrated_mutation_expression_matrix.tsv", sep="\t")
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
fig.savefig("mutation_vs_expression.png", dpi=150)
print("Saved mutation_vs_expression.png")
print("Top mutated genes:")
print(top)
