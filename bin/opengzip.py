import gzip

with gzip.open(r"C:\Users\lyssa\Documents\GitHub\project130_lung_cancer\cohortMAF.2026-07-20.maf.gz","rt") as f:
    for line in f:
        if not line.startswith("#"):
            print(line.strip())