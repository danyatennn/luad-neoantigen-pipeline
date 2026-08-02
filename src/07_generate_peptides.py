import gzip
import re
import pandas as pd

import config

# 3-letter -> 1-letter amino acid codes.
AA3 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
    "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
    "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V",
}

# Load Ensembl human proteome (release 116, matches VEP cache).
# Index by unversioned ENSP (e.g. ENSP00000493376).
proteins = {}
with gzip.open(config.PROTEOME, "rt") as fh:
    pid, seq = None, []
    for line in fh:
        if line.startswith(">"):
            if pid:
                proteins[pid] = "".join(seq)
            pid = line[1:].split()[0].split(".")[0]
            seq = []
        else:
            seq.append(line.strip())
    if pid:
        proteins[pid] = "".join(seq)

variants = pd.read_csv(config.VARIANT_ANNOTATION, sep="\t")

# Parse "ENSPxxxxx.N:p.Ile239Met" -> (ref_aa, pos, alt_aa).
change_re = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")

def windows(seq, m, k):
    # All k-mer windows [start..start+k-1] containing 1-based position m, within the protein.
    for start in range(max(1, m - k + 1), min(len(seq) - k + 1, m) + 1):
        yield start, seq[start - 1 : start - 1 + k]

rows = []
skipped_bad_hgvsp = skipped_no_seq = skipped_bad_ref = 0

for _, v in variants.iterrows():
    match = change_re.search(str(v["Protein_Change"]))
    if not match:
        skipped_bad_hgvsp += 1
        continue
    ref3, pos, alt3 = match.group(1), int(match.group(2)), match.group(3)
    ref, alt = AA3.get(ref3), AA3.get(alt3)
    # Strip a possible version suffix from ENSP (defensive; our column has none).
    protein_id = str(v["Protein_ID"]).split(".")[0]
    seq = proteins.get(protein_id)
    if not seq:
        skipped_no_seq += 1
        continue
    # Verify the reference amino acid agrees with the protein sequence.
    if pos > len(seq) or seq[pos - 1] != ref:
        skipped_bad_ref += 1
        continue
    for k in (9, 15):
        for start, wt_pep in windows(seq, pos, k):
            mut_pos = pos - start + 1  # 1-based position of mutation inside the peptide
            mut_pep = wt_pep[: mut_pos - 1] + alt + wt_pep[mut_pos:]
            # Nail down all four task requirements.
            assert len(wt_pep) == k and len(mut_pep) == k
            assert wt_pep[mut_pos - 1] == ref
            assert mut_pep[mut_pos - 1] == alt
            assert sum(a != b for a, b in zip(wt_pep, mut_pep)) == 1
            common = (v["Genomic_Variant"], v["Gene_Symbol"], v["Transcript_ID"],
                      protein_id, match.group(0), pos, k, mut_pos)
            rows.append((*common, wt_pep, "WildType"))
            rows.append((*common, mut_pep, "Mutant"))

peptides = pd.DataFrame(rows, columns=[
    "GenomicVariant", "GeneName", "TranscriptID", "ProteinID",
    "ProteinChange", "ProteinPosition", "Length", "MutPos", "Peptide", "Type",
])
peptides.to_csv(config.PEPTIDES, sep="\t", index=False)

print(f"Variants processed: {len(variants)}")
print(f"  skipped (unparseable HGVSp):    {skipped_bad_hgvsp}")
print(f"  skipped (protein not in FASTA): {skipped_no_seq}")
print(f"  skipped (ref aa mismatch):      {skipped_bad_ref}")
print(f"Peptides written: {len(peptides)}")
print(peptides.head(4))
