from Bio import SeqIO
from Bio.Align import PairwiseAligner
import csv

#Configuration
ref_path = "sequence.fasta"        # rCRS
sample_path = "sequence-2.fasta"     #Sample mtDNA sequence
output_csv = "mutations_report.csv"

# Load sequences
ref = SeqIO.read(ref_path, "fasta")
sample = SeqIO.read(sample_path, "fasta")

# Align Sequences 
aligner = PairwiseAligner()
aligner.mode = "global" # Use global alignment
alignment = aligner.align(ref.seq, sample.seq)[0]
# Extract aligned sequences from alignment object
aligned_ref = alignment.target
aligned_sample = alignment.query

# Detect mutations
mutations = []

pos_ref = pos_sample = 0
for i, (base_r, base_s) in enumerate(zip(aligned_ref, aligned_sample)):
    if base_r != '-' and base_s != '-':
        pos_ref += 1
        pos_sample += 1
        if base_r != base_s:
            mutations.append(('SNP', pos_ref, base_r, base_s))
    elif base_r == '-' and base_s != '-':
        pos_sample += 1
        mutations.append(('insertion', pos_ref + 1, '-', base_s))
    elif base_r != '-' and base_s == '-':
        pos_ref += 1
        mutations.append(('deletion', pos_ref, base_r, '-'))

# Print results
print(f"\n🔍 Found {len(mutations)} mutations:\n")
for m in mutations[:10]: 
    print(m)


