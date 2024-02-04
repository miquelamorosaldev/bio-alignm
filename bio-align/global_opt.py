# Needleman-Wunsch Algorithm
import numpy        as np
from pathlib        import Path
from util           import print_matrix
from Bio            import SeqIO
from Bio.SeqRecord  import SeqRecord

gap = -1
match = 1
mismatch = -2

#####
# xs = "GCGTAACACGTGCGGATGATAGATGAATCGCTCAGCATACGCTCCAGTAACGTACGACATCATCAGACT"
# ys = "ACAACCCGTGCGACCACTACGACTTAATCGCTCAGTACTACGTCAGTCATTGCAAGACAGTACTACGTT"

xs = "GCGTAACACGTGCGGATGATAGATGAATCGCTCAGCATACGCTCCAGTAACGTACGACATCATCAGACT"
ys = "ACAACCCGTGCGACCACTACGACTTAATCGCTCAGTACTACGTCAGTCATTGCAAGACAGTACTACGTT"

#####
# Agafem un parell de cadenes de tamany, prèviament descarregades de l'NCBI.
current_dir: Path = Path(__file__).parent
print(current_dir)

# Seqüència de la proteïna humana "Hemoglobina beta" (HBB):
# HBB a NCBI

# Seqüència de la proteïna de ratolí "Hemoglobina beta" (Hbb-b1):
# Hbb-b1 a NCBI
# Aquestes seqüències tenen longitud similar i representen les hemoglobines beta en humans i ratolins, respectivament. 

hbbb1: Path = current_dir/'biopython'/'data'/'hbb-b1-protein.faa'
HBB:   Path = current_dir/'biopython'/'data'/'HBB-protein.faa'
hbbb1_record: SeqRecord = SeqIO.read(hbbb1, 'fasta')
HBB_record:   SeqRecord = SeqIO.read(HBB,   'fasta')
xs = hbbb1_record.seq
ys = HBB_record.seq

print(xs)

# xs = hbb_mus
# ys = hbb_hom


M = np.zeros((len(xs) + 1, len(ys) + 1), np.dtype("int8"))
M[:, 0] = np.linspace(0, len(xs) * gap, len(xs) + 1)
M[0, :] = np.linspace(0, len(ys) * gap, len(ys) + 1)

max_value = float('-inf')  # Variable per optimitzar

for i in range(1, len(xs)+1):
    for j in range(1, len(ys)+1):
        M[i, j] = max(
            M[i - 1][j - 1] + (match if xs[i-1] == ys[j-1] else mismatch),
            M[i-1][j] + gap,
            M[i][j-1] + gap)
        max_value: int = max(max_value, M[i, j])

print_matrix(xs,ys,M)

alignment_symbols = []  # Llista per acumular els símbols "|"

# Trace back
i, j = len(xs), len(ys)
rxs, rys = [], []
initial_matches = 0 
gap_count = 0 

while i > 0 or j > 0:
    v = max(M[i-1, j-1], M[i-1][j], M[i][j-1])

    if M[i-1, j-1] == v:
        rxs.append(xs[i - 1])
        rys.append(ys[j - 1])
        i -= 1
        j -= 1
    elif M[i-1, j] == v:
        rxs.append(xs[i - 1])
        rys.append("-")
        gap_count += 1
        i -= 1
    else:
        rxs.append("-")
        rys.append(ys[j - 1])
        gap_count += 1
        j -= 1

# Construeix la llista d'alignment_symbols a partir de les seqüències rxs i rys
algs: list[str] = ["|" if rxs[k] == rys[k] and rxs[k] != "-" else " " for k in range(len(rxs))]

# Reverse the strings.
rxs: list[str] = "".join(rxs)[::-1]
rys: list[str] = "".join(rys)[::-1]
algs: list[str] = "".join(algs)[::-1]

print("\n".join([rxs, algs, rys]))

# Calculate percentages
total_symbols = len(algs)
gap_count = algs.count(" ")
final_matches = algs[:min(len(rxs), len(rys))].count("|")

initial_percentage = (initial_matches / total_symbols) * 100
final_percentage = (final_matches / total_symbols) * 100
gap_percentage = (gap_count / total_symbols) * 100

print("Initial Percentage of Matches:", round(initial_percentage,4), "%")
print("Final Percentage of Matches:", round(final_percentage,4), "%")
print("Percentage of Gaps Introduced:", round(gap_percentage,4), "%")
