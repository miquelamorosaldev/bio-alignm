# Needleman-Wunsch Algorithm
import numpy as np
from util import print_matrix

# xs = "GCGTAACACGTGCG"
# ys = "ACAACCCGTGCGAC"

xs = "GCGTAACACGTGCGGATGATAGATGAATCGCTCAGCATACGCTCCAGTAACGTACGACATCATCAGACT"
ys = "ACAACCCGTGCGACCACTACGACTTAATCGCTCAGTACTACGTCAGTCATTGCAAGACAGTACTACGTT"

gap = -1
match = 1
mismatch = -2

#####

initial_matches = sum(1 for x, y in zip(xs, ys) if x == y)

M = np.zeros((len(xs) + 1, len(ys) + 1), np.dtype("int8"))
M[:, 0] = np.linspace(0, len(xs) * gap, len(xs) + 1)
M[0, :] = np.linspace(0, len(ys) * gap, len(ys) + 1)

for i in range(1, len(xs)+1):
    for j in range(1, len(ys)+1):
        M[i, j] = max(
            M[i - 1][j - 1] + (match if xs[i-1] == ys[j-1] else mismatch),
            M[i-1][j] + gap,
            M[i][j-1] + gap)

#print_matrix(xs,ys,M)

# Trace back
i, j = len(xs), len(ys)
rxs, rys = [], []

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
        i -= 1
    else:
        rxs.append("-")
        rys.append(ys[j - 1])
        j -= 1



# Reverse the strings.
rxs = "".join(rxs)[::-1]
rys = "".join(rys)[::-1]

print("\n".join([rxs, rys]))