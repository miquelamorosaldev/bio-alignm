# Levenshtein Distance Algorithm
import numpy as np
from util import print_matrix

xs = "CGA"
ys = "AGAT"

M = np.zeros((len(xs) + 1, len(ys) + 1), dtype=np.dtype("int8"))
M[:, 0] = np.linspace(0, len(xs), len(xs)+1)
M[0, :] = np.linspace(0, len(ys), len(ys)+1)

for i in range(1, len(xs) + 1):
    for j in range(1, len(ys) + 1):
        M[i][j] = min(M[i-1][j-1], M[i-1][j], M[i][j-1])
        if (xs[i-1] != ys[j-1]):
            M[i][j] += 1

# Print
print_matrix(xs, ys, M)
print(f"\nDistance: {M[len(xs)][len(ys)]}")
