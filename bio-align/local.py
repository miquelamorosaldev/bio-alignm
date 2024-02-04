# Smith Waterman local alignment
import numpy as np
from util import print_matrix

xs = "AAUGCCAUUGACGG"
ys = "CAGCCUCGCUUAG"

gap = -4
match = 3
mismatch = -1

#####

M = np.zeros((len(xs)+1, len(ys) + 1), np.dtype("int8"))

for i in range(1, len(xs) + 1):
    for j in range(1, len(ys) + 1):
        M[i, j] = max(
            0,
            M[i-1][j-1] + (match if xs[i-1] == ys[j-1] else mismatch),
            M[i-1][j] + gap,
            M[i][j-1] + gap
        )

# print_matrix(xs, ys, M)

# Find the highest scoring cell
max_score = 0
max_index = (1, 1)

for i in range(1, len(xs) + 1):
    for j in range(1, len(ys) + 1):
        score = M[i, j]
        if score > max_score:
            max_score = score
            max_index = (i, j)


# Trace back
i, j = max_index

rxs, rys = [], []
rxs.append(xs[i - 1])
rys.append(ys[j - 1])

while True:

    v = max(M[i-1, j-1], M[i-1][j], M[i][j-1])
    if v < match:
        break

    if M[i-1, j-1] == v:
        rxs.append(xs[i - 2])
        rys.append(ys[j - 2])
        i -= 1
        j -= 1
    elif M[i-1, j] == v:
        rxs.append(xs[i - 2])
        rys.append("-")
        i -= 1
    else:
        rxs.append("-")
        rys.append(ys[j - 2])
        j -= 1

# Reverse the strings.
rxs = "".join(rxs)[::-1]
rys = "".join(rys)[::-1]

# Print
print(rxs)
print(rys)

