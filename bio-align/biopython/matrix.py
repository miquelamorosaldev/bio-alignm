from Bio.Align import PairwiseAligner
import numpy as np
from util import sequence

aligner = PairwiseAligner(mode="local", scoring="blastn")
alignments = aligner.align(sequence("NR_024570.1"), sequence("NR_112116.2"))

substitutions = alignments[0].substitutions
observed_frequencies = substitutions / substitutions.sum()
observed_frequencies = (observed_frequencies + observed_frequencies.transpose()) / 2.0

#print(format(observed_frequencies,"%.4f"))

background = observed_frequencies.sum(0)
expected_frequencies = background[:, None].dot(background[None, :])

scoring_matrix = observed_frequencies / expected_frequencies
scoring_matrix[scoring_matrix == 0] = 1
scoring_matrix = np.log2(scoring_matrix)

aligner = PairwiseAligner(substitution_matrix = scoring_matrix, gap_score = -3)
alignments = aligner.align("ANCCG","TCCGA")
print(alignments[0])


