from Bio import Align
from Bio.Align import substitution_matrices
from util import sequence

aligner = Align.PairwiseAligner(
  mode ="local",
  substitution_matrix = substitution_matrices.load("BLOSUM62"),
  open_gap_score = -10, extend_gap_score = -0.5)

score = aligner.score(sequence("NP_000549"), sequence("NP_000509"))
assert score == 293.5 , f"score: {score}"

alignments = aligner.align(sequence("NP_000549"), sequence("NP_000509"))
print(alignments[0])