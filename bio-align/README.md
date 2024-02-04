# bio-alignment

## Getting started

Install virtual environment

```sh
$ sudo apt install python3-venv
```

Activate `venv` and install dependencies

```sh
$ python3 -m venv .venv
$ source .venv/bin/activate
(.venv) $ python3 -m pip install --upgrade pip
(.venv) $ pip install -r requirements-dev.txt
```

Or execute init.sh

```sh
chmod +x init.sh
./init.sh
```

## Code

### Step 1. Understand Levenshtein distance algorithm.

distance.py

util.py 

The Levenshtein distance is a string metric for measuring difference between two sequences. Informally, the Levenshtein distance between two words is the minimum number of single-character edits (i.e. insertions, deletions or substitutions) required to change one word into the other.

For that reason, is very useful to compare DNA, RNA or aminoacid sequences.

<em>For example, the Levenshtein distance between “FLOMAX” and “VOLMAX” is 3, since the following three edits change one into the other, and there is no way to do it with fewer than three edits:</em>

### Step 2. Test Global and Local Pairwise alignments in Python

global.py

local.py

A global alignment is defined as the end-to-end alignment of two strings s and t. 
A local alignment of string s and t is an alignment of substrings of s with substrings of t. alignments because we normally do not know the boundaries of genes and only a small domain of the gene may be conserved.

### Step 3. Test BioPython PairwiseAligner from real sequences.

util.py

Code to obtain sequences from NCBI-Entrez in .fasta format.

aligner.py

matrix.py

Alingment examples.

Scoring matrix from the alignment of the 16S ribosomal RNA gene sequences of Escherichia coli [NR_024570.1] and Bacillus subtilis [NR_112116.2].
