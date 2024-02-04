# Needleman-Wunsch Algorithm
import numpy as np
from pathlib import Path
from util import print_matrix
from Bio            import SeqIO
from Bio.SeqRecord  import SeqRecord

#####
# Agafem un parell de cadenes de bon tamany, prèviament descarregades de l'NCBI.
current_dir: Path = Path(__file__).parent

reference_file: Path = current_dir/'data'/'NC_045512.2.genbank'
variant_file:   Path = current_dir/'data'/'OL466363.1.genbank'

reference_record: SeqRecord = SeqIO.read(reference_file, 'genbank')
variant_record:   SeqRecord = SeqIO.read(variant_file,   'genbank')

print(reference_record.seq)
print(variant_record.seq)