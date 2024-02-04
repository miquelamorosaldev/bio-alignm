import re
import json
from   pathlib import Path

from Bio            import SeqIO
from Bio.SeqRecord  import SeqRecord
from Bio.SeqFeature import SeqFeature
# Observació: Hi havia una classe anomenada PairwiseAlignment que ja no la 
# reconeix la documentació de Biopython; però no és imprescindible pel funcionament
# del codi.
from Bio.Align      import PairwiseAligner, PairwiseAlignments, substitution_matrices

from Bio.Align.substitution_matrices import Array


# ---------------------------------------------------------------------
def align_2seq(reference_prot: str, variant_prot: str) -> PairwiseAlignments:

    # Create Global Aligner
    aligner: PairwiseAligner = PairwiseAligner()

    # Get BLOSUM62 matrix
    blosum62_matrix: Array = substitution_matrices.load('BLOSUM62')

    # Put the matrix in the aligner.
    aligner.substitution_matrix = blosum62_matrix
    
    # Get first global alignment
    alignments: PairwiseAlignments = aligner.align(reference_prot, variant_prot)
    alignment:  PairwiseAlignments  = alignments[0]

    return alignment


# - Read the same two files as in Q2:
#   - SARS-CoV-2 reference: NC_045512.2.genbank
#   - SARS-CoV-2 variant:   OL466363.1.genbank
# - Return a dict with the following contents:
#   - keys:   protein identifier
#   - values: protein alignment (reference vs variant)
# - Align proteins from the genes in Q2 whose length is less than 80 AA.
# ---------------------------------------------------------------------
def align_sarscov_var(reference_file: Path, variant_file: Path) -> dict[str, PairwiseAlignments]:

    # Read files
    reference_record: SeqRecord = SeqIO.read(reference_file, 'genbank')
    variant_record:   SeqRecord = SeqIO.read(variant_file,   'genbank')

    # Filter features by type
    ref_cds_feat_list: list[SeqFeature] = [ feature for feature in reference_record.features
                                            if (feature.type == 'CDS')
                                            and (len(feature.qualifiers['translation'][0]) < 80)]

    var_cds_feat_list: list[SeqFeature] = [ feature for feature in variant_record.features
                                            if feature.type == 'CDS'
                                            and (len(feature.qualifiers['translation'][0]) < 80)]

    # Extract protein names
    prot_names: list[str] = [   feature.qualifiers['protein_id'][0]
                                for feature
                                in  ref_cds_feat_list ]

    # Extract proteins
    ref_prot_list: list[str] = [feature.qualifiers['translation'][0]
                                for feature
                                in  ref_cds_feat_list ]

    var_prot_list: list[str] = [feature.qualifiers['translation'][0]
                                for feature
                                in var_cds_feat_list ]

    # Creat dict
    result: dict[str, PairwiseAlignments] = {prot_name: align_2seq(ref_prot, var_prot)
                                            for prot_name, ref_prot, var_prot
                                            in  zip(prot_names, ref_prot_list, var_prot_list)
    }

    return result


# Main
# ---------------------------------------------------------------------
def main():

    # Constants
    current_dir: Path = Path(__file__).parent

    reference_file: Path = current_dir/'data'/'NC_045512.2.genbank'
    variant_file:   Path = current_dir/'data'/'OL466363.1.genbank'

    q3_result: dict[str, PairwiseAlignments] = align_sarscov_var(reference_file, variant_file)

    alignments:     str  = ''.join( f"{prot_name}:\n{alignment}\n\n" for prot_name, alignment in q3_result.items() )
    q3_output_file: Path = current_dir/'answers'/'q3.alignments'
    q3_output_file.write_text(alignments)


# ---------------------------------------------------------------------
this_module: str = __name__
main_module: str = "__main__"

if (this_module == main_module): main()
# ---------------------------------------------------------------------
