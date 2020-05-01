import sys
import argparse
from pkg_resources import get_distribution
from orf_annotate.annotate_bed import annotate_bed

__version__ = get_distribution('orf_annotate').version

def main(args=None):
    if args is None:
        args = sys.argv[1:]
    parser = argparse.ArgumentParser(
            description="""Annotate Open Reading Frames in a BED file.""") 
    parser.add_argument('--version', '-v', action='version', version='v'+__version__)
    parser.add_argument('--bedfile', '-b', type=str, 
                        help="Path to the BED file.", required=True)
    parser.add_argument('--fasta', '-f', type=str, 
                        help="Path to a fasta file with the spliced, stranded sequence of the BED file entries.", required=True)
    parser.add_argument('--min-len', '-m', type=int, 
                        help="Minimum length of ORF (in aa) for considering it", default=50)
    parser.add_argument('--no-stop', '-n', action="store_true", 
                        help="If set, also consider ORFs that don't terminate with a canonical stop codon")
    parser.add_argument('--longest', '-l', action="store_true", 
                        help="""By default, the first ORF is annotated. If this option is used, we
                                annotate the longest intead. In case of ties the first wins.""")
    prot_fasta = parser.add_argument_group("Protein DB")
    prot_fasta.add_argument('--proteins-fasta', type=str, 
                            help="""Path to a fasta file of protein sequences. ORFs identied will be
                                    aligned against records in this file, the the ID of the best hit 
                                    will be appended to the name field of the BED record.""")
    prot_fasta.add_argument('--min-aln-score', type=int, 
                            help="Minimum alignment score to consider a protein hit valid", default=-1)
    args = parser.parse_args()

    annotator = annotate_bed(bedfile=args.bedfile, fasta_file=args.fasta, stop_required=(not args.no_stop), min_len=args.min_len, longest=args.longest, protein_fasta=args.proteins_fasta, min_aln_score_thr=args.min_aln_score)
    for bl in annotator:
        bl.print()


