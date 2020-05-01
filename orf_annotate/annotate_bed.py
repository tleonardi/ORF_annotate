from bedparse.bedline import bedline
from Bio import SeqIO
from Bio.Seq import Seq
from orf_annotate.orf import orf
from orf_annotate.sequence import sequence
from orf_annotate import utils

def annotate_bed(bedfile, fasta_file, stop_required=True, min_len=50, longest=False, protein_fasta=None, min_aln_score_thr=-1):
    """ Given a BED file and corresponding fasta file,
    returns an iterator of bedline objects with ORF annotations.
    Args:
        bedfile (string): path to the bedfile
        fasta_files (string): path to the fasta file
        protein_fasta (str): Fasta file of know proteins against which ORFs are aligned
        min_aln_score_thr (int): Minimum aligmement score to consider hits againts protein_fasta
    Returns:
        bedline 
    """
    seq_dict = utils.load_fasta(fasta_file)
    if(protein_fasta is not None):
        protein_seq_dict = utils.load_fasta(protein_fasta)
    with open(bedfile, "r") as bed:
        for line in bed:
            bl = bedline(line.split('\t'))
            seq = sequence(utils.get_sequence_from_seq_dict(bl, seq_dict), orfs_req_stop=stop_required, min_orf_len=min_len)
            if longest:
                new_orf = seq.longest_orf
            else:
                new_orf = seq.first_orf
            if new_orf is not None:
                orf_start = bl.tx2genome(new_orf.start, stranded=True)
                orf_end = bl.tx2genome(new_orf.stop-1, stranded=True)+1
                if(protein_fasta is not None):
                    orf_id = new_orf.find_orf_best_match(protein_seq_dict, min_score_thr=min_aln_score_thr)
                    bl.name = bl.name + "#" + orf_id
            else:
                orf_start = bl.end
                orf_end = bl.end
            bl.cdsStart = orf_start
            bl.cdsEnd = orf_end
            yield bl


