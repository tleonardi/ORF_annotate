from bedparse.bedline import bedline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from orf_annotate.orf import orf
from orf_annotate.sequence import sequence
from orf_annotate import utils

def annotate_bed(bedfile, fasta_file, stop_required=True, min_len=50, longest=False):
    """ Given a BED file and corresponding fasta file,
    returns an iterator of bedline objects with ORF annotations.
    Args:
        bedfile (string): path to the bedfile
        fasta_files (string): path to the fasta file
    Returns:
        bedline 
    """
    seq_dict = utils.load_fasta(fasta_file)
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
            else:
                orf_start = bl.end
                orf_end = bl.end
            bl.cdsStart = orf_start
            bl.cdsEnd = orf_end
            yield bl


