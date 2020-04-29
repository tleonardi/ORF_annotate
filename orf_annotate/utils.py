from bedparse.bedline import bedline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def load_fasta(fasta):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    return(seq_dict)

def get_sequence_from_seq_dict(bedline, seq_dict, bedtools_format=True):
    """ Extract the sequence of a specific bedline object from a seq_dict.
    The seq dict must contain RNA sequences, i.e. stranded and spliced.
    Args:
        bedline (bedline): bedline for which we have to get the seqeunce
        seq_dict (dict): dictionary of SeqIO objects
        bedtools_format (bool): sequence names are produced by bedtools -nameOnly
    Returns:
        Seq: sequence object
    """
    if bedtools_format:
        name = bedline.name + "(" + bedline.strand + ")"
    else:
        name = bedline.name
    try:
        seq = seq_dict[name].seq
    except KeyError as e:
        raise Exception('The record for '+name+' was not found in the fasta file') from e 

    return(seq)

def starts_in_seq(seq):
    """ Find the indexes of ATG codons in a sequence.
    Args:
        seq (Seq): sequence to scan
    Return:
        iterator: indexes of ATG occurences in sequence
    """
    start_codon="ATG"
    i = seq.find(start_codon)
    if(i == -1): return None
    while i != -1:
        yield i
        i = seq.find(start_codon, i+3)

def trim_sequence(seq):
    """ Trims the end of a sequence so that its length is 
    a multiple of three

    Args:
        seq (Seq): Sequence to be trimmed
    """
    remainder = len(seq)%3
    if remainder != 0:
        seq = seq[:-remainder]
    assert len(seq)%3 == 0
    return(seq)
