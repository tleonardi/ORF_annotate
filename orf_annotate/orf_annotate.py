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

def pad_sequence(seq):
    """ Pads a sequence so that its length is 
    a multiple of three

    Args:
        seq (Seq): Sequence to be padded
    """
    remainder = len(seq)%3
    if remainder: seq = seq + "N"*(3-remainder)
    return(seq)

def find_orfs_in_seq(seq, stop_required=True, min_len=50):
    """ Returns a dictionary with information for every ORF in a sequence.
    Args:
        seq (Seq): sequence
        stop_required (bool): only report ORFs that end with a stop codon
        min_len (int): minimum size (in aa) for reporting ORF
    Returns:
        dict: dictionary of dictionaries with ORF start position as key
    """
    orfs=dict()
    for start_idx in starts_in_seq(str(seq)):
        subseq = pad_sequence(seq[start_idx:])
        tx = subseq.translate(to_stop=False)
        stop_idx = tx.find("*")
        if stop_idx != -1:
            tx = tx[:stop_idx+1]
            # convert the stop idx from aa coord to nt coord
            stop_idx = start_idx + len(tx[:stop_idx+1])*3
            subseq = seq[start_idx:stop_idx]
        if len(tx)>=min_len and (stop_idx != -1 or stop_required is False):
            if stop_idx == -1:
                stop_idx = len(subseq)+1
            orfs[int(start_idx)] = {"start":int(start_idx), "len": len(tx), "tx": str(tx), "stop": stop_idx, "seq": str(subseq)}
    return(orfs)

def choose_first_orf(orfs):
    """ Returns the first ORF from a dict of ORFs
    Args:
        orfs (dict): ORFs dict indexed by (int) ORF start coordinate
    Returns:
        dict:   value of orfs at first key or None if orfs is empty
    """
    if(len(orfs)==0):
        return None
    min_key = min([k for k in orfs.keys()])
    return(orfs[min_key])

def choose_longest_orf(orfs):
    """ Returns the longest ORF from a dict of ORFs.
    In case of a tie, the first one (i.e. lowest index) is returned.
    Args:
        orfs (dict): ORFs dict of dicts indexed by (int) ORF start coordinate
                     and containing a len key.
    Returns:
        dict:   value of longest entry of orfs or None if orfs is empty
    """
    if(len(orfs)==0):
        return None
    max_length = max([v["len"] for v in orfs.values()]) 
    longest_orfs = {k:v for k,v in orfs.items() if v["len"]==max_length}
    if(len(longest_orfs)>1):
        return choose_first_orf(longest_orfs)
    else:
        return longest_orfs.popitem()[1]

def annotate_bed(bedfile, fasta_file, stop_required=True, min_len=50, longest=False):
    """ Given a BED file and corresponding fasta file,
    returns an iterator of bedline objects with ORF annotations.
    Args:
        bedfile (string): path to the bedfile
        fasta_files (string): path to the fasta file
    Returns:
        bedline 
    """
    seq_dict = load_fasta(fasta_file)
    with open(bedfile, "r") as bed:
        for line in bed:
            bl = bedline(line.split('\t'))
            seq = get_sequence_from_seq_dict(bl, seq_dict)
            orfs = find_orfs_in_seq(seq, stop_required=stop_required, min_len=min_len)
            if longest:
                new_orf = choose_longest_orf(orfs)
            else:
                new_orf = choose_first_orf(orfs)
            if new_orf is not None:
                orf_start = bl.tx2genome(new_orf["start"], stranded=True)
                orf_end = bl.tx2genome(new_orf["stop"]-1, stranded=True)+1
            else:
                orf_start = bl.end
                orf_end = bl.end
            bl.cdsStart = orf_start
            bl.cdsEnd = orf_end
            yield bl


