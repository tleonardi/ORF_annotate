from Bio import Align
from orf_annotate import utils

class orf(object):
    """ Class that represents and Open Reading Frame
    Attributes:
        start (int): Start codon coordinate of the ORF in the parent sequence
        stop (int): Stop codon coordinate of the ORF in the parent sequence
        length (int): Number of aminoacids in the ORF
        translation (str): Translation of the ORF
        seq (str): ORF sequence
        orf_id (str): ID of the best-match ORF from a database of sequences
    """
    def __init__(self, start, stop, length, translation=None, sequence=None):
        """ Init an object of class orf.
            If protein_seq_dict is provided, align the ORF translation against the sequences in the
            dict and set the orf_id attribute to the id of sequence with best score.
        Attributes:
            start (int): Start codon coordinate of the ORF in the parent sequence
            stop (int): Stop codon coordinate of the ORF in the parent sequence
            length (int): Number of aminoacids in the ORF
            translation (str): Translation of the ORF
            sequence (str): ORF sequence
        """
        if(type(start) is not int):
            raise Exception("ORF start is not an int")
        if(type(stop) is not int):
            raise Exception("ORF stop is not an int")
        if(start<0 or stop<0):
            raise Exception("ORF start and stop coordinates must be >=0")
        if(type(length) is not int):
            raise Exception("ORF length is not an int")
        if(length<1):
            raise Exception("ORF length is less than 1")
        if(type(length) is not int):
            raise Exception("ORF length is not an int")
        if(translation and type(translation) is not str):
            raise Exception("ORF translation is not a string")
        if(sequence and type(sequence) is not str):
            raise Exception("ORF sequence is not a string")
        
        assert (stop-start) % 3 == 0
        assert (stop-start) // 3 == length     

        self.start = start
        self.stop = stop
        self.length = length
        self.translation = translation
        self.sequence = sequence 

    def find_orf_best_match(self, protein_seq_dict=None, min_score_thr=-1):
        """  Align the ORF translation against the sequences in the protein_seq_dict
            dict and return the id of the sequence with best score.
        Attributes:
            protein_seq_dict (dict): Dictionary of id:sequence of known proteins
            min_score_thr (int): Minimum alignment score to consider a hit
        """
        if(type(protein_seq_dict) is not dict):
            raise Exception("protein_seq_dict must be a dictionary of if:sequence")
        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -10
        aligner.target_end_gap_score = -0.1
        aligner.query_end_gap_score = -0.1
        seq_to_match = self.translation
        aligner.match_score = 2
        aligner.mismatch_score = -1
        # Remove stop codon if present
        if(seq_to_match[-1] == "*"):
            seq_to_match = seq_to_match[:-1]
        best_hit = ["NA", -1]
        for id,seq in protein_seq_dict.items():
            aln = aligner.align(seq_to_match, seq.seq)
            aln = sorted(aln)[0]
            score=aln.score
            identity=utils.percent_identity(aln)
            id="{}({:.2f}%/{:.2f}%)".format(id, identity[0], identity[1])
            if score> best_hit[1]:
                best_hit[0] = id
                best_hit[1] = score
        if best_hit[1]>min_score_thr:
            return best_hit[0]
        else:
            return "Unknown"

