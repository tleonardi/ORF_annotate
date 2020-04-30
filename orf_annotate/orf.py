class orf(object):
    """ Class that represents and Open Reading Frame
    Attributes:
        start (int): Start codon coordinate of the ORF in the parent sequence
        stop (int): Stop codon coordinate of the ORF in the parent sequence
        length (int): Number of aminoacids in the ORF
        translation (str): Translation of the ORF
        seq (str): ORF sequence
    """
    def __init__(self, start, stop, length, translation=None, sequence=None):
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

   def match_orf_to_db(self, fastadb, min_score_thr=-1):
        aligner = Align.PairwiseAligner()
        seq_to_match = self.translation
        seq_dict = SeqIO.to_dict(SeqIO.parse(fastadb, "fasta"))
        best_hit = ("NA", -1)
        for id,seq in seq_dict.items():
            score = aligner.score(seq_to_match, seq)
            if score> best_hit[1]:
                best_hit[0] = id
                best_hit[1] = score
        if best_hist[1]>min_score_thr:
            self.orf_id = best_hit[0]

