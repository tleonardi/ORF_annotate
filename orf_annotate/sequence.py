from Bio.Seq import Seq
from orf_annotate.orf import orf
from orf_annotate import utils

class sequence(object):
    """ Class that represents a sequence
    Attributes:
        sequence (Bio.Seq): Sequence
        orfs (list of orf): list of orf objects
    """
    def __init__(self, sequence, orfs_req_stop=True, min_orf_len=50):
        """ Initialise seq instance
        Args:
            sequence (Bio.Seq): Sequence
            orfs_req_stop (bool): only report ORFs that end with a stop codon
            min_orf_len (int): minimum size (in aa) for reporting ORF
        """
        if(not isinstance(sequence, Seq)):
            raise Exception("Sequence is not an instance of Bio.Seq")
        
        self.sequence = sequence
        self.__find_orfs_in_seq(stop_required=orfs_req_stop, min_len=min_orf_len)
        if(len(self.orfs)==0):
            self.has_orfs = False
            self.longest_orf_len = 0
        else:
            self.has_orfs = True
            self.longest_orf_len = max([orf.length for orf in self.orfs])


        
    def __find_orfs_in_seq(self, stop_required=True, min_len=50):
        """ Returns a dictionary with information for every ORF in a sequence.
        Args:
            seq (Seq): sequence
            stop_required (bool): only report ORFs that end with a stop codon
            min_len (int): minimum size (in aa) for reporting ORF
        Returns:
            list: list of orf objects ordered by start position 
        """
        self.orfs=list()
        seq=self.sequence
        for start_idx in utils.starts_in_seq(str(seq)):
            subseq = utils.trim_sequence(seq[start_idx:])
            tx = subseq.translate(to_stop=False)
            stop_idx = tx.find("*")
            if stop_idx != -1:
                tx = tx[:stop_idx+1]
                # convert the stop idx from aa coord to nt coord
                stop_idx = start_idx + len(tx[:stop_idx+1])*3
                #subseq = seq[start_idx:stop_idx]
            if stop_idx == -1 and stop_required is False:
                # We can't user len(seq) because it wouldn't take into account
                # the trimming applied to subseq
                stop_idx = start_idx+len(subseq)
            if len(tx)>=min_len and stop_idx != -1:
                self.orfs.append(orf(start=int(start_idx), stop=int(stop_idx), length=len(tx), translation=str(tx), sequence=str(subseq)))
        
     
    @property
    def first_orf(self):
        """ Getter of the first ORFi
        Args:
            orfs (dict): ORFs dict indexed by (int) ORF start coordinate
        Returns:
            orf:  orf object representing first ORF in seq or None if there are no ORFs 
        """
        if(not self.has_orfs):
            return None
        else:
            return(self.orfs[0])
    
    @property
    def longest_orf(orfs):
        """ Getter of the longest ORF
        In case of a tie, the first one is returned.
        Args:
            orfs (dict): ORFs dict of dicts indexed by (int) ORF start coordinate
                         and containing a len key.
        Returns:
            orf:  orf object representing longest ORF in seq or None if there are no ORFs 
        """
        if(not self.has_orfs):
            return None
        else:
            longest_orfs = [orf for orf in self.orfs if orf.length==self.longest_orf_len]
            return longest_orfs[0]
