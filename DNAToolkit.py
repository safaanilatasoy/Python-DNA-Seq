# DNA Tolkit file
import collections
from structures import *


# Validating sequence
def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq


def countNucFrequency(seq):
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T":0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict
# return dict(collections.Counter(seq))


def transcription(seq):
    # DNA -> RNA Transctiption
    return seq.replace("T","U")


#DNA Reverse Compliment
def reverse_compliment(seq):
    return ''.join([DNA_ReverseCompliment[nuc] for nuc in seq])[::-1]

# GC Content in a DNA-RNA sequence
def gc_content(seq):
    return round((seq.count('C') + seq.count('G')) / len(seq)*100)

# GC Content calculations with K WİNDOW SİZE ON dna
def gc_content_subsec(seq, k=20):
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        res.append(gc_content(subseq))
    return res 