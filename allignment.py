from Bio import SeqIO
from Bio import pairwise2

def read_fasta(file):
    sequences = []
    with open(file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(str(record.seq))
    return sequences[0]

def count_c_and_g(sequence):
    count_c = sequence.count('C')
    count_g = sequence.count('G')
    
    return count_c, count_g

def perform_sequence_alignment(seq1, seq2, gap_penalty, match_score, mismatch_penalty, start, end, printRes=False):
    sequences = [seq1[start:end], seq2[start:end]]
    
    alignments = pairwise2.align.globalms(sequences[0], sequences[1], match_score, mismatch_penalty, gap_penalty, gap_penalty)
    
    if printRes:
        print(pairwise2.format_alignment(*alignments[0]))
        
    return alignments[0].score