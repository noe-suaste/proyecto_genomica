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

def needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty, start, end):
    len1 = len(seq1[start:end])
    len2 = len(seq2[start:end])
    
    seq1_aux, seq2_aux = seq1[start:end], seq2[start:end]

    score_matrix = [[0] * (len2 + 1) for _ in range(len1 + 1)]

    traceback_matrix = [[0] * (len2 + 1) for _ in range(len1 + 1)]

    for i in range(1, len1 + 1):
        score_matrix[i][0] = gap_penalty * i
        traceback_matrix[i][0] = 'U' 

    for j in range(1, len2 + 1):
        score_matrix[0][j] = gap_penalty * j
        traceback_matrix[0][j] = 'L' 

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if seq1_aux[i - 1] == seq2_aux[j - 1]:
                match = match_score
            else:
                match = mismatch_penalty

            diagonal_score = score_matrix[i - 1][j - 1] + match
            up_score = score_matrix[i - 1][j] + gap_penalty
            left_score = score_matrix[i][j - 1] + gap_penalty

            max_score = max(diagonal_score, up_score, left_score)
            score_matrix[i][j] = max_score

            if max_score == diagonal_score:
                traceback_matrix[i][j] = 'D'  
            elif max_score == up_score:
                traceback_matrix[i][j] = 'U'  
            else:
                traceback_matrix[i][j] = 'L'  

    aligned_seq1 = ""
    aligned_seq2 = ""
    i = len1
    j = len2
    score = 0

    while i > 0 or j > 0:
        if traceback_matrix[i][j] == 'D':
            aligned_seq1 = seq1_aux[i - 1] + aligned_seq1
            aligned_seq2 = seq1_aux[j - 1] + aligned_seq2
            if seq1_aux[i - 1] == seq2_aux[j - 1]:
                score += match_score
            else:
                score += mismatch_penalty
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 'U':
            aligned_seq1 = seq1_aux[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            score += gap_penalty
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2_aux[j - 1] + aligned_seq2
            score += gap_penalty
            j -= 1

    return score