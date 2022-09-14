from Bio import SeqIO
import stumpy
import numpy as np
from copy import copy

CODES = {'R':1,'H':2,'K':3,
              'D':4,'E':5,
              'S':6,'T':7,'N':8,'Q':9,
              'C':10,'G':11,'P':12,
              'A':13,'V':14,'I':15,'L':16,'M':17,'F':18,'Y':19,'W':20,
              'X':21}

def get_sequences_from_fasta(input_file):
    seqs = SeqIO.parse(open(input_file),'fasta')
    return [s.seq for s in seqs]

def proteins_to_series(seqs):
    return [[float(CODES[a]) for a in s ] for s in seqs]

def series_to_proteins(seqs):
    inv_codes = {v: k for k, v in CODES.items()}
    return [[inv_codes[a] for a in s ] for s in seqs]

def get_n_motifs(seqs, m ,n):
    
    seqs = copy(seqs)
    motifs = []
    pos_in_seq = []
    nans = [float('nan') for i in range(m)]

    for j in range(n): 
        # get consensus
        radius , Ts_idx, subseq_idx = stumpy.ostinato(copy(seqs), m)
        seed_motif = seqs[Ts_idx][subseq_idx : subseq_idx + m]

        motifs.append(copy(seed_motif))

        # remove consensus
        position = []
        for i, e in enumerate(seqs):
            motif_idx = np.argmin(stumpy.core.mass(seed_motif, e))
            seqs[i][motif_idx:(motif_idx+m)] = nans
            position.append(motif_idx)

        pos_in_seq.append(position)

    return motifs , pos_in_seq
        

if __name__ == '__main__':
    seqs = get_sequences_from_fasta('samplefiles/test_small.fasta')
    nseqs = proteins_to_series(seqs)
    nseqs = nseqs[12:18]
    
    m = 4
    n = 3
    motifs, pos_in_seq = get_n_motifs(nseqs, m , n)

    print(motifs)
    print(series_to_proteins(motifs))

    print(pos_in_seq)



