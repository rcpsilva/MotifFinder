from Bio import SeqIO

def get_sequences_from_fasta(input_file):
    seqs = SeqIO.parse(open(input_file),'fasta')
    return [s.seq for s in seqs]

def protein_to_series(seqs):
    codes = {'R':1,'H':2,'K':3,
              'D':4,'E':5,
              'S':6,'T':7,'N':8,'Q':9,
              'C':10,'G':11,'P':12,
              'A':13,'V':14,'I':15,'L':16,'M':17,'F':18,'Y':19,'W':20}
    
    return [[codes[a] for a in s ] for s in seqs]

if __name__ == '__main__':
    seqs = get_sequences_from_fasta('test.fasta')
    seqs_num = protein_to_series(seqs)
    print(seqs[0][1])
    print(seqs_num[0][1])

