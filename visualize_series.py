import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import stumpy

def visualize_sequences(nseqs):
    fig, ax = plt.subplots(len(nseqs), sharex=True, sharey=True)
    for i in range(len(nseqs)):
        ax[i].plot(nseqs[i])
    plt.show()

def visualize_consensus_motifs_in_sequence(nseqs, m, Ts_idx, subseq_idx):
    Ts = nseqs    
    fig, ax = plt.subplots(len(nseqs), sharex=True, sharey=True)
    for i in range(len(nseqs)):
        ax[i].plot(nseqs[i])
    
    nn = np.zeros(len(Ts), dtype=np.int64)
    nn[Ts_idx] = subseq_idx
    seed_motif = Ts[Ts_idx][subseq_idx : subseq_idx + m]

    for i, e in enumerate(Ts):
        nn[i] = np.argmin(stumpy.core.mass(seed_motif, e))
    
        r = Rectangle((nn[i], 0), m, 20, alpha=0.3)    
        ax[i].add_patch(r)
        ax[i].plot(np.arange(nn[i],nn[i]+m+1),Ts[i][nn[i]:nn[i]+m+1],'r')

    plt.show()

def visualize_consensus_motifs(nseqs, m, Ts_idx, subseq_idx):
    Ts = nseqs
    seed_motif = Ts[Ts_idx][subseq_idx : subseq_idx + m]
    x = np.linspace(0,1,m)
    nn = np.zeros(len(Ts), dtype=np.int64)
    nn[Ts_idx] = subseq_idx
    for i, e in enumerate(Ts):
        if i != Ts_idx:
            nn[i] = np.argmin(stumpy.core.mass(seed_motif, e))
            lw = 1
            label = None
        else:
            lw = 4
            label = 'Seed Motif'
        plt.plot(e[nn[i]:nn[i]+m], lw=lw, label=label)
    plt.title('The Consensus Motif')
    plt.legend()
    plt.show()
