import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from string import printable as alphabet
from sw import sw
from nw import nw
from bndm import align

TRIALS = 10
DEFAULT_ALPH_SIZE = 10

# Generate length-n string from first s characters of alpha
def generate(n, s, alpha=alphabet):
    assert s <= 90
    return ''.join(np.random.choice(list(alpha[:s]), n))

# Generate random letter from first s characters of alpha
def rand_letter(s, alpha=alphabet):
    return generate(1, s, alpha)

# Generate random length-m substring of txt and add errors
# with probability 1-p (so should be roughly p-similar)
def get_similar(m, txt, p, s):
    i = np.random.randint(0, len(txt)-m+1)
    pat = txt[i:i+m]
    for j in range(m,0,-1):
        if np.random.uniform() < p:
            continue
        err = np.random.choice(['sub', 'ins', 'del'])
        if      err == 'sub':
            pat = pat[:j-1] + rand_letter(s) + pat[j:]
        elif    err == 'ins' and len(pat) < m:
            pat = pat[:j-1] + rand_letter(s) + pat[j-1:]
        else: # err == 'del'
            pat = pat[:j-1] +                  pat[j:]
    return pat

# Graph runtime of chosen algorithm against patlen m
def analyze(algorithm='bndm', alphabet_size=DEFAULT_ALPH_SIZE, \
            sim=0, trials=TRIALS, show_approx=True, verbose=True):
    N0 = 10 # txt size starts at 2^N0
    N  = 6  # txt size ends at 2^(N0+N)
    M  = 8  # pat size goes from 8*1-1 to 8*M-1, never pick M > 8!
    runtimes  = np.zeros((N,M))
    txt_sizes = [1<<b for b in range(N0,N0+N)]
    pat_sizes = [8*b+7 for b in range(M)]
    for k in range(trials):
        if verbose:
            print('>>> Trial ' + str(k+1) + ':')
        for i,n in enumerate(txt_sizes):
            if verbose:
                print('>>  i = ' + str(i+1) + ' out of ' + str(N))
            txt = generate(n, alphabet_size)
            for j,m in enumerate(pat_sizes):
                if verbose:
                    print('>   j = ' + str(j+1) + ' out of ' + str(M))
                if sim:
                    pat = get_similar(m, txt, sim, alphabet_size)
                else:
                    pat = generate(m, alphabet_size)
                if   algorithm == 'sw':
                    runtimes[i,j] += sw(txt, pat)[0]
                elif algorithm == 'nw':
                    runtimes[i,j] += nw(txt, pat)[0]
                else: #           'bndm':
                    runtimes[i,j] += align(txt, pat)
    runtimes /= trials  # average over trials
    if verbose:
        print('Plotting...')
    for i,n in enumerate(txt_sizes):
        plt.figure()    # plot runtime against m for each txt size
        plt.title(algorithm + ' runtime against pattern length at text length ' + str(n))
        plt.xlabel('Pattern length m')
        plt.ylabel('Algorithm runtime in ns')
        plt.plot(pat_sizes, runtimes[i,:], label=algorithm)
        if show_approx: # graph shape of asymptotic data, normalized to last m value
            linear_factor = runtimes[i,-1]/pat_sizes[-1]
            plt.plot(pat_sizes, linear_factor*np.array(pat_sizes), label='O(m)')
            if algorithm == 'bndm':
                log_factor = runtimes[i,-1]/np.log2(pat_sizes[-1])
                plt.plot(pat_sizes, log_factor*np.log2(pat_sizes), label='O(logm)')
                polylog_factor = runtimes[i,-1]/(pat_sizes[-1]*np.log2(pat_sizes[-1]))
                plt.plot(pat_sizes, polylog_factor*np.multiply(pat_sizes, np.log2(pat_sizes)), label='O(mlogm)')
                square_factor = runtimes[i,-1]/(pat_sizes[-1]**2)
                plt.plot(pat_sizes, square_factor*np.multiply(pat_sizes, pat_sizes), label='O(m^2)')
            plt.legend()
        plt.show()

# Load dna.txt into numpy array, strip extra data
def load():
    seqs = np.array(pd.read_csv('dna.txt', lineterminator=':', comment='>', header=None))[1:,0]
    for i in range(len(seqs)):
        seqs[i] = seqs[i].translate(str.maketrans('','','0123456789\n'))
    return seqs

# Pick a random length-m substring from a random sequence
def rand_seq(seqs, m):
    seq_num = np.random.choice(seqs.size)
    start = np.random.randint(len(seqs[seq_num])-m)
    seq = seqs[seq_num][start:start+m]
    return seq,seq_num

# Graph runtime of chosen algorithm against sequence length
def search_dna(algorithm='bndm', alphabet_size=DEFAULT_ALPH_SIZE, \
            sim=0, trials=TRIALS, show_approx=True, verbose=True):
    N = 100 # number of sequences to search through
    M = 8   # pat size goes from 8*1-1 to 8*M-1, never pick M > 8!
    runtimes  = np.zeros((N,M))
    pat_sizes = [8*b+7 for b in range(M)]
    for k in range(trials):
        if verbose:
            print('>>> Trial ' + str(k+1) + ':')
        seqs = np.random.choice(load(), N)
        for j,m in enumerate(pat_sizes):
            if verbose:
                print('>>  j = ' + str(j+1) + ' out of ' + str(M))
            pat,_ = rand_seq(seqs, m)
            # pat = generate(m, 4, 'acgt')
            for i in range(N):
                if verbose:
                    print('>   i = ' + str(i+1) + ' out of ' + str(N))            
                txt = seqs[i]
                if   algorithm == 'sw':
                    runtimes[i,j] += sw(txt, pat)[0]
                elif algorithm == 'nw':
                    runtimes[i,j] += nw(txt, pat)[0]
                else: #           'bndm':
                    runtimes[i,j] += align(txt, pat)
    runtimes = np.average(runtimes, axis=0)/trials  # average over trials
    if verbose:
        print('Plotting...')
    plt.figure()    # plot runtime against m averaged over each sequence searched
    plt.title(algorithm + ' runtime against pattern length in DNA sequence')
    plt.xlabel('Pattern length m')
    plt.ylabel('Algorithm runtime in ns')
    plt.plot(pat_sizes, runtimes, label=algorithm)
    if show_approx: # graph shape of asymptotic data, normalized to last m value
        linear_factor = runtimes[-1]/pat_sizes[-1]
        plt.plot(pat_sizes, linear_factor*np.array(pat_sizes), label='O(m)')
        if algorithm == 'bndm':
            log_factor = runtimes[-1]/np.log2(pat_sizes[-1])
            plt.plot(pat_sizes, log_factor*np.log2(pat_sizes), label='O(logm)')
            polylog_factor = runtimes[-1]/(pat_sizes[-1]*np.log2(pat_sizes[-1]))
            plt.plot(pat_sizes, polylog_factor*np.multiply(pat_sizes, np.log2(pat_sizes)), label='O(mlogm)')
            square_factor = runtimes[-1]/(pat_sizes[-1]**2)
            plt.plot(pat_sizes, square_factor*np.multiply(pat_sizes, pat_sizes), label='O(m^2)')
        plt.legend()
    plt.show()