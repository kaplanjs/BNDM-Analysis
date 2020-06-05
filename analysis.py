import numpy as np
from matplotlib import pyplot as plt
from string import printable as alphabet
from sw import sw
from nw import nw
from bndm import align

TRIALS = 10
DEFAULT_ALPH_SIZE = 10

def generate(n, s):
    assert s <= 90
    return ''.join(np.random.choice(list(alphabet[:s]), n))

def rand_letter(s):
    return np.random.choice(list(alphabet[:s]))

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

def analyze(algorithm='bndm', alphabet_size=DEFAULT_ALPH_SIZE, \
            sim=0, trials=TRIALS, show_approx=True, verbose=True):
    N0 = 10
    N  = 8
    M  = 8  # never > 8!
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
                    runtimes[i,j] += sw(txt, pat)
                elif algorithm == 'nw':
                    runtimes[i,j] += nw(txt, pat)
                else: #           'bndm':
                    runtimes[i,j] += align(txt, pat)
    runtimes /= trials
    if verbose:
        print('Plotting...')
    for i,n in enumerate(txt_sizes):
        plt.figure()
        plt.title(algorithm + ' runtime against pattern length at text length ' + str(n))
        plt.xlabel('Pattern length m')
        plt.ylabel('Algorithm runtime in ns')
        plt.plot(pat_sizes, runtimes[i,:], label=algorithm)
        if show_approx:
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
