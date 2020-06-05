import itertools
import numpy as np
from time import time_ns as time

# https://tiefenauer.github.io/blog/smith-waterman/
def matrix(a, b, match_score=1, gap_cost=1):
    H = np.zeros((len(a) + 1, len(b) + 1), np.int)
    H[0,:] = [-gap_cost*j for j in range(len(b)+1)]
    H[:,0] = 0

    for i, j in itertools.product(range(1, len(a)+1), range(1, len(b)+1)):
        match  = H[i - 1, j - 1] - (0 if a[i - 1] == b[j - 1] else match_score)
        delete = H[i - 1, j]     - gap_cost
        insert = H[i, j - 1]     - gap_cost
        H[i, j] = max(match, delete, insert)
    return H

# https://gist.github.com/radaniba/11019717
def traceback(H, a, b, verbose=0):
    j = len(b)
    t = i = np.argmax(H[:,j])
    a_ = b_ = ''
    while i and j:
        curr = H[i, j]
        diag = H[i-1, j-1]
        up   = H[i-1, j]
        left = H[i, j-1]
        i -= 1
        j -= 1
        if curr >= max(diag, up, left): # match
            a_ += a[i]
            b_ += b[j]
        elif left >= max(diag, up):     # deletion
            a_ += '-'
            b_ += b[j]
            i  += 1
        elif up >= diag:                # insertion
            a_ += a[i]
            b_ += '-'
            j  += 1
        else:                           # substitution
            a_ += '*'
            b_ += '*'
    while j:    # fill in remaining insertions
        j  -= 1
        a_ += '-'
        b_ += b[j]
    if verbose:
        print(b_[::-1])
        print(a_[::-1])
        print('>>> ...' + a[i:t] + '|' + a[t:])

def nw(a, b, match_score=1, gap_cost=1, verbose=0):
    start = time()
    a, b = a.lower(), b.lower()
    H = matrix(a, b, match_score, gap_cost)
    traceback(H, a, b, verbose)
    end = time()
    if verbose:
        print(str(end-start) + ' ns elapsed')
        print()
    return end-start