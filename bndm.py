import numpy as np
from time import sleep
from time import time_ns as time
from nw import nw

# def set_bit(mask, bitnum):
#     mask[bitnum >> 5] |= 1 << (bitnum & 31)
#     return mask

DELAY = 0.2

ULINE  = '\033[4m'
HLIGHT = '\033[0;30;47m'
PLAIN  = '\033[0m'

def uline(string):
    return ULINE + string + PLAIN

def hlight(string):
    return HLIGHT + string + PLAIN

def reconstruct(states, masks, txt, pat, t, verbose):
    align_pat = ''
    align_txt = ''
    m = len(pat)
    k = len(states[:,0])-1  # available errors
    pos = 1 << (m-1)        # state position, align with pat (i)
    i = 0                   # pat position
    j = 0                   # txt position
    while j < t or i < m:
        if k > 0:
            if pos & states[k-1,t-j-1]:         # substitution
                align_pat += '*'
                align_txt += '*'
                pos >>= 1   # align state with pat (i)
                i += 1      # pat char consumed
                j += 1      # txt char consumed
                k -= 1      # error used
                continue
            if pos & (states[k-1,t-j-1] >> 1):  # insertion
                align_pat += '-'
                align_txt += txt[j]
                j += 1      # txt char consumed
                k -= 1      # error used
                continue
            if pos & states[k-1,t-j]:           # deletion
                align_pat += pat[i]
                align_txt += '-'
                pos >>= 1   # align state with pat (i)
                i += 1      # pat char consumed
                k -= 1      # error used
                continue
        if pos & states[k,t-j-1] & masks[ord(txt[j])]:  # match
            align_pat += pat[i]
            align_txt += txt[j]
            pos >>= 1   # align state with pat (i)
            i += 1      # pat char consumed
            j += 1      # txt char consumed
            continue
    if verbose:
        print(align_pat)
        print(align_txt)

def bndm(txt, pat, k=0, verbose=0):
    if verbose > 1 and k:
        print('Simulation currently not supported ' \
            + 'for k-approximate matching.')
        if(input('Press \'c\' to continue with ' \
            + 'simulation anyways, or any other ' \
            + 'key to continue without. ') != 'c'):
            verbose = 1

    n = len(txt)
    m = len(pat)
    assert m <= 63
    states = np.zeros((k+1, m+k+1), dtype=np.int64)   
    masks = np.zeros(256, dtype=np.int64)
    for i,c in enumerate(pat):
        masks[ord(c)] |= 1 << (m-i-1)

    i = 0
    while i <= n-m+k:   # window size m-k
        # initially all 1s, matches prefixes
        all_states = np.full(k+1, (1<<m)-1, dtype=np.int64)
        # only starts with 1, for complete matches
        states[:,0] = [(1<<b)-1 for b in range(1,k+2)]
        # j0 marks end of window
        # shift_to tracks start of next window
        j0 = shift_to = i+m-k
        # j marks window cursor (scan end to start0)
        j = shift_to-1
        if verbose > 1:
            matchstr = np.array([' ']*m)

        while all_states[k] and j >= max(0,i-2*k):

            # update all_states

            prev = all_states[0]                # track state pre-updates
            all_states[0] &= masks[ord(txt[j])] # screen matches
            all_states[0] &= (1<<m)-1           # truncate
            if verbose > 1 and k == 0:
                # prefixes matched to pattern so far
                matchstr[np.array(list(format(all_states[0],'0'+str(m)+'b')))=='1'] = txt[j]
                print(' '*i + pat)
                print(' '*i + ''.join(matchstr))
                print(' '*k + txt[:j] + hlight(txt[j:j0]) + txt[j0:])
                sleep(DELAY)
            all_states[0] <<= 1                 # shift for next round

            for d in range(k):                  # trickle down
                curr = all_states[d+1]
                all_states[d+1] &= masks[ord(txt[j])]
                all_states[d+1] |= (prev >> 1) | prev | all_states[d]
                all_states[d+1] &= (1<<m)-1
                all_states[d+1] <<= 1
                prev = curr

            if all_states[k] >> m and j < shift_to:
                # shift to matching prefix,
                # but always move at least 1
                shift_to = max(i+1,j)

            # repeat updates on states, but
            # track each step for alignment

            # prev = states[0,0]
            states[0,j0-j] = states[0,j0-j-1] & masks[ord(txt[j])]
            states[0,j0-j] <<= 1

            for d in range(k):
                # curr = states[d+1,j0-j]
                states[d+1,j0-j] = states[d+1,j0-j-1] & masks[ord(txt[j])]
                states[d+1,j0-j] |= (states[d,j0-j-1] >> 1) | states[d,j0-j-1] | states[d,j0-j]
                if verbose > 1 and d == k-1:
                    # state matched to pattern so far
                    matchstr[np.array(list(format(states[k,j0-j],'0'+str(m)+'b')))=='1'] = txt[j]
                    print(' '*i + pat)
                    print(' '*i + ''.join(matchstr))
                    print(' '*k + txt[:j] + hlight(txt[j:j0]) + txt[j0:])
                    sleep(DELAY)
                states[d+1,j0-j] <<= 1
                # prev = curr

            if states[k,j0-j] >> m:
                # construct alignment from match,
                # return text starting from match
                # and use '|' to mark match end
                reconstruct(states, masks, txt[j:j0], pat, j0-j, verbose)
                return txt[j:j0] + '|' + txt[j0:]
            j -= 1

        if verbose > 1:
            # prefix alignment for shift
            print(' '*shift_to + hlight(pat[:j0-shift_to]) + pat[j0-shift_to:])
            print()
            print(' '*k + txt[:shift_to] + hlight(txt[shift_to:j0]) + txt[j0:])
            sleep(1.5*DELAY)
        i = shift_to
    return ''

def align(txt, pat, max_err=-1, verbose=0):
    start = time()
    n = len(txt)
    m = len(pat)
    s = np.unique(list(txt+pat)).size
    if m > n:
        pat,txt = txt,pat
        m,n = n,m
    high = (int(m/(2+np.log2(m)/np.log2(s))) if max_err == -1 else max_err)
    low = 0

    # binary search on edit distance
    # for optimal alignment
    while high > low:
        avg = int((high+low)/2)
        if bndm(txt, pat, avg) == '':
            if avg == low:
                high,low = high,high
            else:
                high,low = high,avg
        else:
            high,low = avg,low
    if bndm(txt, pat, high, 0) == '':
        if nw(txt, pat)[1] < m/3:
            print('Alignment not found.')
    else:
        print('Alignment found!.')
    if high > m/(2+np.log2(m)/np.log2(s)):
        print('Runtime bound exceeded.')
    end = time()
    if verbose:
        print('>>> ...' + \
            bndm(txt, pat, high, verbose))
        print(str(end-start) + ' ns elapsed')
        print()
    return end-start

def demo():
    print('>>> ...' + bndm('acggtacga' \
        + 'tcatgcgagtcgtagcaggttgatgt' \
        + 'accggactgcacgctgatgctcgtag' \
        + 'tgctcgatagtcgtagctgatcgatg' \
        + 'ctcgata', 'tagct', 0, 2))

def demo2():
    print('>>> ...' + bndm('ctgatgcat' \
        + 'gctagcgctatgctgtcgatcgatcg' \
        + 'cgcagagatgtcgcgaatcatgcaag' \
        + 'tcgcagcgtatacgcgcggatccgga' \
        + 'tcacgc', 'gactgatacgtatcga' \
        + 'tcgata', 6, 2))