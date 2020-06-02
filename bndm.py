import numpy as np
from time import sleep

# def set_bit(mask, bitnum):
#     mask[bitnum >> 5] |= 1 << (bitnum & 31)
#     return mask

ULINE = '\033[4m'
PLAIN = '\033[0m'

def uline(string):
    return ULINE + string + PLAIN

def bndm(txt, pat, k=0, verbose=False):
    if verbose and k:
        print('Simulation currently not supported' \
            + 'for k-approximate matching.')
        verbose = False

    n = len(txt)
    m = len(pat)
    assert m <= 63
    masks = np.zeros(256, dtype=np.int64)
    for i,c in enumerate(pat):
        masks[ord(c)] |= 1 << (m-i-1)
    
    i = 0
    while i <= n-m+k:
        all_states = np.full(k+1, (1<<m)-1, dtype=np.int64)
        states = np.array([(1<<b)-1 for b in range(1,k+2)], dtype=np.int64)
        shift_to = i+m-k
        j = shift_to-1
        if verbose:
            matchstr = np.array([' ']*m)

        while all_states[k] and j >= max(0,i-2*k):
            prev = all_states[0]
            all_states[0] &= masks[ord(txt[j])]
            if verbose:
                matchstr[np.array(list(format(all_states[0],'0'+str(m)+'b')))=='1'] = txt[j]
                print(' '*i + pat)
                print(' '*i + ''.join(matchstr))
                print(' '*k + txt[:j] + uline(txt[j:i+m-k]) + txt[i+m-k:])
                sleep(1)
            all_states[0] &= (1<<m)-1
            all_states[0] <<= 1

            for d in range(k):
                curr = all_states[d+1]
                all_states[d+1] &= masks[ord(txt[j])]
                all_states[d+1] |= (prev >> 1) | prev | all_states[d]
                all_states[d+1] &= (1<<m)-1
                all_states[d+1] <<= 1
                prev = curr

            if all_states[k] >> m and j < shift_to:
                    shift_to = max(i+1,j)

            prev = states[0]
            states[0] &= masks[ord(txt[j])]
            # if verbose:
            #     matchstr[np.array(list(format(states[0],'0'+str(m)+'b')))=='1'] = txt[j]
            #     print(' '*i + pat)
            #     print(' '*i + ''.join(matchstr))
            #     print(' '*k + txt[:j] + uline(txt[j:i+m-k]) + txt[i+m-k:])
            #     sleep(1)
            states[0] <<= 1

            for d in range(k):
                curr = states[d+1]
                states[d+1] &= masks[ord(txt[j])]
                states[d+1] |= (prev >> 1) | prev | states[d]
                states[d+1] <<= 1
                prev = curr

            if states[k] >> m:
                return txt[j:]
            j -= 1

        if verbose:
            print(' '*shift_to + uline(pat[:i+m-k-shift_to]) + pat[i+m-k-shift_to:])
            print()
            print(' '*k + txt[:shift_to] + uline(txt[shift_to:i+m-k]) + txt[i+m-k:])
            sleep(1)
        i = shift_to
    return ''

def demo():
    print('>>> ...' + bndm('acggtacga' \
        + 'tcatgcgagtcgtagcaggttgatgt' \
        + 'accggactgcacgctgatgctcgtag' \
        + 'tgctcgatagtcgtagctgatcgatg' \
        + 'ctcgata', 'tagct', 0, True))