import numpy as np
from time import sleep

# def set_bit(mask, bitnum):
#     mask[bitnum >> 5] |= 1 << (bitnum & 31)
#     return mask

ULINE = '\033[4m'
PLAIN = '\033[0m'

def uline(string):
    return ULINE + string + PLAIN

# Standard Approx Shift-AND Algorithm
# def verify(txt, pat, k):
#     states = np.array([(1<<b)-1 for b in range(k+1)], dtype=np.int64)


def bndm(txt, pat, k=0, verbose=False):
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
        # print(all_states)
        # print(states)
        # state = (1<<m)-1
        shift_to = i+m-k
        j = shift_to-1
        if verbose:
            matchstr = np.array([' ']*m)

        while states[k] != 0 and j >= 0:
            if j < i:
                if states[k] >> m:
                    return txt[j+1:]
            else:
                prev = all_states[0]
                all_states[0] &= masks[ord(txt[j])]
                if verbose:
                    matchstr[np.array(list(format(all_states[0],'0'+str(m)+'b')))=='1'] = txt[j]
                    print(' '*i + pat)
                    print(' '*i + ''.join(matchstr))
                    print(txt[:j] + uline(txt[j:i+m]) + txt[i+m:])
                    sleep(1)
                all_states[0] <<= 1
                for d in range(k):
                    curr = all_states[d+1]
                    all_states[d+1] &= masks[ord(txt[j])]
                    all_states[d+1] |= (prev >> 1) | prev | all_states[d]
                    all_states[d+1] <<= 1
                    prev = curr
            prev = states[0]
            states[0] &= masks[ord(txt[j])]
            states[0] <<= 1
            for d in range(k):
                curr = states[d+1]
                states[d+1] &= masks[ord(txt[j])]
                states[d+1] |= (prev >> 1) | prev | states[d]
                states[d+1] <<= 1
                prev = curr
            # state &= masks[ord(txt[j])]
            # state <<= 1
            if all_states[k] >> m:
                shift_to = j
            j -= 1

        if verbose:
            print(' '*shift_to + uline(pat[:i+m-shift_to]) + pat[i+m-shift_to:])
            print()
            print(txt[:shift_to] + uline(txt[shift_to:i+m]) + txt[i+m:])
            sleep(1)
        i = shift_to
    return ''

def demo():
    print('>>> ...' + bndm('acggtacga' \
        + 'tcatgcgagtcgtagcaggttgatgt' \
        + 'accggactgcacgctgatgctcgtag' \
        + 'tgctcgatagtcgtagctgatcgatg' \
        + 'ctcgata', 'tagct', 0, True))