import numpy as np
from time import sleep

# def set_bit(mask, bitnum):
#     mask[bitnum >> 5] |= 1 << (bitnum & 31)
#     return mask

ULINE = '\033[4m'
PLAIN = '\033[0m'

def uline(string):
    return ULINE + string + PLAIN

def bndm(txt, pat, verbose=False):
    n = len(txt)
    m = len(pat)
    assert m <= 64
    masks = np.zeros(256, dtype=np.int64)
    for i,c in enumerate(pat):
        masks[ord(c)] |= 1 << (m-i-1)
    
    i = 0
    while i <= n-m:
        state = (1<<m)-1
        shift_to = i+m
        j = i+m-1
        if verbose:
            matchstr = np.array([' ']*m)

        while state != 0:
            if j < i:
                return txt[i:]
            state &= masks[ord(txt[j])]
            if verbose:
                matchstr[np.array(list(format(state,'0'+str(m)+'b')))=='1'] = txt[j]
                print(' '*i + pat)
                print(' '*i + ''.join(matchstr))
                print(txt[:j] + uline(txt[j:i+m]) + txt[i+m:])
                sleep(1)
            state <<= 1
            if state >> m:
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
        + 'ctcgata', 'tagct', True))