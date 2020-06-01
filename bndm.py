import numpy as np

# def set_bit(mask, bitnum):
#     mask[bitnum >> 5] |= 1 << (bitnum & 31)
#     return mask

# def get_mask(masks, c):
#     if c in masks:
#         return masks[c]
#     return 0

def bndm(txt, pat):
    n = len(txt)
    m = len(pat)
    assert m <= 64
    masks = np.zeros(256, dtype=np.int64)
    for i,c in enumerate(pat):
        masks[ord(c)] |= 1 << (m-i-1)
        # print(bin(masks[ord(c)]))
    
    i = 0
    while i <= n-m:
        state = (1<<m)-1
        # print(bin(state))
        shift_to = i+m
        j = i+m-1
        while state != 0:
            if j < i:
                return txt[i:]
            state &= masks[ord(txt[j])] # get_mask(masks, txt[j])
            state <<= 1
            # print(bin(state))            
            if state >> m:
                shift_to = j
            j -= 1
        i = shift_to
    return ''