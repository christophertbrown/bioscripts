#!/usr/bin/env python3

"""
script for converting fastq file to fasta file
"""

import sys
import os
from itertools import cycle

def fq2fa(fq):
    """
    convert fq to fa
    """
    c = cycle([1, 2, 3, 4])
    for line in fq:
        n = next(c)
        if n == 1:
            seq = ['>%s' % (line.strip().split('@', 1)[1])]
        if n == 2:
            seq.append(line.strip())
            yield seq

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('specify fastq file')
        exit()
    fq = sys.argv[1]
    if fq == '-':
        fq = sys.stdin
    else:
        fq = open(fq)
    for seq in fq2fa(fq):
        print('\n'.join(seq))
