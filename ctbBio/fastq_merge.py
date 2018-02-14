#!/usr/bin/env python3

"""
script for merging separate fastq files into
an interleaved fastq file
"""

import sys
import os
import itertools
import gzip

def fq_merge(R1, R2):
    """
    merge separate fastq files
    """
    c = itertools.cycle([1, 2, 3, 4])
    for r1, r2 in zip(R1, R2):
        n = next(c)
        if n == 1:
            pair = [[], []]
        pair[0].append(r1.strip())
        pair[1].append(r2.strip())
        if n == 4:
            yield pair

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('usage: fastq_merge.py <R1.fastq> <R2.fastq>')
        exit()
    R1, R2 = sys.argv[1], sys.argv[2]
    if R1.rsplit('.', 1)[1] == 'gz':
        R1, R2 = gzip.open(R1, 'rt'), gzip.open(R2, 'rt')
    else:
        R1, R2 = open(R1), open(R2)
    for pair in fq_merge(R1, R2):
        print('\n'.join(itertools.chain(*pair)))
