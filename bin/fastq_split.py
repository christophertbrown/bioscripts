#!/usr/bin/env python3

'''
script for taking an interleaved fastq file and printing the left and
right reads in separate fasta files
'''

import sys
import os
from itertools import cycle

def split(fastq, prefix):
    f1 = open('%s.R1.fastq' % (prefix), 'w')
    f2 = open('%s.R2.fastq' % (prefix), 'w')
    c = cycle([1, 1, 1, 1, 2, 2, 2, 2])
    for line in fastq:
        n = next(c)
        if n == 1:
            f1.write(line)
        else:
            f2.write(line)
    f1.close()
    f2.close()
    return [f1.name, f2.name]

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('specify fastq file and file prefix')
        exit()
    fastq, prefix = sys.argv[1], sys.argv[2]
    if fastq == '-':
        fastq = sys.stdin
    else:
        fastq = open(fastq)
    split(fastq, prefix)
