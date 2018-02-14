#!/usr/bin/env python3

"""
script for reporting the length of sequences in a fasta file
"""

import sys
import os
from fasta import iterate_fasta as parse_fasta
import fasta as fasta_parser

def get_length(sequence):
    return [sequence[0], len(sequence[1].replace('-', '').replace('.', ''))]

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print('usage: fasta_length.py <seqs.fa or -> <length threshold or 0>', \
                file=sys.stderr)
        exit()
    fasta = sys.argv[1]
    threshold = float(sys.argv[2])
    if fasta == '-':
        fasta = sys.stdin
    for sequence in parse_fasta(fasta):
        length = get_length(sequence)
        if threshold == 0:
            length = [length[0].split('>')[1].split()[0], length[0].split('>')[1], str(length[1])]
            print('>%s' % ('\t'.join(length[1:])))
            print(sequence[1])
        elif length[1] >= threshold:
            length = [length[0].split('>')[1].split()[0], length[0].split('>')[1], str(length[1])]
            print('\n'.join(sequence))
