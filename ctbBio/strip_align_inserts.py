#!/usr/bin/env python3

"""
script for removing insertion columns ('.' and lower case bases) from 
alignment fasta file
"""

import sys
import os
from ctbBio.fasta import iterate_fasta as parse_fasta

def strip_inserts(fasta):
    """
    remove insertion columns from aligned fasta file
    """
    for seq in parse_fasta(fasta):
        seq[1] = ''.join([b for b in seq[1] if b == '-' or b.isupper()])
        yield seq

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('specify aligned fasta file')
        exit()
    fasta = sys.argv[1]
    if fasta == '-':
        fasta = sys.stdin
    else:
        fasta = open(fasta)
    for seq in strip_inserts(fasta):
        print('\n'.join(seq))
