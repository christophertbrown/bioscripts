#!/usr/bin/env python3

"""
script for concatenating alignments
"""

import sys
import os
from ctbBio.fasta import iterate_fasta as parse_fasta

def concat_align(fastas):
    """
    concatenate alignments
    """
    # read in sequences
    fa2len = {}
    seqs = {}
    IDs = []
    for fasta in fastas:
        seqs[fasta] = {}
        for seq in parse_fasta(fasta):
            ID = seq[0].split('>')[1].split()[0]
            IDs.append(ID)
            seqs[fasta][ID] = seq[1]
        fa2len[fasta] = len(seq[1])
    # concat sequences
    IDs = set(IDs)
    concat = {}
    for fasta in fastas:
        for ID in IDs:
            if ID not in concat:
                concat[ID] = []
            if ID not in seqs[fasta]:
                concat[ID].append('-'*fa2len[fasta])
            else:
                concat[ID].append(seqs[fasta][ID])
    return concat

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('usage: concat_align.py <fastas*.afa>')
        exit()
    fastas = sys.argv[1:]
    for id, c in list(concat_align(fastas).items()):
        print('\n'.join(['>%s' % (id), ''.join(c)]))
