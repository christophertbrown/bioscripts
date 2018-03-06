#!/usr/bin/env python3

"""
script for dereplicating fasta files based on sequence names (not sequence similarity)
"""

import sys
import os
from ctbBio.fasta import iterate_fasta as parse_fasta

def append_index_id(id, ids):
    """
    add index to id to make it unique wrt ids
    """
    index = 1
    mod = '%s_%s' % (id, index)
    while mod in ids:
        index += 1
        mod = '%s_%s' % (id, index)
    ids.append(mod)
    return mod, ids

def de_rep(fastas, append_index, return_original = False):
    """
    de-replicate fastas based on sequence names
    """
    ids = []
    for fasta in fastas:
        for seq in parse_fasta(fasta):
            header = seq[0].split('>')[1].split()
            id = header[0]
            if id not in ids:
                ids.append(id)
                if return_original is True:
                    yield [header, seq]
                else:
                    yield seq
            elif append_index == True:
                new, ids = append_index_id(id, ids) 
                if return_original is True:
                    yield [header, ['>%s %s' % (new, ' '.join(header[1::])), seq[1]]]
                else:
                    yield ['>%s %s' % (new, ' '.join(header[1::])), seq[1]]

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('usage: nr_fasta.py <rename or exclude> <sequence.fasta>')
        exit()
    option, fastas = sys.argv[1], sys.argv[2:]
    for fasta in fastas:
        if fasta == '-':
            fastas = [sys.stdin]
    if option == 'rename':
        append_index = True
    elif option == 'exclude':
        append_index = False
    else:
        print('specify rename or exclude for redundant sequences')
        exit()
    for seq in de_rep(fastas, append_index):
        print('\n'.join(seq))
