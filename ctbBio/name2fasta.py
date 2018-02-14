#!/usr/bin/env python3

import sys
import os
from fasta import iterate_fasta as parse_fasta

def split_fasta(f, id2f):
    """
    split fasta file into separate fasta files based on list of scaffolds
    that belong to each separate file
    """
    opened = {}
    for seq in parse_fasta(f):
        id = seq[0].split('>')[1].split()[0]
        if id not in id2f:
            continue
        fasta = id2f[id]
        if fasta not in opened:
            opened[fasta] = '%s.fa' % fasta
        seq[1] += '\n'
        with open(opened[fasta], 'a+') as f_out:
            f_out.write('\n'.join(seq))

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('usage: name2fasta.py <sequences.fa> <sequence-ID-to-fasta-name.tsv>')
        exit()
    fasta, id2fasta = sys.argv[1:]
    if fasta == '-':
        fasta = sys.stdin
    else:
        fasta = open(fasta)
    if id2fasta == '-':
        id2fasta = sys.stdin
    else:
        id2fasta = open(id2fasta)
    id2fasta = {s.strip().split()[0]:s.strip().replace(' ', '_').split()[1] for s in id2fasta}
    split_fasta(fasta, id2fasta)
