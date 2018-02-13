#!/usr/bin/env python3

"""
script for converting a stockholm formatted alignment to fasta
"""

import os
import re
import sys

def stock2fa(stock):
    """
    convert stockholm to fasta
    """
    seqs = {}
    for line in stock:
        if line.startswith('#') is False and line.startswith(' ') is False and len(line) > 3:
            id, seq = line.strip().split()
            id = id.rsplit('/', 1)[0]
            id = re.split('[0-9]\|', id, 1)[-1]
            if id not in seqs:
                seqs[id] = []
            seqs[id].append(seq)
        if line.startswith('//'):
            break
    return seqs

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('specify stockholm formatted alignment')
        exit()
    stock = sys.argv[1]
    if stock == '-':
        stock = sys.stdin
    else:
        stock = open(stock)
    for id, seq in list(stock2fa(stock).items()):
        print('\n'.join(['>%s' % (id), ''.join(seq)]))
