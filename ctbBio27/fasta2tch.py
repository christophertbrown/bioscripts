#!/usr/bin/python2.

"""
script for turning a fasta file into a tch file
"""

import sys
import os
import fasta as fasta_parser
from tokyocabinet import hash

def tch(fasta):
    tch = '%s.tch' % (fasta)
    db = hash.Hash()
    db.open(tch)
    for sequence in fasta_parser.iterate_fasta(fasta):
        db[sequence[0].split('>')[1].split()[0]] = '\n'.join(sequence)
    db.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print 'usage: fasta2tch.py <fasta.fa>'
        exit()
    fasta = sys.argv[1]
    tch(fasta)
