#!/usr/bin/env python3

"""
fix fasta headers so that they don't have weird characters
"""

import sys
import os
from ctbBio.fasta import iterate_fasta as parse_fasta

remove_characters = ['/', '\\', ':', ',', '(', ')', ' ', '|', ';',] # characters to remove from headers

def remove_char(header):
    for character in remove_characters:
        header = header.replace(character, '_')
    return header

def fix_fasta(fasta):
    """
    remove pesky characters from fasta file header
    """
    for seq in parse_fasta(fasta):
        seq[0] = remove_char(seq[0])
        if len(seq[1]) > 0:
            yield seq

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('specify fasta file')
        exit()
    fasta = sys.argv[1]
    if fasta == '-':
        fasta = sys.stdin
    else:
        fasta = open(fasta)
    for seq in fix_fasta(fasta):
        print('\n'.join(seq))
