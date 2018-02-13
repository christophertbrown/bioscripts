#!/usr/bin/python2.7

"""
script for stripping out masked sequences if the masked sequence
is above a specified length
"""

import sys
import os
from fasta import iterate_fasta as parse_fasta
import argparse

def parse_masked(seq, min_len):
    """
    parse masked sequence into non-masked and masked regions
    """
    nm, masked = [], [[]]
    prev = None
    for base in seq[1]:
        if base.isupper():
            nm.append(base)
            if masked != [[]] and len(masked[-1]) < min_len:
                nm.extend(masked[-1])
                del masked[-1]
            prev = False 
        elif base.islower():
            if prev is False:
                masked.append([])
            masked[-1].append(base)
            prev = True
    return nm, masked

def strip_masked(fasta, min_len, print_masked):
    """
    remove masked regions from fasta file as long as
    they are longer than min_len
    """
    for seq in parse_fasta(fasta):
        nm, masked = parse_masked(seq, min_len)
        nm = ['%s removed_masked >=%s' % (seq[0], min_len), ''.join(nm)]
        yield [0, nm]
        if print_masked is True:
            for i, m in enumerate([i for i in masked if i != []], 1):
                m = ['%s insertion:%s' % (seq[0], i), ''.join(m)]
                yield [1, m]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# remove masked portion of sequences in fasta file')
    parser.add_argument(\
            '-f', required = True, \
            help = 'fasta file')
    parser.add_argument(\
            '-l', default = 0, \
            type = int, \
            help = 'minimum length of masked region required for removal')
    parser.add_argument(\
            '--print-masked', action = 'store_true', \
            help = 'print masked sequences to stderr')
    args = vars(parser.parse_args())
    fasta, min_len, print_masked = \
            args['f'], args['l'], args['print_masked']
    if fasta == '-':
        fasta = sys.stdin
    else:
        fasta = open(fasta)
    for i in strip_masked(fasta, min_len, print_masked):
        if i[0] == 0:
            print '\n'.join(i[1])
        else:
            print >> sys.stderr, '\n'.join(i[1])
