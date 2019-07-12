#!/usr/bin/env python3

"""
script for randomly subsetting reads
in fastq file
"""

import os
import sys
import gzip
import random
import argparse
from itertools import cycle
from subprocess import Popen, PIPE

def parse_fq(fq):
    """
    parse fastq file
    """
    c = cycle([1, 2, 3, 4])
    read = []
    for line in fq:
        n = next(c)
        read.append(line.strip())
        if n == 4:
            yield read
            read = []

def sub_fq(R1, R2, percent):
    """
    randomly subset fastq files
    """
    pool = [1 for i in range(0, percent)] + [0 for i in range(0, 100 - percent)]
    for r1, r2 in zip(parse_fq(R1), parse_fq(R2)):
        if random.choice(pool) == 1:
            yield r1
            yield r2

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# randomly subset fastq files')
    parser.add_argument(\
            '-1', required = True, type = str, \
            help = 'path to forward reads')
    parser.add_argument(\
            '-2', required = True, type = str, \
            help = 'path to reverse reads')
    parser.add_argument(\
            '-p', required = True, type = int,\
            help = 'percent of reads to report, e.g. 50 (approximate)')
    args = vars(parser.parse_args())
    R1, R2, percent = args['1'], args['2'], args['p']
    if R1.rsplit('.', 1)[1] == 'gz':
        R1, R2 = gzip.open(R1, 'rt'), gzip.open(R2, 'rt')
    else:
        R1, R2 = open(R1), open(R2)
    for read in sub_fq(R1, R2, percent):
        print('\n'.join(read))
