#!/usr/bin/env python3

"""
script for randomly subsetting a sam file
"""

import sys
import os
from itertools import cycle
from subprocess import Popen, PIPE
import argparse
import random

def sort_sam(sam, sort):
    """
    sort sam file
    """
    tempdir = '%s/' % (os.path.abspath(sam).rsplit('/', 1)[0])
    if sort is True:
        mapping = '%s.sorted.sam' % (sam.rsplit('.', 1)[0])
        if sam != '-':
            if os.path.exists(mapping) is False:
                os.system("\
                    sort -k1 --buffer-size=%sG -T %s -o %s %s\
                    " % (sbuffer, tempdir, mapping, sam)) 
        else:
            mapping = 'stdin-sam.sorted.sam'
            p = Popen("sort -k1 --buffer-size=%sG -T %s -o %s" \
                    % (sbuffer, tempdir, mapping), stdin = sys.stdin, shell = True) 
            p.communicate()
        mapping = open(mapping)
    else:
        if sam == '-':
            mapping = sys.stdin
        else:
            mapping = open(sam)
    return mapping

def sub_sam(sam, percent, sort = True, sbuffer = False):
    """
    randomly subset sam file
    """
    mapping = sort_sam(sam, sort)
    pool = [1 for i in range(0, percent)] + [0 for i in range(0, 100 - percent)]
    c = cycle([1, 2])
    for line in mapping:
        line = line.strip().split()
        if line[0].startswith('@'): # get the sam header
            yield line
            continue
        if int(line[1]) <= 20: # is this from a single read?
            if random.choice(pool) == 1:
                yield line
        else:
            n = next(c)
            if n == 1:
                prev = line
            if n == 2 and random.choice(pool) == 1:
                yield prev
                yield line

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# randomly subset sam file')
    parser.add_argument(\
            '-s', required = True, help = 'path to sorted sam file (- for stdin)')
    parser.add_argument(\
            '-p', required = True, type = int,\
            help = 'percent of reads to report, e.g. 50 (approximate)')
    parser.add_argument(\
            '--sort', action = 'store_true', help = 'sort the sam file')
    parser.add_argument(\
            '-b', default = "100", help = 'buffer size (GB) to use when sorting sam file (default = 100)')
    args = vars(parser.parse_args())
    sam, percent, sort, buff = args['s'], args['p'], args['sort'], args['b']
    for line in sub_sam(sam, percent, sort, buff):
        print('\t'.join(line))
