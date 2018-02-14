#!/usr/bin/env python3

"""
get unmapped reads
"""

import sys
import os
import argparse

def unmapped(sam, mates):
    """
    get unmapped reads
    """
    for read in sam:
        if read.startswith('@') is True:
            continue
        read = read.strip().split()
        if read[2] == '*' and read[6] == '*':
                yield read
        elif mates is True:
            if read[2] == '*' or read[6] == '*':
                yield read
            for i in read:
                if i == 'YT:Z:UP':
                    yield read

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# unmapped reads from sam file')
    parser.add_argument(\
            '-s', required = True, help = 'path to sam file (- for stdin)')
    parser.add_argument(\
            '--mates', action = 'store_true', help = 'return both mates if one did not map (default: return neither mate if one mapped)')
    args = vars(parser.parse_args())
    sam, mates = args['s'], args['mates']
    if sam == '-':
        sam = sys.stdin
    else:
        sam = open(sam)
    for read in unmapped(sam, mates):
        print('\t'.join(read))
