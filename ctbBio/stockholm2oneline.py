#!/usr/bin/env python3

"""
script for converting a stockholm formatted alignment to single line
"""

import os
import re
import sys

def print_line(l):
    """
    print line if starts with ...
    """
    print_lines = ['# STOCKHOLM', '#=GF', '#=GS', ' ']
    if len(l.split()) == 0:
        return True
    for start in print_lines:
        if l.startswith(start):
            return True
    return False

def stock2one(stock):
    """
    convert stockholm to single line format
    """
    lines = {}
    for line in stock:
        line = line.strip()
        if print_line(line) is True:
            yield line
            continue
        if line.startswith('//'):
            continue
        ID, seq = line.rsplit(' ', 1)
        if ID not in lines:
            lines[ID] = ''
        else:
            # remove preceding white space
            seq = seq.strip()
        lines[ID] += seq
    for ID, line in lines.items():
        yield '\t'.join([ID, line])
    yield '\n//'

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('convert to single line stockholm formatted alignment')
        exit()
    stock = sys.argv[1]
    if stock == '-':
        stock = sys.stdin
    else:
        stock = open(stock)
    for line in stock2one(stock):
        print(line)
