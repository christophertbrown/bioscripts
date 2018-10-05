#!/usr/bin/env python3

"""
script for getting a specified number of hits from a blast tsv file
under a specific evalue or bit score threshold, if provided

"""

import sys, os
from operator import itemgetter

def top_hits(hits, max):
    hits = sorted(hits, key = itemgetter(-1, 11), reverse = True)
    for hit in hits[0:max]:
        yield [str(i) for i in hit[0:-1]]

def best(blast, max, evalue = False, bit = False):
    prev, hits = None, []
    for line in blast:
        line = line.strip().split('\t')
        if line[10] == '*':
            line[10] = line[11] = float(line[2])
        else:
            line[10], line[11] = float(line[10]), float(line[11])
        id = line[0]
        line.append(float(line[10]) / -1)
        if id != prev:
            if len(hits) > 0:
                for hit in top_hits(hits, max):
                    yield hit
            hits = []
        if evalue == False and bit == False:
            hits.append(line)
        elif line[10] <= evalue and bit == False:
            hits.append(line)
        elif line[10] <= evalue and line[11] >= bit:
            hits.append(line)
        elif evalue == False and line[11] >= bit:
            hits.append(line)
        prev = id
    for hit in top_hits(hits, max):
        yield hit

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print('usage: numblast.py <blast tsv file> <number of hits> <evalue threshold> <bit score threshold>')
        print('# use - if from stdin or in place of a threshold')
        exit()
    blast, max, thresholds = sys.argv[1], int(sys.argv[2]), sys.argv[3:6]
    if blast == '-':
        blast = sys.stdin
    else:
        blast = open(blast)
    for i, t in enumerate(thresholds):
        if t == '-':
            t = False
        else:
            t = float(t)
        thresholds[i] = t
    e, bit = thresholds
    for hit in best(blast, max, e, bit):
        print('\t'.join(hit))
