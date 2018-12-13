#!/usr/bin/env python3

"""
filter specified number of hits from a blast or HMM search
that pass evalue and bit score thresholds

note: file must be sorted by query ID (and domain # for domtblout),
but does not have to be sorted by significance
"""

import os
import sys
import argparse
import pandas as pd
from operator import itemgetter

def top_hits(hits, num, column, reverse):
    """
    get top hits after sorting by column number
    """
    hits = sorted(hits, key = itemgetter(-1, column), reverse = reverse)
    for hit in hits[0:num]:
        yield hit

def numBlast(blast, numHits, evalueT = False, bitT = False):
    header = ['#query', 'target', 'pident', 'alen', 'mismatch', 'gapopen',
              'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bitscore']
    yield header
    prev, hits = None, []
    for line in blast:
        line = line.strip().split('\t')
        ID = line[0]
        evalue, bit = float(line[10]), float(line[11])
        line[10], line[11] = evalue, bit
        if id != prev:
            if len(hits) > 0:
                for hit in top_hits(hits, numHits, 10, False):
                    yield hit
            hits = []
        if evalueT == False and bitT == False:
            hits.append(line)
        elif evalue <= evalueT and bitT == False:
            hits.append(line)
        elif evalue <= evalueT and bit >= bitT:
            hits.append(line)
        elif evalueT == False and bit >= bitT:
            hits.append(line)
        prev = ID
    for hit in top_hits(hits, numHits, 10, False):
        yield hit

def numDomtblout(domtblout, numHits, evalueT, bitT):
    """
    parse hmm domain table output
    """
    header = ['#target name', 'target accession', 'tlen',
              'query name', 'query accession', 'qlen',
              'full E-value', 'full score', 'full bias',
              'domain #', '# domains',
              'domain c-Evalue', 'domain i-Evalue', 'domain score', 'domain bias',
              'hmm from', 'hmm to', 'seq from', 'seq to', 'env from', 'env to',
              'acc', 'target description']
    yield header
    prev, hits = None, []
    for line in domtblout:
        if line.startswith('#'):
            continue
        # parse line and get description
        line = line.strip().split()
        desc = ' '.join(line[18:])
        line = line[0:18]
        line.append(desc)
        # create ID based on query name and domain number
        ID = line[0] + line[9]
        # domain c-Evalue and domain score thresholds
        evalue, bitscore = float(line[11]), float(line[13])
        line[11], line[13] = evalue, bitscore
        if ID != prev:
            if len(hits) > 0:
                for hit in top_hits(hits, numHits, 11, False):
                    yield hit
            hits = []
        if evalueT == False and bitT == False:
            hits.append(line)
        elif evalue <= evalueT and bitT == False:
            hits.append(line)
        elif evalue <= evalueT and bit >= bitT:
            hits.append(line)
        elif evalueT == False and bit >= bitT:
            hits.append(line)
        prev = ID
    for hit in top_hits(hits, numHits, 11, False):
        yield hit

def numTblout(tblout, numHits, evalueT, bitT):
    """
    parse hmm table output
    """
    header = ['#target name', 'target accession',
              'query name', 'query accession',
              'full E-value', 'full score',  'full bias',
              'best E-value', 'best score',  'best bias',
              'exp', 'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc',
              'description of target']
    yield header
    prev, hits = None, []
    for line in tblout:
        if line.startswith('#'):
            continue
        # parse line and get description
        line = line.strip().split()
        desc = ' '.join(line[18:])
        line = line[0:18]
        line.append(desc)
        # ID and scores
        ID = line[0]
        evalue, bitscore = float(line[4]), float(line[5])
        line[4], line[5] = evalue, bitscore
        if ID != prev:
            if len(hits) > 0:
                for hit in top_hits(hits, numHits, 4, False):
                    yield hit
            hits = []
        if evalueT == False and bitT == False:
            hits.append(line)
        elif evalue <= evalueT and bitT == False:
            hits.append(line)
        elif evalue <= evalueT and bit >= bitT:
            hits.append(line)
        elif evalueT == False and bit >= bitT:
            hits.append(line)
        prev = ID
    for hit in top_hits(hits, numHits, 4, False):
        yield hit

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# filter blast or HMM tab output')
    parser.add_argument(\
            '-i', default = '-', \
            help = 'path to search results (sorted by query; default = stdin)')
    parser.add_argument(\
            '-n', default = 1, type = int, \
            help = 'number of hits (default = 1)')
    parser.add_argument(\
            '-e', default = False, type = float, \
            help = 'e-value threshold (default = None)')
    parser.add_argument(\
            '-b', default = False, type = float, \
            help = 'bit score threshold (default = None)')
    parser.add_argument(\
            '-f', default = 'b6', type = str,\
            help = 'format (default = b6, options: b6, domtblout, tblout)')
    args = vars(parser.parse_args())
    # check if file is from stdin
    if args['i'] == '-':
        args['i'] = sys.stdin
    else:
        args['i'] = open(args['i'])
    if args['f'] == 'b6':
        for hit in numBlast(args['i'], args['n'], args['e'], args['b']):
            print('\t'.join([str(i) for i in hit]))
    elif args['f'] == 'domtblout':
        for hit in numDomtblout(args['i'], args['n'], args['e'], args['b']):
            print('\t'.join([str(i) for i in hit]))
    elif args['f'] == 'tblout':
        for hit in numTblout(args['i'], args['n'], args['e'], args['b']):
            print('\t'.join([str(i) for i in hit]))
    else:
        print('unsupported format:', args['f'])
