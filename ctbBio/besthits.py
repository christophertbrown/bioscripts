#!/usr/bin/env python3

"""
filter specified number of hits from a blast or HMM search
that pass evalue and bit score thresholds

note: file must be sorted by query ID, but does not have to be
sorted by significance
"""

import os
import sys
import argparse
import pandas as pd
from operator import itemgetter

def top_hits(hits, num, column):
    hits = sorted(hits, key = itemgetter(-1, column), reverse = True)
    for hit in hits[0:num]:
        yield [str(i) for i in hit[0:-1]]

def numBlast(blast, numHits, evalue = False, bit = False):
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
                for hit in top_hits(hits, numHits, 11):
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
    for hit in top_hits(hits, numHits, 11):
        yield hit

def numDomtblout(domtblout, numHits, evalue, bit):
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
    hmm = {h:[] for h in header}
    for line in domtblout:
        if line.startswith('#'):
            continue
        line = line.strip().split()
        if evalue is not False and float(line[11]) > evalue:
            continue
        if bit is not False and float(line[13]) < bit:
            continue
        line, desc = line[0:22], ' '.join(line[22:])
        line.append(desc)
        for i, h in zip(line, header):
            hmm[h].append(i)
    hmm = pd.DataFrame(hmm)
    hmm['domain score'] = [float(i) for i in hmm['domain score']]
    for query, df in hmm.groupby(by = ['#target name', 'domain #']):
        df = df.sort_values(by = ['domain score'], ascending = False)
        for hit in df[header].values[0:numHits]:
            yield [str(i) for i in hit]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# filter blast or HMM tabulat output')
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
            help = 'format (default = b6, options: b6, domtblout)')
    #    parser.add_argument(\
    #            '-t', default = 6, type = int,\
    #            help = 'number of threads (default = 6)')
    args = vars(parser.parse_args())
    # check if file is from stdin
    if args['i'] == '-':
        args['i'] = sys.stdin
    else:
        args['i'] = open(args['i'])
    if args['f'] == 'b6':
        for hit in numBlast(args['i'], args['n'], args['e'], args['b']):
            print('\t'.join(hit))
    elif args['f'] == 'domtblout':
        for hit in numDomtblout(args['i'], args['n'], args['e'], args['b']):
            print('\t'.join(hit))
    else:
        print('unsupported format:', args['f'])
