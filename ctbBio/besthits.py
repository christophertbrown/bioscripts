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
    hits.sort(key = itemgetter(column), reverse = reverse)
    for hit in hits[0:num]:
        yield hit

def numBlast_sort(blast, numHits, evalueT, bitT):
    """
    parse b6 output with sorting
    """
    header = ['#query', 'target', 'pident', 'alen', 'mismatch', 'gapopen',
              'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bitscore']
    yield header
    hmm = {h:[] for h in header}
    for line in blast:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        # Evalue and Bitscore thresholds
        line[10], line[11] = float(line[10]), float(line[11])
        evalue, bit = line[10], line[11]
        if evalueT is not False and evalue > evalueT:
            continue
        if bitT is not False and bit < bitT:
            continue
        for i, h in zip(line, header):
            hmm[h].append(i)
    hmm = pd.DataFrame(hmm)
    for query, df in hmm.groupby(by = ['#query']):
        df = df.sort_values(by = ['bitscore'], ascending = False)
        for hit in df[header].values[0:numHits]:
            yield hit

def numBlast(blast, numHits, evalueT = False, bitT = False, sort = False):
    """
    parse b6 output
    """
    if sort is True:
        for hit in numBlast_sort(blast, numHits, evalueT, bitT):
            yield hit
        return
    header = ['#query', 'target', 'pident', 'alen', 'mismatch', 'gapopen',
              'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bitscore']
    yield header
    prev, hits = None, []
    for line in blast:
        line = line.strip().split('\t')
        ID = line[0]
        line[10], line[11] = float(line[10]), float(line[11])
        evalue, bit = line[10], line[11]
        if ID != prev:
            if len(hits) > 0:
                # column is 1 + line index
                for hit in top_hits(hits, numHits, 11, True):
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
    for hit in top_hits(hits, numHits, 11, True):
        yield hit

def numDomtblout_sort(domtblout, numHits, evalueT, bitT):
    """
    parse hmm domain table output with sorting
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
        # domain c-Evalue and domain score thresholds
        line[11], line[13] = float(line[11]), float(line[13])
        evalue, bit = line[11], line[13]
        if evalueT is not False and evalue > evalueT:
            continue
        if bitT is not False and bit < bitT:
            continue
        line, desc = line[0:22], ' '.join(line[22:])
        line.append(desc)
        for i, h in zip(line, header):
            hmm[h].append(i)
    hmm = pd.DataFrame(hmm)
    for query, df in hmm.groupby(by = ['#target name', 'domain #']):
        df = df.sort_values(by = ['domain score'], ascending = False)
        for hit in df[header].values[0:numHits]:
            yield hit

def numDomtblout(domtblout, numHits, evalueT, bitT, sort):
    """
    parse hmm domain table output
    this version is faster but does not work unless the table is sorted
    """
    if sort is True:
        for hit in numDomtblout_sort(domtblout, numHits, evalueT, bitT):
            yield hit
        return
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
        line[11], line[13] = float(line[11]), float(line[13])
        evalue, bitscore = line[11], line[13]
        line[11], line[13] = evalue, bitscore
        if ID != prev:
            if len(hits) > 0:
                for hit in top_hits(hits, numHits, 13, True):
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
    for hit in top_hits(hits, numHits, 13, True):
        yield hit

def numTblout_sort(tblout, numHits, evalueT, bitT):
    """
    parse hmm hmm table output with sorting
    """
    header = ['#target name', 'target accession',
              'query name', 'query accession',
              'full E-value', 'full score',  'full bias',
              'best E-value', 'best score',  'best bias',
              'exp', 'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc',
              'description of target']
    yield header
    hmm = {h:[] for h in header}
    for line in tblout:
        if line.startswith('#'):
            continue
        line = line.strip().split()
        # domain full Evalue and full score thresholds
        line[4], line[5] = float(line[4]), float(line[5])
        evalue, bit = line[4], line[5]
        if evalueT is not False and evalue > evalueT:
            continue
        if bitT is not False and bit < bitT:
            continue
        line, desc = line[0:18], ' '.join(line[18:])
        line.append(desc)
        for i, h in zip(line, header):
            hmm[h].append(i)
    hmm = pd.DataFrame(hmm)
    for query, df in hmm.groupby(by = ['#target name']):
        df = df.sort_values(by = ['full score'], ascending = False)
        for hit in df[header].values[0:numHits]:
            yield hit

def numTblout(tblout, numHits, evalueT, bitT, sort):
    """
    parse hmm table output
    this version is faster but does not work unless the table is sorted
    """
    if sort is True:
        for hit in numTblout_sort(tblout, numHits, evalueT, bitT):
            yield hit
        return
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
        line[4], line[5] = float(line[4]), float(line[5])
        evalue, bitscore = line[4], line[5]
        line[4], line[5] = evalue, bitscore
        if ID != prev:
            if len(hits) > 0:
                for hit in top_hits(hits, numHits, 5, True):
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
    for hit in top_hits(hits, numHits, 5, True):
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
    parser.add_argument(\
            '--sort', action = 'store_true', \
            help = 'sort hits by query name (evalues are sorted by default)')
    args = vars(parser.parse_args())
    # check if file is from stdin
    if args['i'] == '-':
        args['i'] = sys.stdin
    else:
        args['i'] = open(args['i'])
    if args['f'] == 'b6':
        for hit in numBlast(args['i'], args['n'], args['e'], args['b'],
                args['sort']):
            print('\t'.join([str(i) for i in hit]))
    elif args['f'] == 'domtblout':
        for hit in numDomtblout(args['i'], args['n'], args['e'], args['b'],
                args['sort']):
            print('\t'.join([str(i) for i in hit]))
    elif args['f'] == 'tblout':
        for hit in numTblout(args['i'], args['n'], args['e'], args['b'],
                args['sort']):
            print('\t'.join([str(i) for i in hit]))
    else:
        print('unsupported format:', args['f'])
