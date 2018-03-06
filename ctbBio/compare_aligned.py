#!/usr/bin/env python3

import sys
import os
from multiprocessing import Pool as multithread
import argparse
import itertools
import numpy as np
#from Levenshtein import ratio as lr

from ctbBio.nr_fasta import de_rep as nr_fasta

def calc_pident(a, b):
    """
    calculate percent identity
    """
    m = 0 # matches
    mm = 0 # mismatches
    for A, B in zip(list(a), list(b)):
        if A == '-' or A == '.' or B == '-' or B == '.':
            continue
        if A == B:
            m += 1
        else:
            mm += 1
    return float(float(m)/float((m + mm))) * 100

def remove_gaps(A, B):
    """
    skip column if either is a gap
    """
    a_seq, b_seq = [], []
    for a, b in zip(list(A), list(B)):
        if a == '-' or a == '.' or b == '-' or b == '.':
            continue
        a_seq.append(a)
        b_seq.append(b)
    return ''.join(a_seq), ''.join(b_seq)

def compare_seqs(seqs):
    """
    compare pairs of sequences
    """
    A, B = seqs
    a, b = A[1], B[1] # actual sequences
    if len(a) != len(b):
        print('# reads are not the same length', file=sys.stderr)
        exit()
    pident = calc_pident(a, b)
    return A[0], B[0], pident

def compare_seqs_leven(seqs):
    """
    calculate Levenshtein ratio of sequences
    """
    A, B = seqs
    a, b = remove_gaps(A[1], B[1]) # actual sequences
    if len(a) != len(b):
        print('# reads are not the same length', file=sys.stderr)
        exit()
    pident = lr(a, b) * 100
    return A[0], B[0], pident

def pairwise_compare(afa, leven, threads, print_list):
    """
    make pairwise sequence comparisons between aligned sequences
    """
    # load sequences into dictionary
    seqs = {seq[0]: seq for seq in nr_fasta([afa], append_index = True)}
    # define all pairs
    pairs = itertools.combinations(list(seqs.values()), 2)
    pool = multithread(threads)
    # calc percent identity between all pairs - parallelize
    if leven is True:
        pident = pool.map(compare_seqs_leven, pairs)
    else:
        pident = pool.map(compare_seqs, pairs)
    pool.close()
    pool.join()
    return to_dictionary(pident, print_list)

def to_dictionary(pw, print_list):
    """
    - convert list of comparisons to dictionary
    - print list of pidents (if requested) to stderr
    """
    pairs = {}
    for p in pw:
        a, b, pident = p
        if a not in pairs:
            pairs[a] = {a: '-'}
        if b not in pairs:
            pairs[b] = {b: '-'}
        pairs[a][b] = pident
        pairs[b][a] = pident
        if print_list is True:
            A, B = a.split('>')[1], b.split('>')[1]
            print('\t'.join([str(i) for i in [A, B, pident]]), file=sys.stderr)
            print('\t'.join([str(i) for i in [B, A, pident]]), file=sys.stderr)
    return pairs

def print_pairwise(pw, median = False):
    """
    print matrix of pidents to stdout
    """
    names = sorted(set([i for i in pw]))
    if len(names) != 0:
        if '>' in names[0]:
            yield ['#'] + [i.split('>')[1] for i in names if '>' in i]
        else:
            yield ['#'] + names
        for a in names:
            if '>' in a:
                yield [a.split('>')[1]] + [pw[a][b] for b in names]
            else:
                out = []
                for b in names:
                    if b in pw[a]:
                        if median is False:
                            out.append(max(pw[a][b]))
                        else:
                            out.append(np.median(pw[a][b]))
                    else:
                        out.append('-')
                yield [a] + out

def print_comps(comps):
    """
    print stats for comparisons
    """
    if comps == []:
        print('n/a')
    else:
        print('# min: %s, max: %s, mean: %s' % \
            (min(comps), max(comps), np.mean(comps)))

def compare_clades(pw):
    """
    print min. pident within each clade and then matrix of between-clade max.
    """
    names = sorted(set([i for i in pw]))
    for i in range(0, 4):
        wi, bt = {}, {}
        for a in names:
            for b in pw[a]:
                if ';' not in a or ';' not in b:
                    continue
                pident = pw[a][b]
                cA, cB = a.split(';')[i], b.split(';')[i]
                if i == 0 and '_' in cA and '_' in cB:
                    cA = cA.rsplit('_', 1)[1]
                    cB = cB.rsplit('_', 1)[1]
                elif '>' in cA or '>' in cB:
                    cA = cA.split('>')[1]
                    cB = cB.split('>')[1]
                if cA == cB:
                    if cA not in wi:
                        wi[cA] = []
                    wi[cA].append(pident)
                else:
                    if cA not in bt:
                        bt[cA] = {}
                    if cB not in bt[cA]:
                        bt[cA][cB] = []
                    bt[cA][cB].append(pident)
        print('\n# min. within')
        for clade, pidents in list(wi.items()):
            print('\t'.join(['wi:%s' % str(i), clade, str(min(pidents))]))
        # print matrix of maximum between groups
        comps = []
        print('\n# max. between')
        for comp in print_pairwise(bt):
            if comp is not None:
                print('\t'.join(['bt:%s' % str(i)] + [str(j) for j in comp]))
                if comp[0] != '#':
                    comps.extend([j for j in comp[1:] if j != '-'])
        print_comps(comps)
        # print matrix of median between groups
        comps = []
        print('\n# median between')
        for comp in print_pairwise(bt, median = True):
            if comp is not None:
                print('\t'.join(['bt:%s' % str(i)] + [str(j) for j in comp]))
                if comp[0] != '#':
                    comps.extend([j for j in comp[1:] if j != '-'])
        print_comps(comps)

def matrix2dictionary(matrix):
    """
    convert matrix to dictionary of comparisons
    """
    pw = {}
    for line in matrix:
        line = line.strip().split('\t')
        if line[0].startswith('#'):
            names = line[1:]
            continue
        a = line[0]
        for i, pident in enumerate(line[1:]):
            b = names[i]
            if a not in pw:
                pw[a] = {}
            if b not in pw:
                pw[b] = {}
            if pident != '-':
                pident = float(pident)
            pw[a][b] = pident
            pw[b][a] = pident
    return pw

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# calculate percent identity of aligned reads')
    parser.add_argument(\
            '-a', default = False, \
            help = 'aligned fasta file')
    parser.add_argument(\
            '-m', default = False, \
            help = 'matrix of comparisons (for clade calcs.)')
    parser.add_argument(\
            '--list', action = 'store_true', \
            help = 'print list of pair-wise identities to stderr')
    parser.add_argument(\
            '--no-matrix', action = 'store_false', \
            help = 'do not print matrix')
    parser.add_argument(\
            '--clades', action = 'store_true', \
            help = 'compare clades based on header, e.g. >[0]Bacteria;[1]OD1;[2]unknown or >Bacteria;OD1;unknown')
    parser.add_argument(\
            '--leven', action = 'store_true', \
            help = 'calculate Levenshtein ratio')
    parser.add_argument(\
            '-t', default = 6, type = int,\
            help = 'number of threads (default: 6)')
    args = vars(parser.parse_args())
    afa, matrix, print_list, print_matrix, clades, leven, threads = \
            args['a'], args['m'], args['list'], args['no_matrix'], \
            args['clades'], args['leven'], args['t']
    if (afa is False and matrix is False) or (afa is not False and matrix is not False):
        print('# use -a or -m; -h for help', file=sys.stderr)
        exit()
    if afa is not False:
        if afa == '-':
            afa = sys.stdin
        else:
            afa = open(afa)
        pairwise = pairwise_compare(afa, leven, threads, print_list)
        if print_matrix is True:
            for i in print_pairwise(pairwise):
                print('\t'.join([str(j) for j in i]))
        if clades is True:
            compare_clades(pairwise)
    if matrix is not False:
        pairwise = matrix2dictionary(open(matrix))
        compare_clades(pairwise)
