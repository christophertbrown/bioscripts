#!/usr/bin/env python3

"""
script for determining 16S copy number from read mapping
"""

import sys
import os
from ctbBio.fasta import iterate_fasta as parse_fasta
from ctbBio.mapped import count_mismatches as count_mismatches

def get_overlap(a, b):
    """
    get overlap between ranges
    """
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def rna_bases(rna_cov, scaffold, bases, line):
    """
    determine if read overlaps with rna, if so count bases
    """
    start = int(line[3])
    stop = start + bases - 1
    if scaffold not in rna_cov:
        return rna_cov
    for pos in rna_cov[scaffold][2]:
        ol = get_overlap([start, stop], pos)
        rna_cov[scaffold][0] += ol
    return rna_cov

def parse_s2bins(s2bins):
    """
    parse ggKbase scaffold-to-bin mapping
        - scaffolds-to-bins and bins-to-scaffolds
    """
    s2b = {}
    b2s = {}
    for line in s2bins:
        line = line.strip().split()
        s, b = line[0], line[1]
        if 'UNK' in b:
            continue
        if len(line) > 2:
            g = ' '.join(line[2:])
        else:
            g = 'n/a'
        b = '%s\t%s' % (b, g)
        s2b[s] = b 
        if b not in b2s:
           b2s[b] = []
        b2s[b].append(s)
    return s2b, b2s

def parse_rna(rna, s2bins, min_rna):
    """
    parse [16,23]SfromHMM.py output
    - rna_cov[scaffold] = [0, 0, []] # [bases, length, ranges]
    """
    rna_cov = {}
    for seq in parse_fasta(rna):
        # check that length passes threshold
        length = len(seq[1])
        if length < min_rna:
            continue
        # check if sequence is binnned
        s = seq[0].split('>')[1].split()[0]
        if s not in s2bins:
            continue
        if s not in rna_cov:
            rna_cov[s] = [0, 0, []]
        position = [int(i) for i in seq[0].rsplit('pos=', 1)[1].split()[0].split('-')]
        rna_cov[s][2].append(position)
        rna_cov[s][1] += length
    return rna_cov

def filter_missing_rna(s2bins, bins2s, rna_cov):
    """
    remove any bins that don't have 16S
    """
    for bin, scaffolds in list(bins2s.items()):
        c = 0
        for s in scaffolds:
            if s in rna_cov:
                c += 1
        if c == 0:
            del bins2s[bin]
    for scaffold, bin in list(s2bins.items()):
        if bin not in bins2s:
            del s2bins[scaffold]
    return s2bins, bins2s

def calc_bin_cov(scaffolds, cov):
    """
    calculate bin coverage
    """
    bases = sum([cov[i][0] for i in scaffolds if i in cov])
    length = sum([cov[i][1] for i in scaffolds if i in cov])
    if length == 0:
        return 0
    return float(float(bases)/float(length))

def copies(mapping, s2bins, rna, min_rna = 800, mismatches = 0):
    """
    1. determine bin coverage
    2. determine rRNA gene coverage
    3. compare
    """
    cov = {} # cov[scaffold] = [bases, length]
    s2bins, bins2s = parse_s2bins(s2bins)
    rna_cov = parse_rna(rna, s2bins, min_rna)
    s2bins, bins2s = filter_missing_rna(s2bins, bins2s, rna_cov)
    # count bases mapped to scaffolds and rRNA gene regions
    for line in mapping:
        line = line.strip().split()
        # get scaffold lengths
        if line[0].startswith('@'):
            if line[0].startswith('@SQ') is False:
                continue
            s = line[1].split(':')[1]
            l = int(line[2].split(':')[1])
            # check if scaffold is binned
            if s not in s2bins:
                continue
            if s not in cov:
                cov[s] = [0, l]
        # check mismatch threshold
        mm = count_mismatches(line)
        if mm is False or mm > mismatches:
            continue
        # check that scaffold is in bin
        s, bases = line[2], len(line[9])
        if s not in cov:
            continue
        cov[s][0] += bases 
        rna_cov = rna_bases(rna_cov, s, bases, line) 
    print('# mismatches threshold: %s' % (mismatches))
    header = ['#rRNA scaffold', 'rRNA genes >=%sbp on scaffold' % (min_rna), \
            'rRNA coverage', \
            'bin', 'bin info', 'bin coverage', \
            'rRNAs >=%sbp in bin' % (min_rna), \
            'rRNA coverage/bin coverage', \
            'estimated number of copies']
    print('\t'.join(header))
    for bin, scaffolds in list(bins2s.items()):
        rna_count = sum([len(rna_cov[s][2]) for s in scaffolds if s in rna_cov])
        for s in scaffolds:
            if s not in rna_cov:
                continue
            out = []
            counts = rna_cov[s]
            bin_cov = calc_bin_cov(bins2s[bin], cov)
            num_genes = len(counts[2])
            rna_coverage = float(float(counts[0])/float(counts[1]))
            if bin_cov == 0:
                rna_div_bin = 0
            else:
                rna_div_bin = float(rna_coverage/bin_cov)
            est = int(max([rna_count, counts, rna_div_bin]))
            out = [s, num_genes, rna_coverage, bin, bin_cov, rna_count, rna_div_bin, est]
            print('\t'.join([str(i) for i in out]))

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('usage: rRNA_copies.py <mapping.sam> <scaffolds2bins.tsv> <[16S,23s]fromHMM.py sequences>')
        exit()
    mapping = sys.argv[1]
    if mapping == '-':
        mapping = sys.stdin
    else:
        mapping = sys.argv[1]
    s2bins, rna = [open(i) for i in sys.argv[2:]]
    copies(mapping, s2bins, rna)
