#!/usr/bin/env python3

"""
script for calculating genome coverage
"""

import os
import sys
import argparse
import pandas as pd
from ctbBio.fasta import iterate_fasta as parse_fasta

def parse_cov(cov_table, scaffold2genome):
    """
    calculate genome coverage from scaffold coverage table
    """
    size   = {} # size[genome] = genome size
    mapped = {} # mapped[genome][sample] = mapped bases
    # parse coverage files
    for line in open(cov_table):
        line = line.strip().split('\t')
        if line[0].startswith('#'):
            samples = line[1:]
            samples = [i.rsplit('/', 1)[-1].split('.', 1)[0] for i in samples]
            continue
        scaffold, length = line[0].split(': ')
        length = float(length)
        covs  = [float(i) for i in line[1:]]
        bases = [c * length for c in covs]
        if scaffold not in scaffold2genome:
            continue
        genome = scaffold2genome[scaffold]
        if genome not in size:
            size[genome] = 0
            mapped[genome] = {sample:0 for sample in samples}
        # keep track of genome size
        size[genome] += length
        # keep track of number of mapped bases
        for sample, count in zip(samples, bases):
            mapped[genome][sample] += count
    # calculate coverage from base counts and genome size
    coverage = {'genome':[], 'genome size (bp)':[], 'sample':[], 'coverage':[]}
    for genome, length in size.items():
        for sample in samples:
            cov = mapped[genome][sample] / length
            coverage['genome'].append(genome)
            coverage['genome size (bp)'].append(length)
            coverage['sample'].append(sample)
            coverage['coverage'].append(cov)
    return pd.DataFrame(coverage)

def genome_coverage(covs, s2b):
    """
    calculate genome coverage from scaffold coverage
    """
    COV = []
    for cov in covs:
        COV.append(parse_cov(cov, s2b))
    return pd.concat(COV)

def parse_s2bs(s2bs):
    """
    convert s2b files to dictionary
    """
    s2b = {}
    for s in s2bs:
        for line in open(s):
            line = line.strip().split('\t')
            s, b = line[0], line[1]
            s2b[s] = b
    return s2b

def fa2s2b(fastas):
    """
    convert fastas to s2b dictionary
    """
    s2b = {}
    for fa in fastas:
        for seq in parse_fasta(fa):
            s = seq[0].split('>', 1)[1].split()[0]
            s2b[s] = fa.rsplit('/', 1)[-1].rsplit('.', 1)[0]
    return s2b

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# calculate genome coverage from scaffold coverage')
    parser.add_argument(\
            '-c', required = True, nargs = '*', \
            help = 'calculate_coverage.py scaffold coverage file(s) - required')
    parser.add_argument(\
            '-s', required = False, nargs = '*', \
            help = 'scaffold to bin files(s)')
    parser.add_argument(\
            '-f', required = False, nargs = '*', \
            help = 'fasta file(s) for each genome - use instead of -s')
    args = vars(parser.parse_args())
    s2bs, fastas, coverages = args['s'], args['f'], args['c']
    if s2bs is None and fastas is None:
        print('-s or -f is required')
        exit()
    if s2bs is not None:
        s2b = parse_s2bs(s2bs)
    else:
        s2b = fa2s2b(fastas)
    df = genome_coverage(coverages, s2b)
    df['genome: length (bp)'] = ['%s: %s' % (g, l) for g, l in zip(df['genome'], df['genome size (bp)'])]
    print(df.pivot('genome: length (bp)', 'sample', 'coverage').to_csv(sep = '\t'))
