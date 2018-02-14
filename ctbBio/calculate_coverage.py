#!/usr/bin/env python3

import sys
import os
import argparse

def length_and_bases(coverage, sam):
    for line in sam:
        line = line.strip()
        if line.startswith('@SQ'):
            line = line.strip().split()
            scaffold, length = line[1].split(':', 1)[1], float(line[-1].split(':', 1)[1])
            if scaffold not in coverage:
                coverage[scaffold] = [length, {}]
        elif line.startswith('@') is False:
            line = line.split('\t')
            scaffold, bases = line[2], float(len(line[9]))
            if scaffold not in coverage:
                coverage[scaffold] = [length, {}]
            map = int(line[1])
            if map != 4 and map != 8:
                if scaffold != '*':
                    if sam not in coverage[scaffold][1]:
                        coverage[scaffold][1][sam] = 0
                    coverage[scaffold][1][sam] += bases
    return coverage

def combine_by_sample(coverage):
    combined_coverage = {}
    sams = []
    for scaffold in coverage:
        length, combined = coverage[scaffold][0], {}
        for sam in coverage[scaffold][1]:
            comb = '.'.join(sam.name.split('.')[0:2])
            if comb not in sams:
                sams.append(comb)
            if comb not in combined:
                combined[comb] = 0
            combined[comb] += coverage[scaffold][1][sam]
        combined_coverage[scaffold] = [length, combined]
    return combined_coverage, sams

def calculate_coverage(bases):
    coverage = {}
    for scaffold in bases:
        length, counts = bases[scaffold][0], bases[scaffold][1]
        scaffold_coverage = {}
        for count in counts:
            scaffold_coverage[count] = float(counts[count] / length)
        coverage[scaffold] = [length, scaffold_coverage]
    return coverage

def print_coverage(coverage, sams):
    out = ['# scaffold: length']
    for sam in sams:
        out.append(sam.name)
    yield out
    for scaffold in coverage:
        length, cov = coverage[scaffold][0], coverage[scaffold][1]
        out = ['%s: %s' % (scaffold, length)]
        for sam in sams:
            if sam in cov:
                out.append(cov[sam])
            else:
                out.append(0)
        yield out

def iterate_sams(sams, combine = False):
    coverage = {}
    for sam in sams:
        coverage = length_and_bases(coverage, sam)
    if combine is True:
        coverage, sams = combine_by_sample(coverage)
    coverage = calculate_coverage(coverage)
    return coverage, sams

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# calculate coverage from sam file')
    parser.add_argument(\
            '-s', nargs = '*', action = 'store', required = True, \
            help = 'sam(s)')
    args = vars(parser.parse_args())
    sams = []
    for sam in sorted(args['s']):
        if sam == '-':
            sam = sys.stdin
        else:
            sam = open(sam)
        sams.append(sam)
    coverage, sams = iterate_sams(sams)
    for i in print_coverage(coverage, sams):
        print('\t'.join([str(j) for j in i]))
