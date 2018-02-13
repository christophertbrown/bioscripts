#!/usr/bin/env python3

"""
script for calculating n50 and other genome stats
"""

import sys
import os
import argparse

# ctb scripts
from fasta import iterate_fasta as parse_fasta

def gc(sequence):
    count = 0
    for base in sequence:
        base = base.lower()
        if base == 'g' or base == 'c':
            count += 1
    return float(float(count) / float(len(sequence))) * float(100)

def n50(fasta):
    length_list = []
    sequences = []
    for sequence in parse_fasta(fasta):
        length_list.append(float(len(sequence[1])))
        sequences.append(sequence[1])
    length_list.sort()
    length_list.reverse()
    total = float(sum(length_list))
    n = total * float(0.50)
    n50_value = running = length_list[0]
    for length in length_list:
        if running >= n:
            return n50_value, total, \
                                len(length_list), gc(''.join(sequences))
        else:
            n50_value = length
            running += n50_value

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# calculate n50 and other genome stats')
    parser.add_argument(\
            '-i', nargs = '*', action = 'store', required = True, \
            help = 'fasta(s)')
    args = vars(parser.parse_args())
    genomes = []
    for genome in args['i']:
        if genome == '-':
            genome = sys.stdin
        else:
            genome = open(genome)
        genomes.append(genome)
    print('\t'.join(['#genome', 'contigs', 'bases', 'gc', 'n50']))
    for genome in genomes:
        n50_value, total_bases, num_contigs, gc_cont = n50(genome)
        total_bases = '{:,}'.format(int(total_bases))
        gc_cont = round(gc_cont, 2)
        n50_value = '{:,}'.format(int(n50_value))
        print('\t'.join([str(i) for i in [genome.name, num_contigs, total_bases, gc_cont, n50_value]]))
