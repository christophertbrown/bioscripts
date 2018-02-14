#!/usr/bin/env python3

"""
script for randomly messing up a genome
"""

import os
import sys
import random
import argparse
import numpy as np

# ctb
from fasta import iterate_fasta as parse_fasta

def plot_dist_normal(s, mu, sigma):
    """
    plot distribution
    """
    import matplotlib.pyplot as plt
    count, bins, ignored = plt.hist(s, 30, normed=True)
    plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) \
            * np.exp( - (bins - mu)**2 / (2 * sigma**2) ), \
            linewidth = 2, color = 'r')
    plt.show()

def rev_c(read):
    """
    return reverse completment of read
    """
    rc = []
    rc_nucs = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
    for base in read:
        rc.extend(rc_nucs[base.upper()])
    return rc[::-1]

def shuffle_genome(genome, cat, fraction = float(100), plot = True, \
        alpha = 0.1, beta = 100000, \
        min_length = 1000, max_length = 200000):
    """
    randomly shuffle genome
    """
    header = '>randomized_%s' % (genome.name)
    sequence = list(''.join([i[1] for i in parse_fasta(genome)]))
    length = len(sequence)
    shuffled = []
    # break genome into pieces
    while sequence is not False:
        s = int(random.gammavariate(alpha, beta))
        if s <= min_length or s >= max_length:
            continue
        if len(sequence) < s:
            seq = sequence[0:]
        else:
            seq = sequence[0:s]
        sequence = sequence[s:]
#        if bool(random.getrandbits(1)) is True:
#            seq = rev_c(seq)
#            print('fragment length: %s reverse complement: True' % ('{:,}'.format(s)), file=sys.stderr)
#        else:
#            print('fragment length: %s reverse complement: False' % ('{:,}'.format(s)), file=sys.stderr)
        shuffled.append(''.join(seq))
        if sequence == []:
            break
    # shuffle pieces
    random.shuffle(shuffled)
    # subset fragments
    if fraction == float(100):
        subset = shuffled
    else:
        max_pieces = int(length * fraction/100)
        subset, total = [], 0
        for fragment in shuffled:
            length = len(fragment)
            if total + length <= max_pieces:
                subset.append(fragment)
                total += length
            else:
                diff = max_pieces - total
                subset.append(fragment[0:diff])
                break
    # combine sequences, if requested
    if cat is True:
        yield [header, ''.join(subset)]
    else:
        for i, seq in enumerate(subset):
            yield ['%s fragment:%s' % (header, i), seq]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# randomly re-arrange genome')
    parser.add_argument(\
            '-f', nargs = '*', action = 'store', required = True, \
            help = 'fasta(s)')
    parser.add_argument(\
            '-p', type = float, default = 100, 
            help = 'percent of genome to return (default = 100)')
    parser.add_argument(\
            '--cat', action = 'store_true', \
            help = 'concatenate random fragments')
    args = vars(parser.parse_args())
    for genome in args['f']:
        if genome == '-':
            genome = sys.stdin
        else:
            genome = open(genome)
        for seq in shuffle_genome(genome, args['cat'], fraction = args['p']):
            print('\n'.join(seq))
