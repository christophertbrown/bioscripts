#!/usr/bin/env python3

"""
script for 'stripping' out columns in a MSA that represent gaps for X percent of sequences
"""

import sys
import os
from fasta import iterate_fasta as parse_fasta

def plot_gaps(plot, columns):
	"""
	plot % of gaps at each position
	"""
	from plot_window import window_plot_convolve as plot_window
#	plot_window([columns], len(columns)*.01, plot)
	plot_window([[100 - i for i in columns]], len(columns)*.01, plot)

def strip_msa_100(msa, threshold, plot = False):
	"""
	strip out columns of a MSA that represent gaps for X percent (threshold) of sequences
	"""
	msa = [seq for seq in parse_fasta(msa)]
	columns = [[0, 0] for pos in msa[0][1]] # [[#bases, #gaps], [#bases, #gaps], ...]
	for seq in msa:
		for position, base in enumerate(seq[1]):
			if base == '-' or base == '.':
				columns[position][1] += 1
			else:
				columns[position][0] += 1
	columns = [float(float(g)/float(g+b)*100) for b, g in columns] # convert to percent gaps
	for seq in msa:
		stripped = []
		for position, base in enumerate(seq[1]):
			if columns[position] < threshold:
				stripped.append(base)
		yield [seq[0], ''.join(stripped)]
	if plot is not False:
		plot_gaps(plot, columns)

def strip_msa(msa, threshold, plot = False):
	"""
	strip out columns of a MSA that represent gaps for X percent (threshold) of sequences
	"""
	msa = [seq for seq in parse_fasta(msa)]
	columns = [[0, 0] for pos in msa[0][1]] # [[#bases, #gaps], [#bases, #gaps], ...]
	for seq in msa:
		for position, base in enumerate(seq[1]):
			if base == '-' or base == '.':
				columns[position][1] += 1
			else:
				columns[position][0] += 1
	columns = [float(float(g)/float(g+b)*100) for b, g in columns] # convert to percent gaps
	for seq in msa:
		stripped = []
		for position, base in enumerate(seq[1]):
			if columns[position] <= threshold:
				stripped.append(base)
		yield [seq[0], ''.join(stripped)]
	if plot is not False:
		plot_gaps(plot, columns)


if __name__ == '__main__':
	if len(sys.argv) != 4:
		print('specify MSA, threshold, and file name for pdf or False')
		exit()
	msa, threshold, plot = sys.argv[1], float(sys.argv[2]), sys.argv[3]
	if msa == '-':
		msa = sys.stdin
	else:
		msa = open(msa)
	if plot == 'False':
		plot = False
	if threshold == 100:
		for seq in strip_msa_100(msa, threshold, plot):
			print('\n'.join(seq))
	else:
		for seq in strip_msa(msa, threshold, plot):
			print('\n'.join(seq))

