#!/usr/bin/env python3

"""
script for getting a specific region of a fasta file
"""

import sys
import os
import fasta as fasta_parser

def positions(region):
	region = region.split('-')
	start = int(region[0])
	if region[1] == '':
		stop = len(sequence[1])
	else:
		stop = int(region[1])
	return start, stop

def extract(sequence, region, length=0):
	start, stop = positions(region)
	header = '%s %s' % (sequence[0], region)
	cut = sequence[1][start:stop]
	yield fasta_parser.format_print([header, cut], length)

if __name__ == "__main__":
	if len(sys.argv) == 2:
		region = sys.argv[1]
		for sequence in fasta_parser.iterate_fasta(sys.stdin):
			for split in extract(sequence, region):
				print('\n'.join(split))
	elif len(sys.argv) == 3:
		region, length = sys.argv[1], int(sys.argv[2])
		for sequence in fasta_parser.iterate_fasta(sys.stdin):
			for split in extract(sequence, region, length):
				print('\n'.join(split))
	elif len(sys.argv) != 4:
		print('please specify the fasta file, the region that you would like to extract, \
'			'and the number of characters that you want per line, or 0 for one line \
'			'eg: 0-500 or 0- or 100-500')
		exit()
	else:
		fasta, region, length = sys.argv[1], sys.argv[2], int(sys.argv[3])
		if fasta == '-':
			fasta = sys.stdin
		else:
			fasta = open(fasta)
		for sequence in fasta_parser.iterate_fasta(fasta):
			for split in extract(sequence, region, length):
				print('\n'.join(split))
