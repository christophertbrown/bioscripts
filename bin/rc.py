#!/usr/bin/env python3

"""
script for getting the reverse complement of a nucleotide sequence
"""

import sys
import os
from fasta import iterate_fasta as parse_fasta

rc = {'A': 'T', \
      'T': 'A', \
      'G': 'C', \
      'C': 'G', \
      'N': 'N', \
      'a': 't', \
      't': 'a', \
      'g': 'c', \
      'c': 'g', \
      'n': 'n'}

def complement(seq):
	rev_c = []
	for base in seq[1]:
		rev_c.append(rc[base])
	return [seq[0], ''.join(rev_c)]

def reverse_complement(seq):
	rev_c = []
	for base in seq[1][::-1]:
		if base not in rc:
			rev_c.append('N')
		else:
			rev_c.append(rc[base])
	return [seq[0], ''.join(rev_c)]
	
if __name__ == '__main__':
	if len(sys.argv) != 3:
		print('specify fasta or - if from stdin and c (for complement) or rc (for reverse complement)')
		exit()
	fasta, option = sys.argv[1], sys.argv[2]
	if fasta == '-':
		fasta = sys.stdin
	else:
		fasta = open(fasta)
	if option == 'c':
		for seq in parse_fasta(fasta):
			print('\n'.join(complement(seq)))
	elif option == 'rc':
		for seq in parse_fasta(fasta):
			print('\n'.join(reverse_complement(seq)))
	else:
		print('specify fasta or - if from stdin \
				and c (for complement) or rc (for reverse complement)')
		exit()
