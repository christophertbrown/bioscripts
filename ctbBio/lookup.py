#!/usr/bin/env python3

"""
script for finding and replacing elements in a file
"""

import sys
import os

def f2lookup(f, lookup):
	"""
	find and replace elements in lookup within file
	"""
	lookup = {i: r for i, r in [l.strip().split('\t')[0:2] for l in lookup]}
	for line in f:
		line = line.strip()
		for find, replace in list(lookup.items()):
			line = line.replace(find, replace)
		yield line

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print('specify file and lookup')
		exit()
	for c, i in enumerate(sys.argv[1:], 1):
		if i == '-':
			i = sys.stdin
		else:
			i = open(i)
		sys.argv[c] = i
	f, lookup = sys.argv[1:]
	for line in f2lookup(f, lookup):
		print(line)
