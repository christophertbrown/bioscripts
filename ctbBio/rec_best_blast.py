#!/usr/bin/env python3

"""
script for getting rec. best blast hits
"""

import sys
import os
from ctbBio.numblast import best as bestblast

def rec_hits(blast, evalue = float(0.01), bit = False):
	genes = {}
	rec_hits = []
	for out in blast:
		out = open(out)
		for hit in bestblast(out, 1, evalue, bit):
			query, match = hit[0].split()[0], hit[1]
			if query not in genes:
				genes[query] = {}
			genes[query][match] = 0
	for out in blast:
		out = open(out)
		for hit in bestblast(out, 1, evalue, bit):
			query, match = hit[0].split()[0], hit[1]
			if match in genes and query in genes[match]:
				genes[match][query] = 1
	for out in blast:
		out = open(out)
		for hit in bestblast(out, 1, evalue, bit):
			query, match = hit[0].split()[0], hit[1]
			if genes[query][match] == 1:
				rec_hits.append('\t'.join(hit))
	return set(rec_hits)

if __name__ == '__main__':
	blast = sys.argv[1:]
	for rec in rec_hits(blast):
		print(rec)
