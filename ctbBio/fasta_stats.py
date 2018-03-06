#!/usr/bin/python

import sys
import os
from ctbBio.fasta import iterate_fasta as fasta_parser



if __name__ == '__main__':
	fasta = sys.argv[1]
	stats(fasta)
