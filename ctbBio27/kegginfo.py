#!/usr/bin/env python2.7

"""
script for getting kegg functional analysis from fasta file with KOs in the header or a list of KOs
"""

import sys
import os
from tokyocabinet import hash
from ctbBio27.fasta import iterate_fasta as parse_fasta

keggpath = '%s/kegg' % (os.environ['databases'])

def option2kegg(option):
	if option == 'class':
		tch = '%s/ko2class.tch' % (keggpath)
	elif option == 'module':
		tch = '%s/ko2module.tch' % (keggpath)
	elif option == 'pathway':
		tch = '%s/ko2pathway.tch' % (keggpath)
	elif option == 'definition':
		tch = '%s/ko2definition.tch' % (keggpath)
	else:
		print('please specify pathway, class, module, or definition for kegg analysis')
		exit()
	return tch

def find_ko(line):
	ks = []
	remove = [';', ',', ':', '.', '-']
	nums = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
	for i in line:
		if i.startswith('K'):
			for c in remove:
				i = i.replace(c, '')
			if len(i) == 6:
				switch = 0
				for c in i[1:]:
					if c not in nums:
						switch = 1
						break
				if switch == 0:
					ks.append(i)
	return ks

def ko2kegg(file, option, file_type):
	tch = option2kegg(option)
	kegg = hash.Hash()
	kegg.open(tch)
	if file_type == 'fasta':
		for sequence in parse_fasta(file):
			header = sequence[0].split('>')[1]
			id = header.split()[0]
			yield header
			ks = set(find_ko(header.split()))
			for k in ks:
				if k in kegg:
					for function in kegg[k].split('|'):
						# - id - k - function
						yield '\t%s\t%s\t%s' % (id, k, function)
				else:
						yield '\t%s\t%s\tn/a' % (id, k)

	elif file_type == 'list':
		for line in file:
			line = line.strip()
			if len(line.split()) != 0:
				id = line.split()[0]
				yield line
				ks = set(find_ko(line.split()))
				for k in ks:
					if k in kegg:
						for function in kegg[k].split('|'):
							# - id - k - function
							yield '\t%s\t%s\t%s' % (id, k, function)
					else:
						yield '\t%s\t%s\tn/a' % (id, k)

	else:
		ks = set(find_ko(file))
		for k in ks:
			if k in kegg:
				for function in kegg[k].split('|'):
					yield [k, function]
			else:
				yield [k, 'n/a']
	kegg.close()

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print('please specify file (or \'-\' if from stdin), option (kegg: pathway, class, module, or definition), and file type (list or fasta)')
		exit()
	file = sys.argv[1]
	if file == '-':
		file = sys.stdin
	else:
		file = open(file)
	option = sys.argv[2]
	file_type = sys.argv[3]
	for annotation in ko2kegg(file, option, file_type):
		print(annotation)
