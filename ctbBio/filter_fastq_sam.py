#!/usr/bin/env python3

import sys, os
from glob import glob as glob
from itertools import cycle as cycle

def sam_list(sam):
	"""
	get a list of mapped reads
	"""
	list = []
	for file in sam:
		for line in file:
			if line.startswith('@') is False:
				line = line.strip().split()
				id, map = line[0], int(line[1])
				if map != 4 and map != 8:
					list.append(id)
	return set(list)

def sam_list_paired(sam):
	"""
	get a list of mapped reads
	require that both pairs are mapped in the sam file in order to remove the reads
	"""
	list = []
	pair = ['1', '2']
	prev = ''
	for file in sam:
		for line in file:
			if line.startswith('@') is False:
				line = line.strip().split()
				id, map = line[0], int(line[1])
				if map != 4 and map != 8:
					read = id.rsplit('/')[0]
					if read == prev:
						list.append(read)
					prev = read
	return set(list)

def filter_paired(list):
	"""
	require that both pairs are mapped in the sam file in order to remove the reads
	"""
	pairs = {}
	filtered = []
	for id in list:
		read = id.rsplit('/')[0]
		if read not in pairs:
			pairs[read] = []
		pairs[read].append(id)
	for read in pairs:
		ids = pairs[read]
		if len(ids) == 2:
			filtered.extend(ids)
	return set(filtered)

def filter_fastq(fastq, list):
	c = cycle([1, 2, 3, 4])
	for file in fastq:
		if '/' in file:
			filtered = '%s.filtered.fastq' % (file.rsplit('.', 1)[0].rsplit('/', 1)[1])
			matched = '%s.matched.fastq' % (file.rsplit('.', 1)[0].rsplit('/', 1)[1])
		else:
			filtered = '%s.filtered.fastq' % (file.rsplit('.', 1)[0])
			matched = '%s.matched.fastq' % (file.rsplit('.', 1)[0])
		filtered, matched = open(filtered, 'w'), open(matched, 'w')
		switch = 1
		for line in open(file):
			n = next(c)
			if n == 1:
				id = line.strip().split()[0].split('@')[1]
				if id in list:
					switch = 0
				else:
					switch = 1
			if switch == 1:
				print(line.strip(), file=filtered)
			else:
				print(line.strip(), file=matched)
		filtered.close()
		matched.close()

def filter_fastq_paired(fastq, list):
	c = cycle([1, 2, 3, 4])
	for file in fastq:
		if '/' in file:
			filtered = '%s.filtered.fastq' % (file.rsplit('.', 1)[0].rsplit('/', 1)[1])
			matched = '%s.matched.fastq' % (file.rsplit('.', 1)[0].rsplit('/', 1)[1])
		else:
			filtered = '%s.filtered.fastq' % (file.rsplit('.', 1)[0])
			matched = '%s.matched.fastq' % (file.rsplit('.', 1)[0])
		filtered, matched = open(filtered, 'w'), open(matched, 'w')
		switch = 1
		for line in open(file):
			n = next(c)
			if n == 1:
				read = line.strip().split()[0].split('@')[1].rsplit('/', 1)[0]
				if read in list:
					switch = 0
				else:
					switch = 1
			if switch == 1:
				print(line.strip(), file=filtered)
			else:
				print(line.strip(), file=matched)
		filtered.close()
		matched.close()


def filter(fastq, sam, paired = False):
	""" 
	filter sequences that are shown to be mapped in the sam file
	reads not in sam file are in *.filtered.fastq
	reads that are in the sam file are in *.matched.fastq
	"""
	if paired is False:
		list = sam_list(sam)
	else:
		list = sam_list_paired(sam)
	if paired is False:
		filter_fastq(fastq, list)
	else:
		filter_fastq_paired(fastq, list)

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print('specify fastq, sam, and paired or unpaired')
		exit()
	fastq, sam, paired = glob(sys.argv[1]), glob(sys.argv[2]), sys.argv[3]
	if paired == 'paired':
		paired = True
	else:
		paired = False
	sam = [open(i) for i in sam]
	list(filter(fastq, sam, paired))
