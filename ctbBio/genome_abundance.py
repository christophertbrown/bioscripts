#!/usr/bin/env python3

"""
calculate genome coverage, relative abundance, and absolute abundance
"""

import sys
import os
from fasta import iterate_fasta as parse_fasta
from glob import glob as glob
import numpy

def calc_custom(custom, genome, scaffold, sequence, scaffold_coverage, total_bases):
	"""
	custom = {(reads mapped to scaffold)/(total reads for sample)}/(length of scaffold)
	"""
	index = 0
	if scaffold in scaffold_coverage: # what if the scaffold does not have bases mapped back to it? (this *should* not happen)
		if genome not in custom:
			custom[genome] = [[] for i in scaffold_coverage[scaffold]]
		for cov in scaffold_coverage[scaffold]:
			length = float(len(sequence[1]))
			bases = cov * length
			custom_value = ((bases) / (total_bases[index])) / length
			custom[genome][index].append(custom_value)
			index += 1
	return custom

def sum_coverage(coverage, std, genome, scaffold, sequence, scaffold_coverage):
	index = 0
	if scaffold in scaffold_coverage: # what if the scaffold does not have bases mapped back to it? (this *should* not happen)
		if genome not in std:
			std[genome] = [[] for i in scaffold_coverage[scaffold]]
		if genome not in coverage:
			coverage[genome] = [[0, 0] for i in scaffold_coverage[scaffold]]
		for cov in scaffold_coverage[scaffold]:
			length = float(len(sequence[1]))
			bases = cov * length
			coverage[genome][index][0] += bases
			coverage[genome][index][1] += length
			std[genome][index].append(cov)
			index += 1
	return coverage, std

def absolute_abundance(coverage, total_bases):
	"""
	absolute abundance = (number of bases mapped to genome / total number of bases in sample) * 100
	"""
	absolute = {}
	for genome in coverage:
		absolute[genome] = []
		index = 0
		for calc in coverage[genome]:
			bases = calc[0]
			total = total_bases[index]
			absolute[genome].append((bases / total) * float(100))
			index += 1
	total_assembled = [0 for i in absolute[genome]]
	for genome in absolute:
		index = 0
		for cov in absolute[genome]:
			total_assembled[index] += cov
			index += 1
	absolute['Unassembled'] = [(100 - i) for i in total_assembled]
	return absolute

def relative_abundance(coverage):
	"""
	cov = number of bases / length of genome
	relative abundance = [(cov) / sum(cov for all genomes)] * 100
	"""
	relative = {}
	sums = []
	for genome in coverage:
		for cov in coverage[genome]:
			sums.append(0)
		break
	for genome in coverage:
		index = 0
		for cov in coverage[genome]:
			sums[index] += cov
			index += 1
	for genome in coverage:
		index = 0
		relative[genome] = []
		for cov in coverage[genome]:
			if sums[index] == 0:
				relative[genome].append(0)
			else:
				relative[genome].append((cov / sums[index]) * float(100))
			index += 1
	return relative

def calc_total_mapped_bases(coverage):
	total = []
	for genome in coverage:
		for cov in coverage[genome]:
			total.append(0)
		break
	for genome in coverage:
		index = 0
		for cov in coverage[genome]:
			bases = cov[0]
			total[index] += bases
			index += 1
	return total

def calc_std(cov):
	std = {}
	for genome in cov:
		std[genome] = []
		for sample in cov[genome]:
			std[genome].append(numpy.std(sample))
	return std

def genome_coverage(genomes, scaffold_coverage, total_bases):
	"""
	coverage = (number of bases / length of genome) * 100
	"""
	coverage = {}
	custom = {}
	std = {}
	for genome in genomes:
		for sequence in parse_fasta(genome):
			scaffold = sequence[0].split('>')[1].split()[0]
			coverage, std = sum_coverage(coverage, std, genome, scaffold, sequence, scaffold_coverage)
			custom = calc_custom(custom, genome, scaffold, sequence, scaffold_coverage, total_bases)
	std = calc_std(std)
	custom_std = calc_std(custom)
	custom_av = {}
	for genome in custom:
		custom_av[genome] = []
		for sample in custom[genome]:
			custom_av[genome].append(numpy.mean(sample))
	for genome in coverage:
		print('%s\t%s' % (genome, coverage[genome][0][1]))
	if total_bases is True:
		total_bases = calc_total_mapped_bases(coverage)
	absolute = absolute_abundance(coverage, total_bases)
	for genome in coverage:
		calculated = []
		for calc in coverage[genome]:
			calculated.append(calc[0] / calc[1])
		coverage[genome] = calculated
	relative = relative_abundance(coverage)
	return coverage, std, absolute, relative, custom_av, custom_std

def print_genome_calculations(genomes, scaffold_coverage, total_bases, samples):
	coverage, std, absolute, relative, custom, custom_std = genome_coverage(genomes, scaffold_coverage, total_bases)
	for i in ['coverage', coverage], ['coverage_std', std], ['absolute abundance', absolute], ['relative abundance', relative], ['custom', custom], ['custom_std', custom_std]:
		out = open('%s.tsv' % (i[0]).replace(' ', '_'), 'w')
		i = i[1]
		print('#\t%s\ttotal' % ('\t'.join(samples)), file=out)
		for genome in i:
			print('%s\t%s' % (genome, '\t'.join([str(j) for j in i[genome]])), file=out)
		out.close()
			
def parse_scaffold_coverage(coverage):
	cov = {}
	for line in coverage:
		if line.startswith('#') is False:
			line = line.strip().split('\t')
			scaffold, coverage = line[0].split(':')[0], [float(i) for i in line[1:]]
			coverage.append(sum(coverage))
			cov[scaffold] = coverage
	return cov

def parse_total_bases(total_bases, samples):
	total = {}
	for line in total_bases:
		if line.startswith('#') is False:
			line = line.strip().split()
			sample, bases = '.'.join(line[0].split('.')[0:2]), float(line[2])
			if sample not in total:
				total[sample] = 0
			total[sample] += bases
	totals = [float(total[i]) for i in samples if 'total' not in i]
	totals.append(sum(totals))
	return totals

def sample_names(file):
	for line in file:
		if line.startswith('#') is True:
			header = line.strip().split('\t')[1:]
			return header
			break

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print('please provide coverage file, total reads file, and fasta files for each genome for calculating coverage')
		exit()
	genomes = glob(sys.argv[3])
	scaffold_coverage, total_bases = sys.argv[1], sys.argv[2]
#	scaffold_coverage = '/Users/ctb/banfield_lab/dora/5.assembly/coverage_ctb.tsv'
	print('using %s for coverage information' % (scaffold_coverage))
	scaffold_coverage = open(scaffold_coverage)
#	total_bases = '/Users/ctb/banfield_lab/dora/5.assembly/dora_base_counts.tsv'
	print('using %s for total bases' % (total_bases))
	print('coverage = (number of bases mapped to genome / length of genome) * 100')
	print('absolute abundance = (number of bases mapped to genome / total number of bases in sample) * 100')
	print('relative abundance = (coverage / sum(coverage for all genomes in sample)) * 100')
	print('custom = average({(reads mapped to scaffold)/(total reads for sample)} / (length of scaffold))')
	total_bases = open(total_bases)
	samples = sample_names(scaffold_coverage)
	scaffold_coverage = parse_scaffold_coverage(scaffold_coverage)
	total_bases = parse_total_bases(total_bases, samples)
	print_genome_calculations(genomes, scaffold_coverage, total_bases, samples)
