#!/usr/bin/env python3

"""
summarize orthologer output
"""

import sys, os
import numpy as np

def count_genes(counts, genes):
    index = 0
    for gene in genes:
        if gene != '-':
            counts[index] += 1
        index += 1
    return counts

def count_orthologs(counts, genes):
    index = 1
    if genes[0] != '-':
        for gene in genes[1:]:
            if gene != '-':
                counts[index] += 1
            index += 1
    return counts

def append_scores(scores, current):
    current = current[0].split(' ** ')
    index = 0 
    for score in current:
        if score == '!':
            scores[index].append(100)
        elif score == '-' or score[0] == '-':
            continue
        else:
            pident = float(score.split()[2])
            scores[index].append(pident)
        index += 1
    return scores

def print_summary(genomes, genes, orthologs, pident):
    if len(genomes) > 2:
        header = ['# genome', 'genes', 'orthologs', 'average percent identity of orthologs']
        print('\t'.join(header))
        index = 0
        for genome in genomes:
            out = [genome]
            out.append(genes[index])
            out.append(orthologs[index])
            out.append(pident[index])
            out = [str(i) for i in out]
            print('\t'.join(out))
            index += 1
    else:
        header = ['# query genome', 'genes in query', 'reference genome', 'genes in reference', 'number of orthologs', 'average percent identity of orthologs']
        print('\t'.join(header))
        out = [genomes[0], genes[0], genomes[1], genes[1], orthologs[1], pident[1]]
        out = [str(i) for i in out]
        print('\t'.join(out))

def summarize(file):
    switch = 0
    for line in file:
        line = line.strip().split('\t')
        if line[0].startswith('### output'):
            switch = 1
            continue
        if switch == 0:
            continue
        if len(line) == 1:
            continue
        if line[0].startswith('#'):
            line[0] = line[0].split('# ')[1]
            genomes = line[::3]
            gene_counts = [0 for i in genomes]
            ortholog_counts = [0 for i in genomes]
            scores = [[] for i in genomes]
            continue
        gene_counts = count_genes(gene_counts, line[::3])
        ortholog_counts = count_orthologs(ortholog_counts, line[::3])
        scores = append_scores(scores, line[1::3])
    average_pident = [np.average(i) for i in scores]
    return genomes, gene_counts, ortholog_counts, average_pident

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('usage: orthologer_summary.py <orthologer_output.tsv or - if from stdin>')
        exit()
    file = sys.argv[1]
    if file == '-':
        file = sys.stdin
    else:
        file = open(file)
    genomes, gene_counts, ortholog_counts, average_pident = summarize(file)
    print_summary(genomes, gene_counts, ortholog_counts, average_pident)
