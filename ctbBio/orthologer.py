#!/usr/bin/env python3

"""
chris brown
banfield lab
ctb@berkeley.edu

draft script

script for finding orthologs and making a table for comparing
fasta files or genomes (.faa or .fna ORF predictions) by ordering the output

global or reference modes
blast / usearch every genome against every other genome
make a dictionary using every gene ID as a key:
    genes[ID] = [sort id, fasta, description], match
    match[ID] = [rec. hit?, [alignment length, pident, evalue, bitscore]]
output = [id, scores, description, id2, scores2, description2, id3, scores3, description3, ...]
scores = [! ** alignment lenght2 | pidnet2 | evalue2 | bitscore2 ** alignment length3 | pidnet3 | evalue3 | bitscore3 ** ...]
"""

# improve sorting to use a list (priority sort the list)
# improve sorting module
# documentation
# for global - must examine every combination of possible orthologs
# - do this in a non-redundant manner (eg only half of total gene list)
# - will not need to convert to a set
# - will need to evaluate the list, and the complete list for redundnacies

import os
import sys
from operator import itemgetter
from datetime import datetime as datetime
from multiprocessing import Pool as multithread
from itertools import permutations as permutations

#ctbBio
from ctbBio.search import search as search
from ctbBio.numblast import best as gethits
from ctbBio.fasta import iterate_fasta as parse_fasta
from ctbBio.rec_best_blast import rec_hits as rec_best_blast

# global variables
threads = 4
e = float(0.01)
bit = float(40)
length = float(.65)
algorithm = 'usearch'
genes = {} # info and search results for each gene and if match is a rec. best match (under threshold)
g2index = {} # index value for each gene ID

def log(mode, fastas):
    print('### running: %s' % (sys.argv[0]))
    print('### start time: %s' % (datetime.now()))
    print('### mode: %s' % (mode))
    print('### algorithm: %s' % (algorithm))
    print('### evalue threshold: %s' % (e))
    print('### bit score threshold: %s' % (bit))
    print('### alignment length threshold (currently only for global mode): %s' % (length))
    print('### comparing files: \n ... %s \n' % ('\n ... '.join(fastas)))

def find_genes(fastas):
    index = 0
    for fasta in fastas:
        previous = 0
        for sequence in parse_fasta(fasta):
            header = sequence[0].split('>')[1].split(' ', 1)
            id = header[0]
            if len(header) > 1:
                description = header[1]
            else:
                description = 'n/a'
            genes[id] = [[id, fasta, description, previous], {}]
            g2index[id] = index
            index += 1
            previous = id

def run_search(comparison):
    query, database = comparison[0], comparison[1]
    return search(query, database, algorithm, max_hits = 5)

def search_fastas(fastas, mode):
    if mode == 'global':
        comparisons = [comparison for comparison in permutations(fastas, 2)]
    elif mode == 'reference':
        comparisons = [[fastas[0], fasta] for fasta in fastas if fasta != fastas[0]]
        reverse = [[i[1], i[0]] for i in comparisons]
        comparisons.extend(reverse)
    search_output = [] # list of blast or usearch output files
    pool = multithread(threads)
    search_output.append(pool.map(run_search, comparisons))
    pool.close()
    pool.join()
    return search_output

def rec_hits(outfile):
    for hit in rec_best_blast(outfile, evalue = e, bit = bit):
        hit = hit.split('\t')
        query, match = hit[0].split()[0], hit[1]
        length, pidnet, evalue, bitscore = hit[3], hit[2], hit[10], hit[11]
        info = [length, pidnet, evalue, bitscore]
        genes[query][1][match] = [1, info]

def global_orthologs(fastas):
    from neto import neto as neto
    import networkx as nx
    graph = neto(fastas, algorithm = algorithm, e = e, bit = bit, length = length, norm_bit = False)
    results = []
    for group in nx.connected_components(graph):
        f2g = {}
        r = []
        for g in group:
            f2g[genes[g][0][1]] = g
            for other in group:
                if other != g:
                    if other in graph[g]:
                        s = graph[g][other]
                        scores = [s['length_fraction'], s['percent_id'], s['e_value'], s['bit_score']]
                        scores = [str(i) for i in scores]
                        rbh = 1
                    else:
                        rbh = 0
                        scores = ['-', '-', '-', '-']
                    if other not in genes[g][1]:
                        genes[g][1][other] = [rbh, []]
                    genes[g][1][other][1] = scores
        for fasta in fastas:
            if fasta not in f2g:
                r.append('-')
            else:
                r.append(f2g[fasta])
        results.append('*'.join(r))
    return set(results)

def reference_orthologs(fastas):
    results = []
    for gene in genes:
        gene_info = genes[gene][0]
        gene_matches = genes[gene][1]
        fasta_index = fastas.index(gene_info[1])
        query = [gene, fasta_index]
        orthologs = ['-' for i in fastas]
        if len(gene_matches) == 0:
            orthologs[fasta_index] = gene
        elif fasta_index == 0:
            for match in gene_matches:
                if gene_matches[match][0] == 1:
                    gene_info = genes[match][0]
                    fasta_index = fastas.index(gene_info[1])
                    orthologs[query[1]] = query[0]
                    orthologs[fasta_index] = match
        if orthologs != ['-' for i in fastas]:
            results.append('*'.join(orthologs))
    return set(results)

def get_scores(match, matches):
    scores = []
    for hit in matches:
        if hit == match:
            scores.append('!')
        elif hit == '-':
            scores.append('-')
        elif hit in genes[match][1]:
            scores.append(' | '.join(genes[match][1][hit][1]))
        else:
            scores.append('-')
    scores = ' ** '.join(scores)
    return scores

def find_previous(matches, fasta):
    matches = [match for match in matches if match != '-']
    candidates = []
    for match in matches:
        previous = '-'
        while previous == '-':
            other_previous = genes[match][0][3]
            if other_previous == 0:
                previous = [100, match]
            else:
                for hit in genes[other_previous][1]:
                    if other_previous in genes[hit][1]:
                        if genes[hit][0][1] == fasta and genes[hit][1][other_previous][0] == 1:
                            score = genes[hit][1][other_previous][1][2]
                            previous = [float(score), hit]
            match = other_previous
        candidates.append(previous)
    previous = sorted(candidates)[0][1]
    return previous

def find_sort_id(matches):
    id = []
    index = 0
    for gene in matches:
        if gene == '-':
            fasta = fastas[index]
            id.append(g2index[find_previous(matches, fasta)])
        else:
            id.append(g2index[gene])
        index += 1
    return id
#    return ''.join(id)

def format_results(results):
    formatted = []
    for matches in results:
        matches = matches.split('*')
        format = []
        sort_id = find_sort_id(matches)
        for match in matches:
            if match != '-':
                description = genes[match][0][2]
                scores = get_scores(match, matches)
                format.append('\t'.join([match, scores, description]))
            else:
                format.append('\t'.join(['-', '-', '-']))
        formatted.append(sort_id + ['\t'.join(format)])
    return sorted(formatted)

def print_results(results, fastas):
    header = ['%s\tscores\tdescription' % fasta for fasta in fastas]
    results = format_results(results)
    print('### end time: %s \n' % (datetime.now()))
    print('### output: \n')
    print('# %s' % ('\t'.join(header)))
    for result in results:
        print(result[-1])

def orthologer(mode, fastas):
    find_genes(fastas) # find all genes in all fasta files
    if mode == 'global':
        results = global_orthologs(fastas)
    if mode == 'reference':
        search_output = search_fastas(fastas, mode)[0] 
            # search fasta files against one another and get a list of output files
        rec_hits(search_output) 
            # find rec. best blast hits under e-value threshold
        results = reference_orthologs(fastas)
    print_results(results, fastas)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print('please specify \'global\' or \'reference\' mode and the fasta files to compare')
        exit()
    mode = sys.argv[1]
    fastas = sys.argv[2:]
    if mode != 'global' and mode != 'reference':
        print('please specify \'global\' or \'reference\' mode and the fasta files to compare')
        exit()
    log(mode, fastas)
    orthologer(mode, fastas)
