#!/usr/bin/env python3

"""
* script for finding set of 16 syntenous ribosomal proteins
    in ORFs predicted for scaffolds (Prodigal)
* evaluates each scaffold independently
    (does not expect that the input as that from a single genome)
"""

import sys
import os
import argparse
from random import choice
from operator import itemgetter

# ctb
from ctbBio.search import search as search
from ctbBio.numblast import best as numblast
from ctbBio.fasta import iterate_fasta as parse_fasta

def find_databases(databases):
    """
    define ribosomal proteins and location of curated databases
    """
    # 16 ribosomal proteins in their expected order
    proteins = ['L15', 'L18', 'L6', 'S8', 'L5', 'L24', 'L14',
            'S17', 'L16', 'S3', 'L22', 'S19', 'L2', 'L4', 'L3', 'S10']
    # curated databases
    protein_databases = {
                'L14': 'rpL14_JGI_MDM.filtered.faa',
                'L15': 'rpL15_JGI_MDM.filtered.faa',
                'L16': 'rpL16_JGI_MDM.filtered.faa',
                'L18': 'rpL18_JGI_MDM.filtered.faa',
                'L22': 'rpL22_JGI_MDM.filtered.faa',
                'L24': 'rpL24_JGI_MDM.filtered.faa',
                'L2': 'rpL2_JGI_MDM.filtered.faa',
                'L3': 'rpL3_JGI_MDM.filtered.faa',
                'L4': 'rpL4_JGI_MDM.filtered.faa',
                'L5': 'rpL5_JGI_MDM.filtered.faa',
                'L6': 'rpL6_JGI_MDM.filtered.faa',
                'S10': 'rpS10_JGI_MDM.filtered.faa',
                'S17': 'rpS17_JGI_MDM.filtered.faa',
                'S19': 'rpS19_JGI_MDM.filtered.faa',
                'S3': 'rpS3_JGI_MDM.filtered.faa',
                'S8': 'rpS8_JGI_MDM.filtered.faa'}
    protein_databases = {key: '%s/%s' % (databases, database) \
            for key, database in list(protein_databases.items())}
    return proteins, protein_databases

def scaffold_hits(searches, fasta, max_hits):
    """
    get hits from each search against each RP
    scaffolds[scaffold] = # ORfs
    s2rp[scaffold] = {rp:[hits]}
    """
    # initialize
    ## scaffolds[scaffold] = # ORFs
    scaffolds = {}
    for seq in parse_fasta(fasta):
        scaffold = seq[0].split()[0].split('>', 1)[1].rsplit('_', 1)[0]
        if scaffold not in scaffolds:
            scaffolds[scaffold] = 0
        scaffolds[scaffold] += 1
    s2rp = {s: {r[0]: []
        for r in searches}
        for s in scaffolds}
    # get hits from blast
    for search in searches:
        rp, blast = search
        hits = [i for i in numblast(open(blast), max_hits, evalue_thresh, bit_thresh)]
        for hit in hits:
            s = hit[0].split()[0].rsplit('_', 1)[0]
            hit[10], hit[11] = float(hit[10]), float(hit[11])
            s2rp[s][rp].append(hit)
    return scaffolds, s2rp

def find_next(start, stop, i2hits):
    """
    which protein has the best hit, the one to the 'right' or to the 'left?'
    """
    if start not in i2hits and stop in i2hits:
        index = stop
    elif stop not in i2hits and start in i2hits:
        index = start
    elif start not in i2hits and stop not in i2hits:
        index = choice([start, stop])
        i2hits[index] = [[False]]
    else:
        A, B = i2hits[start][0], i2hits[stop][0]
        if B[10] <= A[10]:
            index = stop
        else:
            index = start
    if index == start:
        nstart = start - 1
        nstop = stop
    else:
        nstop = stop + 1
        nstart = start
    match = i2hits[index][0]
    rp = match[-1]
    return index, nstart, nstop, rp, match

def find_block(rps, num_orfs, hits, best, max_errors):
    best_index, best_rp = int(best[0].split()[0].rsplit('_', 1)[1]), best[-1]
    index2hits = {i: [] for i in range(1, num_orfs + 1)}
    errors = 0
    found = [best_rp]
    block = {best_index: best}
    for rp, matches in list(hits.items()):
        for match in matches:
            i = int(match[0].split()[0].rsplit('_', 1)[1])
            index2hits[i].append(match + [rp])
    for i, matches in list(index2hits.items()):
        if matches == []:
            del index2hits[i]
            continue
        index2hits[i] = sorted(matches, key = itemgetter(10))
    index, start, stop, rp = best_index, best_index - 1, best_index + 1, best_rp
    while errors < max_errors and len(set(found)) < 17: # and index <= num_orfs:
        new_index, start, stop, rp, match = find_next(start, stop, index2hits)
        if rp is False: # count as error if there is no hit
            errors += 1
        index = new_index
        found.append(rp)
        block[index] = match
    trimmed_block = {}
    for i in sorted(block):
        if block[i] is not False:
            out = [str(j) for j in block[i]]
            trimmed_block[out[-1]] = out
    return trimmed_block

def find_ribosomal(rps, scaffolds, s2rp, min_hits, max_hits_rp, max_errors):
    """
    determine which hits represent real ribosomal proteins, identify each in syntenic block
    max_hits_rp = maximum number of hits to consider per ribosomal protein per scaffold
    """
    for scaffold, proteins in list(s2rp.items()):
        # for each scaffold, get best hits for each rp
        hits = {p: [i for i in sorted(hits, key = itemgetter(10))][0:max_hits_rp]
            for p, hits in list(proteins.items()) if len(hits) > 0}
        # skip if fewer than min_hits RPs are identified
        if len(hits) < min_hits:
            continue
        best = sorted([hit[0] + [p]
            for p, hit in list(hits.items())], key = itemgetter(10))[0]
        block = find_block(rps, scaffolds[scaffold], hits, best, max_errors)
        if (len(block) - 1) >= min_hits:
            yield scaffold, block

def ribosomal(scaffolds, DBdir, min_hits, evalue_thresh, bit_thresh, \
                method = 'usearch', threads = 6, \
                max_hits = 1, max_hits_rp = 1, max_errors = 35):
    """
    find ribosomal proteins
    max_hits = maximum number of blast hits to consider for an orf
                  if 1, only consider best blast hit for each ORF
    max_hits_rp = maximum number of hits to consider per ribosomal protein per scaffold
                     if 1, only consider best RP match to contig
    max_errors = maximum number of errors when looking for block of proteins (e.g. out of order or gap)
    """
    # rps = list (in syntenic order) of ribosomal proteins
    # rp_db = dictionary to find the database files
    rps, rp_db = find_databases(DBdir)
    searches = [[rp, search(scaffolds, rp_db[rp], method = method, threads = str(threads), max_hits = 10)]
            for rp in rp_db]
    scaffolds, scaffold2rp = scaffold_hits(searches, scaffolds, max_hits)
    print('# scaffold\t%s' % ('\t'.join(rps)))
    for scaffold, block in \
            find_ribosomal(rps, scaffolds, scaffold2rp, min_hits, max_hits_rp, max_errors):
        id_rps = []
        for rp in rps:
            if rp in block:
                id_rps.append(block[rp][0].split()[0])
            else:
                id_rps.append('-')
        print('%s\t%s' % (scaffold, '\t'.join(id_rps)))


if __name__ == '__main__':
    desc = '# find syntenic group of 16 ribosomal proteins'
    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument(\
            '-f', required = True, type = str, \
            help = 'ORF predictions (Prodigal-format fasta)')
    parser.add_argument(\
            '-d', required = False, default = False, type = str, \
            help = 'directory with ribosomal protein databases \
            (default = check for databases env. variable)')
    parser.add_argument(\
            '-m', required = False, default = 3, type = int, \
            help = 'min. # of RPs to include as match (default = 3)')
    parser.add_argument(\
            '-e', required = False, default = float(1e-6), type = float, \
            help = 'maximum evalue to consider as hit (default = 1e-6)')
    parser.add_argument(\
            '-a', required = False, default = 'usearch', type = str, \
            help = 'algorithm: usearch (default), usearch-cluster, blast')
    parser.add_argument(\
            '-t', required = False, default = 6, type = int, \
            help = 'threads')
    parser.add_argument(\
            '-b', required = False, default = float(40), type = float, \
            help = 'minimum bit score to consider as hit (default = 40)')
    args = vars(parser.parse_args())
    evalue_thresh, bit_thresh = args['e'], args['b']
    scaffolds, min_hits = args['f'], args['m']
    method = args['a']
    threads = args['t']
    DBdir = args['d']
    if 'databases' not in os.environ and DBdir is False:
        print('# specify databases directory', file = sys.stderr)
        exit()
    if DBdir is False:
        DBdir = '%s/rp16/Laura/' % (os.environ['databases'])
    ribosomal(scaffolds, DBdir, min_hits, evalue_thresh, bit_thresh, method = method, threads = threads)
