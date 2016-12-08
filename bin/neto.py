#!/usr/bin/env python3

"""
chris brown
banfield lab
ctb@berkeley.edu

draft script

script for building reciprocal best blast hit networks
"""

import sys, os
import itertools
import networkx as nx
from glob import glob as glob
import numpy
#from dh_hist import histogram as histogram # this is a modified version of the histogram module
from operator import itemgetter

# my scripts
from numblast import best as best_blast
from search import search as search 
from fasta import iterate_fasta as parse_fasta

def normalize_bit(A, B, bit, id2desc):
    """
    normalize the bit score:
    normalization factor = average max bit score for the two ORFs
    normalized = bit score / normalization factor
    """
    Amax, Bmax = id2desc[A][-1], id2desc[B][-1]
    norm_factor = float(numpy.average([Amax, Bmax]))
    normalized = bit / norm_factor
    return normalized

def get_descriptions(fastas):
    """
    get the description for each ORF 
    """
    id2desc = {}
    for fasta in fastas:
        for seq in parse_fasta(fasta):
            header = seq[0].split('>')[1].split(' ')
            id = header[0]
            if len(header) > 1:
                desc = ' '.join(header[1:])
            else:
                desc = 'n/a'
            length = float(len([i for i in seq[1].strip() if i != '*']))
            id2desc[id] = [fasta, desc, length]
    return id2desc

def gene_summary(cluster, gene, g, id2desc):
    bits = []
    norm_bits = []
    for hit in g[gene]:
        scores = g[gene][hit]
        bit, norm_bit = scores['bit_score'], scores['norm_bit']
        bits.append(bit)
        norm_bits.append(norm_bit)
    if len(bits) == 1:
        average, min_bit, max_bit = bits[0], bits[0], bits[0]
        std = 0
        average_norm = norm_bits[0]
    elif len(bits) == 0:
        average, min_bit, max_bit = 0, 0, 0
        std = 0
        average_norm = 0
    else:
        average = numpy.average(bits)
        std = numpy.std(bits)
        min_bit = min(bits)
        max_bit = max(bits)
        average_norm = numpy.average(norm_bits)
    fasta, desc = id2desc[gene][0:2]
    out = [cluster, fasta, gene, desc, len(bits), average, min_bit, max_bit, std, average_norm]
    out = [str(i) for i in out]
    return out

def cluster_stats(cluster, group, g, id2desc):
    bits = []
    for gene in group:
        for match in g[gene]:
            bit = [cluster, id2desc[gene][0], g[gene][match]['bit_score']]
            bits.append(bit)
    for bit in sorted(bits, key = itemgetter(2), reverse = True):
        print('\t'.join([str(i) for i in bit]))

def print_summary(g, fastas, id2desc, file_name):
    header = ['#cluster', 'fasta', 'gene id', 'gene', \
            'number of connections', 'average bit', 'min bit', 'max bit', \
            'std bit', 'average norm bit']
    out = open(file_name, 'w')
    print('\t'.join(header), file=out)
    cluster = 0
    for group in nx.connected_components(g):
        cluster += 1
        fasta2genes = {}
        for fasta in fastas:
            fasta2genes[fasta] = []
        for gene in group:
            fasta = id2desc[gene][0]
            fasta2genes[fasta].append(gene)
        for fasta in fasta2genes:
            genes = set(fasta2genes[fasta])
            for gene in genes:
                print('\t'.join(gene_summary(cluster, gene, g, id2desc)), file=out)
        print('', file=out)
    out.close()
#            if len(genes) > 1:
#                cluster_stats(cluster, group, g, id2desc)
#                break
#                for gene in genes:
#                    print '\t'.join(gene_summary(cluster, gene, g, id2desc))
#        print ''

def print_network_matrix(g, fastas, id2desc, file_name):
    out = open(file_name, 'w')
    fastas = sorted(fastas)
    print('# group\t%s' % ('\t'.join(fastas)), file=out)
    cluster = 1
    for group in nx.connected_components(g):
        line = ['%s: %s' % (cluster, ' | '.join(group))]
        hit_fastas = {}
        for gene in group:
            fasta = id2desc[gene][0]
            if fasta not in hit_fastas:
                hit_fastas[fasta] = 0
            hit_fastas[fasta] += 1
        for genome in fastas:
            if genome in hit_fastas:
                line.append(str(hit_fastas[genome]))
            else:
                line.append('0')
        cluster += 1
        print('\t'.join(line), file=out)
    out.close()

def print_genome_matrix(hits, fastas, id2desc, file_name):
    """
    optimize later? slow ...
    should combine with calculate_threshold module
    """
    out = open(file_name, 'w')
    fastas = sorted(fastas)
    print('## percent identity between genomes', file=out)
    print('# - \t %s' % ('\t'.join(fastas)), file=out)
    for fasta in fastas:
        line = [fasta]
        for other in fastas:
            if other == fasta:
                average = '-'
            else:
                average = numpy.average([hits[fasta][other][i][3] for i in hits[fasta][other]])
            line.append(str(average))
        print('\t'.join(line), file=out)
    print('', file=out)
    print('## percent of orfs that are orthologous between genomes', file=out)
    print('# - \t %s' % ('\t'.join(fastas)), file=out)
    for fasta in fastas:
        line = [fasta]
        for other in fastas:
            if other == fasta:
                percent = '-'
            else:
                orthologs = float(len(hits[fasta][other]))
                orfs = float(len([i for i in id2desc if id2desc[i][0] == fasta]))
                percent = float(orthologs / orfs) * 100
            line.append(str(percent))
        print('\t'.join(line), file=out)

def extra_connections(group, id2desc):
    connected = []
    for gene in group:
        fasta = id2desc[gene][0]
        if fasta in connected:
            return True
        else:
            connected.append(fasta)
    return False

def break_group(g, group, id2desc):
    # find most connected gene in connections (for one fasta)
    # remove all other genes from that fasta from that group
    original = [i for i in group]
    extra = True
    while extra is True:
        fastas = {}
        for gene in group:
            fa = id2desc[gene][0]
            if fa not in fastas:
                fastas[fa] = []
            fastas[fa].append(gene)
        excess = []
        for genes in list(fastas.values()):
            if len(genes) > 1:
                excess.extend(genes)
        counts = []
        for gene in excess:
            count = len([match for match in g[gene] if match in group])
            counts.append([gene, count])
        counts = sorted(counts, key = itemgetter(1), reverse = True)
        most_connected, remove = counts[0][0], [i[0] for i in counts[1:]]
#        if len(remove) == 0:
#            print counts
#            remove = most_connected
        new_group = []
        for gene in group:
            if gene not in remove:
                new_group.append(gene)
        group = [i for i in new_group]
        extra = extra_connections(group, id2desc)
    return group, [i for i in original if i not in group]

def update_graph(g, old, new):
    remove = []
    for gene in new:
        for match in g[gene]:
            if match not in new:
                remove.append(' '.join(sorted([gene, match])))
    for r in set(remove):
        a, b = r.split()
        g.remove_edge(a, b)
    return g

def split_network(g, id2desc, file_name):
    for group in nx.connected_components(g):
        extra = extra_connections(group, id2desc)
        if extra is False:
            continue
        while extra is not False:
            new, removed = break_group(g, group, id2desc)
            g = update_graph(g, group, new)
            group = removed
            extra = extra_connections(group, id2desc)
    nx.write_edgelist(g, file_name, delimiter = '\t', data = ['match_type', 'length_fraction', 'percent_id', 'e_value', 'bit_score', 'norm_bit'])
    return g

def self_compare(fastas, id2desc, algorithm):
    """
    compare genome to self to get the best possible bit score for each ORF
    """
    for fasta in fastas:
        blast = open(search(fasta, fasta, method = algorithm, alignment = 'local'))
        for hit in best_blast(blast, 1):
            id, bit = hit[0].split()[0], float(hit[-1])
            id2desc[id].append(bit)
    return id2desc

def compare_genomes(fastas, id2desc, algorithm):
    """
    make all pairwise genome comparisons to get rbhs
    rbh = {genome: {other_genome: {rbh_id: [rbh_match, pident, e, bit, norm_bit]]}}
    """
    hits = {}
    for pair in itertools.permutations(fastas, 2):
        if pair[0] not in hits:
            hits[pair[0]] = {}
        hits[pair[0]][pair[1]] = {}
    for pair in itertools.combinations(fastas, 2):
        A, B = pair
        blastF = open(search(A, B, method = algorithm, alignment = 'local'))
        blastR = open(search(B, A, method = algorithm, alignment = 'local'))
        for compare in [A, B, blastF], [B, A, blastR]:
            query, ref, blast = compare
            for hit in best_blast(blast, 1):
                a, b, pident, qstart, qend, e, bit = \
                        hit[0].split()[0], hit[1].split()[0], float(hit[2]), float(hit[6]), float(hit[7]), float(hit[-2]), float(hit[-1])
                alignment_length = abs(qstart - qend + 1)
                norm_bit = normalize_bit(a, b, bit, id2desc)
                length_fraction = float((alignment_length) / id2desc[a][2])
                hits[query][ref][a] = ['fbh', a, b, pident, length_fraction, e, bit, norm_bit]
    return hits

def find_rbh(hits, id2desc):
    rbh = {}
    for a in hits:
        rbh[a] = {}
        for b in hits[a]:
            rbh[a][b] = {}
    for genome in hits:
        for compare in hits[genome]:
            for orf in hits[genome][compare]:
                type = hits[genome][compare][orf][0]
                match = hits[genome][compare][orf][2]
                if type == 'fbh' and match in hits[compare][genome]:
                    if hits[compare][genome][match][2] == orf:
                        hits[genome][compare][orf][0] = 'rbh'
                        hits[compare][genome][match][0] = 'rbh'
#                        # for rbh, use the scores from the longest orf to represent the rbh
#                        if id2desc[orf][2] > id2desc[match][2]:
#                            scores = hits[genome][compare][orf]
#                        elif id2desc[match][2] > id2desc[orf][2]:
#                            scores = hits[compare][genome][match]
#                        else:
#                            scores = sorted([hits[genome][compare][orf], hits[compare][genome][match]], key = itemgetter(-1), reverse = True)[0]
#                        rbh[genome][compare][orf] = rbh[compare][genome][match] = scores
                        rbh[genome][compare][orf] = hits[genome][compare][orf]
                        rbh[compare][genome][match] = hits[compare][genome][match]
    return hits, rbh

def calc_thresholds(rbh, file_name, thresholds = [False, False, False, False], stdevs = 2):
    """
    if thresholds are not specififed, calculate based on the distribution of normalized bit scores
    """
    calc_threshold = thresholds[-1]
    norm_threshold = {}
    for pair in itertools.permutations([i for i in rbh], 2):
        if pair[0] not in norm_threshold:
            norm_threshold[pair[0]] = {}
        norm_threshold[pair[0]][pair[1]] = {}
    out = open(file_name, 'w')
    print('#### summary of rbh comparisons\n', file=out)
    comparisons = []
    for genome in rbh:
        for compare in rbh[genome]:
            pair = ''.join(sorted([genome, compare]))
            if pair in comparisons:
                continue
            comparisons.append(pair)
            scores = {'percent identity': [], 'e-value': [], 'bit score': [], 'normalized bit score': [], 'alignment length fraction': []}
            print('### blast between %s and %s\n' % (genome, compare), file=out)
            for id in rbh[genome][compare]:
                pident, length_fraction, e, bit, norm_bit = rbh[genome][compare][id][3:]
                scores['percent identity'].append(pident)
                scores['alignment length fraction'].append(length_fraction)
                scores['e-value'].append(e)
                scores['bit score'].append(bit)
                scores['normalized bit score'].append(norm_bit)
            if calc_threshold is True:
                norms = scores['normalized bit score']
                average = numpy.average(norms) 
                std = numpy.std(norms)
                normal_thresh = average - (std * stdevs)
                print('## average normalized bit score: %s' % average, file=out)
                print('## standard deviation of normalized bit scores: %s' % std, file=out)
                print('## normalized bit score threshold set to: %s\n' % (normal_thresh), file=out)
                norm_threshold[genome][compare], norm_threshold[compare][genome] = normal_thresh, normal_thresh
            for score in scores:
                print('## %s' % (score), file=out)
                if len(scores[score]) > 0:
                    print('## average: %s' % numpy.average(scores[score]), file=out)
#                    hist = histogram(scores[score], [])
#                    for line in hist:
#                        print >> out, line
                print('', file=out)
    out.close()
    if calc_threshold is True:
        return thresholds[0:-1] + [norm_threshold]
    else:
        return thresholds

def rbh_network(id2desc, rbh, file_name, thresholds = [False, False, False, False]):
    """
    make the network based on rbhs and score thresholds
    """
    g = nx.Graph() # network graph for storing rbhs
    filtered = {}
    e_thresh, bit_thresh, length_thresh, norm_thresh = thresholds
    for genome in rbh:
        filtered[genome] = {}
        for other in rbh:
            if other != genome:
                filtered[genome][other] = {}
    comparisons = []
    for genome in rbh:
        for compare in rbh[genome]:
            pair = ''.join(sorted([genome, compare]))
            if pair in comparisons: # make sure you only have to make rbh comparison once
                continue
            comparisons.append(pair)
            for orf in rbh[genome][compare]:
                scoresA = rbh[genome][compare][orf]
                match = scoresA[2]
                if match in rbh[compare][genome]:
                    scoresB = rbh[compare][genome][match]
                else:
                    scoresB = scoresA
                typeA, AA, BA, pidentA, lengthA, eA, bitA, norm_bitA = scoresA
                typeB, AB, BB, pidentB, lengthB, eB, bitB, norm_bitB = scoresB
                if norm_thresh is not False:
                    if norm_bitA < norm_thresh[genome][compare] \
                        or norm_bitB < norm_thresh[compare][genome]:
                        continue
                if e_thresh is not False:
                    if eA > e_thresh \
                        or eB > e_thresh:
                        continue
                if bit_thresh is not False:
                    if bitA < bit_thresh \
                        or bitB < bit_thresh:
                        continue
                if length_thresh is not False:
                    if lengthA < length_thresh \
                        or lengthB < length_thresh:
                        continue
                if id2desc[orf][2] > id2desc[match][2]:
                    scores = scoresA
                elif id2desc[orf][2] < id2desc[match][2]:
                    scores = scoresB
                else:
                    scores = sorted([scoresA, scoresB], key = itemgetter(-1), reverse = True)[0]
                type, A, B, pident, length, e, bit, norm_bit = scores
                g.add_edge(A, B, match_type = type, length_fraction = length, \
                        percent_id = pident, e_value = e, bit_score = bit, norm_bit = norm_bit)
                filtered[genome][compare][orf] = scoresA
                filtered[compare][genome][match] = scoresB
    missing = set([i for i in id2desc]).difference(set([i for i in g]))
    for orf in missing:
        g.add_edge(orf, orf, percent_id = 0, e_value = 0, bit_score = 0, norm_bit = 0, \
                    length_fraction = 0)
    nx.write_edgelist(g, file_name, delimiter = '\t', data = ['match_type', 'length_fraction', 'percent_id', 'e_value', 'bit_score', 'norm_bit'])
    return g, filtered

def neto(fastas, algorithm = 'usearch', e = 0.01, bit = 40, length = .65, norm_bit = False):
    """
    make and split a rbh network
    """
    thresholds = [e, bit, length, norm_bit]
    id2desc = get_descriptions(fastas)
            # get [fasta, description, length] for ORF id
    id2desc = self_compare(fastas, id2desc, algorithm)
            # get best possible bit score for each ORF 
            # (comparing with itself) [fasta, description, length, bestbit]
    hits = compare_genomes(fastas, id2desc, algorithm)
            # pair wise genome comparisons {genome: {id: [match_type = 'rbh' or 'fbh', scores]}}
    calc_thresholds(hits, file_name = 'fbh.scores.summary.txt')
    rbh_network(id2desc, hits, file_name = 'fbh.network.edges.txt')
    hits, rbh = find_rbh(hits, id2desc)
            # remove hits that are not reciprocal best blast hits
    thresholds = calc_thresholds(rbh, 'rbh.scores.summary.txt', thresholds)
            # print rbh score summary to rbh_score_summary.txt and
            # calculate normalized bit score cutoff for each pair of
            # genomes, if desired
    g = rbh_network(id2desc, rbh, file_name = 'rbh.network.edges.txt')
    filtered_g, filtered_rbh = rbh_network(id2desc, rbh, 'rbh.filtered.network.edges.txt', thresholds)
    calc_thresholds(filtered_rbh, file_name = 'rbh.filtered.scores.summary.txt')
    print_summary(filtered_g, fastas, id2desc, file_name = 'rbh.filtered.network.nodes.txt')
    print_network_matrix(filtered_g, fastas, id2desc, file_name = 'rbh.filtered.network.matrix.txt')
    print_genome_matrix(filtered_rbh, fastas, id2desc, file_name = 'rbh.filtered.network.genome_matrix.txt')
    split_g = split_network(filtered_g, id2desc, file_name = 'rbh.filtered.split.network.edges.txt')
    print_summary(split_g, fastas, id2desc, file_name = 'rbh.filtered.split.network.nodes.txt')
    print_network_matrix(split_g, fastas, id2desc, file_name = 'rbh.filtered.split.network.matrix.txt')
    return split_g

def network(thresholds, fastas, algorithm = 'usearch'):
    """
    make a rbh network for all pair-wise genome comparisons
    - filter network based on normalized bit score (default, automatic) or specified e-value / bit score
    - evaluate the scores for each genome pair compared
    - build second network that is filtered so that clusters have only one ORF per genome
    - evaluate the scores in this network
        - consider including fbhs under threshold if they don't violate cluster
    - compare 'connected-ness' for each genome
    """
    id2desc = get_descriptions(fastas)
            # get [fasta, description, length] for ORF id
    id2desc = self_compare(fastas, id2desc, algorithm)
            # get best possible bit score for each ORF 
            # (comparing with itself) [fasta, description, length, bestbit]
    hits = compare_genomes(fastas, id2desc, algorithm)
            # pair wise genome comparisons {genome: {id: [match_type = 'rbh' or 'fbh', scores]}}
    calc_thresholds(hits, file_name = 'fbh.scores.summary.txt')
    rbh_network(id2desc, hits, file_name = 'fbh.network.edges.txt')
    hits, rbh = find_rbh(hits, id2desc)
            # remove hits that are not reciprocal best blast hits
    thresholds = calc_thresholds(rbh, 'rbh.scores.summary.txt', thresholds)
            # print rbh score summary to rbh_score_summary.txt and
            # calculate normalized bit score cutoff for each pair of
            # genomes, if desired
    g = rbh_network(id2desc, rbh, file_name = 'rbh.network.edges.txt')
    filtered_g, filtered_rbh = rbh_network(id2desc, rbh, 'rbh.filtered.network.edges.txt', thresholds)
    calc_thresholds(filtered_rbh, file_name = 'rbh.filtered.scores.summary.txt')
    print_summary(filtered_g, fastas, id2desc, file_name = 'rbh.filtered.network.nodes.txt')
    print_network_matrix(filtered_g, fastas, id2desc, file_name = 'rbh.filtered.network.matrix.txt')
    print_genome_matrix(filtered_rbh, fastas, id2desc, file_name = 'rbh.filtered.network.genome_matrix.txt')
    split_g = split_network(filtered_g, id2desc, file_name = 'rbh.filtered.split.network.edges.txt')
    print_summary(split_g, fastas, id2desc, file_name = 'rbh.filtered.split.network.nodes.txt')
    print_network_matrix(split_g, fastas, id2desc, file_name = 'rbh.filtered.split.network.matrix.txt')
    return split_g

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print('* specify evalue threshold (or -), bit score threshold (or -), length fraction threshold (or -), normalized bit score threshold (True or False), and fastas (list or glob)')
        print('* every ORF must have a unique identifier') 
        exit()
    evalue, bit, length, norm_bit, fastas = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:]
    if norm_bit == 'True':
        norm_bit = True
    else:
        norm_bit = False
    thresholds = []
    for i in evalue, bit, length:
        if i == '-':
            thresholds.append(False)
        else:
            thresholds.append(float(i))
    thresholds = thresholds + [norm_bit]
    if len(fastas) == 1:
        fastas = glob(fastas[0])
    network(thresholds, fastas)
