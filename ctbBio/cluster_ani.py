#!/usr/bin/env python3

"""
script for clustering genomes based on average nucleotide
identity 
"""

import os
import sys
import argparse
import subprocess
from networkx import Graph as Graph
from networkx import connected_components as connected_components

# ctb
from fasta import iterate_fasta as parse_fasta

def make_mashes(fastas, mash_file, threads, kmer = 21, force = False):
    """
    Create mash files for multiple fasta files
    Input:
        fastas <list[str]>  -- paths to fasta files
        mash_file <str>     -- path to output mash file
        threads <int>       -- # threads for parallelization
        kmer <int>          -- kmer size for mash sketching
        force <boolean>     -- force overwrite of all mash files
    """
    mash_processes = set()
    sketches = [fasta + '.msh' for fasta in fastas]
    devnull = open(os.devnull, 'w')
    # Perform the sketching
    for fasta, sketch in zip(fastas, sketches):
        if os.path.isfile(sketch):
            continue
        mash_cmd = ['/opt/bin/bio/mash', 'sketch', '-o', fasta, '-k', str(kmer), fasta]
        mash_processes.add(subprocess.Popen(mash_cmd, stderr=devnull))
        if len(mash_processes) >= threads:
            os.wait()
            mash_processes.difference_update([mp for mp in mash_processes if mp.poll() is not None])
    # Collect stragglers
    for mp in mash_processes:
        if mp.poll() is None:
            mp.wait()
    # Paste sketches into single mash
    paste_mashes(sketches, mash_file, force = force)
    return

def paste_mashes(sketches, pasted_mash, force = False):
    """
    Combine mash files into single sketch
    Input:
        sketches <list[str]>  -- paths to sketch files
        pasted_mash <str>     -- path to output mash file
        force <boolean>     -- force overwrite of all mash file
    """
    if os.path.isfile(pasted_mash):
        if force:
            subprocess.Popen(['rm', pasted_mash]).wait()
        else:
            return
    pasted_mash = pasted_mash.rsplit('.msh')[0]
    mash_cmd = ['/opt/bin/bio/mash', 'paste', pasted_mash]
    mash_cmd.extend(sketches)
    process = subprocess.Popen(mash_cmd)
    process.wait()
    return

def ani(fastas, mash_file, sim_threshold, threads):
    """
    Use mash to estimate ANI of genomes
    Input:
        fastas <list[str]>      -- paths to fasta files
        sim_threshold <float>   -- fractional cutoff % identity for cluster joining
        mash_file <str>         -- pasted sketch file of all fastas being compared
        threads <int>           -- number threads for distance estimation
    """
    ANI = Graph()
    # use Mash to estimate ANI
    for fasta in fastas:
        indiv_mash = fasta + '.msh'
        if os.path.isfile(indiv_mash):
            cmp_file = indiv_mash
        else:
            cmp_file = fasta
        mash_cmd = ['/opt/bin/bio/mash', 'dist', cmp_file, mash_file]
        process = subprocess.Popen(mash_cmd, stdout = subprocess.PIPE)
        for pair in process.communicate()[0].splitlines():
            a, b, dist, p, shared = pair.decode().strip().split()
            a = a.rsplit('.', 1)[0].rsplit('/', 1)[-1].rsplit('.contigs')[0]
            b = b.rsplit('.', 1)[0].rsplit('/', 1)[-1].rsplit('.contigs')[0]
            p = float(p)
            similarity = (1 - float(dist)) * 100
            if similarity >= sim_threshold:
                ANI.add_edge(a, b, si = similarity, pval = p, sharedK = shared)
        process.wait()
    return ANI

def genome_info(genome, info):
    """
    return genome info for choosing representative
    
    if ggKbase table provided - choose rep based on SCGs and genome length
        - priority for most SCGs - extra SCGs, then largest genome

    otherwise, based on largest genome
    """
    try:
        scg       = info['#SCGs']
        dups      = info['#SCG duplicates']
        length    = info['genome size (bp)']
        return [scg - dups, length, genome]
    except:
        return [False, False, info['genome size (bp)'], genome]

def print_clusters(fastas, info, ANI):
    """
    choose represenative genome and 
    print cluster information

    *if ggKbase table is provided, use SCG info to choose best genome
    """
    header = ['#cluster', 'num. genomes', 'rep.', 'genome', '#SCGs', '#SCG duplicates', \
            'genome size (bp)', 'fragments', 'list']
    yield header
    in_cluster = []
    for cluster_num, cluster in enumerate(connected_components(ANI)):
        cluster = sorted([genome_info(genome, info[genome]) \
                            for genome in cluster], \
                            key = lambda x: x[0:], reverse = True)
        rep = cluster[0][-1]
        cluster = [i[-1] for i in cluster]
        size = len(cluster)
        for genome in cluster:
            in_cluster.append(genome)
            try:
                stats = [size, rep, genome, \
                            info[genome]['#SCGs'], info[genome]['#SCG duplicates'], \
                            info[genome]['genome size (bp)'], info[genome]['# contigs'], cluster]
            except:
                stats = [size, rep, genome, \
                            'n/a', 'n/a', \
                            info[genome]['genome size (bp)'], info[genome]['# contigs'], cluster]
            if rep == genome:
                stats = ['*%s' % (cluster_num)] + stats
            else:
                stats = [cluster_num] + stats
            yield stats 
    # print singletons
    try:
        start = cluster_num + 1
    except:
        start = 0
    fastas = set([i.rsplit('.', 1)[0].rsplit('/', 1)[-1].rsplit('.contigs')[0] for i in fastas])
    for cluster_num, genome in \
            enumerate(fastas.difference(set(in_cluster)), start):
        try:
            stats = ['*%s' % (cluster_num), 1, genome, genome, \
                        info[genome]['#SCGs'], info[genome]['#SCG duplicates'], \
                        info[genome]['genome size (bp)'], info[genome]['# contigs'], [genome]]
        except:
            stats = ['*%s' % (cluster_num), 1, genome, genome, \
                        'n/a', 'n/a', \
                        info[genome]['genome size (bp)'], info[genome]['# contigs'], [genome]]
        yield stats 

def to_int(i):
    """
    convert to integer, if possible
    """
    try:
        return int(i)
    except:
        return i

def parse_ggKbase_tables(tables, id_type):
    """
    convert ggKbase genome info tables to dictionary
    """
    g2info = {}
    for table in tables:
        for line in open(table):
            line = line.strip().split('\t')
            if line[0].startswith('name'):
                header = line
                header[4] = 'genome size (bp)'
                header[12] = '#SCGs'
                header[13] = '#SCG duplicates'
                continue
            name, code, info = line[0], line[1], line
            info = [to_int(i) for i in info]
            if id_type is False: # try to use name and code ID
                if 'UNK' in code or 'unknown' in code:
                    code = name
                if (name != code) and (name and code in g2info):
                    print('# duplicate name or code in table(s)', file=sys.stderr)
                    print('# %s and/or %s' % (name, code), file=sys.stderr)
                    exit()
                if name not in g2info:
                    g2info[name] = {item:stat for item, stat in zip(header, info)}
                if code not in g2info:
                    g2info[code] = {item:stat for item, stat in zip(header, info)}
            else:
                if id_type == 'name':
                    ID = name
                elif id_type == 'code':
                    ID = code
                else:
                    print('# specify name or code column using -id', file=sys.stderr)
                    exit()
                ID = ID.replace(' ', '')
                g2info[ID] = {item:stat for item, stat in zip(header, info)}
                if g2info[ID]['genome size (bp)'] == '':
                    g2info[ID]['genome size (bp)'] = 0
    return g2info

def parse_checkM_tables(tables):
    """
    convert checkM genome info tables to dictionary
    """
    g2info = {}
    for table in tables:
        for line in open(table):
            line = line.strip().split('\t')
            if line[0].startswith('Bin Id'):
                header = line
                header[8] = 'genome size (bp)'
                header[5] = '#SCGs'
                header[6] = '#SCG duplicates'
                continue
            ID, info = line[0], line
            info = [to_int(i) for i in info]
            ID = ID.replace(' ', '')
            g2info[ID] = {item:stat for item, stat in zip(header, info)}
            if g2info[ID]['genome size (bp)'] == '':
                g2info[ID]['genome size (bp)'] = 0
    return g2info

def genome_lengths(fastas, info):
    """
    get genome lengths
    """
    if info is False:
        info = {}
    for genome in fastas:
        name = genome.rsplit('.', 1)[0].rsplit('/', 1)[-1].rsplit('.contigs')[0]
        if name in info:
            continue
        length = 0
        fragments = 0
        for seq in parse_fasta(genome):
            length += len(seq[1])
            fragments += 1
        info[name] = {'genome size (bp)':length, '# contigs':fragments}
    return info 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = \
            '# cluster genomes based on average nucleotide identity (ani)')
    parser.add_argument(\
            '-f', nargs = '*', action = 'store', required = True, \
                help = 'fastas')
    parser.add_argument(\
            '-m', action = 'store', required = True, type = str, \
            help = 'mash file (will be created if it does not exist)')
    parser.add_argument(\
            '-s', default = 98, type = float, required = False, \
                help = 'percent similarity (default = 98)')
    parser.add_argument(\
            '-g', nargs = '*', action = 'store', required = False, \
                default = False, \
                help = 'ggKbase genome table for selecting representative (optional)')
    parser.add_argument(\
            '-c', nargs = '*', action = 'store', required = False, \
                default = False, \
                help = 'checkM genome table for selecting representative (optional)')
    parser.add_argument(\
            '-id', default = False, \
            help = 'use name or code column in ggKbase table (default: try both)')
    parser.add_argument(\
            '-t', required = False, default = 6, type = int, \
            help = 'threads (default = 6)')
    args = vars(parser.parse_args())
    fastas, similarity, id_type, threads, mash_file = \
            args['f'], args['s'], args['id'], args['t'], args['m']
    gg, cm = args['g'], args['c']
    if '.msh' not in mash_file:
        mash_file = '%s.msh' % (mash_file)
    info = False # assume no marker gene file is given (either ggKbase or checkM)
    if gg is not False:
        info = parse_ggKbase_tables(gg, id_type)
    elif cm is not False:
        info = parse_checkM_tables(cm)
    info = genome_lengths(fastas, info)
    make_mashes(fastas, mash_file, threads)
    ANI = ani(fastas, mash_file, similarity, threads)
    for genome in print_clusters(fastas, info, ANI):
        print('\t'.join([str(i) for i in genome]))
