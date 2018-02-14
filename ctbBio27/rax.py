#!/usr/bin/python2.7

"""
script for quickly running FastTree and raxml
"""

import sys
import os
import random
import argparse
from subprocess import Popen
import string as string_gen

from fasta import iterate_fasta as parse_fasta

def check_type(fasta):
    nucl = ['A', 'T', 'G', 'C']
    junk = ['N', 'U', '.', '-', ' ']
    type = 'nucl'
    for seq in parse_fasta(fasta):
        seq = seq[1].upper()
        for residue in seq:
            if residue in junk:
                continue
            if residue not in nucl:
                type = 'prot'
                break
        break
    return type

def check(file):
    """
    if a file exists, return 'True,' else, 'False'
    """
    try:
        open(file)
        return True
    except (OSError, IOError), e:
        return False

def remove_bad(string):
    """
    remove problem characters from string
    """
    remove = [':', ',', '(', ')', ' ', '|', ';', '\'']
    for c in remove:
        string = string.replace(c, '_')
    return string

def id_generator(size = 10, chars = string_gen.ascii_uppercase):
    return ''.join(random.choice(chars) for _ in range(size))

def get_ids(a):
    """
    make copy of sequences with short identifier
    """
    a_id = '%s.id.fa' % (a.rsplit('.', 1)[0])
    a_id_lookup = '%s.id.lookup' % (a.rsplit('.', 1)[0])
    if check(a_id) is True:
        return a_id, a_id_lookup
    a_id_f = open(a_id, 'w')
    a_id_lookup_f = open(a_id_lookup, 'w')
    ids = []
    for seq in parse_fasta(open(a)):
        id = id_generator() 
        while id in ids:
            id = id_generator() 
        ids.append(id)
        header = seq[0].split('>')[1]
        name = remove_bad(header)
        seq[0] = '>%s %s' % (id, header)
        print >> a_id_f, '\n'.join(seq)
        print >> a_id_lookup_f, '%s\t%s\t%s' % (id, name, header)
    return a_id, a_id_lookup

def convert2phylip(convert):
    """
    convert fasta to phylip because RAxML is ridiculous
    """
    out = '%s.phy' % (convert.rsplit('.', 1)[0])
    if check(out) is False:
        from Bio import AlignIO
        convert = open(convert, 'rU')
        out_f = open(out, 'w')
        alignments = AlignIO.parse(convert, "fasta")
        AlignIO.write(alignments, out, "phylip")
    return out

def run_fast(aligned, threads, cluster, node):
    """
    run FastTree
    """
    tree = '%s.fasttree.nwk' % (aligned.rsplit('.', 1)[0])
    if check(tree) is False:
        if 'FastTreeV' in os.environ:
            ft = os.environ['FastTreeV']
            os.environ['OMP_NUM_THREADS'] = str(threads)
        else:
            ft = 'FastTreeMP'
            os.environ['OMP_NUM_THREADS'] = str(threads)
        if check_type(aligned) == 'nucl':
            type = '-nt -gamma -spr 4 -mlacc 2 -slownni'
        else:
            type = '-spr 4 -mlacc 2 -slownni'
        dir = os.getcwd()
        command = 'cat %s/%s | cut -d \' \' -f 1 | %s -log %s/%s.log %s > %s/%s 2>>%s/%s.log' % \
                (dir, aligned, ft, dir, tree, type, dir, tree, dir, tree)
#        print '\t %s' % (command)
        if cluster is False:
            p = Popen(command, shell = True)
        else:
            if threads > 24:
                ppn = 24
            else:
                ppn = threads
            re_call = 'cd %s; %s --no-rax' % (dir.rsplit('/', 1)[0], ' '.join(sys.argv))
            if node is False:
                node = '1'
            qsub = 'qsub -l nodes=%s:ppn=%s -m e -N FastTree' % (node, ppn)
            p = Popen('echo "%s;%s" | %s' % (command, re_call, qsub), shell = True)
        p.communicate()
    return tree

def run_raxml(rax_out, boot, a_id_phylip, threads, aligned, model, cluster, node):
    """
    run raxml
    """
    if 'raxml' in os.environ:
        raxml = os.environ['raxml']
        threads = '-T %s' % (threads)
    else:
        raxml = 'raxml'
        threads = ''
    rax_tree = 'RAxML_bipartitions.%s' % (rax_out)
    if check(rax_tree) is False:
        seed = random.randint(123456789, 12345678910000000)
        print >> open('seed.txt', 'w'), seed
        if check_type(aligned) == 'nucl' and model is False:
            model = 'GTRCAT'
        elif model is False:
            model = 'PROTCATJTT'
        dir = os.getcwd()
        command = '%s -f a -m %s -n %s -N %s -s %s -x %s -p %s %s > %s.log 2>>%s.log' % \
                    (raxml, model, rax_out, boot, a_id_phylip, seed, seed, threads, rax_out, rax_out)
        if cluster is False:
            p = Popen(command, shell = True)
        else:
            if threads > 24:
                ppn = 24
            else:
                ppn = threads
            if node is False:
               node = '1'
            qsub = 'qsub -l nodes=%s:ppn=%s -m e -N raxml' % (node, ppn)
            command = 'cd /tmp; mkdir raxml_%s; cd raxml_%s; cp %s/%s .; %s; mv * %s/; rm -r ../raxml_%s' \
                    % (seed, seed, dir, a_id_phylip, command, dir, seed)
            re_call = 'cd %s; %s --no-fast' % (dir.rsplit('/', 1)[0], ' '.join(sys.argv))
            p = Popen('echo "%s;%s" | %s' % (command, re_call, qsub), shell = True)
        p.communicate()
    return rax_tree

def fix_tree(tree, a_id_lookup, out):
    """
    get the names for sequences in the raxml tree
    """
    if check(out) is False and check(tree) is True:
        tree = open(tree).read()
        for line in open(a_id_lookup):
            id, name, header = line.strip().split('\t')
            tree = tree.replace(id+':', name+':')
        out_f = open(out, 'w')
        print >> out_f, tree.strip()
    return out

def rax(a, boot, threads, fast = False, run_rax = False, model = False, cluster = False, node = False):
    """
    run raxml on 'a' (alignment) with 'boot' (bootstraps) and 'threads' (threads)
    store all files in raxml_a_b
    1. give every sequence a short identifier
    2. convert fasta to phylip
    3. run raxml
    4. convert ids in raxml tree to original names
    """
    a = os.path.abspath(a)
    a_base = a.rsplit('/', 1)[1]
    out_dir = '%s/%s_rax_boots_%s' % \
                (a.rsplit('/', 1)[0], a_base.rsplit('.', 1)[0], boot)
    os.system('mkdir -p %s' % (out_dir))
    os.system('ln -sf %s %s/%s' % (os.path.abspath(a), out_dir, a.rsplit('/', 1)[1]))
    os.chdir(out_dir)
    a_id, a_id_lookup = get_ids(a_base)
    a_id_phylip = convert2phylip(a_id)
    rax_out = '%s.raxml.txt' % (a_id_phylip)
    if fast is True:
        final_fast = '%s.fasttree.tree' % (a_id_lookup.rsplit('.', 2)[0])
        fast_tree = run_fast(a_id, threads, cluster, node)
        good_fast = fix_tree(fast_tree, a_id_lookup, final_fast)
        yield '%s/%s' % (out_dir, final_fast)
    if run_rax is True:
        final_rax = '%s.raxml.tree' % (a_id_lookup.rsplit('.', 2)[0])
        rax_tree = run_raxml(rax_out, boot, a_id_phylip, threads, a_id, model, cluster, node)
        good_tree = fix_tree(rax_tree, a_id_lookup, final_rax)
        yield '%s/%s' % (out_dir, final_rax)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = \
            '# run raxml and FastTree on aligned fasta file')
    parser.add_argument(\
            '-a', required = True, \
            help = 'aligned fasta file')
    parser.add_argument(\
            '-b', default = 100, required = False, \
            help = 'bootstraps (default: 100)')
    parser.add_argument(\
            '-t', default = 6, required = False, type = int, \
            help = 'threads (default: 6)')
    parser.add_argument(\
            '-m', default = False, required = False, \
            help = 'model (only for raxml, default: GTRCAT/PROTCATJTT)')
    parser.add_argument(\
            '--no-fast', action = 'store_false', required = False, \
            help = 'do not run FastTree')
    parser.add_argument(\
            '--no-rax', action = 'store_false', required = False, \
            help = 'do not run raxml')
    parser.add_argument(\
            '--cluster', action = 'store_true', required = False, \
            help = 'run on cluster')
    parser.add_argument(\
            '-node', default = False, required = False, \
            help = 'name of cluster node (optional: for use with --cluster and -t <threads>)')
    args = vars(parser.parse_args()) 
    alignment, bootstraps, threads, \
            fasttree, run_rax, model, cluster, node = \
            args['a'], args['b'], args['t'], args['no_fast'], \
            args['no_rax'], args['m'], args['cluster'], args['node']
    if cluster is False and node is not False:
        print >> sys.stderr, '# use --cluster with -node'
        exit()
    if cluster is True:
        if node is False:
            threads = 48
        else:
            threads = args['t']
    [i for i in rax(alignment, bootstraps, threads, fasttree, run_rax, model, cluster, node)]
