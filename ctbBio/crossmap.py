#!/usr/bin/env python3

"""
script for mapping reads against scaffolds
"""

import os
import sys
import json
import numpy
import shutil
import random
import argparse
import subprocess

def bowtiedb(fa, keepDB):
    """
    make bowtie db
    """
    btdir = '%s/bt2' % (os.getcwd())
    # make directory for 
    if not os.path.exists(btdir):
        os.mkdir(btdir)
    btdb = '%s/%s' % (btdir, fa.rsplit('/', 1)[-1])
    if keepDB is True:
        if os.path.exists('%s.1.bt2' % (btdb)):
            return btdb
    p = subprocess.Popen('bowtie2-build -q %s %s' \
        % (fa, btdb), shell = True)
    p.communicate()
    return btdb

def bowtie(sam, btd, f, r, u, opt, no_shrink, threads):
    """
    generate bowtie2 command
    """
    bt2 = 'bowtie2 -x %s -p %s ' % (btd, threads)
    if f is not False:
        bt2 += '-1 %s -2 %s ' % (f, r)
    if u is not False:
        bt2 += '-U %s ' % (u)
    bt2 += opt
    if no_shrink is False:
        if f is False:
            bt2 += ' | shrinksam -u -k %s-shrunk.sam ' % (sam)
        else:
            bt2 += ' | shrinksam -k %s-shrunk.sam ' % (sam)
    else:
        bt2 += ' > %s.sam' % (sam)
    return bt2

def chunks(l, n):
    return numpy.array_split(numpy.array(l), n)

def crossmap(fas, reads, options, no_shrink, keepDB, threads, cluster, nodes):
    """
    map all read sets against all fasta files
    """
    if cluster is True:
        threads = '48'
    btc = []
    for fa in fas:
        btd = bowtiedb(fa, keepDB)
        F, R, U = reads
        if F is not False:
            if U is False:
                u = False
            for i, f in enumerate(F):
                r = R[i]
                if U is not False:
                    u = U[i]
                sam = '%s/%s-vs-%s' % (os.getcwd(), \
                        fa.rsplit('/', 1)[-1], f.rsplit('/', 1)[-1].rsplit('.', 3)[0])
                btc.append(bowtie(sam, btd, f, r, u, options, no_shrink, threads))
        else:
            f = False
            r = False
            for u in U:
                sam = '%s/%s-vs-%s' % (os.getcwd(), \
                        fa.rsplit('/', 1)[-1], u.rsplit('/', 1)[-1].rsplit('.', 3)[0])
                btc.append(bowtie(sam, btd, f, r, u, options, no_shrink, threads))
    if cluster is False:
        for i in btc:
            p = subprocess.Popen(i, shell = True)
            p.communicate()
    else:
        ID = ''.join(random.choice([str(i) for i in range(0, 9)]) for _ in range(5))
        for node, commands in enumerate(chunks(btc, nodes), 1):
            bs = open('%s/crossmap-qsub.%s.%s.sh' % (os.getcwd(), ID, node), 'w')
            print('\n'.join(commands), file=bs)
            bs.close()
            p = subprocess.Popen(\
                    'qsub -V -N crossmap %s' \
                        % (bs.name), \
                    shell = True)
            p.communicate()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# cross map using bowtie')
    parser.add_argument(\
            '-f', nargs = '*', action = 'store', \
            required = True, help = 'path to fasta(s)')
    parser.add_argument(\
            '-1', nargs = '*', action = 'store', \
            default = False, help = 'path to forward reads')
    parser.add_argument(\
            '-2', nargs = '*', action = 'store', \
            default = False, help = 'path to reverse reads')
    parser.add_argument(\
            '-U', nargs = '*', action = 'store', \
            default = False, help = 'path to single reads')
    parser.add_argument(\
            '-p', default = '--very-fast --reorder --quiet', \
            help = 'bowtie options (default = "--very-fast --reorder --quiet"')
    parser.add_argument(\
            '--no-shrink', action = 'store_true', help = 'do not use shrinksam')
    parser.add_argument(\
            '--keepDB', action = 'store_true', help ='do not overwrite bowtie database')
    parser.add_argument(\
            '-t', default = '6', help = 'number of cpus (default = 6)')
    parser.add_argument(\
            '--cluster', action = 'store_true', help = 'run on cluster')
    parser.add_argument(\
            '-n', default = 1, type = int, help = 'number of cluster nodes (default = 1)')
    args = vars(parser.parse_args())
    fa = [os.path.abspath(i) for i in args['f']]
    if (args['1'] is False or args['2'] is False) and args['U'] is False:
        print('# specify -1 and -2 and/or -U', file=sys.stderr)
        exit()
    if args['1'] is not False:
        f = [os.path.abspath(i) for i in args['1']]
        r = [os.path.abspath(i) for i in args['2']]
    else: 
        f = r = False
    if args['U'] is not False:
        u = [os.path.abspath(i) for i in args['U']]
    else:
        u = False
    crossmap(fa, [f, r, u], args['p'], args['no_shrink'], args['keepDB'], args['t'], \
            args['cluster'], args['n'])
