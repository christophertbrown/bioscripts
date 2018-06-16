#!/usr/bin/env python3

"""
script for downloading genomes from NCBI

ctb@berkeley.edu

The MIT License

Copyright (c) 2011 Christopher T. Brown

Permission is hereby granted, free of charge,
to any person obtaining a copy of this software and
associated documentation files (the "Software"), to
deal in the Software without restriction, including
without limitation the rights to use, copy, modify,
merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom
the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice
shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

# python modules
import os
import sys
import argparse
from tqdm import tqdm
from subprocess import Popen
from glob import glob as glob
from multiprocessing import Pool

def wget(ftp, f = False, exclude = False, tries = 10):
    """
    use wget to download file
    """
    # file name
    if f is False:
        f = ftp.rsplit('/', 1)[-1]
    # downloaded file if it does not already exist
    t = 0
    while len(glob(f)) == 0:
        t += 1
        if exclude is False:
            command = 'wget -q --random-wait %s' % (ftp)
        else:
            command = 'wget -q --random-wait -R %s %s' % (exclude, ftp)
        p = Popen(command, shell = True)
        p.communicate()
        if t >= tries:
            print('not downloaded:', f, ftp)
            return False
    return f

def check(line, queries):
    """
    check that at least one of
    queries is in list, l
    """
    line = line.strip()
    spLine = line.replace('.', ' ').split()
    if len(set(spLine).intersection(queries)) > 0:
        return line.split('\t')
    return False

def getFTPs(accessions, ftp, search, exclude):
    """
    download genome info from NCBI
    """
    info = wget(ftp)
    for genome in open(info, encoding = 'utf8'):
        genome = str(genome)
        genomeInfo = check(genome, accessions)
        if genomeInfo is not False:
            f = genomeInfo[0] + search
            ftp = genomeInfo[19]
            ftp = ftp + '/' + search
            yield (ftp, f, exclude)

def wgetGenome(pars):
    """
    """
    ftp, f, exclude = pars
    return wget(ftp, f, exclude)

def download(args):
    """
    download genomes from NCBI
    """
    accessions, infoFTP = set(args['g']), args['i']
    search, exclude = args['s'], args['e']
    FTPs = getFTPs(accessions, infoFTP, search, exclude)
    pool = Pool(args['t'])
    pool = pool.imap_unordered(wgetGenome, FTPs)
    files = []
    for f in tqdm(pool, total = len(accessions)):
        files.append(f)
    return files

if __name__ == '__main__':
    ftp = 'ftp://ftp.ncbi.nih.gov/genomes/genbank/assembly_summary_genbank.txt'
    parser = argparse.ArgumentParser(description='# download genomes from NCBI')
    parser.add_argument(\
            '-g', nargs = '*', action = 'store',
            required = True, help = 'list of genome accession numbers (- for stdin)')
    parser.add_argument(\
            '-s', default = '*.fna.gz',
            required = False, help = 'search term for download (default = "*.fna.gz")')
    parser.add_argument(\
            '-e', default = '*from_genomic*',
            required = False, help = 'search exclusion term (default = "*from_genomic*")')
    parser.add_argument(\
            '-i', default = ftp,
            required = False, help = 'genome info FTP (default: %s)' % (ftp))
    parser.add_argument(\
            '-t', default = 3, type = int,
            required = False, help = 'threads (default = 3)')
    args = vars(parser.parse_args())
    if args['g'][0] == '-':
        args['g'] = [i.strip() for i in sys.stdin]
    download(args)
