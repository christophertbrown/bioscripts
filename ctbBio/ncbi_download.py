#!/usr/bin/env python3

"""
script for downloading genomes from NCBI
"""

# python modules
import os
import sys
import argparse
from tqdm import tqdm
from subprocess import Popen
from glob import glob as glob
from multiprocessing import Pool
from subprocess import Popen, PIPE

def wget(ftp, f = False, exclude = False, name = False, tries = 10):
    """
    use wget to download file
    """
    # file name
    if f is False:
        f = ftp.rsplit('/', 1)[-1]
    # downloaded file if it does not already exist
    t = 0
    if name is not False:
        print('# downloading:', name, f)
    while len(glob(f)) == 0:
        t += 1
        if exclude is False:
            command = 'wget -q --random-wait %s' % (ftp)
        else:
            command = 'wget -q --random-wait -R %s %s' % (exclude, ftp)
        p = Popen(command, shell = True)
        p.communicate()
        if t >= tries:
            print('not downloaded:', name, f)
            return [f, False]
    return [f, True]

def check(line, queries):
    """
    check that at least one of
    queries is in list, l
    """
    line = line.strip()
    spLine = line.replace('.', ' ').split()
    matches = set(spLine).intersection(queries)
    if len(matches) > 0:
        return matches, line.split('\t')
    return matches, False

def entrez(db, acc):
    """
    search entrez using specified database
    and accession
    """
    c1 = ['esearch', '-db', db, '-query', acc]
    c2 = ['efetch', '-db', 'BioSample', '-format', 'docsum']
    p1 = Popen(c1, stdout = PIPE, stderr = PIPE)
    p2 = Popen(c2, stdin = p1.stdout, stdout = PIPE, stderr = PIPE)
    return p2.communicate()

def searchAccession(acc):
    """
    attempt to use NCBI Entrez to get
    BioSample ID
    """
    # try genbank file
    # genome database
    out, error = entrez('genome', acc)
    for line in out.splitlines():
        line = line.decode('ascii').strip()
        if 'Assembly_Accession' in line or 'BioSample' in line:
            newAcc = line.split('>')[1].split('<')[0].split('.')[0].split(',')[0]
            if len(newAcc) > 0:
                return (True, acc, newAcc)
    # nucleotide database
    out, error = entrez('nucleotide', acc)
    for line in out.splitlines():
        line = line.decode('ascii').strip()
        if 'Assembly_Accession' in line or 'BioSample' in line:
            newAcc = line.split('>')[1].split('<')[0].split('.')[0].split(',')[0]
            if len(newAcc) > 0:
                return (True, acc, newAcc)
    # assembly database
    out, error = entrez('assembly', acc)
    for line in out.splitlines():
        line = line.decode('ascii').strip()
        if 'Assembly_Accession' in line or 'BioSample' in line:
            newAcc = line.split('>')[1].split('<')[0].split('.')[0].split(',')[0]
            if len(newAcc) > 0:
                return (True, acc, newAcc)
    for error in error.splitlines():
        error = error.decode('ascii').strip()
        if '500 Can' in error:
            return (False, acc, 'no network')
    return (False, acc, 'efetch failed')

def getFTPs(accessions, ftp, search, exclude, convert = False, threads = 1, attempt = 1,
            max_attempts = 2):
    """
    download genome info from NCBI
    """
    info = wget(ftp)[0]
    allMatches = []
    for genome in open(info, encoding = 'utf8'):
        genome = str(genome)
        matches, genomeInfo = check(genome, accessions)
        if genomeInfo is not False:
            f = genomeInfo[0] + search
            Gftp = genomeInfo[19]
            Gftp = Gftp + '/' + search
            allMatches.extend(matches)
            yield (Gftp, f, exclude, matches)
    # print accessions that could not be matched
    # and whether or not they could be converted (optional)
    newAccs = []
    missing = accessions.difference(set(allMatches))
    if convert is True:
        pool = Pool(threads)
        pool = pool.imap_unordered(searchAccession, missing)
        for newAcc in tqdm(pool, total = len(missing)):
            status, accession, newAcc = newAcc
            if status is True:
                newAccs.append(newAcc)
            print('not found:', accession, '->', newAcc)
    else:
        for accession in missing:
            print('not found:', accession)
    # re-try after converting accessions (optional)
    if len(newAccs) > 0 and attempt <= max_attempts:
        print('convert accession attempt', attempt)
        attempt += 1
        for hit in getFTPs(set(newAccs), ftp, search, exclude, convert,
                threads = 1, attempt = attempt):
            yield hit

def wgetGenome(pars):
    """
    """
    ftp, f, exclude, matches = pars
    name = ';'.join(list(matches))
    return wget(ftp, f, exclude, name)

def download(args):
    """
    download genomes from NCBI
    """
    accessions, infoFTP = set(args['g']), args['i']
    search, exclude = args['s'], args['e']
    FTPs = getFTPs(accessions, infoFTP, search, exclude, threads = args['t'],
            convert = args['convert'])
    if args['test'] is True:
        for genome in FTPs:
            print('found:', ';'.join(genome[-1]), genome[0])
        return FTPs
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
            required = False,
            help = 'search exclusion term, or False (default = "*from_genomic*")')
    parser.add_argument(\
            '-i', default = ftp,
            required = False, help = 'genome info FTP (default: %s)' % (ftp))
    parser.add_argument(\
            '-t', default = 3, type = int,
            required = False, help = 'threads (default = 3)')
    parser.add_argument(\
            '--convert', action = 'store_true', required = False,
            help = 'convert missing accessions using Entrez Direct (slow; requires `esearch` and `efetch`)')
    parser.add_argument(\
            '--test', action = 'store_true', required = False,
            help = 'look for genomes, but do not download them')
    args = vars(parser.parse_args())
    if args['e'] == 'False' or args['e'] == 'FALSE':
        args['e'] = False
    if args['g'][0] == '-':
        args['g'] = [i.strip() for i in sys.stdin]
    print('# downloading genome info:', args['i'])
    download(args)

