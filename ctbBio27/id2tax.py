#!/usr/bin/python2.7

'''
script for getting full taxonomy given a gene id (GI), TaxID, or partial taxonomic description
'''

import os
import sys
import argparse
from tokyocabinet import hash
from ctbBio27.check import check as check

# ncbi files
ncbi = {'gi_taxid_prot.dmp': ['gi2taxid.tch'], \
        'gi_taxid_nucl.dmp': ['gi2taxid.tch'], \
        'names.dmp': ['names2ids.tch', 'ids2names.tch'], \
        'nodes.dmp': ['nodes.tch'] \
        }
if 'databases' not in os.environ:
    directory = '/home/cbrown/databases/ncbi'
else:
    directory = '%s/ncbi' % (os.environ['databases'])
taxdir = directory + '/taxonomy'

# hierarchy
group = {'domain':'[0]', 'superkingdom':'[0]', \
    'phylum':'[1]', 'class':'[2]', \
    'order':'[3]', 'family':'[4]', \
    'genus':'[5]', 'species':'[6]', \
    'species group':'[6]', 'subspecies':'[7]', \
    'species subgroup':'[7]', 'query':'[8]' \
    }
reject = ['cellular organisms', 'root']

def get_tchs():
    [download_file(file) for file in ncbi]
    [make_tch(file) for file in ncbi]
    # database directory
    tchs = {} # dictionary of tch files with their open files as keys
    for tch in set([f for i in ncbi.values() for f in i]):
            full_tch = '%s/%s' % (taxdir, tch)
            db = hash.Hash()
            db.open(full_tch)
            tchs[tch] = db
    return tchs

def download_file(file):
    tch = [ncbi[file][0], '%s/%s' % (taxdir, ncbi[file][0])]
    if check(tch[1]) is False:
        os.system('mkdir -p %s' % (taxdir))
        location = '%s/%s' % (taxdir, file)
        if check(location) is False:
            if tch[0] == 'gi2taxid.tch':
                os.system('wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/%s.gz -P %s' % (file, taxdir))
                os.system('gunzip %s/%s.gz' % (taxdir, file))
            else:
                os.system('wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -P %s' % (taxdir))
                os.system('tar -zxf %s/taxdump.tar.gz -C %s' % (taxdir, taxdir))
                os.system('rm %s/taxdump.tar.gz' % (taxdir))
        
def make_tch(file):
    tch = [ncbi[file][0], '%s/%s' % (taxdir, ncbi[file][0])]
    if check(tch[1]) is False:
        if tch[0] == 'gi2taxid.tch':
            db = hash.Hash()
            db.open(tch[1])
            gi2taxid = ['%s/%s' % (taxdir, file) for file in 'gi_taxid_prot.dmp', 'gi_taxid_nucl.dmp']
            for file in gi2taxid:
                for line in open(file):
                    gi, taxid = line.strip().split()
                    db[gi] = taxid
            db.close()
        elif file == 'names.dmp':
            names2ids, ids2names = ['%s/%s' % (taxdir, i) for i in ncbi[file]]
            names2idsdb = hash.Hash()
            names2idsdb.open(names2ids)
            ids2namesdb = hash.Hash()
            ids2namesdb.open(ids2names)
            for line in open('%s/%s' % (taxdir, file)):
                if 'scientific name' in line:
                    id, sciname = [i.strip() for i in line.strip().split('|')[0:2]]
                    names2idsdb[sciname] = id
                    ids2namesdb[id] = sciname
            names2idsdb.close()
            ids2namesdb.close()
        elif file == 'nodes.dmp':
            db = hash.Hash()
            db.open(tch[1])
            for line in open('%s/%s' % (taxdir, file)):
                id, parent_id, parent = [i.strip() for i in line.strip().split('|')[0:3]]
                db[id] = '%s %s' % (parent_id, parent)
            db.close()

def level_id(query, tchs):
    current = ''
    if query in tchs['names2ids.tch']:
        id = tchs['names2ids.tch'][query]
        level = tchs['nodes.tch'][id].split()[1]
        if level in group:
            current = group[level]
    return current

def get_level(id, level, lineage, tchs):
    parent = tchs['nodes.tch'][id].split(' ', 1)
    query = tchs['ids2names.tch'][parent[0]]
    if query not in reject:
        lineage.append('%s%s' % (level_id(query, tchs), query))
    id = parent[0]
    level = parent[1]
    return id, level, lineage

def above(query, lineage, tchs):
    if query.strip() in tchs['names2ids.tch']:
        id = tchs['names2ids.tch'][query.strip()]
        level = tchs['nodes.tch'][id].split(' ', 1)[1]
        while id != '1':
            id, level, lineage = get_level(id, level, lineage, tchs)
        lineage.reverse()
    return lineage

def find_hierarchy(query, tchs):
    if query == 'n/a': # or if query not in db
        return 'n/a;%s' % (query)
#    lineage = ['%s%s' % (group['query'], query)]
    lineage = ['%s%s' % (level_id(query, tchs), query)]
    lineage = above(query, lineage, tchs)
    if len(lineage) == 1:
            query = query.split(' ')[0]
            query = query.replace('[', '')
            query = query.replace(']', '')
            lineage = above(query, lineage, tchs)
    return ';'.join(lineage)

def id2name(query, tchs, gi):
    """
    convert id number to name
    """
    if query.isdigit() is True:
        if gi is True:
            if query in tchs['gi2taxid.tch']:
                taxid = tchs['gi2taxid.tch'][query]
            else:
                query = 'n/a'
        else:
            taxid = query
        if taxid in tchs['ids2names.tch']:
            query = tchs['ids2names.tch'][taxid]
        else:
            query = 'n/a'
    return query

def get_taxa(args):
    tchs = get_tchs()
    search = args['i']
    for term in search:
        query = term
        query = id2name(query, tchs, args['gi'])
        yield [term, find_hierarchy(query, tchs)]
    for tch in tchs:
        tchs[tch].close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='# get NCBI taxonomy hierarchy')
    parser.add_argument(\
            '-i', nargs = '*', action = 'store',
            required = True, help = 'list of IDs or names (GI, NCBI TaxID, or name)')
    parser.add_argument(\
            '--gi', action = 'store_true', \
            help = 'IDs are GI numbers')
    args = vars(parser.parse_args())
    if args['i'][0] == '-':
        args['i'] = [i.strip() for i in sys.stdin]
    args['i'] = set(args['i'])
    for hierarchy in get_taxa(args):
        print '\t'.join(hierarchy)
