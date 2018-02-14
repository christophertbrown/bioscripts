#!/usr/bin/python2.7

"""
script for searching a query fasta against a database using either 
usearch or blast
"""

import sys
import os
import argparse
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

def blastdb(fasta, maxfile = 10000000):
    """
    make blast db
    """
    db = fasta.rsplit('.', 1)[0]
    type = check_type(fasta)
    if type == 'nucl':
        type = ['nhr', type]
    else:
        type = ['phr', type]
    if os.path.exists('%s.%s' % (db, type[0])) is False \
            and os.path.exists('%s.00.%s' % (db, type[0])) is False:
        print >> sys.stderr, '# ... making blastdb for: %s' % (fasta)
        os.system('makeblastdb \
                -in %s -out %s -dbtype %s -max_file_sz %s >> log.txt' \
                % (fasta, db, type[1], maxfile))
    else:
        print >> sys.stderr, '# ... database found for: %s' % (fasta)
    return db

def usearchdb5(fasta):
    """
    make usearch db
    """
    if '.udb' in fasta:
        print >> sys.stderr, '# ... database found: %s' % (fasta)
        return db
    type = check_type(fasta)
    if type == 'nucl':
        type = ['wdb', type]
    else:
        type = ['udb', type]
    db = '%s.%s' % (fasta.rsplit('.', 1)[0], type[0])
    if os.path.exists(db) is False:
        print >> sys.stderr, '# ... making usearch db for: %s' % (fasta)
        os.system('usearch -make%s %s -output %s >> log.txt' % (type[0], fasta, db))
    else:
        print >> sys.stderr, '# ... database found for: %s' % (fasta)
    return db

def usearchdb(fasta, alignment = 'local', usearch_loc = 'usearch'):
    """
    make usearch db
    """
    if '.udb' in fasta:
        print >> sys.stderr, '# ... database found: %s' % (fasta)
        return fasta
    type = check_type(fasta)
    db = '%s.%s.udb' % (fasta.rsplit('.', 1)[0], type)
    if os.path.exists(db) is False:
        print >> sys.stderr, '# ... making usearch db for: %s' % (fasta)
        if alignment == 'local':
            os.system('%s -makeudb_ublast %s -output %s >> log.txt' % (usearch_loc, fasta, db))
        elif alignment == 'global':
            os.system('%s -makeudb_usearch %s -output %s >> log.txt' % (usearch_loc, fasta, db))
    else:
        print >> sys.stderr, '# ... database found for: %s' % (fasta)
    return db

def phmmer2blast(phmmer, out):
    out = open(out, 'w')
    na = 'n/a'
    for line in open(phmmer):
        if line.startswith('#'):
            continue
        line = line.strip().split()
        na = 'n/a'
        if len(line) >= 6:
            blast = [line[2], line[0], na, na, na, na, na, na, na, na, line[4], line[5]]
            print >> out, '\t'.join(blast)
    out.close()

def phmmer(query, db, type, out, threads = '4', evalue = '0.01'):
    """
    run phmmer
    """
    if os.path.exists(out) is False:
        print '# ... running phmmer with %s as query and %s as database' % (query, db)
        os.system('phmmer -o %s.ph1 --tblout %s.ph2 --acc --noali --notextw -E %s --cpu %s %s %s' % (out, out, evalue, threads, query, db))
    else:
        print '# ... phmmer output found for %s as query and %s as database' % (query, db)
    phmmer2blast('%s.ph2' % out, out)

def blast(query, db, type, out, threads = '4', maxtargets = '100', megablast = False):
    """
    run blast
    """
    if os.path.exists(out) is False:
        db = blastdb(db) # make the database file, if necessary 
        print '# ... running blast with %s as query and %s as database' % (query, db)
        if type == 'nucl':
            blast = 'blastn'
            if megablast == True:
                blast = 'blastn -task megablast'
        else:
            blast = 'blastp'
        os.system('%s \
                -query %s -db %s -out %s -outfmt 6 \
                -max_target_seqs %s -num_threads %s >> log.txt' \
                % (blast, query, db, out, maxtargets, threads))
    else:
        print '# ... blast output found for %s as query and %s as database' % (query, db)

def usearch5(query, db, type, out, threads = '4', evalue = '100', alignment = 'local'):
    """
    run usearch
    """
    if os.path.exists(out) is False:
        print '# ... running usearch with %s as query and %s as database' % (query, db)
        if type[1] == 'nucl':
            threads = ''
        else:
            threads = '-threads %s' % (threads)
        os.system('usearch \
                -query %s -%s %s -blast6out %s \
                -evalue %s %s -%s >> log.txt' \
                % (query, type[0], db, out, evalue, threads, alignment))
    else:
        print '# ... usearch output found for %s as query and %s as database' % (query, db)

def usearch(query, db, type, out, threads = '6', evalue = '100', alignment = 'local', max_hits = 100, cluster = False):
    """
    run usearch
    """
    if 'usearch64' in os.environ:
        usearch_loc = os.environ['usearch64']
    else:
        usearch_loc = 'usearch'
    if os.path.exists(out) is False:
        db = usearchdb(db, alignment, usearch_loc) # make the database file, if neceesary
        print >> sys.stderr, '# ... running usearch with %s as query and %s as database' % (query, db)
        if type == 'nucl':
            strand = '-strand both'
        else:
            strand = ''
        if alignment == 'local' and cluster is False:
            os.system('%s \
                    -ublast %s -db %s -blast6out %s \
                    -evalue %s -threads %s %s -maxhits %s >> log.txt' \
                    % (usearch_loc, query, db, out, evalue, threads, strand, max_hits))
        elif alignment == 'global' and cluster is False:
            os.system('%s \
                    -usearch_global %s -db %s -blast6out %s \
                    -id 0.10 -threads %s %s >> log.txt' \
                    % (usearch_loc, query, db, out, threads, strand))
        elif alignment == 'local' and cluster is True:
            qsub = 'qsub -l nodes=1:ppn=24 -m e -N usearch'
            os.system('echo "%s -ublast `pwd`/%s -db %s -blast6out `pwd`/%s -evalue %s -threads %s %s -maxhits %s >> `pwd`/log.txt" | %s' \
                    % (usearch_loc, query, db, out, evalue, threads, strand, max_hits, qsub))
        else:
            print >> sys.stderr, 'specify local or global alignment'
            exit()
    else:
        print >> sys.stderr, '# ... usearch output found for %s as query and %s as database' % (query, db)

def outfile(query, method, database, prefix):
    type = check_type(query)
    query = query.rsplit('.', 1)[0]
    database = database.rsplit('.', 1)[0]
    if '/' in query:
        query = query.rsplit('/', 1)[1]
    if '/' in database:
        database = database.rsplit('/', 1)[1]
    out = '%s-%s_%s-%s.b6' % (query, method.split('-')[0], type, database)
    if prefix is not False:
        out = '%s%s' % (prefix, out)
    return out, type

def search(query, database, method = 'usearch', alignment = 'local', max_hits = 100, threads = '6', prefix = False):
    out, type = outfile(query, method, database, prefix)
    if method == 'usearch':
        usearch(query, database, type, out, alignment = alignment, max_hits = max_hits, threads = threads)
    elif method == 'usearch-cluster':
        usearch(query, database, type, out, alignment = alignment, max_hits = max_hits, threads = threads, cluster = True)
    elif method == 'blast':
        blast(query, database, type, out, threads = threads)
    elif method == 'phmmer':
        phmmer(query, database, type, out, threads = threads)
    return out

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = \
            '# search sequences against database')
    parser.add_argument(\
            '-q', required = True, \
            help = 'query')
    parser.add_argument(\
            '-d', required = True, \
            help = 'database')
    parser.add_argument(\
            '-a', default = 'usearch', \
            help = 'algorithm: usearch (default), usearch-cluster, blast, phmmer')
    parser.add_argument(\
            '-t', default = "6", \
            help = 'threads (default = 6)')
    args = vars(parser.parse_args())
    threads, query, database, method = args['t'], args['q'], args['d'], args['a']
    if method != 'usearch-cluster':
        os.system('cat %s' % (search(query, database, method, threads = threads)))
    else:
        search(query, database, method, threads = '48')
