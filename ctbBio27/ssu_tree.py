#!/usr/bin/env python2.7

"""
script for making a quick 16S/18S rRNA gene tree
"""

# python modules
import sys
import os
import argparse
import subprocess
from glob import glob as glob
from tokyocabinet import hash

# ctb modules
from ctbBio27.nr_fasta import de_rep as fix_fasta
from ctbBio27.search import search as search_seqs
from ctbBio27.numblast import best as numblast
from ctbBio27.fasta import iterate_fasta as parse_fasta
from ctbBio27.ssufromHMM import find_16S as find_16S
from ctbBio27.ssufromHMM import run_cmsearch as run_cmsearch
from ctbBio27.strip_masked import parse_masked as strip_masked
from ctbBio27.stockholm2fa import stock2fa as stock2fa
from ctbBio27.rax import rax as rax
from ctbBio27.fasta2tch import tch as fasta2tch

def remove_char(header):
    remove_characters = ['\'', '/', '\\', ':', ',', '(', ')', ' ', '|', ';',] \
            # characters to remove from headers
    for character in remove_characters:
        header = header.replace(character, '_')
    return header

def clean_seqs(fasta, out):
    """
    clean sequence files
    """
    clean = open('%s/%s.clean.fa' % (out, fasta.split('/')[-1].rsplit('.', 1)[0]), 'w')
    lookup = open('%s/%s.clean.lookup' % (out, fasta.split('/')[-1].rsplit('.', 1)[0]), 'w')
    for seq in fix_fasta([fasta], append_index = True, return_original = True):
        original_header, seq = seq
        header = remove_char(seq[0].split()[0])
        print >> clean, '\n'.join([header, seq[1].upper()])
        print >> lookup, '\t'.join(['>%s' % (' '.join(original_header)), header])
    clean.close()
    lookup.close()
    return clean.name

def trim_sequences(seqs, trim_seqs, ssucmdb, threads):
    """
    identifiy 16S rRNA genes
    """
    if trim_seqs is False:
        return seqs
    print >> sys.stderr, '# identifying 16S rRNA genes'
    trimmed = open('%s.trimmed.fa' % (seqs.rsplit('.', 1)[0]), 'w')
    for seq in find_16S([seqs], \
            run_cmsearch([open(seqs)], threads, ssucmdb), masking = True):
        print >> trimmed, '\n'.join(seq)
    os.system('rm -f cmsearch.log')
    return trimmed.name


def combine_with_hits(clean, s_db, search_out, hits):
    """
    combine sequences with best hits from search
    """
    best = set([hit[1].split()[0] for hit in numblast(open(search_out), hits, False, False)])
    combo = '%s.best%srefs.fa' % (clean.rsplit('.', 1)[0], hits)
    if os.path.exists(combo) is True:
        return combo
    combo = open(combo, 'w')
    for seq in parse_fasta(clean):
        print >> combo, '\n'.join(seq)
    # create/open tch for search database
    s_tch = '%s.tch' % (s_db)
    if os.path.exists(s_tch) is False:
        fasta2tch(s_db)
    id2seq = hash.Hash()
    id2seq.open(s_tch)
    # get sequences for best hits from tch
    for hit in best:
        seq = id2seq[hit].split('\n')
        header = remove_char(seq[0].split()[0]).replace('>', '>best-hit_')
        print >> combo, '\n'.join([header, seq[1].upper()])
    combo.close()
    id2seq.close()
    return combo.name

def remove_insertions(seqs, ignore_inserts, max, ssucmdb, threads):
    """
    remove insertions
    """
    if ignore_inserts is True:
        return seqs
    print >> sys.stderr, '# removing insertion sequences'
    masked = open('%s.masked.fa' % (seqs.rsplit('.', 1)[0]), 'w')
    stripped = open('%s.i%s-stripped.fa' % (seqs.rsplit('.', 1)[0], max), 'w')
    nr_list = []
    for seq in find_16S([seqs], \
            run_cmsearch([open(seqs)], threads, ssucmdb), masking = True):
        id = seq[0].split()[0]
        if id in nr_list:
            continue
        nr_list.append(id)
        print >> masked, '\n'.join(seq)
        nm, m = strip_masked(seq, max)
        print >> stripped, '\n'.join(['%s removed masked >=%s' % (seq[0], max), ''.join(nm)])
    os.system('rm -f cmsearch.log')
    return stripped.name

def stk2afa(stk, min_aligned, out):
    """
    1. convert stk to afa
    2. remove insertion columns
    3. remove seqs with <= min_aligned bases
    """
    name = '.'.join(stk.rsplit('/', 1)[1].rsplit('.', 2)[0:2])
    afa = open('%s/%s.afa' % (out, name), 'w')
    # convert to fasta
    for id, seq in stock2fa(open(stk)).items():
        # strip insertion columns
        seq = ''.join([b for b in seq[0] if b == '-' or b.isupper()])
        if len([i for i in seq if i != '-']) >= min_aligned:
            print >> afa, '\n'.join(['>%s' % (id), ''.join(seq)])
    return [name.rsplit('.', 1)[1], afa.name]

def align_seqs(seqs, archaea, bacteria, eukarya, min_aligned, out, threads):
    """
    align sequences using ssu-align
    """
    if archaea is True:
        model = "-n archaea --no-search"
    elif bacteria is True:
        model = "-n bacteria --no-search"
    elif eukarya is True:
        model = "-n eukarya --no-search"
    else:
        model = ""
    ssu_dir = '%s.ssu-align' % (seqs.rsplit('.', 1)[0].rsplit('/', 1)[-1])
    # check size of fasta file - if large use ssu-align in multi-threaded mode
    if int(os.path.getsize(seqs)) > 150000: # this is a large fasta file
        ssu_sh = '%s.ssu-align.sh' % (ssu_dir)
        p = subprocess.Popen('\
                cd %s; \
                ssu-prep -x %s %s %s %s >>log.txt 2>>log.txt; \
                ./%s >>log.txt 2>>log.txt'
                % (out, \
                model, os.path.abspath(seqs), ssu_dir, threads, \
                ssu_sh), \
                shell = True)
        p.communicate()
    else:
        p = subprocess.Popen('\
                cd %s; \
                ssu-align %s %s %s >>log.txt 2>>log.txt' \
                % (out, \
                model, seqs, ssu_dir), \
                shell = True)
        p.communicate()
    stks = glob('%s/%s/*stk' % (out, ssu_dir))
    afas = [stk2afa(stk, min_aligned, out) for stk in stks]
    return {i[0]:i[1] for i in afas}

def run_tree(seqs, refs, threads, raxml, cluster = False, no_fast = False):
    """
    combine alignments and run tree
    """
    aligned = int(os.path.getsize(seqs))
    if aligned == 0:
        print >> sys.stderr, '# sequences could not be aligned'
        exit()
    combo = open('%s.refs.afa' % (seqs.rsplit('.', 1)[0]), 'w')
    to_combine = [i for i in seqs, refs if i is not False]
    for seq in fix_fasta(to_combine, append_index = False, return_original = False):
        print >> combo, '\n'.join(seq)
    combo.close()
    if no_fast == False:
        fast = True
    else:
        fast = False
    return rax(combo.name, 100, threads, fast = fast, \
            run_rax = raxml, model = False, cluster = False)

def final_lookup(pLook, tLook):
    """
    convert names in tree to original (safe version)
    """
    tree = {}
    lookup = {}
    for line in open(tLook):
        ID, safe_short, short = line.replace('>', '').strip().split('\t')
        tree[short] = safe_short
    for line in open(pLook):
        name, short = line.replace('>', '').strip().split('\t')
        safe_name = remove_char(name)
        if short in tree:
            safe_short = tree[short]
            lookup[safe_short] = safe_name
    return lookup

def tree(seqs, trim_seqs, s_db, hits, rb, ra, re, ru, \
        no_search, no_ref, archaea, bacteria, eukarya, \
        length, raxml, ignore_inserts, ssucmdb, \
        max_ins, out, threads, no_fast):
    """
    generate 16S rRNA gene tree:
    1. clean up sequence format
    2. search against database and collect best hits
    3. identify and remove large insertion sequences
    4. align sequences and hits using ssu-align
        - if reference set is un-aligned, align the references
    5. remove insertion columns from alignment and combine with reference alignment
    6. run tree
    """
    out = os.path.abspath(out)
    os.system('mkdir -p %s' % (out))
    clean = clean_seqs(seqs, out)
    trimmed = trim_sequences(clean, trim_seqs, ssucmdb, threads)
    if no_search is False:
        search_out = search_seqs(trimmed, s_db, method = 'usearch', \
                alignment = 'local', max_hits = hits, threads = threads, \
                prefix = '%s/' % (out))
        print >> sys.stderr, '# collecting reference sequences'
        with_hits = combine_with_hits(trimmed, s_db, search_out, hits)
        os.system('rm -f log.txt')
    else:
        with_hits = trimmed
    no_inserts = remove_insertions(with_hits, ignore_inserts, max_ins, ssucmdb, threads)
    print >> sys.stderr, '# aligning sequences'
    aligned = align_seqs(no_inserts, archaea, bacteria, eukarya, length, out, threads)
    if ru is not False:
        refs = align_seqs(ru, archaea, bacteria, eukarya, length, out, threads)
    else:
        refs = {'bacteria':rb, 'archaea':ra, 'eukarya':re}
    if no_fast is False or (no_fast is False and raxml is False):
        print >> sys.stderr, '# running tree inference'
    projectLookup = glob('%s/*lookup' % (out))[0]
    if raxml is False and no_fast is True:
            print >> sys.stderr, '# skipping tree inference'
    for domain, alignment in aligned.items():
        print >> sys.stderr, '# %s alignment: %s' % (domain, alignment)
        if domain not in refs:
            ref = False
        elif no_ref is True and ru is False:
            ref = False
        else:
            ref = refs[domain]
        if raxml is False and no_fast is True:
            continue
        tree = [tree for tree in \
                run_tree(alignment, ref, threads, raxml, cluster, no_fast)][0]
        treeLookup = tree.rsplit('.', 2)[0] + '.id.lookup'
        lookup = final_lookup(projectLookup, treeLookup)
        final = open('%s/%s-%s-FINAL.tree' \
                    % (out, out.rsplit('/', 1)[-1], domain), 'w')
        print >> sys.stderr, 'tree: %s' % (final.name)
        for line in open(tree):
            line = line.strip()
            for find, replace in lookup.items():
                line = line.replace(find, replace)
            print >> final, line
        final.close()

def check_options(s_db, no_search, ssucmdb, ignore_inserts, \
        rb, ra, re, ru, no_ref, prok, archaea, bacteria, eukarya):
    """
    check that options are used correctly
    """
    # use one: --archaea, --bacteria, --eukarya
    if (archaea is not False and bacteria is not False) \
            or (eukarya is not False and bacteria is not False) \
            or (eukarya is not False and archaea is not False):
        print >> sys.stderr, '# use --archaea, --bacteria, or --eukarya'
        exit()

    # blast search database (not if --no-search is true)
    if (s_db is False and 'ssuref' not in os.environ) and no_search is False:
        print >> sys.stderr, '# specify search database with -d or ssuref env. variable'
        exit()

    # cmsearch database (not if --ignore-inserts is true)
    if (ssucmdb is False and 'ssucmdb' not in os.environ) and ignore_inserts is False:
        print >> sys.stderr, \
                '# specify CM for 16S HMM search with --ssu-cmdb or ssucmdb env. variable'
        exit()

    # don't use ru with rb, ra, or re 
    if ru is not False:
        if ra is not False or rb is not False or re is not False:
            print >> sys.stderr, '# do not use -rb, -ra, or -re with -ru'
            exit()

    # when using --archaea, --bacteria, or --eukarya, use paired reference set
    if archaea is True:
        if rb is not False or re is not False:
            print >> sys.stderr, \
                    '# use -ra with --archaea, -rb with --bacteria, and -re with --eukarya'
            exit()
    if bacteria is True:
        if ra is not False or re is not False:
            print >> sys.stderr, \
                    '# use -ra with --archaea, -rb with --bacteria, and -re with --eukarya'
            exit()
    if eukarya is True:
        if rb is not False or ra is not False:
            print >> sys.stderr, \
                    '# use -ra with --archaea, -rb with --bacteria, and -re with --eukarya'
            exit()

    # use bacteria and archaea reference sets with --prok option
    if prok is True:
        if ra is not False or rb is not False or ra is not False:
            print >> sys.stderr, \
                    '# do not use -ra, -rb, or -ra with --prok'
            exit()
        if archaea is False and bacteria is False:
            print >> sys.stderr, \
                    '# use --prok with --archaea or --bacteria'
            exit()
        if archaea is True:
            if ra is False and 'ssual_ref_prok_archaea' not in os.environ:
                print >> sys.stderr, \
                   '# specify archaea refs with -ra or ssual_ref_prok_archaea env. variable'
                exit()
        if bacteria is True:
            if rb is False and 'ssual_ref_prok_bacteria' not in os.environ:
                print >> sys.stderr, \
                    '# specify bacteria refs with -rb or ssual_ref_prok_bacteria env. variable'
                exit()
        if eukarya is True:
            print >> sys.stderr, \
                    '# do not use --prok with --eukarya'
            exit()
    else: # make sure bacteria, archaea, and eukarya reference sets are provided
    # ra
        if ra is False and 'ssual_ref_archaea' not in os.environ:
          print >> sys.stderr, \
                  '# specify archaea refs with -ra or ssual_ref_archaea env. variable'
          exit()
        # rb
        if rb is False and 'ssual_ref_bacteria' not in os.environ:
            print >> sys.stderr, \
                    '# specify bacteria refs with -rb or ssual_ref_bacteria env. variable'
            exit()
        # re
        if re is False and 'ssual_ref_eukarya' not in os.environ:
            print >> sys.stderr, \
                    '# specify eukarya refs with -re or ssual_ref_eukarya env. variable'
            exit()

def define_dbs(s_db, no_search, ssucmdb, ignore_inserts, \
        rb, ra, re, ru, no_ref, prok, archaea, bacteria, eukarya, out):
    """
    define database to search and use for reference sets
    """
    # check options
    check_options(s_db, no_search, ssucmdb, ignore_inserts, \
            rb, ra, re, ru, no_ref, prok, archaea, bacteria, eukarya)

    # output directory
    if out is False:
        if archaea is True:
            out = '%s.aligned2archaea.tree' % (seqs.rsplit('.', 1)[0].rsplit('/', 1)[-1])
        elif bacteria is True:
            out = '%s.aligned2bacteria.tree' % (seqs.rsplit('.', 1)[0].rsplit('/', 1)[-1])
        elif eukarya is True:
            out = '%s.aligned2eukarya.tree' % (seqs.rsplit('.', 1)[0].rsplit('/', 1)[-1])
        else:
            out = '%s.tree' % (seqs.rsplit('.', 1)[0].rsplit('/', 1)[-1])

    # blast search database (not if --no-search is true)
    if no_search is True:
        s_db = False
    elif s_db is False:
        s_db = os.environ['ssuref']
    else:
        s_db = os.path.abspath(s_db)

    # cmsearch database (not if --ignore-inserts is true)
    if ignore_inserts is True:
        ssucmdb = False
    elif ssucmdb is False:
        ssucmdb = os.environ['ssucmdb']
    else:
        ssucmdb = os.path.abspath(ssucmdb)

    # unaligned reference set provided
    if ru is not False:
        return s_db, ssucmdb, False, False, False, os.path.abspath(ru), out

    # default reference sets

    # force archaea
    if archaea is True:
        if ra is not False:
            return s_db, ssucmdb, False, ra, False, ru, out
        if prok is True:
            ra = os.environ['ssual_ref_prok_archaea']
        else:
            ra = os.environ['ssual_ref_archaea']
        return s_db, ssucmdb, False, ra, False, ru, out

    # force bacteria
    if bacteria is True:
        if rb is not False:
            return s_db, ssucmdb, rb, False, False, ru, out
        if prok is True:
            rb = os.environ['ssual_ref_prok_bacteria']
        else:
            rb = os.environ['ssual_ref_bacteria']
        return s_db, ssucmdb, rb, False, False, ru, out

    # force eukarya
    if eukarya is True:
        if re is not False:
            return s_db, ssucmdb, False, False, re, ru, out
        re = os.environ['ssual_ref_eukarya']
        return s_db, ssucmdb, False, False, re, ru, out

    # no force
    if ra is False:
        ra = os.environ['ssual_ref_archaea']
    if rb is False:
        rb = os.environ['ssual_ref_bacteria']
    if re is False:
        re = os.environ['ssual_ref_eukarya']
    return s_db, ssucmdb, rb, ra, re, ru, out

def print_databases(s_db, ssucmdb, rb, ra, re, ru, no_ref):
    """
    print which databases are being used to stdout
    """
    if no_ref is True:
        rb = ra = re = False
    print >> sys.stderr, "# database for finding best-hits: %s" % (s_db)
    print >> sys.stderr, "# CM for ssu analysis: %s" % (ssucmdb)
    print >> sys.stderr, "# aligned bacteria reference set: %s" % (rb)
    print >> sys.stderr, "# aligned archaea reference set: %s" % (ra)
    print >> sys.stderr, "# aligned eukarya reference set: %s" % (re)
    print >> sys.stderr, "# unaligned reference set: %s" % (ru)

if __name__ == '__main__':
    # desc.
    parser = argparse.ArgumentParser(\
            description = '# make tree from ssu (16S/18S) rRNA gene sequences')
    # input
    parser.add_argument(\
            '-i', required = True, help = 'input fasta file')
    # options
    parser.add_argument(\
            '--trim', action = 'store_true', \
            help = 'identify 16S/18S rRNA genes (use with contigs)')
    parser.add_argument(\
            '--no-search', action = 'store_true', \
            help = 'do not search for best hits')
    parser.add_argument(\
            '--no-ref', action = 'store_true', \
            help = 'do not include reference sequences')
    parser.add_argument(\
            '--archaea', action = 'store_true', \
            help = 'force alignment to archaea model')
    parser.add_argument(\
            '--bacteria', action = 'store_true', \
            help = 'force alignment to bacteria model')
    parser.add_argument(\
            '--eukarya', action = 'store_true', \
            help = 'force alignment to eukarya model')
    parser.add_argument(\
            '--prok', action = 'store_true', \
            help = 'incldue bacteria and archaea refs. (use with --archaea or --bacteria)')
    parser.add_argument(\
            '-l', default = 800, \
            required = False, type = int, \
            help = 'minimum aligned length required for inclusion (default = 800)')
    parser.add_argument(\
            '--max-ins', default = 10, \
            required = False, type = int, \
            help = 'maximum insertion length allowed prior to alignment (default = 10)')
    parser.add_argument(\
            '--raxml', action = 'store_true', \
            help = 'run RAxML in addition to FastTree2')
    parser.add_argument(\
            '--no-fast', action = 'store_true', \
            help = 'do not run FastTree2')
    parser.add_argument(\
            '--ignore-insertions', action = 'store_true', \
            help = 'ignore insertions in genes')
    # searching databases
    parser.add_argument(\
            '-d', required = False, \
            default = False, \
            help = 'database to search for closely related sequences')
    parser.add_argument(\
            '--hits', required = False, type = int, \
            default = 3, \
            help = 'number of best hits to include (default = 3)')
    # reference sets
    parser.add_argument(\
            '-rb', required = False, \
            default = False, \
            help = 'references aligned to ssu-align bacteria model')
    parser.add_argument(\
            '-ra', required = False, \
            default = False, \
            help = 'references aligned to ssu-align archaea model')
    parser.add_argument(\
            '-re', required = False, \
            default = False, \
            help = 'references aligned to ssu-align eukarya model')
    parser.add_argument(\
            '-ru', default = False, required = False, \
            help = 'unaligned reference set')
    parser.add_argument(\
            '--ssu-cmdb', default = False, \
            help = 'CM for 16S/18S HMM search')
    parser.add_argument(\
            '-t', required = False, default = 6, \
            help = 'number of threads (default = 6)')
    parser.add_argument(\
            '--cluster', action = 'store_true', \
            help = 'run on cluster')
    parser.add_argument(\
            '-o', default = False, required = False, \
            help = 'name for output directory')
    args = vars(parser.parse_args())
    seqs, s_db, hits = args['i'], args['d'], args['hits']
    rb, ra, re, ru, prok = args['rb'], args['ra'], args['re'], args['ru'], args['prok']
    no_search, no_ref = args['no_search'], args['no_ref']
    archaea, bacteria, eukarya = args['archaea'], args['bacteria'], args['eukarya']
    max_ins, length, raxml, out = args['max_ins'], args['l'], args['raxml'], args['o']
    ignore_inserts, ssucmdb, threads = args['ignore_insertions'], args['ssu_cmdb'], args['t']
    trim_seqs, cluster = args['trim'], args['cluster']
    no_fast = args['no_fast']
    if cluster is True:
        command = [i for i in sys.argv if i != "--cluster"]
        command.extend(['-t', '48'])
        command = 'echo \"cd %s; %s\" | qsub -m e -l nodes=1:ppn=24' % (os.getcwd(), ' '.join(command))
        p = subprocess.Popen(command, shell = True)
        p.communicate()
        exit()
    s_db, ssucmdb, rb, ra, re, ru, out = \
            define_dbs(s_db, no_search, ssucmdb, \
            ignore_inserts, rb, ra, re, ru, no_ref, prok, \
            archaea, bacteria, eukarya, out)
    print_databases(s_db, ssucmdb, rb, ra, re, ru, no_ref)
    tree(seqs, trim_seqs, s_db, hits, rb, ra, re, ru, no_search, no_ref, \
            archaea, bacteria, eukarya, length, raxml, ignore_inserts, \
            ssucmdb, max_ins, out, threads, no_fast)
