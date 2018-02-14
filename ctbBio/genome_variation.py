#!/usr/bin/env python3

"""
script for quantifying genommic variation from
metagenomic sequencing read mapping
"""

# python
import os
import sys
import random
import argparse
from itertools import cycle

# biopython
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC

# ctb
from fasta import iterate_fasta as parse_fasta

def parse_sam(sam, qual):
    """
    parse sam file and check mapping quality
    """
    for line in sam:
        if line.startswith('@'):
            continue
        line = line.strip().split()
        if int(line[4]) == 0 or int(line[4]) < qual:
            continue
        yield line

def save_refs(sam, fastas, genomes, s2b):
    """
    genomes = {} # genomes[genome][contig][sample] = {'bp_stats':[]}
    """
    if s2b is False:
        s2b = {}
    for fasta in fastas:
        # if no reference sequence supplied, get length from sam file
        if fasta is False:
            for line in sam:
                if line.startswith('@PG'):
                    break
                if line.startswith('@SQ') is False:
                    continue
                line = line.strip().split()
                contig, length = line[1].split(':', 1)[1], int(line[2].split(':', 1)[1])
                if contig not in s2b:
                    genome = 'n/a'
                    s2b[contig] = genome
                else:
                    genome = s2b[contig]
                if genome not in genomes:
                    genomes[genome] = {}
                if contig not in genomes[genome]:
                    genomes[genome][contig] = {}
                bp_stats = []
                for i in range(0, length):
                    bp_stats.append({'A':0, 'T':0, 'G':0, 'C':0, 'N':0, \
                                        'In':[], 'Del':[], 'ref':'N'})
                genomes[genome][contig][sam.name] = {'bp_stats':bp_stats}
                return genomes, s2b
        # save reference sequences, if available
        genome = fasta.rsplit('.', 1)[0]
        if genome not in genomes:
            genomes[genome] = {}
        bp_stats = []
        for seq in parse_fasta(fasta):
            contig = seq[0].split('>')[1].split()[0]
            s2b[contig] = genome
            if contig not in genomes[genome]:
                genomes[genome][contig] = {}
            bp_stats = []
            for base in seq[1]:
                bp_stats.append({'A':0, 'T':0, 'G':0, 'C':0, 'N':0, \
                                    'In':[], 'Del':[], 'ref':base.upper()})
            genomes[genome][contig][sam.name] = {'bp_stats':bp_stats}
        return genomes, s2b

def parse_cigar(cigar):
    """
    parse cigar string into list of operations
      e.g.: 28M1I29M2I6M1I46M ->
        [['28', 'M'], ['1', 'I'], ['29', 'M'], ['2', 'I'], ['6', 'M'], ['1', 'I'], ['46', 'M']]
    """
    cigar = cigar.replace('M', 'M ').replace('I', 'I ').replace('D', 'D ').split()
    cigar = [c.replace('M', ' M').replace('I', ' I').replace('D', ' D').split() for c in cigar]
    return [(int(c[0]), c[1]) for c in cigar]

def get_bp_stats(sam, genomes, s2b, qual):
    """
    save base pair frequency statistics for each base from sam file
    genomes = {} # genomes[genome][contig][sample] = {'bp_stats':{}}
    """
    for read in parse_sam(sam, qual):
        ref = read[2] # scaffold that read mapped to
        if ref not in s2b:
            continue
        genome = s2b[ref] # genome that scaffold belongs to
        refs = genomes[genome][ref][sam.name]['bp_stats']
        ref_start = int(read[3]) - 1 # position of start of alignment on reference
        ref_pos = int(read[3]) - 1 # for keeping track of reference region
        sequence = list(read[9]) # read sequence
        cigar = parse_cigar(read[5]) # parsed cigar string
        cigar_start = 0 # for keeping track of start of cigar regions
        bases = [] # bases to compare with reference
        for cigar_pos, status in cigar:
            if status == 'D': # deletion compared to reference
                refs[ref_pos - 1]['Del'].append(cigar_pos)
                for b in range(0, cigar_pos):
                    bases.append(False)
                ref_pos += cigar_pos
            else:
                cigar_stop = cigar_start + cigar_pos
                if status == 'M': # aligned to reference
                    for b in sequence[cigar_start:cigar_stop]: # bases for cigar region
                        bases.append(b)
                    ref_pos += cigar_pos
                elif status == 'I': # insertion compared to reference
                    refs[ref_pos - 1]['In'].append(cigar_pos)
                else:
                    print('# unrecognized cigar character: %s' % (status), file=sys.stderr)
                    exit()
                cigar_start += cigar_pos
        # add base to frequency at each position
        for base, position in zip(bases, list(range(ref_start, ref_start + len(bases)))):
            if base is False:
                continue
            try:
                refs[position][base.upper()] += 1
            except IndexError:
                continue
    return genomes

def rc_stats(stats):
    """
    reverse completement stats
    """
    rc_nucs = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
    rcs = []
    for pos in reversed(stats):
        rc = {}
        rc['reference frequencey'] = pos['reference frequency']
        rc['consensus frequencey'] = pos['consensus frequency']
        rc['In'] = pos['In']
        rc['Del'] = pos['Del']
        rc['ref'] = rc_nucs[pos['ref']]
        rc['consensus'] = (rc_nucs[pos['consensus'][0]], pos['consensus'][1])
        for base, stat in list(pos.items()):
            if base in rc_nucs:
                rc[rc_nucs[base]] = stat
        rcs.append(rc)
    return rcs

def parse_codons(ref, start, end, strand):
    """
    parse codon nucleotide positions in range start -> end, wrt strand
    """
    codon = []
    c = cycle([1, 2, 3])
    ref = ref[start - 1:end]
    if strand == -1:
        ref = rc_stats(ref)
    for pos in ref:
        n = next(c)
        codon.append(pos)
        if n == 3:
            yield codon
            codon = []

def calc_coverage(ref, start, end, length, nucs):
    """
    calculate coverage for positions in range start -> end
    """
    ref = ref[start - 1:end]
    bases = 0
    for pos in ref:
        for base, count in list(pos.items()):
            if base in nucs:
                bases += count
    return float(bases)/float(length)

def parse_gbk(gbks):
    """
    parse gbk file
    """
    for gbk in gbks:
        for record in SeqIO.parse(open(gbk), 'genbank'):
            for feature in record.features:
                if feature.type == 'gene':
                    try:
                        locus = feature.qualifiers['locus_tag'][0]
                    except:
                        continue
                if feature.type == 'CDS':
                    try:
                        locus = feature.qualifiers['locus_tag'][0]
                    except:
                        pass
                    start = int(feature.location.start) + int(feature.qualifiers['codon_start'][0])
                    end, strand = int(feature.location.end), feature.location.strand
                    if strand is None:
                        strand = 1
                    else:
                        strand = -1
                    contig = record.id
#                    contig = record.id.rsplit('.', 1)[0]
                    yield contig, [locus, \
                            [start, end, strand], \
                            feature.qualifiers]

def parse_fasta_annotations(fastas, annot_tables, trans_table):
    """
    parse gene call information from Prodigal fasta output
    """
    if annot_tables is not False:
        annots = {}
        for table in annot_tables:
            for cds in open(table):
                ID, start, end, strand = cds.strip().split()
                annots[ID] = [start, end, int(strand)]
    for fasta in fastas:
        for seq in parse_fasta(fasta):
            if ('# ;gc_cont' not in seq[0] and '# ID=' not in seq[0]) and annot_tables is False:
                print('# specify fasta from Prodigal or annotations table (-t)', file=sys.stderr)
                exit()
            if 'ID=' in seq[0]:
                ID = seq[0].rsplit('ID=', 1)[1].split(';', 1)[0]
                contig = seq[0].split()[0].split('>')[1].rsplit('_%s' % (ID), 1)[0]
            else:
                contig = seq[0].split()[0].split('>')[1].rsplit('_', 1)[0]
            locus = seq[0].split()[0].split('>')[1]
            # annotation info from Prodigal
            if ('# ;gc_cont' in seq[0] or '# ID=' in seq[0]):
                info = seq[0].split(' # ')
                start, end, strand = int(info[1]), int(info[2]), info[3]
                if strand == '1':
                    strand = 1
                else:
                    strand = -1
                product = [''.join(info[4].split()[1:])]
            # annotation info from table
            else:
                start, end, strand = annots[locus]
                product = seq[0].split(' ', 1)[1]
            info = {'transl_table':[trans_table], \
                    'translation':[seq[1]], \
                    'product':product}
            yield contig, [locus, [start, end, strand], info]

def parse_annotations(annots, fmt, annot_tables, trans_table):
    """
    parse annotations in either gbk or Prodigal fasta format
    """
    annotations = {} # annotations[contig] = [features]
    # gbk format
    if fmt is False:
        for contig, feature in parse_gbk(annots):
            if contig not in annotations:
                annotations[contig] = []
            annotations[contig].append(feature)
    # fasta format
    else:
        for contig, feature in parse_fasta_annotations(annots, annot_tables, trans_table):
            if contig not in annotations:
                annotations[contig] = []
            annotations[contig].append(feature)
    return annotations

def count_mutations(codon, AA, alleles, counts, nucs, trans_table):
    """
    count types of mutations in codon
    counts = {'obs syn':#, 'pos syn':#, 'obs non-syn':#, 'pos non-syn':#}
    """
    # find alternative codons based on SNPs
    obs_codons = [] # codons observed from SNPs
    for pos, pos_mutations in enumerate(alleles):
        mutations = [b for b in pos_mutations if b in nucs]
        for mutation in mutations:
            alt_codon = codon[:]
            alt_codon[pos] = mutation
            obs_codons.append(alt_codon)
    obs_codons = [i for i in obs_codons if i != codon]
    obs_AAs = [codon2aa(i, trans_table) for i in obs_codons]
    obs_AAs = [i for i in obs_AAs if i != AA]
    num_mutations = len(obs_codons)
    num_nonSyn    = len(obs_AAs)
    num_syn       = num_mutations - num_nonSyn
    counts['obs syn']     += num_syn
    counts['obs non-syn'] += num_nonSyn
    # find all possible alternative codons based on single base changes
    pos_codons = [] # codons inferred from all possible SNPs
    for pos in range(0, 3):
        for nuc in nucs:
            pos_codon = codon[:]
            pos_codon[pos] = nuc
            pos_codons.append(pos_codon)
    pos_codons = [i for i in pos_codons if i != codon]
    pos_AAs = [codon2aa(i, trans_table) for i in pos_codons]
    pos_AAs = [i for i in pos_AAs if i != AA]
    num_mutations = len(pos_codons)
    num_nonSyn    = len(pos_AAs)
    num_syn       = num_mutations - num_nonSyn
    counts['pos syn']     += num_syn
    counts['pos non-syn'] += num_nonSyn
    return counts

def codon2aa(codon, trans_table):
    """
    convert codon to amino acid
    """
    return Seq(''.join(codon), IUPAC.ambiguous_dna).translate(table = trans_table)[0]

def calc_PnPs(counts):
    """
    counts = {'obs syn':#, 'pos syn':#, 'obs non-syn':#, 'pos non-syn':#}
    """
    if counts['pos non-syn'] == 0 or counts['pos syn'] == 0:
        return 'n/a'
    NonSyn = float(counts['obs non-syn'])/float(counts['pos non-syn'])
    Syn    = float(counts['obs syn'])    /float(counts['pos syn'])
    if Syn == 0:
        return 'n/a'
    return NonSyn/Syn

def sub_rates(stats, annots):
    """
    calculate syn and non-syn sub. rates for each locus/gene
    """
    nucs = ['A', 'T', 'G', 'C']
    subs = {} # subs[locus]{pn/ps reference, pn/ps consensus, snp density}
    for cds in annots:
        locus, location, info = cds
        subs[locus] = {}
        start, end, strand = location
        trans_table = info['transl_table'][0]
        wrt_reference = {'obs syn':0, 'pos syn':0, 'obs non-syn':0, 'pos non-syn':0}
        wrt_consensus = {'obs syn':0, 'pos syn':0, 'obs non-syn':0, 'pos non-syn':0} 
        for codon in parse_codons(stats, start, end, strand):
            # reference and consensus translations
            ref_codon = [i['ref'] for i in codon]
            con_codon = [i['consensus'][0] for i in codon]
            ref_aa = codon2aa(ref_codon, trans_table)
            con_aa = codon2aa(con_codon, trans_table)
            # count types of mutations with respect to reference
            wrt_reference = \
                    count_mutations(ref_codon, ref_aa, codon, wrt_reference, nucs, trans_table) 
            # count types of mutations with respect to consensus 
            wrt_consensus = \
                    count_mutations(con_codon, con_aa, codon, wrt_consensus, nucs, trans_table)
        # calculate pn/ps
        subs[locus]['ref PnPs']       = calc_PnPs(wrt_reference)
        subs[locus]['consensus PnPs'] = calc_PnPs(wrt_consensus)
        # SNP density
        locus_length = end - start + 1
        subs[locus]['ref SNP density'] \
                = float(wrt_reference['obs syn'] + wrt_reference['obs non-syn']) / locus_length
        subs[locus]['consensus SNP density'] \
                = float(wrt_consensus['obs syn'] + wrt_consensus['obs non-syn']) / locus_length
        # sequencing coverage
        subs[locus]['cov'] = calc_coverage(stats, start, end, locus_length, nucs)
        info['length'] = locus_length
        info['position'] = ((start, end), strand)
        # point to cds info
        subs[locus]['info'] = info
    return subs

def find_consensus(bases):
    """
    find consensus base based on nucleotide
    frequencies
    """
    nucs = ['A', 'T', 'G', 'C', 'N']
    total = sum([bases[nuc] for nuc in nucs if nuc in bases])
    # save most common base as consensus (random nuc if there is a tie)
    try:
        top = max([bases[nuc] for nuc in nucs if nuc in bases])
    except:
        bases['consensus'] = ('N', 'n/a')
        bases['consensus frequency'] = 'n/a'
        bases['reference frequency'] = 'n/a'
        return bases
    top = [(nuc, bases[nuc]) for nuc in bases if bases[nuc] == top]
    if top[0][1] == 0:
        bases['consensus'] = ('n/a', 0)
    else:
        bases['consensus'] = random.choice(top)
    if total == 0:
        c_freq = 'n/a'
        ref_freq = 'n/a'
    else:
        c_freq = float(bases['consensus'][1]) / float(total)
        if bases['ref'] not in bases:
            ref_freq = 0
        else:
            ref_freq = float(bases[bases['ref']]) / float(total)
    bases['consensus frequency'] = c_freq
    bases['reference frequency'] = ref_freq
    return bases

def calc_frequencies(genomes, bp_table, min_cov, min_per):
    """
    print bp frequencies to table
    genomes = {} # genomes[genome][contig][sample] = {'bp_stats':{}}
    """
    nucs = ['A', 'T', 'G', 'C', 'N']
    if bp_table is not False:
        bp_table = open(bp_table, 'w')
        header = ['#genome', 'contig', 'sample', 'position', \
                    'reference', 'ref. frequency', \
                    'consensus', 'con. frequency', \
                    'A', 'T', 'G', 'C', 'N', '# insertions', '# deletions']
        print('\t'.join(header), file=bp_table)
    for genome, contigs in list(genomes.items()):
        for contig, samples in list(contigs.items()):
            for sample, stats in list(samples.items()):
                for pos, ps in enumerate(stats['bp_stats'], 1):
                    coverage = sum([ps[nuc] for nuc in nucs])
                    for nuc in nucs:
                        # make sure support for base passes thresholds
                        nuc_cov = ps[nuc]
                        if coverage == 0:
                            nuc_per = 0
                        else:
                            nuc_per = (float(nuc_cov)/coverage)*100
                        if nuc_cov < min_cov or nuc_per < min_per:
                            del ps[nuc]
                    ps = find_consensus(ps)
                    genomes[genome][contig][sample][pos] = ps
                    if bp_table is not False:
                        out = [genome, contig, sample, pos]
                        for i in ['ref', 'reference frequency', \
                                  'consensus', 'consensus frequency', \
                                  'A', 'T', 'G', 'C', 'N', \
                                  'In', 'Del']:
                            try:
                                if i == 'consensus':
                                    out.append(ps[i][0])
                                elif i in ['In', 'Del']:
                                    out.append(len(ps[i]))
                                else:
                                    out.append(ps[i])
                            except:
                                out.append('n/a')
                        print('\t'.join([str(i) for i in out]), file=bp_table)
    return genomes

def print_consensus(genomes):
    """
    print consensensus sequences for each genome and sample
    """
    # generate consensus sequences
    cons = {} # cons[genome][sample][contig] = consensus
    for genome, contigs in list(genomes.items()):
        cons[genome] = {}
        for contig, samples in list(contigs.items()):
            for sample, stats in list(samples.items()):
                if sample not in cons[genome]:
                    cons[genome][sample] = {}
                seq = cons[genome][sample][contig] = []
                for pos, ps in enumerate(stats['bp_stats'], 1):
                    ref, consensus = ps['ref'], ps['consensus'][0]
                    if consensus == 'n/a':
                        consensus = ref.lower()
                    seq.append(consensus)
    # print consensus sequences
    for genome, samples in cons.items():
        for sample, contigs in samples.items():
            fn = '%s.%s.consensus.fa' % (genome, sample)
            f = open(fn, 'w')
            for contig, seq in contigs.items():
                print('>%s' % (contig), file = f)
                print(''.join(seq), file = f)
            f.close()
    return cons

def print_stats(genomes):
    """
    print substitution rate data to table
    genomes[genome][contig][sample] = \
            {'bp_stats':{}, 'sub_rates'[locus] = {ref PnPs, consensus PnPs}}
    """
    header = ['#genome', 'contig', 'locus', 'position', 'strand', 'length', \
                'sample', 'coverage', \
                'ref. Pn/Ps', 'ref SNP density', \
                'consensus Pn/Ps', 'consensus SNP density', \
                'product']
    print('\t'.join(header))
    for genome, contigs in list(genomes.items()):
        for contig, samples in list(contigs.items()):
            for sample, stats in list(samples.items()):
                for locus, rates in list(stats['sub_rates'].items()):
                    length = rates['info']['length']
                    position, strand = rates['info']['position']
                    position = '%s-%s' % position
                    out = [genome, contig, locus, position, strand, length, \
                            sample, '%.2f' % (rates['cov']), \
                            rates['ref PnPs'], rates['ref SNP density'], \
                            rates['consensus PnPs'], rates['consensus SNP density'], \
                            rates['info']['product'][0]]
                    print('\t'.join([str(i) for i in out]))

def variation(fastas, sams, s2b, annots, qual, min_cov, min_per, no_sub, bp_table, print_con):
    """
    quantify genomic variation from mapping
    genomes[genome][contig][sample] = \
            {'bp_stats':[], 'sub_rates'[locus] = {ref PnPs, consensus PnPs, ...}}
    """
    genomes = {}

    # save base-pair data structures for each genome-sample pair
    for sam in sams:
        genomes, s2b = save_refs(sam, fastas, genomes, s2b)

    # get base-pair frequencies for each genome-sample pair
    for sam in sams:
        genomes = get_bp_stats(sam, genomes, s2b, qual)

    # calculate base-pair frequencies
    genomes = calc_frequencies(genomes, bp_table, min_cov, min_per)

    # print consensus genome
    if print_con is True:
        print('# saving consensus sequences', file = sys.stderr)
        genomes = print_consensus(genomes)

    if no_sub is True:
        return genomes

    # calculate substitution rates
    for genome, contigs in list(genomes.items()):
        for contig, samples in list(contigs.items()):
            for sample in samples:
#                print(genomes.keys())
#                print(genomes[genome].keys())
#                print(genomes[genome][contig].keys())
#                print(annots.keys())
                genomes[genome][contig][sample]['sub_rates'] = \
                    sub_rates(genomes[genome][contig][sample]['bp_stats'], annots[contig])

    return genomes

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = \
            '# quantify genomic variation from sam mapping file')
    parser.add_argument(\
            '-s', nargs = '*', action = 'store', required = True, \
            help = 'path to sam file(s) (required)')
    parser.add_argument(\
            '-f', nargs = '*', action = 'store', default = [False], \
            help = 'path to reference genome fasta(s) (required for calculating substitutions)')
    parser.add_argument(\
            '-b', default = False, \
            help = 'use genome information from scaffold2bin.tsv')
    parser.add_argument(\
            '-a', nargs = '*', action = 'store', default = False, \
            help = 'annotation file(s) (default = GenBank; required for Pn/Ps)')
    parser.add_argument(\
            '--faa', action = 'store_true', \
            help = 'annotations file(s) are in fasta format (protein sequences)')
    parser.add_argument(\
            '-t', nargs = '*', action = 'store', default = False, \
            help = 'orf table: <id> <start> <stop> <strand> (required if fasta not from Prodigal')
    parser.add_argument(\
            '-tt', default = '11', \
            help = 'translation table (use with --faa; default = 11)')
    parser.add_argument(\
            '-c', default = 3, type = int, \
            help = 'minimum coverage required to support base call (default = 3)')
    parser.add_argument(\
            '-p', default = 1, type = int, \
            help = 'minimum percent coverage to support base call (default = 1)')
    parser.add_argument(\
            '-q', default = 10, type = int, \
            help = 'minimum mapping quality score (default = 10)')
    parser.add_argument(\
            '--no-sub', action = 'store_true', \
            help = 'do not calculate non-syn/syn substitutions')
    parser.add_argument(\
            '--con', action = 'store_true', \
            help = 'save consensus sequence(s) (optional)')
    parser.add_argument(\
            '-bp', default = False, \
            help = 'filename for saving per-base nucleotide frequencies (optional)')
    args = vars(parser.parse_args())
    sams, fastas, s2b = args['s'], args['f'], args['b']
    sams = [open(s) for s in sams]
    annots, annot_fa, annot_tables, trans_table = args['a'], args['faa'], args['t'], args['tt']
    min_cov, min_per, qual = args['c'], args['p'], args['q']
    no_sub, bp_file, con = args['no_sub'], args['bp'], args['con']
    if annots is False and no_sub is False:
        print('# annotations are requred for calculating substition rates', file=sys.stderr)
        print('# use --no-sub to calculate only bp frequencies', file=sys.stderr)
        exit()
    if s2b is False and fastas[0] is False and no_sub is False:
        print('# scaffold2bin or fasta file(s) required for calculating substitution rates')
        exit()
    if s2b is not False:
        s2b = \
            {i.strip().split()[0]:i.strip().split()[1] \
                for i in open(s2b)}
    if no_sub is not True:
        annots = parse_annotations(annots, annot_fa, annot_tables, trans_table)
    else:
        annots = False
        if bp_file is False and con is False:
            print('# use -bp to save base-pair frequencies to file', file=sys.stderr)
            print('# and/or use --con to save consensus sequence(s)', file=sys.stderr)
            exit()
    genomes = variation(fastas, sams, s2b, annots, qual, min_cov, min_per, no_sub, bp_file, con)
    if no_sub is not True:
        print_stats(genomes)
