#!/usr/bin/env python3

"""
convert rRNA_insertions.py iTable to gff
"""

import os
import sys
import argparse
import pandas as pd
from ctbBio.fasta import iterate_fasta as parse_fasta

def parse_catalytic(insertion, gff):
    """
    parse catalytic RNAs to gff format
    """
    offset = insertion['offset']
    GeneStrand = insertion['strand']
    if type(insertion['intron']) is not str:
        return gff
    for intron in parse_fasta(insertion['intron'].split('|')):
        ID, annot, strand, pos = intron[0].split('>')[1].split()
        Start, End = [int(i) for i in pos.split('-')]
        if strand != GeneStrand:
            if strand == '+':
                strand = '-'
            else:
                strand = '+'
            Start, End = End - 2, Start - 2
        Start, End = abs(Start + offset) - 1, abs(End + offset) - 1
        gff['#seqname'].append(insertion['ID'])
        gff['source'].append('Rfam')
        gff['feature'].append('Catalytic RNA')
        gff['start'].append(Start)
        gff['end'].append(End)
        gff['score'].append('.')
        gff['strand'].append(strand)
        gff['frame'].append('.')
        gff['attribute'].append('ID=%s; Name=%s' % (ID, annot))
    return gff

def parse_orf(insertion, gff):
    """
    parse ORF to gff format
    """
    offset = insertion['offset']
    if type(insertion['orf']) is not str:
        return gff
    for orf in parse_fasta(insertion['orf'].split('|')):
        ID = orf[0].split('>')[1].split()[0]
        Start, End, strand = [int(i) for i in orf[0].split(' # ')[1:4]]
        if strand == 1:
            strand = '+'
        else:
            strand = '-'
        GeneStrand = insertion['strand']
        if strand != GeneStrand:
            if strand == '+':
                strand = '-'
            else:
                strand = '+'
            Start, End = End - 2, Start - 2
        Start, End = abs(Start + offset) - 1, abs(End + offset) - 1
        annot = orf[0].split()[1]
        if annot == 'n/a':
            annot = 'unknown'
        gff['#seqname'].append(insertion['ID'])
        gff['source'].append('Prodigal and Pfam')
        gff['feature'].append('CDS')
        gff['start'].append(Start)
        gff['end'].append(End)
        gff['score'].append('.')
        gff['strand'].append(strand)
        gff['frame'].append('.')
        gff['attribute'].append('ID=%s; Name=%s' % (ID, annot))
    return gff

def parse_insertion(insertion, gff):
    """
    parse insertion to gff format
    """
    offset = insertion['offset']
    for ins in parse_fasta(insertion['insertion sequence'].split('|')):
        strand = insertion['strand']
        ID = ins[0].split('>')[1].split()[0]
        Start, End = [int(i) for i in ins[0].split('gene-pos=', 1)[1].split()[0].split('-')]
        Start, End = abs(Start + offset), abs(End + offset)
        if strand == '-':
            Start, End = End, Start
        gff['#seqname'].append(insertion['ID'])
        gff['source'].append(insertion['source'])
        gff['feature'].append('IVS')
        gff['start'].append(Start)
        gff['end'].append(End)
        gff['score'].append('.')
        gff['strand'].append(strand) # same as rRNA
        gff['frame'].append('.')
        gff['attribute'].append('ID=%s' % (ID))
    return gff

def parse_masked(seq, min_len):
    """
    parse masked sequence into non-masked and masked regions
    """
    nm, masked = [[]], [[]]
    prev = None
    for base in seq[1]:
        if base.isupper():
            nm[-1].append(base)
            if masked != [[]] and len(masked[-1]) < min_len:
                nm.extend(masked[-1])
                del masked[-1]
            prev = False
        elif base.islower():
            if prev is False:
                masked.append([])
                nm.append([])
            masked[-1].append(base)
            prev = True
    return nm, masked

def parse_rRNA(insertion, seq, gff):
    """
    parse rRNA to gff format
    """
    offset = insertion['offset']
    strand = insertion['strand']
    for rRNA in parse_masked(seq, 0)[0]:
        rRNA = ''.join(rRNA)
        Start = seq[1].find(rRNA) + 1
        End = Start + len(rRNA) - 1
        if strand == '-':
            Start, End = End - 2, Start - 2
        pos = (abs(Start + offset) - 1, abs(End + offset) - 1)
        Start, End = min(pos), max(pos)
        source = insertion['source']
        annot = '%s rRNA' % (source.split('from', 1)[0])
        gff['#seqname'].append(insertion['ID'])
        gff['source'].append(source)
        gff['feature'].append('rRNA')
        gff['start'].append(Start)
        gff['end'].append(End)
        gff['score'].append('.')
        gff['strand'].append(strand)
        gff['frame'].append('.')
        gff['attribute'].append('Name=%s' % (annot))
    return gff

def iTable2GFF(iTable, fa, contig = False):
    """
    convert iTable to gff file
    """
    columns = ['#seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    gff = {c:[] for c in columns}
    for insertion in iTable.iterrows():
        insertion = insertion[1]
        if insertion['ID'] not in fa:
            continue
        # rRNA strand
        strand = insertion['sequence'].split('strand=', 1)[1].split()[0]
        # set rRNA positions for reporting features on contig or extracted sequence
        if contig is True:
            gene = [int(i) for i in insertion['sequence'].split('pos=', 1)[1].split()[0].split('-')]
            if strand == '-':
                offset = -1 * (gene[1])
            else:
                offset = gene[0]
        else:
            strand = '+'
            gene = [1, int(insertion['sequence'].split('total-len=', 1)[1].split()[0])]
            offset = gene[0]
        insertion['strand'] = strand
        insertion['offset'] = offset
        # source for prediction
        source = insertion['sequence'].split('::model', 1)[0].rsplit(' ', 1)[-1]
        insertion['source'] = source
        # rRNA gene
        geneAnnot = '%s rRNA gene' % (source.split('from', 1)[0])
        geneNum = insertion['sequence'].split('seq=', 1)[1].split()[0]
        gff['#seqname'].append(insertion['ID'])
        gff['source'].append(source)
        gff['feature'].append('Gene')
        gff['start'].append(gene[0])
        gff['end'].append(gene[1])
        gff['score'].append('.')
        gff['strand'].append(strand)
        gff['frame'].append('.')
        gff['attribute'].append('ID=%s; Name=%s' % (geneNum, geneAnnot))
        # rRNA
        gff = parse_rRNA(insertion, fa[insertion['ID']], gff)
        # insertions
        gff = parse_insertion(insertion, gff)
        # orfs
        gff = parse_orf(insertion, gff)
        # catalytic RNAs
        gff = parse_catalytic(insertion, gff)
    return pd.DataFrame(gff)[columns].drop_duplicates()

def name2id(name):
    """
    convert header to id (check gene #)
    """
    try:
        return '%s_%s' % (name.split()[0], name.split('seq=', 1)[1].split()[0])
    except:
        print('name error:', name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='# convert rRNA_insertions.py iTable to gff file')
    parser.add_argument(\
            '-i', required = True, \
            help = 'path to rRNA iTable')
    parser.add_argument(\
            '-f', required = True, \
            help = 'path to rRNA fasta file')
    parser.add_argument(\
            '--contigs', action = 'store_true', \
            help = 'report positions relative to contigs, not rRNA genes')
    args = vars(parser.parse_args())
    # load iTable
    rn = {'insertion.1':'insertion sequence', '#sequence':'sequence'}
    iTable = pd.read_csv(args['i'], sep = '\t').rename(columns = rn)
    if args['contigs'] == True:
        iTable['ID'] = [name.split()[0] for name in iTable['sequence']]
    else:
        iTable['ID'] = [name2id(name) for name in iTable['sequence']]
    iTable['insertion ID'] = [name.split('>')[1].split()[0] for name in iTable['insertion sequence']]
    # load sequences
    fa = args['f']
    if args['contigs'] == True:
        fa = {seq[0].split('>')[1].split()[0]:seq for seq in parse_fasta(fa)}
    else:
        fa = {name2id(seq[0].split('>')[1]):seq for seq in parse_fasta(fa)}
    if args['contigs'] is False:
        gffFile = '%s.gff' % (args['i'].rsplit('.', 1)[0])
        # print fasta with seq IDs
        faFile = '%s.gff.fa' % (args['i'].rsplit('.', 1)[0])
        faFile = open(faFile, 'w')
        for ID, seq in fa.items():
            print('>%s' % (ID), file = faFile)
            print(seq[1], file = faFile)
        # convert to gff
    else:
        gffFile = '%s.contigs.gff' % (args['i'].rsplit('.', 1)[0])
    gff = iTable2GFF(iTable, fa, contig = args['contigs'])
    gff.to_csv(open(gffFile, 'w'), sep = '\t', index = False)
