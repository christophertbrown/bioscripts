#!/usr/bin/env python2.7

"""
script for 1) finding 16S rRNA gene sequences from contigs based on ssu-cmsearch or cmsearch tabular output and 2)
identifying 16S rRNA gene insertions

ctb@berkeley.edu
"""

# python modules
import sys
import os
from operator import itemgetter
import argparse
from subprocess import Popen

# ctb modules
from fasta import iterate_fasta as parse_fasta
from rc import reverse_complement as rc

def best_model(seq2hmm):
    """
    determine the best model: archaea, bacteria, eukarya (best score) 
    """
    for seq in seq2hmm:
        best = []
        for model in seq2hmm[seq]:
            best.append([model, sorted([i[-1] for i in seq2hmm[seq][model]], reverse = True)[0]])
        best_model = sorted(best, key = itemgetter(1), reverse = True)[0][0]
        seq2hmm[seq] = [best_model] + [seq2hmm[seq][best_model]]
    return seq2hmm

def check_gaps(matches, gap_threshold = 0):
    """
    check for large gaps between alignment windows
    """
    gaps = []
    prev = None
    for match in sorted(matches, key = itemgetter(0)):
        if prev is None:
            prev = match
            continue
        if match[0] - prev[1] >= gap_threshold:
            gaps.append([prev, match])
            prev = match
    return [[i[0][1], i[1][0]] for i in gaps]

def get_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def check_overlap(current, hit, overlap = 200):
    """
    determine if sequence has already hit the same part of the model,
    indicating that this hit is for another 16S rRNA gene
    """
    for prev in current:
        p_coords = prev[2:4]
        coords = hit[2:4]
        if get_overlap(coords, p_coords) >= overlap:
            return True
    return False

def hit_groups(hits):
    """
    * each sequence may have more than one 16S rRNA gene
    * group hits for each gene
    """
    groups = []
    current = False
    for hit in sorted(hits, key = itemgetter(0)):
        if current is False:
            current = [hit]
        elif check_overlap(current, hit) is True:
            groups.append(current)
            current = [hit]
        else:
            current.append(hit)
    groups.append(current)
    return groups

def find_coordinates(hmms, bit_thresh):
    """
    find 16S rRNA gene sequence coordinates
    """
    # get coordinates from cmsearch output
    seq2hmm = parse_hmm(hmms, bit_thresh)
    seq2hmm = best_model(seq2hmm)
    group2hmm = {} # group2hmm[seq][group] = [model, strand, coordinates, matches, gaps]
    for seq, info in list(seq2hmm.items()):
        group2hmm[seq] = {}
        # info = [model, [hit1], [hit2], ...]
        for group_num, group in enumerate(hit_groups(info[1])): 
            # group is a group of hits to a single 16S gene
            # determine matching strand based on best hit
            best = sorted(group, reverse = True, key = itemgetter(-1))[0]
            strand = best[5]
            coordinates = [i[0] for i in group] + [i[1] for i in group]
            coordinates = [min(coordinates), max(coordinates), strand]
            # make sure all hits are to the same strand
            matches = [i for i in group if i[5] == strand]
            # gaps = [[gstart, gend], [gstart2, gend2]]
            gaps = check_gaps(matches)
            group2hmm[seq][group_num] = [info[0], strand, coordinates, matches, gaps]
    return group2hmm

def get_info(line, bit_thresh):
    """
    get info from either ssu-cmsearch or cmsearch output
    """
    if len(line) >= 18: # output is from cmsearch
        id, model, bit, inc = line[0].split()[0], line[2], float(line[14]), line[16]
        sstart, send, strand = int(line[7]), int(line[8]), line[9]
        mstart, mend = int(line[5]), int(line[6])
    elif len(line) == 9: # output is from ssu-cmsearch
        if bit_thresh == 0:
            print >> sys.stderr, '# ssu-cmsearch does not include a model-specific inclusion threshold, '
            print >> sys.stderr, '# please specify a bit score threshold'
            exit()
        id, model, bit = line[1].split()[0], line[0], float(line[6])
        inc = '!' # this is not a feature of ssu-cmsearch
        sstart, send = int(line[2]), int(line[3])
        mstart, mend = int(4), int(5)
        if send >= sstart:
            strand = '+'
        else:
            strand = '-'
    else:
        print >> sys.stderr, '# unsupported hmm format:'
        print >> sys.stderr, '# provide tabular output from ssu-cmsearch and cmsearch supported'
        exit()
    coords = [sstart, send]
    sstart, send = min(coords), max(coords)
    mcoords = [mstart, mend]
    mstart, mend = min(mcoords), max(mcoords)
    return id, model, bit, sstart, send, mstart, mend, strand, inc

def parse_hmm(hmms, bit_thresh):
    seq2hmm = {}
    for hmm in hmms:
        for line in hmm:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            id, model, bit, sstart, send, mstart, mend, strand, inc = get_info(line, bit_thresh)
            if bit >= bit_thresh and inc == '!':
                if id not in seq2hmm:
                    seq2hmm[id] = {}
                if model not in seq2hmm[id]:
                    seq2hmm[id][model] = []
                length = abs(sstart - send) + 1
                seq2hmm[id][model].append([sstart, send, mstart, mend, length, strand, bit])
    return seq2hmm

def mask_sequence(seq, gaps):
    """
    mask (make lower case) regions of sequence found in gaps between model alignments
    """
    seq = [i.upper() for i in seq]
    for gap in gaps:
        for i in range(gap[0] - 1, gap[1]):
            seq[i] = seq[i].lower()
    return ''.join(seq)

def check_buffer(coords, length, buffer):
    """
    check to see how much of the buffer is being used
    """
    s = min(coords[0], buffer)
    e = min(length - coords[1], buffer)
    return [s, e]

def find_16S(fastas, hmms, bit_thresh = float(20), length_thresh = 500, masking = True, buffer = 0):
    """
    1) parse hmm output into dictionary (sequence must pass bit_thresh and inc == '!')
        seq2hmm[seq] = {model: [sstart, ssend, length, strand, score]} 
    2) determine which model (archaea, bacteria, eukarya) the sequence most closely matches
        seq2hmm[seq] = [model, sstart, send, length, strand, score], [model2, sstart2, send2, length2, strand2, score2], ...]
    3) identify regions that match to 16S (for best model)
    4) mask internal regions that do not align to model
    5) length threshold applies to aligned regions of 16S sequence
    5) export 16S sequnece based on complete gene (including masked insertions)
    """
    # identify start/stop positions
    # group2hmm[seq][group] = [model, strand, coordinates, matches, gaps]
    group2hmm = find_coordinates(hmms, bit_thresh)
    # get sequences from fasta file
    for fasta in fastas:
        for seq in parse_fasta(fasta):
            id = seq[0].split('>')[1].split()[0]
            if id not in group2hmm:
                continue
            seq[1] = seq[1].upper()
            count = 0 # how many 16S genes are there on the contig?
            for group, info in list(group2hmm[id].items()):
                model, strand, coords, matches, gaps = info
                # count insertion bases (ib) from gaps
                ib = sum([i[1] - i[0] + 1 for i in gaps])
                # calcualte length of non-insertion regions (don't include buffer)
                tl = coords[1] - coords[0] + 1
                length = tl - ib
                if length < length_thresh:
                    continue 
                # count sequence
                count += 1
                # set retrieval coords based on buffer
                ret_coords = [max([coords[0] - buffer, 1]), \
                        min([coords[1] + buffer, len(seq[1])]), coords[2]]
                buffer_ends = check_buffer(coords, len(seq[1]), buffer)
                # mask insertion sequences
                if masking is True:
                    seq[1] = mask_sequence(seq[1], gaps)
                S = seq[1][(ret_coords[0] - 1):(ret_coords[1])]
                inserts = [gap[1] - gap[0] + 1 for gap in gaps]
                inserts.append('end')
                model_pos = ';'.join(['%s-%s(%s)' % (match[2], match[3], insert) for match, insert in zip(matches, inserts)])
                header = '%s 16SfromHMM::model=%s seq=%s pos=%s-%s strand=%s total-len=%s 16S-len=%s model-pos(ins-len)=%s buffer-len=%s/%s ins-bases=%s' % \
                        (seq[0], model, count, ret_coords[0], ret_coords[1], strand, tl, length, model_pos, buffer_ends[0], buffer_ends[1], ib)
                # reverse complement if strand is reverse
                if strand == '-':
                    S = rc(['', S])[1]
                yield [header, S]
    
def run_cmsearch(fastas, threads, cm):
    """
    run cmsearch: comapre 16S sequences to ssu-align's CM
    """
    out = []
    for fasta in fastas:
        cmsearch = '%s.16S.cmsearch' % (fasta.name.rsplit('.', 1)[0])
        if os.path.exists(cmsearch) is False:
            p = Popen('\
                    cmsearch --cpu %s --hmmonly --acc --noali -T -1 --tblout %s %s %s >> cmsearch.log' \
                    % (threads, cmsearch, cm, fasta.name), shell = True)
            p.communicate()
        else:
            print >> sys.stderr, '# cmsearch output found: %s' % (cmsearch)
        out.append(open(cmsearch))
    return out

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='# find 16S rRNA genes using cmsearch')
    parser.add_argument(\
            '-f', nargs = '*', action = 'store', \
            required = True, help = 'path to fasta file(s) for searching')
    parser.add_argument(\
            '-c', nargs = '*', action = 'store', default = False, \
            help = 'path to cmsearch table(s)')
    parser.add_argument(\
            '-t', default = False, help = 'number of cpus to use')
    parser.add_argument(\
            '-b', default = 0, help = 'size of flanking region to collect, default = 0')
    parser.add_argument(\
            '-m', action = 'store_true', \
            help = 'mask insertion sequences, default = False')
    parser.add_argument(\
            '-d', default = False, help = 'path to cm')
    parser.add_argument(\
            '-bit', default = 0, type = float, help = 'minimum bit score to consider, default = 0')
    parser.add_argument(\
            '-l', default = 500, type = int, \
            help = 'minimum length 16S rRNA gene sequence to report, default = 500')
    args = vars(parser.parse_args())
    # make sure either a cmsearch table or the number of cpus to use for running cmsearch is specified
    if args['c'] is False and args['t'] is False:
        print >> sys.stderr, 'specify: -c <cmsearch table> or -t <number of cpus for running cmsearch>'
        exit()
    # make sure that fasta is not from stdin if cmsearch table is not provided
    if args['c'] is False and args['f'][0] == '-':
        print >> 'specify: -c <cmsearch table> or -f <fasta file, not stdin>'
        exit()
    # check if either fasta or cmsearch is from stdin
    if args['f'][0] == '-' and args['c'] == False:
        print >> 'specify: -f <fasta file, not stdin> or -c <cmsearch table>'
        exit()
    for i in 'f' 'c':
        if args[i] is False:
            continue
        if args[i][0] == '-':
            args[i] = [sys.stdin]
        else:
            args[i] = [open(j) for j in args[i]]
    # if cmsearch output not specified, run cmsearch
    if args['c'] is False:
        if args['d'] is False and 'ssucmdb' not in os.environ:
            print >> 'specify: -d <cmsearch database>'
            exit()
        if args['d'] is False:
            cm = os.environ['ssucmdb']
            print >> '# 16S rRNA CM: %s' % (cm)
        else:
            cm = args['d']
        args['c'] = run_cmsearch(args['f'], args['t'], cm)
    fastas, hmms, bit_thresh, length_thresh, masking, buffer = \
            args['f'], args['c'], float(args['bit']), int(args['l']), args['m'], int(args['b'])
    for seq in find_16S(fastas, hmms, bit_thresh, length_thresh, masking, buffer):
        print '\n'.join(seq)
