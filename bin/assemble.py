#!/usr/bin/env python3

import os
import sys
import subprocess
from glob import glob as glob

"""
script for assembling reads using velvet by exploring different parameters
"""

def run(out, outdir, command, mv = False, special = False, silent = False):
    outdir = outdir.rsplit('/', 1)[0]
    if os.path.exists(out) is False or special is True:
        command.append('>>%s/velvet.log' % (outdir))
        if silent is True:
            command.append('2>>%s/velvet.log' % (outdir))
        command = ' '.join([str(i) for i in command])
        if silent is False:
            print(command)
        os.system(command)
        if mv is not False:
            contigs = '%s/%s/contigs.fa' % (mv[0], mv[1])
            if os.path.exists(contigs) is True:
                subprocess.Popen(['mv', contigs, out]).wait()

def check_type(reads):
    for line in open(reads):
        if line.startswith('>') is True:
            type = 'fasta'
            break
        elif line.startswith('@') is True:
            type = 'fastq'
            break
    else:
        type = False
    return type

def velveth(out, kmer, paired, single, unused_reads, silent):
    out = '%s/%s' % (out, kmer)
    vh = ['velveth', out, kmer]
    for reads in [paired, 'shortPaired'], [single, 'short']:
        for count, file in enumerate(reads[0]):
            if count > 0:
                cat = '%s%s' % (reads[1], count)
            else:
                cat = reads[1]
            type = check_type(file)
            if type is False:
                continue
            vh.extend(['-' + type, '-' + cat, file])
    out_file = '%s/stats.txt' % (out)
    run(out_file, out, vh, special = unused_reads, silent = silent)
    return out

def velvetg(paired, single, out, out_file, vh, kmer, exp_cov, max_coverage, cov_cutoff, unused_reads, read_trkg, min_contig, scaffolding, silent):
    vg = ['velvetg', vh]
    if min_contig is not False:
        vg.append('-min_contig_lgth %s' % (min_contig))
    if max_coverage is not False:
        vg.append('-max_coverage %s' % (max_coverage))
    if cov_cutoff is not False:
        vg.append('-cov_cutoff %s' % (cov_cutoff))
    if unused_reads is not False:
        unused = '%s/%s/UnusedReads.fa' % (out, kmer)
        vg.append('-unused_reads yes')
    if read_trkg is not False:
        vg.append('-read_trkg yes')
    if exp_cov is not False:
        vg.append('-exp_cov %s' % (exp_cov))
    else:
        vg.append('-exp_cov auto')
    if scaffolding is True:
        vg.append('-scaffolding yes')
    else:
        vg.append('-scaffolding no')
    run(out_file, out, vg, mv = [out, kmer], silent = silent)
    if unused_reads is False:
        return paired, single
    else:
        return [], [unused]

def velvet(paired = [], single = [], out = 'velvet', \
            kmer_min = 21, kmer_max = 71, \
            kmer_increase = 10, min_contig = 400, \
            exp_cov = False, cov_diff = 100, \
            max_coverage = False, cov_cutoff = 'auto', \
            unused_reads = False, read_trkg = False, scaffolding = True, \
            silent = False):
    subprocess.Popen(['mkdir',  '-p', out]).wait()
    kmer_max = kmer_max + kmer_increase
    for kmer in range(kmer_min, kmer_max)[::kmer_increase]:
        vh = velveth(out, kmer, paired, single, unused_reads, silent)
        if exp_cov is False:
            out_file = '%s/kmer_%s.fasta' % (out, kmer)
            paired, single = velvetg(paired, single, out, out_file, vh, kmer, \
                                        exp_cov, max_coverage, cov_cutoff, unused_reads, \
                                        read_trkg, min_contig, scaffolding, silent)
        else:
            min, max = exp_cov - cov_diff, exp_cov + 2*cov_diff
            if min == max:
                cov = exp_cov
                out_file = '%s/kmer_%s_cov_%s.fasta' % (out, kmer, cov)
                paired, single = velvetg(paired, single, out, out_file, vh, \
                                            kmer, cov, max_coverage, cov_cutoff, \
                                            unused_reads, read_trkg, min_contig, scaffolding, \
                                            silent)
                if unused_reads is True:
                    vh = velveth(out, kmer, paired, single, unused_reads, silent)
            else:
                for cov in range(min, max)[::cov_diff]:
                    out_file = '%s/kmer_%s_cov_%s.fasta' % (out, kmer, cov)
                    paired, single = velvetg(paired, single, out, out_file, vh, \
                                                kmer, cov, max_coverage, cov_cutoff, \
                                                unused_reads, read_trkg, min_contig, \
                                                scaffolding, silent)
                    if unused_reads is True:
                        vh = velveth(out, kmer, paired, single, unused_reads, silent)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print('specify scaffolding (True or False), paired reads, single reads, and name for output directory')
        exit()
    scaffold, paired, single, out = sys.argv[1], glob(sys.argv[2]), glob(sys.argv[3]), sys.argv[4]
    if scaffold == 'False' or scaffold == 'false' or scaffold == 'FALSE':
        scaffold = False
    else:
        scaffold = True
    velvet(paired = paired, single = single, out = out, scaffolding = scaffold)
