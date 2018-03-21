#!/usr/bin/env python

from distutils.core import setup

version = '0.14'

packages = ['ctbBio', 'ctbBio27']

scripts = ['ctbBio/16SfromHMM.py', 'ctbBio/23SfromHMM.py', 'ctbBio/numblast.py',
           'ctbBio/calculate_coverage.py', 'ctbBio/cluster_ani.py', 'ctbBio/compare_aligned.py',
           'ctbBio/concat_align.py', 'ctbBio/crossmap.py', 'ctbBio/fasta.py', 'ctbBio/fasta_length.py',
           'ctbBio/fasta_region.py', 'ctbBio/fasta_stats.py', 'ctbBio/fastq2fasta.py', 'ctbBio/fastq_merge.py',
           'ctbBio/fastq_split.py', 'ctbBio/filter_fastq_sam.py', 'ctbBio/fix_fasta.py', 'ctbBio/genome_abundance.py',
           'ctbBio/genome_coverage.py', 'ctbBio/genome_variation.py', 'ctbBio/lookup-word.py', 'ctbBio/lookup.py',
           'ctbBio/mapped.py', 'ctbBio/n50.py', 'ctbBio/name2fasta.py', 'ctbBio/neto.py', 'ctbBio/rec_best_blast.py',
           'ctbBio/nr_fasta.py', 'ctbBio/numblast-pident.py', 'ctbBio/orthologer.py', 'ctbBio/orthologer_summary.py',
           'ctbBio/parallel.py', 'ctbBio/rRNA_copies.py', 'ctbBio/rRNA_insertions.py', 'ctbBio/rax.py',
           'ctbBio/rc.py', 'ctbBio/rp16.py', 'ctbBio/rp16_retreive.sh', 'ctbBio/sam2fastq.py', 'ctbBio/search.py',
           'ctbBio/shuffle_genome.py', 'ctbBio/sixframe.py', 'ctbBio/stats.py', 'ctbBio/stockholm2fa.py',
           'ctbBio/stockholm2oneline.py', 'ctbBio/strip_align.py', 'ctbBio/strip_align_inserts.py',
           'ctbBio/strip_masked.py', 'ctbBio/subset_sam.py', 'ctbBio/transform.py',
           'ctbBio/unmapped.py', 'ctbBio27/id2tax.py', 'ctbBio27/kegginfo.py', 'ctbBio27/ssu_tree.py']

classifiers = ['Programming Language :: Python', 'Programming Language :: Python :: 3', 'Programming Language :: Python :: 2.7']

requirements = ['networkx', 'python-Levenshtein', 'numpy', 'pandas', 'biopython'],

setup(name='ctbBio',
      author='Chris Brown',
      author_email='ctb@berkeley.edu',
      packages=packages,
      scripts=scripts,
      version=version,
      license='MIT',
      url='https://github.com/christophertbrown/bioscripts',
      description='scripts for working with sequencing data',
      install_requires=requirements,
      classifiers=classifiers
      )
