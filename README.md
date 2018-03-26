# bioscripts

some useful scripts for working with genomics and sequencing data

see also [bioscripts27](https://github.com/christophertbrown/bioscripts27)

# installation

`pip install ctbBio`

## rRNA identification using 16SfromHMM.py and 23SfromHMM.py

The scripts use `cmsearch` from the Infernal package to do a sequence-based HMM search against 16S and 23S covariance models. The curated model from SSU-Align is used for 16S, and a custom-built model for 23S. The method is similar to what SSU-Align does by default, but accounts for the fact that rRNA genes may contain large insertion sequences. 

The method is described in:

["Unusual biology across a group comprising more than 15% of domain Bacteria"](http://dx.doi.org/10.1038/nature14486) - Christopher T. Brown, Laura A. Hug, Brian C. Thomas, Itai Sharon, Cindy J. Castelle, Andrea Singh, Michael J. Wilkins, Kelly C. Wrighton, Kenneth H. Williams & Jillian F. Banfield (*Nature* 2015).

If using this software, please cite our paper along with Infernal and SSU-Align. 

E. P. Nawrocki, D. L. Kolbe, and S. R. Eddy, "Infernal 1.0: inference of RNA alignments.," Bioinformatics, vol. 25, no. 10, pp. 1335–1337, May 2009.

E. P. Nawrocki, "Structural RNA Homology Search and Alignment using Covariance Models," Washington University in Saint Louis, School of Medicine, 2009.

### requirements

* python3 
* [infernal](http://eddylab.org/infernal/)
* rRNA_insertions.py requires HMMER3 and Pfam (use databases env. variable: $databases/pfam/Pfam-A.hmm). 

### databases

* 16S CM: databases/ssu-align-0p1.1.cm
* 23S CM: databases/23S.cm

### use env. variable to reference databases (optional)

`export ssucmdb="databases/ssu-align-0p1.1.cm"`

`export lsucmdb="databases/23S.cm"`

### example usage for finding and analyzing rRNA insertions

* find 16S rRNA genes and insertions

`$ 16SfromHMM.py -f <seqs.fa> -m -t 6 > <seqs.16S.fa>`

* remove insertions (useful for phylogenetic analyses)

`$ strip_masked.py -f <seqs.16S.fa> -l 10 > <seqs.16S.s10.fa>`

note: -l 10 specifies that insertions >= 10 bp are removed

### analyze insertions

`$ rRNA_insertions.py <seqs.16S.fa> <False or tax_table.tsv>`

## ortholog identification between pairs of genomes using orthologer.py 

* orthologer.py conducts reciprocal usearch similarity searches between pairs of provided genomes to identify reciprocal best hits

* genomes can be supplied as either gene or protein multi-fasta files (one file per genome; each ORF must have a unique identifier)

### requirements

* python3
* `usearch`

### usage

`$ orthologer.py <mode> <genome1> <genome2> <genome...> > orthologer.out`

Mode can be either "reference" or "global." 

In "reference" mode all searches will be conducted against the first genome that is listed. In "global" mode all possible pairwise searches are conducted between the listed genomes (# searches = #genomes^#genomes). 
