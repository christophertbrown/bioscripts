# bioscripts

some useful scripts for working with genomics and sequencing data

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

## (quick) phylogenetic analysis using ssu_tree.py

Script for quickly generating a 16/18S rRNA gene tree. The script searches a database, usually Silva, for similar reference sequences (choosing the top 3 hits, by default), removes insertions from all sequences, aligns sequences using SSU-Align, combines these sequences with a pre-aligned reference set, and then runs a quick tree using FastTree2. By default, SSU-Align will create separate alignments for sequences from bacteria, archaea, and eukaryotes. In this case, you end up with three separate trees if your sequence set includes sequences from all three domains.

There are options for customizing trees with respect to reference sets and rooting (use -h option for a list).

This is a combination of several scripts that can be used independently for phylogenetic analyses.

`16SfromHMM.py` find ssu rRNA gene sequences (see also `23SfromHMM.py`)

`rax.py` run FastTree and/or RAxML

`search.py` run blast and/or usearch

`numblast.py` get top blast/usearch hits (b6 format)

`stockholm2fa.py` convert stockholm to fasta

`nr_fasta.py` remove sequences in fasta file with same name, or re-name them

`strip_masked.py` remove masked sequence in fasta file

### requirements

* python2.7
* biopython
* tokyocabinet
* [usearch](https://www.drive5.com/usearch/)
* [SSU-Align](http://eddylab.org/software/ssu-align/)
* [FastTree2](http://www.microbesonline.org/fasttree/)
* [RAxML (optional)](https://sco.h-its.org/exelixis/web/software/raxml/index.html)

### example usage

$ ssu_tree.py -i 16S_sequences.fasta

The output will be stored in a directory with a name based on the name of the fasta file you provide. In this case it would be called “16S_sequences.tree". Tree(s) generated will be in this directory and have "FINAL" in the name and end with “.tree".

If you have not identified the 16S genes, you can give the script scaffolds with the --trim option, for example:

$ ssu_tree.py -i sequences.fasta --trim

By default this will include reference sequences and best hits. If you already have your reference set together, you can turn those features off, for example with bacteria sequences:

$ ssu_tree.py -i <unaligned_sequences.fasta> --no-search --no-ref --bacteria

This will align the sequences and run a FastTree. You can also have it run a raxml tree with the --raxml option.

### env. variables for specifying database locations

#### CMs

export ssucmdb="databases/ssu-align-0p1.1.cm"

export lsucmdb="/data1/bio_db/ssu/23S.cm"

#### ssu rRNA gene database (optional)

export ssuref="PATH TO DATABASE FASTA e.g. [Silva](https://www.arb-silva.de)"

#### pre-aligned reference sets (optional)

export ssual_ref_archaea="databases/ssu_refs_archaea-archaea.afa"

export ssual_ref_bacteria="databases/ssu_refs_bacteria-bacteria.afa"

export ssual_ref_eukarya="databases/ssu_refs_eukarya-eukarya.afa"

export ssual_ref_prok_archaea="databases/ssu_refs_prok-archaea.afa"

export ssual_ref_prok_bacteria="databases/ssu_refs_prok-bacteria.afa"

**Note:** include full paths.

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

## hierarchical taxonomy using id2tax.py

id2tax.py is a script for getting the NCBI hierarchical taxonomy for a lineage based on the name of the lineage, an NCBI TaxID, or an organism GI number.

### requirements

* python2.7 
* tokyocabinet
* "databases" env. variable that specifies the path to where an "ncbi" directory will be created for storing taxonomy databases

### notes

* taxonomy databases are created the first time the script is run (usually takes ~1 day to generate) 
* taxonomy databases are downloaded from NCBI and will need to be updated periodically (delete or move the old databases and re-run id2tax.py)
* for help see `id2tax.py -h`

### example usage:

`$ id2tax.py -i "Bacillus subtilis"`

### or, if you have a list of names in a text file:

`$ cat <list.txt> | id2tax.py -i -`

### output:

Major taxonomic subdivisions are delineated by numbers, and minor subdivisions by semi-colons (if present). You can provide names from any taxonomic subdivision and it will give you the parent lineages. 

There are some funny issues that I have not been able to work out. For example, the name "Bacillus" is redundant in the database, so you will not get the Bacteria taxonomy unless you specify the species name:

`$ id2tax.py -i Bacillus`

Bacillus	[0]Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Protostomia;Ecdysozoa;Panarthropoda;[1]Arthropoda;Mandibulata;Pancrustacea;Hexapoda;[2]Insecta;Dicondylia;[5]Pterygota;Neoptera;Polyneoptera;[3]Phasmatodea;Verophasmatodea;Areolatae;Bacilloidea;[4]Bacillidae;Bacillinae;Bacillini;[5]Bacillus

`$ id2tax.py -i "Bacillus subtilis"`

Bacillus subtilis	[0]Bacteria;Terrabacteria group;[1]Firmicutes;[2]Bacilli;[3]Bacillales;[4]Bacillaceae;[5]Bacillus;[6]Bacillus subtilis group;[6]Bacillus subtilis

Try to watch out for things like this!
