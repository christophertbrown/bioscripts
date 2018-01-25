# bioscripts

some useful scripts for working with genomics and sequencing data

## rRNA identification using 16SfromHMM.py and 23SfromHMM.py

The scripts use `cmsearch` from the Infernal package to do a sequence-based HMM search against 16S and 23S covariance models. The curated model from SSU-Align is used for 16S, and a custom-built model for 23S. The method is similar to what SSU-Align does by default, but accounts for the fact that rRNA genes may contain large insertion sequences. 

The method is described in:

["Unusual biology across a group comprising more than 15% of domain Bacteria"](http://dx.doi.org/10.1038/nature14486) - Christopher T. Brown, Laura A. Hug, Brian C. Thomas, Itai Sharon, Cindy J. Castelle, Andrea Singh, Michael J. Wilkins, Kelly C. Wrighton, Kenneth H. Williams & Jillian F. Banfield (*Nature* 2015).

If using this software, please cite our paper along with Infernal and SSU-Align. 

E. P. Nawrocki, D. L. Kolbe, and S. R. Eddy, "Infernal 1.0: inference of RNA alignments.," Bioinformatics, vol. 25, no. 10, pp. 1335–1337, May 2009.

E. P. Nawrocki, "Structural RNA Homology Search and Alignment using Covariance Models," Washington University in Saint Louis, School of Medicine, 2009.

### requirements

Scripts require Python 3 and Infernal. rRNA_insertions.py also requires HMMER 3 and Pfam (use databases env. variable: $databases/pfam/Pfam-A.hmm). 

### databases

16S CM: databases/ssu-align-0p1.1.cm

23S CM: databases/23S.cm

* use env. variable to reference databases (optional)

`export ssucmdb="databases/ssu-align-0p1.1.cm"`

`export lsucmdb="databases/23S.cm"`

### example usage for finding and analyzing rRNA insertions

* find 16S rRNA genes and insertions

`$ 16SfromHMM.py -f <seqs.fa> -m -t 6 > <seqs.16S.fa>`

* remove insertions (useful for phylogenetic analyses)

`$ strip_masked.py -f <sequences.16S.fa> -l 10`

note: -l 10 specifies that insertions >= 10 bp are removed (this is what I usually do)

### analyze insertions

`$ rRNA_insertions.py <seqs.16S.fa> <False or tax_table.tsv>`

## id2tax.py

id2tax.py is a script for getting the NCBI hierarchical taxonomy for a lineage based on the name of the lineage, an NCBI TaxID, or an organism GI number.

* requires python 2.7 and that tokyocabinet be installed within your python 2.7 path
* requires "databases" env. variable that specifies the path where a "ncbi" directory will be created for storing taxonomy databases
* taxonomy databases are created the first time the script is run (usually takes ~1 day to generate) 
* Note: taxonomy databases are downloaded from NCBI and will need to be updated periodically
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
