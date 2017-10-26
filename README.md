# bioscripts
some useful scripts for working with genomics and sequencing data

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


## output:

Major taxonomic subdivisions are delineated by numbers, and minor subdivisions by semi-colons (if present). You can provide names from any taxonomic subdivision and it will give you the parent lineages. 

There are some funny issues that I have not been able to work out. For example, the name "Bacillus" is redundant in the database, so you will not get the Bacteria taxonomy unless you specify the species name:

`$ id2tax.py -i Bacillus`

Bacillus	[0]Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Protostomia;Ecdysozoa;Panarthropoda;[1]Arthropoda;Mandibulata;Pancrustacea;Hexapoda;[2]Insecta;Dicondylia;[5]Pterygota;Neoptera;Polyneoptera;[3]Phasmatodea;Verophasmatodea;Areolatae;Bacilloidea;[4]Bacillidae;Bacillinae;Bacillini;[5]Bacillus

`$ id2tax.py -i "Bacillus subtilis"`

Bacillus subtilis	[0]Bacteria;Terrabacteria group;[1]Firmicutes;[2]Bacilli;[3]Bacillales;[4]Bacillaceae;[5]Bacillus;[6]Bacillus subtilis group;[6]Bacillus subtilis

Try to watch out for things like this!
