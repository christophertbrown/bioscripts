#!/bin/bash

for i in "$1";
do
	echo $i | tr '\n' '\t'; 
	n50.py -i $i | tr '\n' '\t'; 
	echo 'scaffolds: ' $(grep -c "^>" $i) | tr '\n' '\t'; 
	echo 'largest scaffold: ' $(fasta_length.py $i 0 | grep "^>" | cut -f 2 | sort -n -r | head -n 1);
done | tr '\t' '!' | tr2 '!' '!|' | tr '!' '\n' | tr '|' '\t'
