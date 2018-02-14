#!/bin/bash

if [ "$#" -lt 1 ]
then
    echo "usage: ribosomal_get_sequences.sh <rp16.py table> <orfs.faa>"
    exit 1
fi

table=$1
fa=$2

IFS=$'\n'

for i in $(cat $table | grep -v database | head -n 1 | cut -f 2- | tr '\t' '\n' | cat -n)
do 
	out=$(echo `basename $fa .faa`.rp$(echo $i | cut -f 2).faa)
	c=$(echo 1+$(echo $i | rev | cut -d ' ' -f 1  | cut -f 2 | rev) | bc)
	cat $table | grep -v "^#" | cut -f $c | pullseq -i $fa -N > $out
done
