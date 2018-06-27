#!/bin/bash

if [ "$#" -lt 1 ]
then
        echo "usage: entrez_genome.sh <accession>"
        exit 1
fi

esearch -db genome -query $1 | efetch -db BioSample -format docsum
