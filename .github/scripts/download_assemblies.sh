#!/bin/bash

mkdir -p .github/data/assemblies

curl -o .github/data/assemblies/CP003200.1.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=CP003200.1&db=nucleotide&rettype=fasta"
