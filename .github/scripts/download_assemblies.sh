#!/bin/bash

mkdir -p .github/data/assemblies

rm -f .github/data/assemblies/GCF_024700185.1.zip
rm -f .github/data/assemblies/GCF024700185.1.fa
rm -f .github/data/assemblies/README.md

curl -o .github/data/assemblies/GCF_024700185.1.zip  "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_024700185.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"

unzip .github/data/assemblies/GCF_024700185.1.zip -d .github/data/assemblies

mv .github/data/assemblies/ncbi_dataset/data/GCF_024700185.1/GCF_024700185.1_ASM2470018v1_genomic.fna .github/data/assemblies/GCF024700185.1.fa

rm -r .github/data/assemblies/ncbi_dataset
rm -f .github/data/assemblies/README.md
