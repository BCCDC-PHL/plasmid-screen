#!/bin/bash

echo 'ID,R1,R2,ASSEMBLY' > .github/data/samplesheet.csv

for i in $(ls ${PWD}/.github/data/fastq/*_R1.fastq.gz); do
  ID=$(basename $i _R1.fastq.gz)
  R1=$i
  R2=${PWD}/.github/data/fastq/${ID}_R2.fastq.gz
  ASSEMBLY=${PWD}/.github/data/assemblies/${ID}.fa
  echo "$ID,$R1,$R2,$ASSEMBLY" >> .github/data/samplesheet.csv
done
