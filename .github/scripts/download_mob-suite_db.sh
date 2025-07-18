#!/bin/bash

source ${HOME}/.bashrc

eval "$(conda shell.bash hook)"

mkdir -p .github/data

rm -rf .github/data/mob-suite-db

pushd .github/data

wget -O data.tar.gz https://zenodo.org/records/10304948/files/data.tar.gz?download=1

tar -xzf data.tar.gz

rm data.tar.gz

mv data mob-suite-db

conda activate mash-and-blast

mash sketch -i mob-suite-db/ncbi_plasmid_full_seqs.fas

makeblastdb -in mob-suite-db/ncbi_plasmid_full_seqs.fas -dbtype nucl

popd
