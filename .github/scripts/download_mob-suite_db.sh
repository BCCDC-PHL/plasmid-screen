#!/bin/bash

source ${HOME}/.bashrc

eval "$(conda shell.bash hook)"

mkdir -p .github/data

pushd .github/data

wget -O data.tar.gz https://zenodo.org/records/10304948/files/data.tar.gz?download=1

tar -xzf data.tar.gz

rm data.tar.gz

mv data mob-suite-db

conda activate plasmid-screen-35d122a137231eda3b8a0039d42f24f6

mash sketch -i mob-suite-db/ncbi_plasmid_full_seqs.fas

makeblastdb -in mob-suite-db/ncbi_plasmid_full_seqs.fas -dbtype nucl

popd
