#!/bin/bash

mkdir -p .github/data

pushd .github/data

wget -O data.tar.gz https://zenodo.org/records/10304948/files/data.tar.gz?download=1

tar -xzf data.tar.gz

rm data.tar.gz

mv data mob-suite-db

popd
