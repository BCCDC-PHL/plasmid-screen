#!/bin/bash

rm -rf .github/data/assemblies/*
rm -rf .github/data/fastq/*
rm -rf .github/data/mob-suite-db
rm -rf .github/data/samplesheet.csv
rm -rf .github/data/test_output

.github/scripts/download_assemblies.sh

.github/scripts/simulate_reads.sh

.github/scripts/download_mob-suite_db.sh

.github/scripts/create_samplesheet.sh

.github/scripts/run_pipeline.sh

