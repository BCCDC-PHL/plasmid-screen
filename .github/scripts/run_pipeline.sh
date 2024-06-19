#!/bin/bash

set -eo pipefail

sed -i 's/cpus = 16/cpus = 4/g' nextflow.config 

nextflow run main.nf \
	 -profile conda \
	 --cache ${HOME}/.conda/envs \
	 --fastq_input .github/data/fastq \
	 --mob_db ${PWD}/.github/data/mob-suite-db \
	 --prokka \
	 --collect_outputs \
	 --collected_outputs_prefix test \
	 --outdir .github/data/test_output \
	 -with-report .github/data/test_output/nextflow_report.html \
 	 -with-trace .github/data/test_output/nextflow_trace.tsv
