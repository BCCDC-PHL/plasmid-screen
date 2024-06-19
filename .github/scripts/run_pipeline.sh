#!/bin/bash

set -eo pipefail

if [ -n "${GITHUB_ACTIONS}" ]; then
    echo "Running in GitHub Actions Environment"
    echo "Adjusting nextflow.config"
    sed -i 's/cpus = 16/cpus = 4/g' nextflow.config 
else
    echo "Not running in GitHub Actions Environment"
fi

nextflow run main.nf \
	 -profile conda \
	 --cache ${HOME}/.conda/envs \
	 --samplesheet_input .github/data/samplesheet.csv \
	 --pre_assembled \
	 --mob_db ${PWD}/.github/data/mob-suite-db \
	 --collect_outputs \
	 --collected_outputs_prefix test \
	 --outdir .github/data/test_output \
	 -with-report .github/data/test_output/nextflow_report.html \
 	 -with-trace .github/data/test_output/nextflow_trace.tsv
