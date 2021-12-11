#!/usr/bin/env python

import argparse
import csv
import json
import sys


def parse_combined_abricate_mobtyper_report(combined_abricate_mobtyper_report_path):
    combined_abricate_mobtyper_report = []
    
    fieldnames = [
        "sample_id",
        "assembly_file",
        "resistance_gene_contig",
        "num_contigs_in_reconstruction",
        "reconstruction_size",
        "resistance_gene",
        "gene_start",
        "gene_end",
        "percent_resistance_gene_coverage",
        "percent_resistance_gene_identity",
        "replicon_types",
        "mob_suite_primary_cluster_id",
        "mob_suite_secondary_cluster_id",
        "mash_nearest_neighbor",
        "mash_neighbor_distance",
        "alignment_ref_plasmid",
        "num_snps_vs_ref_plasmid",
    ]

    int_fields = [
        "num_contigs",
        "size",
        "gene_start",
        "gene_end",
    ]

    float_fields = [
        "percent_resistance_gene_coverage",
        "percent_resistance_gene_identity",
        "mash_neighbor_distance",
    ]

    with open(combined_abricate_mobtyper_report_path, 'r') as f:
        reader = csv.DictReader(f, fieldnames=fieldnames, dialect="excel-tab")
        next(reader, None)  # skip header
        for row in reader:
            for field in int_fields:
                try:
                    row[field] = int(row[field])
                except Exception as e:
                    row[field] = None
            for field in float_fields:
                try:
                    row[field] = float(row[field])
                except Exception as e:
                    row[field] = None

            if row['mash_nearest_neighbor'] != '-':
                combined_abricate_mobtyper_report.append(row)

    return combined_abricate_mobtyper_report


def main(args):

    combined_abricate_mobtyper_report = parse_combined_abricate_mobtyper_report(args.combined_abricate_mobtyper_report)
    
    fieldnames = [
        "sample_id",
        "assembly_file",
        "reference_plasmid_id",
        "mash_neighbor_distance",
    ]
    
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    writer.writeheader()
    
    for input_record in combined_abricate_mobtyper_report:
        output_record = { 'sample_id': input_record['sample_id'],
                          'assembly_file': input_record['assembly_file'],
                          'reference_plasmid_id': input_record['mash_nearest_neighbor'],
                          'mash_neighbor_distance': input_record['mash_neighbor_distance']}
                  
        writer.writerow(output_record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-r', '--combined-abricate-mobtyper-report', help="")
    args = parser.parse_args()
    main(args)
