#!/usr/bin/env python

import argparse
import csv
import json
import re
import sys


def parse_resistance_plasmid_report(resistance_plasmid_report_path):
    resistance_plasmid_report = []
    
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
    ]

    int_fields = [
        "num_contigs_in_reconstruction",
        "reconstruction_size",
        "gene_start",
        "gene_end",
    ]

    float_fields = [
        "percent_resistance_gene_coverage",
        "percent_resistance_gene_identity",
        "mash_neighbor_distance",
    ]

    with open(resistance_plasmid_report_path, 'r') as f:
        reader = csv.DictReader(f, fieldnames=fieldnames, dialect="excel-tab")
        next(reader)
        for row in reader:
            for field in int_fields:
                row[field] = int(row[field])
            for field in float_fields:
                row[field] = float(row[field])
            resistance_plasmid_report.append(row)

    return resistance_plasmid_report


def parse_snp_report(snp_report_path):
    snp_report = {}
    
    fieldnames = [
        "sample_id",
        "alignment_ref_plasmid",
        "num_snps",
    ]

    int_fields = [
        "num_snps",
    ]

    with open(snp_report_path, 'r') as f:
        reader = csv.DictReader(f, fieldnames=fieldnames)
        next(reader)
        for row in reader:
            for field in int_fields:
                row[field] = int(row[field])
        snp_report[row['sample_id']] = (row)

    return snp_report


def main(args):
    resistance_plasmid_report = parse_resistance_plasmid_report(args.resistance_plasmid_report)

    snp_report = parse_snp_report(args.snp_report)

    output_fieldnames = [
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
    csv.register_dialect('unix-tab', delimiter='\t', doublequote=False, lineterminator='\n', quoting=csv.QUOTE_MINIMAL)
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='unix-tab')
    writer.writeheader()


    for record in resistance_plasmid_report:
        alignment_ref_plasmid = snp_report[record['sample_id']]['alignment_ref_plasmid']
        num_snps = snp_report[record['sample_id']]['num_snps']
        record['alignment_ref_plasmid'] = alignment_ref_plasmid
        record['num_snps_vs_ref_plasmid'] = num_snps
        writer.writerow(record)
        
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--sample-id', help="")
    parser.add_argument('--resistance-plasmid-report', help="")
    parser.add_argument('--snp-report', help="")
    args = parser.parse_args()
    main(args)
