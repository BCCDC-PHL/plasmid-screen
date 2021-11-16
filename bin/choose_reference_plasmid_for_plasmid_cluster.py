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
        "num_contigs",
        "size",
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
            if row['mob_suite_primary_cluster_id'] == '-':
                row['mob_suite_primary_cluster_id'] = None
            if row['mob_suite_secondary_cluster_id'] == '-':
                row['mob_suite_secondary_cluster_id'] = None
            if row['mob_suite_primary_cluster_id'] and row['mob_suite_secondary_cluster_id']:
                row['mob_suite_combined_primary_secondary_cluster_id'] = '-'.join([row['mob_suite_primary_cluster_id'], row['mob_suite_secondary_cluster_id']])
            else:
                row['mob_suite_combined_primary_secondary_cluster_id'] = None
            combined_abricate_mobtyper_report.append(row)

    return combined_abricate_mobtyper_report


def main(args):
    combined_abricate_mobtyper_report = parse_combined_abricate_mobtyper_report(args.combined_abricate_mobtyper_report)

    closest_ref_plasmid_by_plasmid_cluster_id = {}

    for record in combined_abricate_mobtyper_report:

        if args.cluster_level == "primary":
            plasmid_cluster_id = record['mob_suite_primary_cluster_id']
        elif args.cluster_level == "secondary":
            plasmid_cluster_id = record['mob_suite_combined_primary_secondary_cluster_id']

        if not plasmid_cluster_id:
            continue
        
        elif plasmid_cluster_id not in closest_ref_plasmid_by_plasmid_cluster_id:
            closest_ref_plasmid_by_plasmid_cluster_id[plasmid_cluster_id] = { 'reference_plasmid_id': record['mash_nearest_neighbor'],
                                                                              'least_mash_neighbor_distance': record['mash_neighbor_distance'] }
        elif record['mash_neighbor_distance'] < closest_ref_plasmid_by_plasmid_cluster_id[plasmid_cluster_id]['least_mash_neighbor_distance']:
            closest_ref_plasmid_by_plasmid_cluster_id[plasmid_cluster_id] = { 'reference_plasmid_id': record['mash_nearest_neighbor'],
                                                                              'least_mash_neighbor_distance': record['mash_neighbor_distance'] }

    
    fieldnames = [
        "plasmid_cluster_id",
        "reference_plasmid_id",
        "least_mash_distance",
    ]
    
    csv.register_dialect('unix-tab', delimiter='\t', doublequote=False, lineterminator='\n', quoting=csv.QUOTE_MINIMAL)
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, dialect='unix-tab')
    writer.writeheader()
    
    for plasmid_cluster_id, reference_plasmid in closest_ref_plasmid_by_plasmid_cluster_id.items():
        record = {'plasmid_cluster_id': plasmid_cluster_id,
                  'reference_plasmid_id': reference_plasmid['reference_plasmid_id'],
                  'least_mash_distance': reference_plasmid['least_mash_neighbor_distance']}
                  
        writer.writerow(record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-r', '--combined-abricate-mobtyper-report', help="")
    parser.add_argument('-c', '--cluster-level', default='primary', help="")
    args = parser.parse_args()
    main(args)
