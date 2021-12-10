#!/usr/bin/env python

import argparse
import csv
import json
import re
import sys


def parse_abricate_report(abricate_report_path):
    abricate_report = []
    
    fieldnames = [
        "file",
        "sequence",
        "start",
        "end",
        "strand",
        "gene",
        "coverage",
        "coverage_map",
        "gaps",
        "percent_coverage",
        "percent_identity",
        "database",
        "accession",
        "product",
        "resistance",
    ]

    int_fields = [
        "start",
        "end",
    ]

    float_fields = [
        "percent_coverage",
        "percent_identity",
    ]

    with open(abricate_report_path, 'r') as f:
        reader = csv.DictReader(f, fieldnames=fieldnames, dialect="excel-tab")
        for row in reader:
            if row["file"].startswith("#"):
                continue
            else:
                for field in int_fields:
                    row[field] = int(row[field])
                for field in float_fields:
                    row[field] = float(row[field])
            abricate_report.append(row)

    return abricate_report

def parse_mob_typer_report(mob_typer_report_path):
    mob_typer_report = []
    
    fieldnames = [
        "sample_id",
        "num_contigs",
	"size",
	"gc",
	"md5",
	"rep_types",
	"rep_type_accessions",
	"relaxase_types",
	"relaxase_type_accessions",
	"mpf_type",
	"mpf_type_accessions",
	"orit_types",
	"orit_accessions",
	"predicted_mobility",
	"mash_nearest_neighbor",
	"mash_neighbor_distance",
	"mash_neighbor_identification",
	"primary_cluster_id",
	"secondary_cluster_id",
	"predicted_host_range_overall_rank",
	"predicted_host_range_overall_name",
	"observed_host_range_ncbi_rank",
	"observed_host_range_ncbi_name",
	"reported_host_range_lit_rank",
	"reported_host_range_lit_name",
	"associated_pmids",
    ]

    int_fields = [
        "num_contigs",
    ]

    float_fields = [
        "mash_neighbor_distance",
        "gc",
    ]

    with open(mob_typer_report_path, 'r') as f:
        reader = csv.DictReader(f, fieldnames=fieldnames, dialect="excel-tab")
        for row in reader:
            if row["sample_id"].startswith("sample_id"):
                continue
            else:
                for field in int_fields:
                    row[field] = int(row[field])
                for field in float_fields:
                    row[field] = float(row[field])
            mob_typer_report.append(row)

    return mob_typer_report


def main(args):
    abricate_report = parse_abricate_report(args.abricate_report)
    # print("[")
    # print(json.dumps(abricate_report))
    # print(",")
    mob_typer_report = parse_mob_typer_report(args.mob_typer_report)
    # print(json.dumps(mob_typer_report))
    # print("]")

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
        "depth_coverage_threshold",
        "percent_ref_plasmid_coverage_above_depth_threshold",
        "num_snps_vs_ref_plasmid",
    ]
    csv.register_dialect('unix-tab', delimiter='\t', doublequote=False, lineterminator='\n', quoting=csv.QUOTE_MINIMAL)
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='unix-tab')
    writer.writeheader()


    for abricate_record in abricate_report:
        if re.search("chromosome", abricate_record['file']) and abricate_record['resistance'] == "CARBAPENEM":
            output_record = {
                "sample_id": args.sample_id,  
                "assembly_file": abricate_record["file"],
                "resistance_gene_contig": abricate_record['sequence'],
                "num_contigs_in_reconstruction": "-",
                "reconstruction_size": "-",
                "resistance_gene": abricate_record["gene"],
                "gene_start": abricate_record["start"],
                "gene_end": abricate_record["end"],
                "percent_resistance_gene_coverage": abricate_record["percent_coverage"],
                "percent_resistance_gene_identity": abricate_record["percent_identity"],
                "replicon_types": "-",
                "mob_suite_primary_cluster_id": "-",
                "mob_suite_secondary_cluster_id": "-",
                "mash_nearest_neighbor": "-",
                "mash_neighbor_distance": "_",
                "alignment_ref_plasmid": "-",
                "depth_coverage_threshold": "-",
                "percent_ref_plasmid_coverage_above_depth_threshold": "-",
                "num_snps_vs_ref_plasmid": "-",
            }
            writer.writerow(output_record)
        else:
            for mob_typer_record in mob_typer_report:
                if re.search(mob_typer_record['primary_cluster_id'], abricate_record['file']) and abricate_record['resistance'] == "CARBAPENEM":
                    output_record = {
                        "sample_id": args.sample_id,  
                        "assembly_file": abricate_record["file"],
                        "resistance_gene_contig": abricate_record['sequence'],
                        "num_contigs_in_reconstruction": mob_typer_record['num_contigs'],
                        "reconstruction_size": mob_typer_record['size'],
                        "resistance_gene": abricate_record["gene"],
                        "gene_start": abricate_record["start"],
                        "gene_end": abricate_record["end"],
                        "percent_resistance_gene_coverage": abricate_record["percent_coverage"],
                        "percent_resistance_gene_identity": abricate_record["percent_identity"],
                        "replicon_types": mob_typer_record["rep_types"],
                        "mob_suite_primary_cluster_id": mob_typer_record["primary_cluster_id"],
                        "mob_suite_secondary_cluster_id": mob_typer_record["secondary_cluster_id"],
                        "mash_nearest_neighbor": mob_typer_record["mash_nearest_neighbor"],
                        "mash_neighbor_distance": mob_typer_record["mash_neighbor_distance"],
                        "alignment_ref_plasmid": "-",
                        "depth_coverage_threshold": "-",
                        "percent_ref_plasmid_coverage_above_depth_threshold": "-",
                        "num_snps_vs_ref_plasmid": "-",
                    }
                    writer.writerow(output_record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--sample-id', help="")
    parser.add_argument('--abricate-report', help="")
    parser.add_argument('--mob-typer-report', help="")
    args = parser.parse_args()
    main(args)
