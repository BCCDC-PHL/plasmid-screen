#!/usr/bin/env python

import argparse
import csv
import json
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


def main(args):
    abricate_report = parse_abricate_report(args.abricate_report)
    # print(json.dumps(abricate_report))

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
    
    csv.register_dialect('unix-tab', delimiter='\t', doublequote=False, lineterminator='\n', quoting=csv.QUOTE_MINIMAL)
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, dialect='unix-tab')
    writer.writeheader()
    for record in abricate_report:
        if record["resistance"] == "CARBAPENEM":
            writer.writerow(record)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--abricate-report', help="")
    parser.add_argument('--sample-id', help="")

    args = parser.parse_args()
    main(args)
