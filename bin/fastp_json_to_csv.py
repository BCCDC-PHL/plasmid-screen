#!/usr/bin/env python

import argparse
import csv
import json
import sys

def main(args):
    with open(args.fastp_json, 'r') as f:
        fastp_report = json.load(f)

    total_reads_before_filtering = fastp_report['summary']['before_filtering']['total_reads']
    total_reads_after_filtering = fastp_report['summary']['after_filtering']['total_reads']
    total_bases_before_filtering = fastp_report['summary']['before_filtering']['total_bases']
    total_bases_after_filtering = fastp_report['summary']['after_filtering']['total_bases']
    read1_mean_length_before_filtering = fastp_report['summary']['before_filtering']['read1_mean_length']
    read2_mean_length_before_filtering = fastp_report['summary']['before_filtering']['read2_mean_length']
    read1_mean_length_after_filtering = fastp_report['summary']['after_filtering']['read1_mean_length']
    read2_mean_length_after_filtering = fastp_report['summary']['after_filtering']['read2_mean_length']
    q20_bases_before_filtering = fastp_report['summary']['before_filtering']['q20_bases']
    q20_bases_after_filtering = fastp_report['summary']['after_filtering']['q20_bases']
    q20_rate_before_filtering = fastp_report['summary']['before_filtering']['q20_rate']
    q20_rate_after_filtering = fastp_report['summary']['after_filtering']['q20_rate']
    q30_bases_before_filtering = fastp_report['summary']['before_filtering']['q30_bases']
    q30_bases_after_filtering = fastp_report['summary']['after_filtering']['q30_bases']
    q30_rate_before_filtering = fastp_report['summary']['before_filtering']['q30_rate']
    q30_rate_after_filtering = fastp_report['summary']['after_filtering']['q30_rate']
    gc_content_before_filtering = fastp_report['summary']['before_filtering']['gc_content']
    gc_content_after_filtering = fastp_report['summary']['after_filtering']['gc_content']
    adapter_trimmed_reads = fastp_report['adapter_cutting']['adapter_trimmed_reads']
    adapter_trimmed_bases = fastp_report['adapter_cutting']['adapter_trimmed_bases']


    output_fields = [
        'total_reads_before_filtering',
        'total_reads_after_filtering',
        'total_bases_before_filtering',
        'total_bases_after_filtering',
        'read1_mean_length_before_filtering',
        'read1_mean_length_after_filtering',
        'read2_mean_length_before_filtering',
        'read2_mean_length_after_filtering',
        'q20_bases_before_filtering',
        'q20_bases_after_filtering',
        'q20_rate_before_filtering',
        'q20_rate_after_filtering',
        'q30_bases_before_filtering',
        'q30_bases_after_filtering',
        'q30_rate_before_filtering',
        'q30_rate_after_filtering',
        'gc_content_before_filtering',
        'gc_content_after_filtering',
        'adapter_trimmed_reads',
        'adapter_trimmed_bases',
    ]

    output_data = []
    if args.sample_id:
        output_fields = ['sample_id'] + output_fields
        output_data = [args.sample_id]

    output_data = output_data + [
        total_reads_before_filtering,
        total_reads_after_filtering,
        total_bases_before_filtering,
        total_bases_after_filtering,
        read1_mean_length_before_filtering,
        read1_mean_length_after_filtering,
        read2_mean_length_before_filtering,
        read2_mean_length_after_filtering,
        q20_bases_before_filtering,
        q20_bases_after_filtering,
        q20_rate_before_filtering,
        q20_rate_after_filtering,
        q30_bases_before_filtering,
        q30_bases_after_filtering,
        q30_rate_before_filtering,
        q30_rate_after_filtering,
        gc_content_before_filtering,
        gc_content_after_filtering,
        adapter_trimmed_reads,
        adapter_trimmed_bases,
    ]
    output_row = dict(zip(output_fields, output_data))
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fields, dialect='unix', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
    writer.writeheader()
    writer.writerow(output_row)
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fastp_json')
    parser.add_argument('-s', '--sample-id')
    args = parser.parse_args()
    main(args)
