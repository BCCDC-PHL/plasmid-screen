#!/usr/bin/env python

import argparse
import sys

def main(args):

    all_positions = 0
    positions_above_threshold = 0
    for line in sys.stdin:
        [chrom, pos, depth] = line.strip('\t').split()
        all_positions += 1
        if int(depth) >= args.min_depth:
            positions_above_threshold += 1

    print(','.join([
        'sample_id',
        'ref_plasmid',
        'ref_plasmid_length',
        'depth_coverage_threshold',
        'num_positions_above_depth_coverage_threshold',
        'percent_positions_above_depth_coverage_threshold',
    ]))

    print(','.join([
        args.sample_id,
        args.plasmid_id,
        str(all_positions),
        str(args.min_depth),
        str(positions_above_threshold),
        str(round(float(positions_above_threshold) / float(all_positions) * 100, 3)),
    ]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample-id')
    parser.add_argument('--plasmid-id')
    parser.add_argument('--min-depth', default=10)
    args = parser.parse_args()
    main(args)
