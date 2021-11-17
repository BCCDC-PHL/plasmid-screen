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
        'ref_plasmid_length',
        'positions_above_coverage_threshold',
        'percent_coverage_above_threshold',
    ])
    print(','.join([
        str(all_positions),
        str(positions_above_threshold),
        str(round(float(positions_above_threshold) / float(all_positions) * 100, 3)),
    ]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--min-depth', default=10)
    args = parser.parse_args()
    main(args)
