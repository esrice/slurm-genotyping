#!/usr/bin/env python3
"""
This is a modified version of the coverage_to_regions.py
script that comes with freebayes, Copyright (c) 2010
Erik Garrison, Gabor Marth.
"""

import sys
import argparse
import pysam

def parse_args():
    parser = argparse.ArgumentParser(description='Generates a list of regions '
            'with even sequencing coverage based on a bam file.')
    parser.add_argument('-n', '--num-regions', type=int, default=10000,
        help='approximate number of regions to divide genome into [10000]')
    parser.add_argument('bamfile', type=lambda f: pysam.AlignmentFile(f, 'rb'),
        help='a bam file to use to calculate coverage')
    return parser.parse_args()

def get_total_coverage_estimate(bamfile, nhead=1000000):
    """
    Estimate the total coverage across all bams and
    positions using the number of mapped reads and
    mean read length.
    """
    num_mapped_reads = bamfile.mapped
    mean_read_length = sum(map(lambda r: r.qlen, bamfile.head(nhead)))/nhead
    return num_mapped_reads * mean_read_length

def main():
    args = parse_args()

    # get a dictionary of reference sequence lengths from the first bam header
    ref_lengths = dict(zip(args.bamfile.references, args.bamfile.lengths))

    total_coverage = get_total_coverage_estimate(args.bamfile)
    bp_per_region = total_coverage / args.num_regions

    lchrom = None
    lpos = 0
    bp_in_region = 0
    for pile in args.bamfile.pileup():
        chrom, pos = pile.reference_name, pile.reference_pos
        depth = pile.nsegments

        if lchrom != chrom:
            if lchrom:
                print('\t'.join(map(str, [lchrom, lpos, ref_lengths[lchrom]])))
                lpos = 0
                bp_in_region = 0
            lchrom = chrom
        bp_in_region += depth
        if bp_in_region > bp_per_region:
            print('\t'.join(map(str, [chrom, lpos, pos])))
            lpos = pos
            bp_in_region = 0

    print('\t'.join(map(str, [lchrom, lpos, ref_lengths[lchrom]])))

if __name__ == '__main__':
    main()
