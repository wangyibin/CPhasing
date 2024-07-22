#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-05-29 17:20

import argparse
import pysam
from portion import closed

pysam.set_verbosity(0)

primary_flags = {0, 16}

def parse_bam(bam, threads):
    
    unmapped = 0
    
    diff_ref_homo = 0
    diff_ref_homo_mapq_gt0 = 0
    diff_ref_homo_mapq_gt1 = 0
    
    diff_ref_nonhomo = 0
    diff_ref_nonhomo_mapq_gt0 = 0
    diff_ref_nonhomo_mapq_gt1 = 0
    
    same_ref_but_diff_loci = 0
    same_ref_but_diff_loci_mapq_gt0 = 0
    same_ref_but_diff_loci_mapq_gt1 = 0
    
    same_locus = 0
    same_locus_mapq_gt0 = 0
    same_locus_mapq_gt1 = 0
    

    with pysam.AlignmentFile(bam, mode='rb', threads=threads) as f:
        for aln in f:
            
            query_name = aln.query_name
            
            if not aln.is_mapped:
                unmapped += 1
                continue
            if aln.is_secondary or aln.is_supplementary:
                continue
            assert aln.flag in primary_flags
            
            fields = query_name.rsplit('_', 6)
            ref_truth, ref_start_truth, ref_end_truth, orient = fields[3], int(fields[4]), int(fields[5]), fields[6]
            ref_name = aln.reference_name
            
            if ref_truth != ref_name:
                if ref_truth[:-1] == ref_name[:-1]:
                    diff_ref_homo += 1
                    if aln.mapping_quality > 1:
                        diff_ref_homo_mapq_gt1 += 1
                    if aln.mapping_quality > 0:
                        diff_ref_homo_mapq_gt0 += 1
                else:
                    diff_ref_nonhomo += 1
                    if aln.mapping_quality > 1:
                        diff_ref_nonhomo_mapq_gt1 += 1
                    if aln.mapping_quality > 0:
                        diff_ref_nonhomo_mapq_gt0 += 1
                continue
            
            ref_start, ref_end = aln.reference_start, aln.reference_end
            overlap = closed(ref_start_truth, ref_end_truth) & closed(ref_start + 1, ref_end)
            
            if overlap:
                ovl_len = overlap.upper - overlap.lower + 1
                min_len = min(ref_end_truth - ref_start_truth + 1, ref_end - ref_start)
                ovl_ratio = ovl_len / min_len
            
            if not overlap or ovl_ratio < 0.8 or (orient == '+' and aln.is_reverse) or (orient == '-' and aln.is_forward):
                same_ref_but_diff_loci += 1
                if aln.mapping_quality > 1:
                    same_ref_but_diff_loci_mapq_gt1 += 1
                if aln.mapping_quality > 0:
                    same_ref_but_diff_loci_mapq_gt0 += 1
            else:
                same_locus += 1
                if aln.mapping_quality > 1:
                    same_locus_mapq_gt1 += 1
                if aln.mapping_quality > 0:
                    same_locus_mapq_gt0 += 1
    

    total_mapped = diff_ref_homo + diff_ref_nonhomo + same_ref_but_diff_loci + same_locus
    total = unmapped + total_mapped 
    
    print('------------ All -------------')

    print('unmapped: {} ({}%)\ndiff_ref_homo: {} ({}%)\ndiff_ref_nonhomo: {} ({}%)\nsame_ref_but_diff_loci: {} ({}%)\nsame_locus: {} ({}%)'.format(
        unmapped, unmapped/total*100, diff_ref_homo, diff_ref_homo/total_mapped*100, diff_ref_nonhomo, diff_ref_nonhomo/total_mapped*100, 
        same_ref_but_diff_loci, same_ref_but_diff_loci/total_mapped*100, same_locus, same_locus/total_mapped*100))
    
    total_mapq_gt0 = diff_ref_homo_mapq_gt0 + diff_ref_nonhomo_mapq_gt0 + same_ref_but_diff_loci_mapq_gt0 + same_locus_mapq_gt0

    print('--------- MAPQ >= 1 -----------')

    print('mapq >= 1: {} ({}%)\ndiff_ref_homo: {} ({}%)\ndiff_ref_nonhomo: {} ({}%)\nsame_ref_but_diff_loci: {} ({}%)\nsame_locus: {} ({}%)'.format(
        total_mapq_gt0, total_mapq_gt0/total_mapped*100, diff_ref_homo_mapq_gt0, diff_ref_homo_mapq_gt0/total_mapq_gt0*100, diff_ref_nonhomo_mapq_gt0,
        diff_ref_nonhomo_mapq_gt0/total_mapq_gt0*100, same_ref_but_diff_loci_mapq_gt0, same_ref_but_diff_loci_mapq_gt0/total_mapq_gt0*100,
        same_locus_mapq_gt0, same_locus_mapq_gt0/total_mapq_gt0*100))

    total_mapq_gt1 = diff_ref_homo_mapq_gt1 + diff_ref_nonhomo_mapq_gt1 + same_ref_but_diff_loci_mapq_gt1 + same_locus_mapq_gt1

    print('--------- MAPQ >= 2 -----------')

    print('mapq >= 2: {} ({}%)\ndiff_ref_homo: {} ({}%)\ndiff_ref_nonhomo: {} ({}%)\nsame_ref_but_diff_loci: {} ({}%)\nsame_locus: {} ({}%)'.format(
        total_mapq_gt1, total_mapq_gt1/total_mapped*100, diff_ref_homo_mapq_gt1, diff_ref_homo_mapq_gt1/total_mapq_gt1*100, diff_ref_nonhomo_mapq_gt1,
        diff_ref_nonhomo_mapq_gt1/total_mapq_gt1*100, same_ref_but_diff_loci_mapq_gt1, same_ref_but_diff_loci_mapq_gt1/total_mapq_gt1*100,
        same_locus_mapq_gt1, same_locus_mapq_gt1/total_mapq_gt1*100))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help='subread alignments in BAM format')
    parser.add_argument('--threads', type=int, default=8, help='threads for bam reading, default: %(default)s')
    args = parser.parse_args()

    parse_bam(args.bam, args.threads)


if __name__ == '__main__':
    main()

