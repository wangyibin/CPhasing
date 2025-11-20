#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

"""

import argparse
import logging
import os
import os.path as op
import sys
import random
import pandas as pd
import numpy as np

from collections import defaultdict
from intervaltree import IntervalTree, Interval


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('contig_bed')
    pReq.add_argument('collapsed_bed', 
            help='')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    # contig_bed = pd.read_csv(args.contig_bed, header=None, sep='\t')
    # contig_bed.columns = ['chrom', 'start', 'end', 'contig']
    collapsed_bed = pd.read_csv(args.collapsed_bed, header=None, sep='\t')
    collapsed_bed.columns = ['chrom', 'start', 'end']

    collapsed_bed['hap'] = collapsed_bed['chrom'].str.slice(0, 1)
    res = {}
    for hap, tmp_df in collapsed_bed.groupby('hap'):
        res[hap] = IntervalTree()
        for row in tmp_df.iterrows():
            res[hap].addi(row[1]['start'], row[1]['end'], row[1]['chrom'])
    
    remove_regions = []
    collapsed_regions = []
    for hap in res:
        if len(res[hap]) <= 1:
            continue 
        points = []
        for interval in res[hap]:
            points.append((interval.begin, 1)) 
            points.append((interval.end, -1)) 
  
        points.sort()
        coverage = 0
        last_pos = points[0][0]
        segments = []
        

        for pos, type in points:
            if pos > last_pos:
                if coverage > 1:
                    segments.append({'start': last_pos, 'end': pos, 'cn': coverage})
            
            coverage += type
            last_pos = pos

        if not segments:
            continue

        merged_regions = []
        current_region = segments[0]
        for i in range(1, len(segments)):
            next_region = segments[i]
            if next_region['start'] == current_region['end'] and next_region['cn'] == current_region['cn']:
        
                current_region['end'] = next_region['end']
            else:

                merged_regions.append(current_region)
                current_region = next_region
        merged_regions.append(current_region) 


        res_regions_by_cn = defaultdict(list)
        for region in merged_regions:
            res_regions_by_cn[region['cn']].append((region['start'], region['end'], region['cn']))
        

        final_regions = list(res_regions_by_cn.values())

        res_regions = []
        for idx, regions in enumerate(final_regions):
            cn = regions[0][2] 
            collapsed_chrom_idx = random.sample(list(range(0, cn)), cn-1)
            chroms = [interval.data for interval in res[hap].items()]

            regions = list(regions[0][:2])
            for idx2, chrom in enumerate(chroms):
                if idx2 in collapsed_chrom_idx:
                    remove_regions.append((chrom, regions[0], regions[1]))
                else:
                    collapsed_regions.append((chrom, regions[0], regions[1]))

            
    # contig_bed['hap'] = contig_bed['chrom'].str.slice(0, 1)
    # res2 = {}
    # for idx, row in contig_bed.iterrows():
    #     hap = row['hap']
    #     if hap not in res:
    #         continue
    
    #     if res[hap].overlaps(Interval(row['start'], row['end'])):
    #         # Find the overlapping intervals
    #         overlapping_intervals = res[hap][row['start']:row['end']]
    #         if not overlapping_intervals:
    #             continue

    #         ## get the max overlap with row 
    #         if len(overlapping_intervals) == 1:
    #             max_overlap = list(overlapping_intervals)[0]
    #         else:
    #             ## get the max overlap by maximum overlap with row
    #             max_overlap = max(overlapping_intervals, key=lambda x: min(x.end, row['end']) - max(x.begin, row['start']))

            
    #         # If the contig is completely covered by the collapsed intervals, we can skip it
    #         if max_overlap.begin <= row['start'] and max_overlap.end >= row['end']:
                
    #             if (max_overlap.begin, max_overlap.end) not in res2:
    #                 res2[(max_overlap.begin, max_overlap.end)] = []
    #             res2[(max_overlap.begin, max_overlap.end)].append(row['contig'])

    remove_regions_df = pd.DataFrame(remove_regions, columns=['chrom', 'start', 'end'])
    collapsed_regions_df = pd.DataFrame(collapsed_regions, columns=['chrom', 'start', 'end'])
    remove_regions_df.to_csv(f"{args.collapsed_bed}.remove_regions.tsv", sep='\t', index=False, header=False)
    collapsed_regions_df.to_csv(f"{args.collapsed_bed}.collapsed_regions.tsv", sep='\t', index=False, header=False)

    cmd = f"bedtools intersect -a {args.contig_bed} -b {args.collapsed_bed}.remove_regions.tsv -v > {args.collapsed_bed}.collapsed_contigs.bed" 

    os.system(cmd)

if __name__ == "__main__":
    main(sys.argv[1:])