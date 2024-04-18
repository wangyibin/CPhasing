#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot the example of hitig correct alignments
"""


import argparse
import logging
import os
import os.path as op
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.collections import LineCollection, PathCollection
from matplotlib.patches import Rectangle

import numpy as np
import pandas as pd 

from collections import defaultdict

def read_paf(paf, minq=0):
    df = pd.read_csv(paf, sep='\t', index_col=None, header=None, usecols=range(13))
    df = df[df[11] >= minq]
    return df 

def parse_region(region):
    chrom, region = region.split(":")
    start, end = region.split("-")
    
    return chrom, int(start), int(end)

def is_overlap(range1, range2):
    start1, end1 = range1
    start2, end2 = range2

    return start1 <= end2 and end1 >= start2

def is_correct(row):
    if row[5] != row.real_chrom:
        return False 

    if not is_overlap((row[7], row[8]), (row.real_start, row.real_end)):
        return False 

    return True
    

    
def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('raw_paf', 
            help='')

    pOpt.add_argument('-r1', help="Plot region", required=True)
    pOpt.add_argument('-r2', help="Plot region", required=True)
    pOpt.add_argument('-q', '--minq', help='minimum mapping quality', 
                      default=0, type=int)
    pOpt.add_argument('-o', '--output', type=str,
            default='output.hitig.alignments.png', help='output file [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    raw_df = read_paf(args.raw_paf, minq=args.minq)
    read_id_df = raw_df[0].str.split("!", expand=True)
    # strand_df = read_id_df[4].str.split("_", expand=True)
    # strand_df.columns = ["strand", "idx"]
    read_id_df.columns = ["read_id", "real_chrom", "real_start", "real_end", 4]
    read_id_df['real_start'] = read_id_df['real_start'].astype(int)
    read_id_df['real_end'] = read_id_df['real_end'].astype(int)

    read_id_df.drop(4, axis=1, inplace=True)
    raw_df = pd.concat([raw_df, read_id_df], axis=1)

    raw_df['is_correct'] = raw_df.apply(is_correct, axis=1)
    

    region1 = parse_region(args.r1)
    region2 = parse_region(args.r2)

    region1_df = raw_df[(raw_df['real_chrom'] == region1[0])]
    region1_df = region1_df[(region1_df['real_end'] >= region1[1])]
    region1_df = region1_df[(region1_df['real_start'] <= region1[2])]
   
    region2_df = raw_df[(raw_df['real_chrom'] == region2[0])]
    region2_df = region2_df[(region2_df['real_end'] >= region2[1])]
    region2_df = region2_df[(region2_df['real_start'] <= region2[2])]

    color1 = "#b02418"
    color2 = "#253761"
    color1_1 = "#cb6e7f"
    color2_1 = "#8896ae"
    
    lines1 = []

    chrom1_start_idx = 100
    
    idx_db = defaultdict(int)
    idx = 2
    colors1 = []
    for i, row in region1_df.sort_values('real_start').iterrows():
        read = row['read_id']
        if read not in idx_db:
            idx += 1
            idx_db[read] = idx
        
        if row.is_correct is True:
           
           colors1.append(color1_1)
        else:
           colors1.append(color2_1)
           
        line = [(row[7] - region1[1], chrom1_start_idx - idx_db[read]), 
                (row[8] - region1[1], chrom1_start_idx - idx_db[read])]
        lines1.append(line)


    chrom2_start_idx = chrom1_start_idx  - max(idx_db.values()) - 10

    lines2 = []
   
    colors2 = []
    idx_db = defaultdict(lambda :5)
    idx = 2
    for i, row in region2_df.sort_values('real_start').iterrows():
        read = row['read_id']
        if read not in idx_db:
            idx += 1
            idx_db[read] = idx
        
        if row.is_correct is True:
            
            colors2.append(color2_1)
        else:
            
            colors2.append(color1_1)

        line = [(row[7] - region2[1], chrom2_start_idx - idx_db[read]), 
                (row[8] - region2[1], chrom2_start_idx - idx_db[read])]
        lines2.append(line)


    plt.rcParams['font.family'] = 'Arial'
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    
    line1 = [(0, chrom1_start_idx), (region1[2] - region1[1], chrom1_start_idx)]
    line2 = [(0, chrom2_start_idx), (region2[2] - region2[1], chrom2_start_idx)]
    line_segments = LineCollection([line1, line2], linewidth=3, color=[color1, color2])
    ax.add_collection(line_segments)

    line_segments = LineCollection(lines1, linewidth=1, color=colors1)
    ax.add_collection(line_segments)

    line_segments = LineCollection(lines2, linewidth=1, color=colors2)
    ax.add_collection(line_segments)

    

    ax.autoscale()
    plt.yticks([])
    
    plt.text(-1, chrom1_start_idx + 5, f"{region1[0]}")
    plt.text(-1, chrom2_start_idx + 5, f"{region2[0]}")
    plt.xlim(0, max(region1[2] - region1[1], region2[2] - region2[1]))
    plt.xticks(np.linspace(0, max(region1[2] - region1[1], region2[2] - region2[1]), 5),
                list(map(lambda x: f"{x:,}", np.linspace(min(region1[1], region2[1]), max(region1[2], region2[2]), 5).astype(int))))

    sns.despine(trim=True, left=True)
    plt.savefig(args.output, dpi=600, bbox_inches='tight')








    
    
    







if __name__ == "__main__":
    main(sys.argv[1:])
