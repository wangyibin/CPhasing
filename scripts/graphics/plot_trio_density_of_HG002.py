#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot bar plot of the phase density 
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
import matplotlib.cm as cm
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

import pandas as pd
import numpy as np

from pathlib import Path
from cphasing.agp import import_agp
from cphasing.utilities import to_humanized, chrom_ticks_convert


def import_phase_result(tsv):
    data = []

    with open(tsv, 'r') as fp:
        for line in fp:
            if line.strip():
                if line.startswith("S"):
                    line_list = line.strip().split()
                    data.append(line_list)
        
    df = pd.DataFrame(data)
    
    df = df.drop([0], axis=1)
    df.columns = ["contig", "pat_kmer", "mat_kmer", "pat_specify", 
                  "pat_shared", "mat_shared", "mat_specify", "seq_len"]

    
    df = df.astype(dict(zip(["pat_kmer", "mat_kmer", "pat_specify", 
                  "pat_shared", "mat_shared", "mat_specify", "seq_len"], [np.int64] * 7)))
    df['total_kmer'] = df['pat_kmer'] + df["mat_kmer"]
    df['p_or_m'] = df[['pat_kmer', 'mat_kmer']].idxmax(axis=1).map(lambda x: 1 if x == 'mat_kmer' else -1)

    df['density'] = (df[['pat_kmer', 'mat_kmer']].max(axis=1) / df['total_kmer']) * df['p_or_m']

    df[df['pat_kmer'] == df['mat_kmer']]['density'] = 0
    
    df[np.abs(df['density']) < 0.999]['density'] = 0

    return df 


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('phase_tsv', 
            help='phase results from `yak trio-eval`')
    pReq.add_argument("agp")
    pOpt.add_argument('-o', '--output', type=str,
            default=None, help='output file [default: %(default)s]')
    
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    output = Path(args.agp).stem if not args.output else args.output

    df = import_phase_result(args.phase_tsv)
    df = df.set_index('contig')
    density_db = df.to_dict()['density']
    agp_df, _ = import_agp(args.agp)
    agp_df['phasing_density'] = agp_df['id'].map(density_db.get)
    agp_df = agp_df.reset_index()
   
    agp_df = agp_df[agp_df['chrom'] != agp_df['id']]
    agp_df['chrom'] = agp_df['chrom'].astype('object')
    
    
    plt.rcParams['font.family'] = 'Arial'
    fig = plt.figure(figsize=(8, 7))

    ax = fig.add_subplot(111)

    patches = []
    colors = []
    cmap = cm.get_cmap('bwr')
  
    norm = plt.Normalize(-1.0, 1.0)
    chroms = []
    for group, tmp_df in agp_df.groupby("chrom"):
        if "Chr" in group:
            group_idx = group.split("g")
       
            if len(group_idx) > 1:
                idx = int(group_idx[0][3:])
                hap = int(group_idx[1])
            else:
                hap = None
                idx = int(group_idx[0])
        elif "group" in group:
            idx = int(group[5:])
            hap = None
        
        elif "chr" in group:
            
            group_idx = group.split("_")
            try:
                idx = int(group_idx[0][3:])
            except ValueError:
                if group_idx[0][3:] == "X":
                    idx = 23
                elif group_idx[0][3:] == "Y":
                    idx = 24

            if group_idx[1] == "PATERNAL":
                hap = 1
            else:
                hap = 2
        if idx not in chroms:
            chroms.append(idx)

        for i, row in tmp_df.iterrows():
            patches.append(Rectangle((row.start, idx + hap/3 if hap is not None else idx), 
                                   row.end - row.start, 0.2, edgecolor='k', linewidth=0))
    
            colors.append(cmap(norm(row['phasing_density'])))

    patch_collection = ax.add_collection(PatchCollection(patches, edgecolor='k', linewidth=0))
    patch_collection.set_facecolor(colors)
    ax.autoscale()
    ax.invert_yaxis()
    _xticks = np.r_[ax.get_xticks()[1:-1]]
    ax.set_xticks(_xticks)
    xticklabels = chrom_ticks_convert(_xticks)
    ax.set_xticklabels(xticklabels, fontsize=10)

    ax.set_yticks(chroms)
    ax.set_yticklabels(chroms)
    



    sns.despine(trim=True, bottom=True, left=True)

    plt.savefig(f"{output}.barplot.png", dpi=600, bbox_inches='tight')
    plt.savefig(f"{output}.barplot.pdf", dpi=300, bbox_inches='tight')



if __name__ == "__main__":
    main(sys.argv[1:])