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
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

import pandas as pd
import numpy as np

from pathlib import Path
from cphasing.agp import import_agp
from cphasing.utilities import to_humanized, chrom_ticks_convert
from cphasing.plot import chrRangeID


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
    df['density'] = df['density'].fillna(0)
    
    # df[np.abs(df['density']) < 0.999]['density'] = 0

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
    pReq.add_argument('contig')
    pOpt.add_argument('-b', '--break_pos', 
                      help="break pos information, e.g. C-Phasing:2000,HapHiC:2000", default=None)
    pOpt.add_argument('-o', '--output', type=str,
            default=None, help='output file [default: %(default)s]')
    
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    output = args.contig if not args.output else args.output 


    df = import_phase_result(args.phase_tsv)
    df2 = pd.DataFrame(df['contig'].map(lambda x: chrRangeID(x, axis=1)).values.tolist())
    df2.columns = ["contig_source", "start", "end"]
    df2 = df2.astype({'start': int, 'end': int})
    df = pd.concat([df, df2], axis=1)
    plt.rcParams['font.family'] = 'Arial'
    fig = plt.figure(figsize=(8, 1))

    ax = fig.add_subplot(111)

    patches = []
    colors = []
    custom_cmap = mcolors.LinearSegmentedColormap.from_list(
        'custom_cmap', ['#0000FF', '#bcbcbc', '#FF0000'], N=256)

    cmap = cm.get_cmap('bwr')
    
    norm = plt.Normalize(-1.0, 1.0)
    chroms = [args.contig]
    chromidx = []
    seq_len = 0
    contig_df = df.loc[df['contig_source'].isin([args.contig])]
    chromidx.append(0.2)
    for idx, row in contig_df.iterrows():
        
        patches.append(Rectangle((row.start, 0), row.end - row.start + 1, 0.2, 
                                 edgecolor='none', linewidth=0.0))
        colors.append(cmap(norm(row['density'])))
        
    
    seq_len = (contig_df['end'] - contig_df['start'] + 1).sum()

 
    patch_collection = ax.add_collection(PatchCollection(patches, edgecolor='k', linewidth=0.15))
    patch_collection.set_facecolor(colors)
    ax.autoscale()
    ax.invert_yaxis()

    _xticks = np.r_[ax.get_xticks()[1:-2], seq_len]
    
    ax.set_xticks(_xticks)
    xticklabels = chrom_ticks_convert(_xticks)
    ax.set_xticklabels(xticklabels, fontsize=10)

    color_db = ["#FBB4AE" ,"#B3CDE3", "#CCEBC5", "#DECBE4"]
    marker_db = ["^", "v", "2", "1"]
    ylim = ax.get_ylim()
    sns.despine(trim=True, bottom=True, left=True)
    if args.break_pos:
        l = args.break_pos.split(",")
        if len(l) > 0:
            l2 = list(map(lambda x: x.split(":"), l))
            l2 = list(filter(lambda x: len(x) == 2, l2))
            if len(l2) > 0:
                for i, (software, pos) in enumerate(l2):
                    
                    pos = int(pos) 
                    y_pos =  -0.5 if i % 2 == 0 else 2
        
                    marker = marker_db[i%3]
                    ax.plot(pos, y_pos, marker=marker, markersize=7.5, 
                    linewidth=1, alpha=1,
                    linestyle='none', color=color_db[i],
                    fillstyle='none', label=software)
                   

        plt.legend(fontsize=8, bbox_to_anchor=(1.1, 1.05))

    
    ax.set_yticks(np.array(chromidx)/2)
    ax.set_yticklabels(chroms)

    plt.ylim(-1, 3)
    
    cax = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), 
                    #    orientation='horizontal',
                       ax=ax, pad=0.1)
    cax.set_ticklabels(['Paternal', '', 'Maternal'], fontsize=8)
    

    plt.savefig(f"{output}.barplot.png", dpi=600, bbox_inches='tight')
    plt.savefig(f"{output}.barplot.pdf", dpi=300, bbox_inches='tight')



if __name__ == "__main__":
    main(sys.argv[1:])