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
from matplotlib.ticker import MaxNLocator


import pandas as pd
import numpy as np

from collections import defaultdict, OrderedDict
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
    pReq.add_argument('break_pos', 
                      help="break pos information file, three columns e.g. utg001 1000 Hitig", default=None)
    pOpt.add_argument('-o', '--output', type=str,
            default="output", help='output file [default: %(default)s]')
    
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    output = args.output 

    df = import_phase_result(args.phase_tsv)
    df2 = pd.DataFrame(df['contig'].map(lambda x: chrRangeID(x, axis=1)).values.tolist())
    df2.columns = ["contig_source", "start", "end"]
    df2 = df2.astype({'start': int, 'end': int})
    df = pd.concat([df, df2], axis=1)
    contigsizes = df2.groupby('contig_source')['end'].max().to_dict()


    break_pos_df = pd.read_csv(args.break_pos, sep='\t', header=None, index_col=None)
    break_pos_df.columns = ['contig', 'pos', 'software']
    

    break_contigs = set(break_pos_df['contig'].values.tolist())
    break_contigs_with_length = dict(zip(break_contigs, map(contigsizes.get, break_contigs)))
    break_contigs = sorted(break_contigs_with_length, key=lambda x: break_contigs_with_length[x], reverse=True)
    

    plt.rcParams['font.family'] = 'Arial'
    fig = plt.figure(figsize=(8, len(break_contigs) / 10))

    ax = fig.add_subplot(111)

    patches = []
    colors = []

    cmap = cm.get_cmap('bwr')
    
    norm = plt.Normalize(-1.0, 1.0)
    chroms = []
    chromidx = []
    seq_len = 0
    contig_y_idx = defaultdict(int)
    idx = 0
    for contig in break_contigs:
        if contig in contig_y_idx:
            y_idx = contig_y_idx[contig]
        else:
            contig_y_idx[contig] = idx
            y_idx = idx
            idx += 1

        contig_y_idx[contig] += idx
        contig_df = df.loc[df['contig_source'].isin([contig])]
 
        chromidx.append(idx)

        for _, row in contig_df.iterrows():
            
            patches.append(Rectangle((row.start, y_idx), row.end - row.start + 1, 0.4, 
                                    edgecolor='none', linewidth=0.0))
            colors.append(cmap(norm(row['density'])))
            
        
        _seq_len = (contig_df['end'] - contig_df['start'] + 1).sum()
        if _seq_len > seq_len:
            seq_len = _seq_len

    
    patch_collection = ax.add_collection(PatchCollection(patches, edgecolor='#bcbcbc', linewidth=0.01))
    patch_collection.set_facecolor(colors)
    ax.autoscale()
    ax.invert_yaxis()

    _xticks = np.r_[ax.get_xticks()[1:-2], seq_len]
    
    ax.set_xticks(_xticks)
    xticklabels = chrom_ticks_convert(_xticks)
    ax.set_xticklabels(xticklabels, fontsize=10)

    color_db = ["#FBB4AE" ,"#B3CDE3", "#CCEBC5", "#DECBE4"]
    softwares = ["Hitig"]

    marker_db = ["^", "v", "2", "1"]
    ylim = ax.get_ylim()
    sns.despine(trim=True, bottom=True, left=True)

    
    for _, row in break_pos_df.iterrows():
        try:
            l = row.pos.split(",")
        except AttributeError:
            l = [row.pos]

        software = row.software
        software_idx = software.index(software)
       
        y_pos =  break_contigs.index(row.contig) + 0.7
        for i, pos in enumerate(l):
            
            pos = int(pos) 
            
            marker = marker_db[software_idx]
            ax.plot(pos, y_pos, marker=marker, markersize=1, 
            linewidth=0.01, alpha=1,
            linestyle='none', color=color_db[software_idx],
            fillstyle='none', label=software)
                

    # plt.legend(fontsize=8, bbox_to_anchor=(1.1, 1.05))

    ax.set_yticks(list(range(len(break_contigs))))
    ax.set_yticklabels(break_contigs, fontsize=4)

    # plt.ylim(-1)
    
    cax = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), 
                    #    orientation='horizontal',
                       ax=ax, pad=0.02, shrink=0.4)
    cax.set_ticks([-1, 0, 1])
    cax.set_ticklabels(['Paternal', '', 'Maternal'], fontsize=8)
    
   

    plt.savefig(f"{output}.barplot.png", dpi=1200, bbox_inches='tight')
    plt.savefig(f"{output}.barplot.pdf", dpi=300, bbox_inches='tight')



if __name__ == "__main__":
    main(sys.argv[1:])