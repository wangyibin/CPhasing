#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot line plot 
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
import patchworklib as pw

import numpy as np 
import pandas as pd 
from copy import deepcopy
from collections import OrderedDict
from functools import reduce
from string import ascii_lowercase

def plot(data, title, x, y, hue, output, merge):
    # colors = {"C-Phasing": "#b02418",
    #           "HapHiC": "#cb6e7f",
    #           "ALLHiC": "#253761",
    #           "ALLHiC_pregroup": "#8896ae"}
    # maker = {"C-Phasing": "s",
    #           "HapHiC": "^",
    #           "ALLHiC": "o",
    #           "ALLHiC_pregroup": "D"}
    colors = {"C-Phasing (Pore-C)": "#b02418",
              "C-Phasing (Hi-C)": "#cb6e7f",
              "HapHiC": "#253761",
              "ALLHiC": "#8896ae"}
    # colors = {"C-Phasing (Pore-C)": "#FBB4AE",
    #           "C-Phasing (Hi-C)": "#B3CDE3",
    #           "HapHiC":  "#CCEBC5",
    #           "ALLHiC": "#DECBE4"}
    maker = {"C-Phasing (Pore-C)": "s",
              "C-Phasing (Hi-C)": "^",
              "HapHiC": "o",
              "ALLHiC": "D"}
    markers = ["D", "^", "s", 'v']
    if merge:
        ax = pw.Brick(figsize=(2.8,4))
    else:
        fig, ax = plt.subplots(figsize=(2.8,4))
    plt.rcParams['font.family'] = 'Arial'
   
    
    if y == "Wall time (s)":
        data[y] = np.array(data[y])

    sns.pointplot(data=data, x="N50", y=y, ax=ax, 
                  hue=hue, palette=colors, 
                  markers=markers, markeredgewidth=0)
    

    plt.setp(ax.collections, alpha=.7) 
    plt.setp(ax.lines, alpha=.5)
    ax.set_xlabel("", fontsize=24)

    if "INH" in y or "IH" in y:
        ax.set_ylabel(y, fontsize=20)
    else:
        ax.set_ylabel(y, fontsize=24)
    
   
    #     ax.set_ylabel("$log_{10}$(Wall time (s))")
    ax.tick_params(axis="x", labelsize=18, rotation=45)
    # ax.set_xticklabels(fontsize=18, rotation=45)
    
    
    if np.max(data[y]) > 1000 and y != "Wall time (s)":
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        formatter = plt.gca().get_yaxis().get_major_formatter()
        ax.yaxis.set_major_formatter(formatter)
        ax.yaxis.get_offset_text().set_fontsize(14)
    
    # ax.set_yticklabels(fontsize=18)
    ax.tick_params(axis="y", labelsize=18)
    if y == "Wall time (s)":
        ax.set_yscale("log")

    ax.set_title(f"Ploidy level = {title}", fontsize=18, fontweight='bold')
    ax.legend(title="Software", bbox_to_anchor=(1.05, 1), loc='upper left' )
    # ax.legend().remove()

    if y == "Group numbers":
        chromosome_number = int(title)*5
        ax.axhline(y=chromosome_number, color='black', linestyle="--")

    ax.spines['bottom'].set_linewidth(2.0)
    ax.spines['left'].set_linewidth(2.0)
    ax.spines['top'].set_linewidth(0)
    ax.spines['right'].set_linewidth(0)
    ax.tick_params(which='both', width=2, length=5)


    if not merge:
        plt.savefig(output, dpi=600, bbox_inches='tight')
        plt.savefig(output.replace("png", "pdf"), dpi=600, bbox_inches='tight')
    
    return ax

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('tsv', 
            help='evalutation results on simulation')
    pOpt.add_argument('--merge',
                      action='store_true',
                      default=False,
                      help="merge plot into one picture")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    df = pd.read_csv(args.tsv, sep='\t', header=0, index_col=None)
    axes = OrderedDict()
    for ploid, tmp_df in df.groupby('Ploidy'):
        tmp_axes = []
        for column in tmp_df.columns:
            if column == "Software" or column == "Ploidy" or \
                column == "N50" or column == "Chromosome numbers":
                continue
            output_middle = column.replace(" ", "_").replace("(", "").replace(")", "")
            ax = plot(tmp_df[["Software", "N50", column]], ploid, "N50", 
                        column, "Software", f"{ploid}.{output_middle}.png", args.merge)
            
            if column not in ["Precision", "Recall"]:
                tmp_axes.append(ax)
        
        axes[ploid] = tmp_axes

 

    if not args.merge:
        return 
    ## merge in to one plot
    res_axes = []
    plt.rcParams['font.family'] = 'Arial'
    for j, ploid in enumerate(axes):
       
        tmp_axes = axes[ploid]
        for i, ax in enumerate(tmp_axes):
            ax.set_title("")
            if i == 0:
                ax.text(-0.5, 0.5, f"Ploidy={ploid}", 
                        ha='right', va='center', 
                        transform=ax.transAxes, fontsize=30,
                        rotation=90)
                ax.text(-0.65, 1.0, ascii_lowercase[j],
                        transform=ax.transAxes,
                        fontweight='bold',
                        fontsize=32)
                
            if i != len(tmp_axes) - 1 and ((j != len(axes) - 1) or j == 0):
                ax.legend().remove()

            

        pw.param['margin'] = 0.2
        tmp_ax = reduce(lambda x,y: x | y, tmp_axes)
        
        res_axes.append(tmp_ax)
    
    pw.param['margin'] = 0.4
    res_ax = reduce(lambda x, y: x/y, res_axes)

    res_ax.savefig("merge.png", dpi=600, bbox_inches='tight')
    res_ax.savefig("merge.pdf", dpi=600, bbox_inches='tight')

        

if __name__ == "__main__":
    main(sys.argv[1:])