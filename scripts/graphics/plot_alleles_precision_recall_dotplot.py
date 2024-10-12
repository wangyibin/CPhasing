#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot retain remove dotplot
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

import pandas as pd
import numpy as np

from matplotlib.ticker import MaxNLocator


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('tsv', 
            help='')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    markers = {"50 kb": "^", 
                "100 kb": "D",
                 "500 kb": "o", 
                 "1 Mb": "s", 
                 "2 Mb": "v"}
    df = pd.read_csv(args.tsv, sep='\t', header=0, index_col=None)

    colors =  ["#cb6e7f", "#253761", "#8896ae"]
    plt.rcParams['font.family'] = 'Arial'
    for ploidy, tmp_df in df.groupby('Ploidy'):

        fig, ax = plt.subplots(figsize=(4, 4))

        sns.scatterplot(data=tmp_df, x="Precision", y="Recall", style=tmp_df['N50'], hue="Software",
                        hue_order=["C-Phasing", "HapHiC", "ALLHiC"],
                        palette=colors, ax=ax, markers=markers, alpha=0.7)
        
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(0)
        ax.spines['right'].set_linewidth(0)
        ax.tick_params(which='both', width=2, length=5)
        ax.set_ylim(-0.1, 1.1)
        ax.set_xlim(-0.1, 1.1)
        ax.set_title(f"Ploidy = {tmp_df['Ploidy'].iloc[0]}", fontsize=16)
        ax.tick_params(axis="both", labelsize=14)
        ax.set_xlabel("Precision", fontsize=16)
        ax.set_ylabel("Recall", fontsize=16)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left' )
        ax.legend().remove()
        plt.savefig(f"{ploidy}.allele.dotplot.png", dpi=600, bbox_inches='tight')
        plt.savefig(f"{ploidy}.allele.dotplot.pdf", dpi=600, bbox_inches='tight')




if __name__ == "__main__":
    main(sys.argv[1:])