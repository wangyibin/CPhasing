#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot the performance of allelic identification
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

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
    

    df = pd.read_csv(args.tsv, sep='\t', header=0, index_col=None)

    colors =  ["#cb6e7f", "#253761", "#8896ae"]
    markers = ["^", "s", 'v']
    figsize = (5, 2)
    plt.rcParams['font.family'] = 'Arial'
    
    
    for item in ["Precision", "Recall", "F1 score"]:

        for ploidy, tmp_df in df.groupby('Ploidy'):
            fig, ax = plt.subplots(figsize=figsize)
            sns.pointplot(data=tmp_df, x="N50", y=item, hue="Software",
                        hue_order=["C-Phasing", "HapHiC", "ALLHiC"],
                        palette=colors, ax=ax, markers=markers,
                        )

            plt.setp(ax.collections, alpha=.7) 
            plt.setp(ax.lines, alpha=.7)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
            ax.spines['bottom'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            ax.spines['top'].set_linewidth(0)
            ax.spines['right'].set_linewidth(0)
            ax.tick_params(which='both', width=2, length=5, labelsize=14)

            plt.ylim(0, 1.05)

            ax.set_title(f"Ploidy = {tmp_df['Ploidy'].iloc[0]}", fontsize=16)
            ax.set_xlabel("N50", fontsize=16)
            ax.set_ylabel(item, fontsize=16)
            # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left' )
            ax.legend().remove()
            plt.savefig(f"{ploidy}.allele.{item.split()[0]}.N50.lineplot.png", dpi=600, bbox_inches='tight')
            plt.savefig(f"{ploidy}.allele.{item.split()[0]}.N50.lineplot.pdf", dpi=600, bbox_inches='tight')


        for N50, tmp_df in df.groupby('N50'):
            fig, ax = plt.subplots(figsize=figsize)
            sns.pointplot(data=tmp_df, x="Ploidy", y=item, hue="Software",
                        hue_order=["C-Phasing", "HapHiC", "ALLHiC"],
                        palette=colors, ax=ax, markers=markers,
                        )

            plt.setp(ax.collections, alpha=.7) 
            plt.setp(ax.lines, alpha=.7)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.spines['bottom'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            ax.spines['top'].set_linewidth(0)
            ax.spines['right'].set_linewidth(0)
            ax.tick_params(which='both', width=2, length=5, labelsize=14)

            plt.ylim(0, 1.05)

            ax.set_title(f"N50 = {tmp_df['N50'].iloc[0]}", fontsize=16)
            ax.set_xlabel("Ploidy", fontsize=16)
            ax.set_ylabel(item, fontsize=16)
            # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left' )
            ax.legend().remove()
            plt.savefig(f"{N50.replace(" ", "_")}.allele.{item.split()[0]}.Ploidy.lineplot.png", dpi=600, bbox_inches='tight')
            plt.savefig(f"{N50.replace(" ", "_")}.allele.{item.split()[0]}.Ploidy.lineplot.pdf", dpi=600, bbox_inches='tight')


if __name__ == "__main__":
    main(sys.argv[1:])