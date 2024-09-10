#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

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
    

    df = pd.read_csv(args.tsv, sep='\s+', header=0, index_col=None)
    
    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(2, 5))

    color = ["#B3CDE3", "#FBB4AE" ]

    ax = sns.barplot(data=df, x="MapQ", y='Rate', hue='Library',
                    palette=color, errorbar=None, edgecolor='k')
    for i in ax.containers:
        ax.bar_label(i, fmt="{:.2f}")
    sns.despine()
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    plt.tick_params(which='both', width=1.5, length=5)
    plt.ylim(0, 1.0)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14,)
    plt.xlabel("Minimum mapq", fontsize=16)
    plt.ylabel("Valid Rate", fontsize=16)
    plt.legend(fontsize=12, bbox_to_anchor=(1, 0.5))
    plt.savefig('validrate.barplot.png', dpi=300, bbox_inches='tight')
    plt.savefig('validrate.barplot.pdf', dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    main(sys.argv[1:])