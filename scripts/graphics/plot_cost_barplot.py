#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot cost/Gb barplot
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


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('table', 
            help='table')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    df = pd.read_csv(args.table, sep='\t', header=0, index_col=None, engine='python')
    
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype'] = 42
    fig, ax = plt.subplots(figsize=(3.5, 5))
    # color = [ "#8896ae", "#cb6e7f"]
    color = ["#B3CDE3", "#FBB4AE" ]
    ax = sns.barplot(data=df, x="Accession", y="Cost/Gb",
                    hue="Method", palette=color,
                    dodge=False, errorbar=None,
                    edgecolor='k')
    
    for i in ax.containers:
        ax.bar_label(i, fmt="{:.2f}", fontsize=12)

    plt.ylim(0, 120)
    sns.despine()
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    plt.tick_params(which='both', width=1.5, length=5)
    
    plt.xticks(fontsize=14, rotation=45)
    plt.yticks(fontsize=14,)
    plt.xlabel("", fontsize=14)
    plt.ylabel("Cost / Gb ($)", fontsize=18)
    plt.legend(fontsize=12)
    plt.savefig('cost.barplot.png', dpi=300, bbox_inches='tight')
    plt.savefig('cost.barplot.pdf', dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    main(sys.argv[1:])
