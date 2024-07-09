#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot the boxplot of each porec library
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

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('table', 
            help='table file')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    df = pd.read_csv(args.table, sep='\t', header=0, index_col=None)
    
    plt.rcParams['font.family'] = 'Arial'
    # plt.rcParams['font.weight'] = 'bold'
    fig, ax = plt.subplots(figsize=(6, 5))
    color = ["#B3CDE3", "#FBB4AE" ]
    # color = [ "#8896ae", "#cb6e7f"]
    sns.boxplot(data=df, x="Accession", y="Yield (Gb)",
                ax=ax, hue='Method', palette=color,
                dodge=False, showfliers=False,
                width=.6)
    s = sns.stripplot(data=df, 
                      x="Accession", y="Yield (Gb)",
                      dodge=False,
                  ax=ax, color=".3")
    
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    plt.tick_params(which='both', width=1.5, length=5)

    sns.despine()
    plt.xticks(fontsize=14, rotation=45)
    plt.yticks(fontsize=14,)
    plt.xlabel("", fontsize=16)
    plt.ylabel("Yield (Gb)", fontsize=18) #, fontdict={'weight': 'bold'})
    plt.legend(fontsize=14)
    plt.savefig('yield.boxplot.png', dpi=300, bbox_inches='tight')
    plt.savefig('yield.boxplot.pdf', dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    main(sys.argv[1:])