#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot the rho boxplot
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

from statannotations.Annotator import Annotator

def import_rho_tsv(tsv):
    df = pd.read_csv(tsv, sep='\t', header=None, index_col=None)

    return np.abs(df[1].values)

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('tsv', nargs="+",
            help='')
    pOpt.add_argument('-o', '--output', type=str,
            default="boxplot.png", help='output picture [default: "boxplot.png"]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    datas = [import_rho_tsv(tsv) for tsv in args.tsv ]
    

    colors = {"C-Phasing (Pore-C)": "#b02418",
              "C-Phasing (Hi-C)": "#cb6e7f",
              "HapHiC": "#253761",
              "ALLHiC": "#8896ae"}

    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(3, 4.5))

    decribes = pd.DataFrame(datas).T.describe()
    decribes.columns = list(colors.keys())
    print(decribes)
    sns.boxplot(datas, width=.6, dodge=False, ax=ax, showfliers=False, palette=colors.values())
    sns.stripplot(datas, dodge=False, color=".3")
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    plt.tick_params(which='both', width=1.5, length=5)
    plt.ylim(0.0, 1.1)
    
    sns.despine()

    plt.xticks(fontsize=16, rotation=90)
    ax.set_xticklabels(list(colors.keys()))
    plt.yticks(fontsize=16)
    plt.ylabel(r'$\rho$', fontsize=20)

    plt.savefig(args.output, dpi=600, bbox_inches='tight')
    plt.savefig(args.output.replace("png", "pdf"), dpi=600, bbox_inches='tight')
    

if __name__ == "__main__":
    main(sys.argv[1:])