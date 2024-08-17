#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot histogram 
"""

import argparse
import logging
import os
import os.path as op
import sys
import gc
import pandas as pd 

import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
import colormaps as cmaps

bluered_12 = list(map(lambda x: mpl.colors.rgb2hex(x.colors),  list(cmaps.bluered_12)))


def read_csv(data):
    df = pd.read_csv(data, sep='\t', header=0, index_col=None, )
    return df


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('data', nargs="+", 
            help='')
    
    pOpt.add_argument('--bins',
            default=50,
            type=int,
            help='number of bins [default: %(default)s]')
    pOpt.add_argument('--stat',
            default='density',
            help='aggregate statistic to compute in each bin. [default: %(default)s]',
            choices=['count', 'frequency', 'probability', 'density'])
    pOpt.add_argument('--x-min', 
            dest="x_min",
            default=None,
            type=float,
            help='minimum value of X axis')
    pOpt.add_argument('--x-max',
            dest="x_max",
            default=None,
            type=float,
            help='maximum value of X axis')
    pOpt.add_argument('--x-label',
            default=None,
            help='x label text [default: %(default)s]')
    pOpt.add_argument('--x-label',
            default=None,
            help='y label text [default: %(default)s]')
    pOpt.add_argument('--kde', default=True, action='store_false', 
                        help="plot kde line [default: %(default)s]")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    fig, ax = plt.subplots(figsize=(5.7, 5))

    df = read_csv(args.data[0])
    res = []
    if args.x_min is not None and args.x_max is not None:
        for group in df.columns:
            tmp_df = df[(df[group] <= args.x_max)]
            res.append(tmp_df)

        df = pd.concat(res, axis=0)
        del res, tmp_df
        gc.collect()

    plt.rcParams['font.family'] = 'Arial'
    
    for i in range(len(df.columns)):
        data = df.iloc[:, i]

        ax = sns.histplot(data, kde=args.kde, color=bluered_12[i], alpha=0.3, stat=args.stat, #linewidth=0,
                      bins=args.bins)
    
    legend_items = [Line2D([], [], color=color, marker=None, linewidth=2.5,
                               markersize=10, label=f'{sample}') for sample, (color, marker) in list(zip(df.columns, zip(bluered_12, range(len(df.columns)))))]
    second_legend = ax.legend(handles=legend_items, fontsize=12, loc='best')
    ax.add_artist(second_legend)

    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    plt.tick_params(which='both', width=1.5, length=5)
    sns.despine()
    
    if args.x_label:
        ax.set_xlabel(args.x_label.capitalize(), fontsize=20)
    ax.set_ylabel(args.stat.capitalize(), fontsize=20)
    if args.x_min is not None and args.x_max is not None:
        plt.xlim(args.x_min, args.x_max)

    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    plt.xticks(fontsize=18, rotation=45, ha='right')

    plt.yticks(fontsize=18)
    
    # plt.legend(fontsize=16)
    


    plt.savefig(f'{args.data[0]}.hist.plot.png', dpi=600, bbox_inches='tight')
    plt.savefig(f'{args.data[0]}.hist.plot.pdf', dpi=600, bbox_inches='tight')

if __name__ == "__main__":
    main(sys.argv[1:])