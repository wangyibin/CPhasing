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

import pandas as pd 

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.ticker import MaxNLocator



def read_csv(data):
    df = pd.read_csv(data, sep='\t', header=0, index_col=None)
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
            help='number of bins [default: %(default)%]')
    pOpt.add_argument('--stat',
            default='density',
            help='aggregate statistic to compute in each bin. [default: %(default)s]',
            choices=['count', 'frequency', 'probability', 'percent', 'density'])
    pOpt.add_argument('--x-min', 
            default=None,
            type=float,
            help='minimum value of X axis')
    pOpt.add_argument('--x-max',
            default=None,
            type=float,
            help='maximum value of X axis')
    pOpt.add_argument('--x-label',
            default=None,
            help='x label text [default: %(default)s]')
    pOpt.add_argument('--x-label',
            default=None,
            help='y label text [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    df = read_csv(args.data[0])

    plt.rcParams['font.family'] = 'Helvetica'
    fig, ax = plt.subplots(figsize=(5.7, 5))
    ax = sns.histplot(df, kde=True, color='r', alpha=0.3, stat=args.stat, #linewidth=0,
                      bins=args.bins)
    
    if args.x_label:
        ax.set_xlabel(args.x_label.capitalize(), fontsize=20)
    ax.set_ylabel(args.stat.capitalize(), fontsize=20)
    if args.x_min and args.x_max:
        plt.xlim(args.x_min, args.x_max)

    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    plt.xticks(fontsize=18)
    
    plt.yticks(fontsize=18)
    # plt.legend(fontsize=16)
    


    plt.savefig(f'{args.data[0]}.hist.plot.png', dpi=600, bbox_inches='tight')
    plt.savefig(f'{args.data[0]}.hist.plot.pdf', dpi=600, bbox_inches='tight')

if __name__ == "__main__":
    main(sys.argv[1:])