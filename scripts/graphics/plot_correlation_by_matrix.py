#!/usr/bin/env python
# -*- coding:utf-8 -*-


import argparse
import logging
import os
import os.path as op
import sys

import cooler 
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np 
import colormaps as cmaps
import scipy

from matplotlib.lines import Line2D

from collections import defaultdict
from itertools import combinations
from math import sqrt


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument("matrix1")
    pReq.add_argument("matrix2")
    pOpt.add_argument('--xlabel', default=0)
    pOpt.add_argument('--ylabel', default=1)
    pOpt.add_argument('-o', '--output', type=str,
            default="regplot.png", help='output file [default: str]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    df1 = pd.read_csv(args.matrix1, sep='\t', header=0, index_col=0).unstack(level=0)
    df2 = pd.read_csv(args.matrix2, sep='\t', header=0, index_col=0).unstack(level=0)

    data = pd.concat([df1.dropna(), df2.dropna()], axis=1)
    
    data = data.reset_index(drop=True)
    # data.columns = [args.xlabel, args.ylabel]

    slope, intercept, r, p, stderr = scipy.stats.linregress(x=data[0], y=data[1])
    label = [r"r = {:.2f}  $\mathit{{P}}$ = {:.2e}".format(r, p)]
    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(5.5, 6))
    
    scatter_params = dict(color='#209093', s=10)
    line_params = dict(color='#032F49', lw=2)
    legend_elements = [Line2D([0], [0], **line_params)]
    sns.regplot(x=data[0], y=data[1], scatter_kws=scatter_params, line_kws=line_params)
    plt.xticks(fontsize=14, rotation=45)
    plt.yticks(fontsize=14)

    if args.xlabel:
        plt.xlabel(args.xlabel, fontsize=18)   
    
    if args.ylabel:
        plt.ylabel(args.ylabel, fontsize=18)   
    plt.legend(legend_elements, label, loc='best', fontsize=14)
    plt.savefig(args.output, bbox_inches='tight', dpi=600)
    plt.savefig(f"{args.output.rsplit('.', 1)[0]}.pdf", bbox_inches='tight', dpi=600)


if __name__ == "__main__":
    main(sys.argv[1:])