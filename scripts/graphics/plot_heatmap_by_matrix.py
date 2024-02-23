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
    pReq.add_argument('matrix', 
            help='matrix format file')
    pOpt.add_argument('-c', '--cmap', default='vik')
    pOpt.add_argument('--vmin', default=None, type=float)
    pOpt.add_argument('--vmax', default=None, type=float)
    pOpt.add_argument('--center', default=None, type=float)
    pOpt.add_argument('--xlabel', default=None)
    pOpt.add_argument('-o', '--output', type=str,
            default="heatmap.png", help='output file [default: str]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    try:
        """
        https://pratiman-91.github.io/colormaps/
        """
        colormap = getattr(cmaps, args.cmap)
    except:
        colormap = args.cmap 
    
    df = pd.read_csv(args.matrix, sep='\t', index_col=0, header=0)
    plt.rcParams['font.family'] = 'Arial'
    heatmap = sns.heatmap(df, cmap=args.cmap, vmax=args.vmax, vmin=args.vmin, cbar=False)

    fig = plt.gcf()
    cbar = fig.colorbar(heatmap.collections[0], ax=heatmap.axes, shrink=0.5)
    
    
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    if args.xlabel:
        plt.xlabel(args.xlabel, fontsize=16)    

    plt.savefig(args.output, bbox_inches='tight', dpi=600)
    plt.savefig(f"{args.output.rsplit('.', 1)[0]}.pdf", bbox_inches='tight', dpi=600)
    

if __name__ == "__main__":
    main(sys.argv[1:])