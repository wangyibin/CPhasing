#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

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

import cooler 
import pandas as pd 
import numpy as np 

from cphasing.plot import plot_heatmap


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('cool', 
            help='')
    pOpt.add_argument('-c', '--chromosome', nargs="+",
                      )
    pOpt.add_argument('--vmax', default=3, type=int, )
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    cool = cooler.Cooler(args.cool)
    binsize = cool.binsize
    bins = cool.bins()[:].reset_index(drop=False)
    bins['chrom'] = bins['chrom'].astype('str')
    bins = bins.set_index('chrom')

    bins = bins.loc[args.chromosome]

    ax = plot_heatmap(args.cool, "test.png", 
                      fontsize=12, scale="none",
                      vmax=args.vmax,
                      triangle=True, chromosomes=args.chromosome)
    pos = 1786500 // binsize
    pos2 = 1703500 // binsize
    pos3 = 1703500 // binsize
#     ax.plot(pos, 10, marker='v', markersize=4, 
#                 linewidth=1, alpha=0.75,
#                 fillstyle='none', color='blue')
#     ax.plot(pos3, 10, marker='v', markersize=4, 
#                 linewidth=1, alpha=0.75,
#                 fillstyle='none', color='green')
    ax.
  
    plt.savefig("test.png", dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    main(sys.argv[1:])
