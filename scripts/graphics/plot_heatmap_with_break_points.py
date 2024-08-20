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
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import seaborn as sns
import colormaps as cmaps

import cooler 
import pandas as pd 
import numpy as np 

from cphasing.plot import plot_heatmap

def darken_color(color, amount=0.9):
    """
    Darken the given color by the specified amount.
    
    Parameters:
    color (str): The color to darken.
    amount (float): The amount to darken the color. Should be between 0 and 1.
    
    Returns:
    str: The darkened color in hex format.
    """
    c = mcolors.hex2color(color)
    c = [max(0, min(1, x * amount)) for x in c]
    return mcolors.to_hex(c)

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('cool', 
            help='')
    pOpt.add_argument('-c', '--chromosome', nargs="+")
    pOpt.add_argument('-b', '--break_pos', 
                      help="break pos information, e.g. C-Phasing:2000,HapHiC:2000", default=None)
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
                      triangle=True, chromosomes=args.chromosome, 
                      rotate_xticks=True)
    
    color_db = ["#FBB4AE" ,"#B3CDE3", "#CCEBC5", "#DECBE4"]
    marker_db = ["^", "v", "2", "1"]
    color_db += []
    ylim = ax.get_ylim()
    if args.break_pos:
        l = args.break_pos.split(",")
        if len(l) > 0:
            l2 = list(map(lambda x: x.split(":"), l))
            l2 = list(filter(lambda x: len(x) == 2, l2))
            if len(l2) > 0:
                for i, (software, pos) in enumerate(l2):
                    
                    pos = int(pos) // binsize
                    y_pos = -(ylim[1]/40) if i % 2 == 0 else (ylim[1]/40)
                    marker = marker_db[i%3]
                    ax.plot(pos, y_pos, marker=marker, markersize=7.5, 
                    linewidth=1, alpha=1,
                    linestyle='none', color=color_db[i],
                    fillstyle='none', label=software)
                   

        plt.legend(fontsize=8)
            
    # ax.plot(pos, -5, marker='^', markersize=5, 
    #             linewidth=1, alpha=0.75,
    #             fillstyle='none', color='blue')
    # ax.plot(pos3, -5, marker='^', markersize=5, 
    #             linewidth=1, alpha=0.75,
    #             fillstyle='none', color='green')

    sns.despine(left=True, bottom=True)
    plt.yticks([])
    plt.xlabel(" ".join(args.chromosome), fontsize=14)
    plt.ylim(-(ylim[1]/15))
    output_prefix = ",".join(args.chromosome)
    plt.savefig(f"{output_prefix}.tria.png", dpi=600, bbox_inches='tight')
    plt.savefig(f"{output_prefix}.tria.pdf", dpi=600, bbox_inches='tight')


if __name__ == "__main__":
    main(sys.argv[1:])
