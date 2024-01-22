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
import colormaps as cmaps
import pandas as pd 
import numpy as np 


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-f', '--hist', nargs="+",
            help='hist files')
    pOpt.add_argument('-o', '--output', type=str,
            default="high.order.distribution.each.png", help='output file [default: high.order.distribution.png]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args) 

    fig, ax = plt.subplots(figsize=(5.5,5))
    plt.rcParams['font.family'] = 'Arial'

    sample_names = list(map(lambda x: x.split(".")[0], args.hist))

  
    hist_list = []
    for hist in args.hist:
        df = pd.read_csv(hist, sep='\t', header=0, index_col=0)
        df = pd.melt(df.reset_index(), id_vars=['order'], var_name='sample')

        hist_list.append(df)
    
    hist_df = pd.concat(hist_list, axis=1)
    hist_df = hist_df.rename(columns={"value": "proportion"})
    print(hist_df)

    colors = ['#a83836', '#253761',  '#df8384', '#8dc0ed',]
    bluered_12 = list(map(lambda x: mpl.colors.rgb2hex(x.colors),  list(cmaps.bluered_12)))
    sns.pointplot(data=hist_df, x='order', y='proportion', palette=bluered_12, hue='sample', ax=ax, 
                  alpha=0.5, errorbar=None)

    plt.setp(ax.collections, alpha=.7) 
    plt.setp(ax.lines, alpha=.5)
    plt.xticks(fontsize=18, rotation=45)
    plt.yticks(fontsize=18)
    plt.xlabel('Contact order', fontsize=24)
    plt.ylabel('Proportion (%)', fontsize=24)


    output = args.output

    plt.savefig(output, dpi=600, bbox_inches='tight')
    plt.savefig(output.replace("png", "pdf"), dpi=600, bbox_inches='tight')



if __name__ == "__main__":
    main(sys.argv[1:])