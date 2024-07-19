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
import matplotlib.lines as mlines
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

    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(6, 5))

    sample_names = list(map(lambda x: x.split(".")[0], args.hist))

  
    hist_list = []
    for hist in args.hist:
        df = pd.read_csv(hist, sep='\t', header=0, index_col=0)
        header_name = df.columns
        df = pd.melt(df.reset_index(), id_vars=['order'], var_name='sample')

        hist_list.append(df)
    
    hist_df = pd.concat(hist_list, axis=1)
    hist_df = hist_df.rename(columns={"value": "proportion"})
    print(hist_df)
    
    # colors = ['#a83836', '#253761',  '#df8384', '#8dc0ed',]
    # colors = ["#B3CDE3", "#B3CDE3", "#FBB4AE", "#FBB4AE", "#FBB4AE"]
    # colors = ['#253761', '#253761',  '#a83836', '#a83836', '#a83836']
    colors = ['#8dc0ed', '#8dc0ed', '#8dc0ed', '#df8384', '#df8384', '#df8384', '#df8384', ]

#     markers = ["D", "^", "p", "s", 'v', '*']
    markers = ["^", "^", '^', 'D', 'D', 'D', 'D']
    # markers = ['D', 'D', 'D', 'D']
    bluered_12 = list(map(lambda x: mpl.colors.rgb2hex(x.colors),  list(cmaps.bluered_12)))
    ax = sns.pointplot(data=hist_df, x='order', y='proportion', markers=markers, palette=bluered_12, hue='sample', ax=ax, 
                  alpha=0.7, errorbar=None, dodge=0.3, color='k')
    ax.get_legend().set_visible(False)
    legend_items = [mlines.Line2D([], [], color=color, marker=None, linewidth=2.5,
                              markersize=10, label=f'{sample}') for sample, (color, marker) in list(zip(header_name, zip(bluered_12, markers)))]
    second_legend = ax.legend(handles=legend_items, fontsize=12, loc='best')
    ax.add_artist(second_legend)
    legend_items = [mlines.Line2D([], [], color="k", label="Pore-C", marker=markers[0], linewidth=0),
                    mlines.Line2D([], [], color="k", label="ePore-C", marker=markers[-1], linewidth=0)]
    second_legend = ax.legend(handles=legend_items, fontsize=12, loc='center right')
    ax.add_artist(second_legend)

    
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    plt.tick_params(which='both', width=1.5, length=5)
    plt.setp(ax.collections, alpha=0.7) 
    plt.setp(ax.lines, alpha=0.7)
    plt.xticks(fontsize=14, rotation=45)
    plt.yticks(fontsize=14)
    plt.xlabel('Contact order', fontsize=18)
    plt.ylabel('Proportion (%)', fontsize=18)


    output = args.output
    sns.despine()
    plt.savefig(output, dpi=600, bbox_inches='tight')
    plt.savefig(output.replace("png", "pdf"), dpi=600, bbox_inches='tight')

    fig, ax = plt.subplots(figsize=(6, 5))
    sns.barplot(data=hist_df, x='order', y='proportion', hue='sample', alpha=0.7, edgecolor='k')
    plt.xticks(fontsize=14, rotation=45)
    plt.yticks(fontsize=14)
    plt.xlabel('Contact order', fontsize=18)
    plt.ylabel('Proportion (%)', fontsize=18)
    sns.despine()
    plt.savefig('barplot.png', dpi=600, bbox_inches='tight')
    plt.savefig('barplot.png'.replace("png", "pdf"), dpi=600, bbox_inches='tight')

if __name__ == "__main__":
    main(sys.argv[1:])