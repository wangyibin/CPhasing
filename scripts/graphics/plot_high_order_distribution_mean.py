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
    pReq.add_argument('-s', '--summary', nargs="+",
            help='porec contact summary from cphasing-rs porec2pairs')
    pOpt.add_argument('-o', '--output', type=str,
            default="high.order.distribution.png", help='output file [default: high.order.distribution.png]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args) 

    fig, ax = plt.subplots(figsize=(5.5,5))
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype'] = 42
    sample_names = list(map(lambda x: x.split(".")[0], args.summary))
  

    df_list = []
    for summary in args.summary:
        df = pd.read_csv(summary, sep='\t', index_col=0, header=None)
    
        df_list.append(df)
    df = pd.concat(df_list, axis=1) if len(df_list) > 1 else df_list[0]

    length_bins = pd.IntervalIndex.from_breaks([1, 2, 3, 4, 6, 11, 
                                                    21, 50, int(1e9)])
    length_bin_labels = {}
    for i in length_bins:
        if i.length == 1:
            label = str(i.right)
        elif i.length >= int(1e8):
            label = f">={i.left}"
        else:
            label = f"{i.left + 1}-{i.right}"
        length_bin_labels[i] = label

    hist_list = []
    for df in df_list:
        df = df.reset_index()

        read_order_hist = (pd.cut(df[0], length_bins, labels=length_bin_labels)
                                    )
        df['index'] = read_order_hist
        
        hist = df.groupby('index').sum().reset_index()
      
        hist.drop(0, axis=1, inplace=True)
        hist.columns = ['order', 'count']
        total_count = hist['count'].sum()
        hist['proportion'] = hist['count'] / total_count * 100

    
        labels = length_bin_labels.values()
        hist['order'] = labels
    
        hist_list.append(hist)
    
    hist_df = pd.concat(hist_list, axis=1) if len(hist_list) > 1 else hist
    hist_df = hist_df.drop(['order', 'count'], axis=1)
    hist_df.columns = sample_names
    hist_df['order'] = labels
    hist_df = hist_df.set_index('order')
    
    
    hist = hist.set_index('order')
    hist = hist.drop('count', axis=1)
    if len(hist_list) > 1: 
        hist['proportion'] = hist_df.mean(axis=1)
        hist = hist.reset_index()
    else:
        hist = hist_df.reset_index()
        hist.columns = ['order', 'proportion']



    colors = ['#a83836', "#253761",  '#df8384', '#8dc0ed',]
    
    sns.pointplot(data=hist, x='order', y='proportion', color=colors[0], ax=ax)

    plt.xticks(fontsize=18, rotation=45)
    plt.yticks(fontsize=18)
    plt.xlabel('Contact order', fontsize=20)
    plt.ylabel('Proportion (%)', fontsize=20)


    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    plt.tick_params(which='both', width=1.5, length=5)
    sns.despine()

    output = args.output
    hist.to_csv(output.replace("png", "mean.hist"), sep='\t', header=True, index=False)
    hist_df.to_csv(output.replace("png", "hist"), sep='\t', header=True, index=True)

    plt.savefig(output, dpi=600, bbox_inches='tight')
    plt.savefig(output.replace("png", "pdf"), dpi=600, bbox_inches='tight')



if __name__ == "__main__":
    main(sys.argv[1:])