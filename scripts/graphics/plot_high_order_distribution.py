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

import pandas as pd 
import numpy as np 


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('summary', 
            help='porec contact summary from cphasing-rs porec2pairs')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args) 

    df = pd.read_csv(args.summary, sep='\t', index_col=0, header=None)
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
    

    fig, ax = plt.subplots(figsize=(5.5,5))
    plt.rcParams['font.family'] = 'Helvetica'
    
    colors = ['#a83836', "#253761",  '#df8384', '#8dc0ed',]
    sns.pointplot(data=hist, x='order', y='proportion', color=colors[0])
    
    plt.xticks(fontsize=18, rotation=45)
    plt.yticks(fontsize=18)
    plt.xlabel('Contact order', fontsize=24)
    plt.ylabel('Proportion (%)', fontsize=24)


    output = args.summary.replace("summary", "png")
    with open(output.replace("png", "hist"), 'w') as fp: 
        print(hist, file=fp)

    plt.savefig(output, dpi=600, bbox_inches='tight')
    plt.savefig(output.replace("png", "pdf"), dpi=600, bbox_inches='tight')



if __name__ == "__main__":
    main(sys.argv[1:])