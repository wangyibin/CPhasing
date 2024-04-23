#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Plot the distribution of contacts in each bins
"""

import argparse
import logging
import os
import os.path as op
import sys

import cooler 
import numpy as np 
import pandas as pd
from scipy import stats

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.ticker import MaxNLocator

def plot(data, output="output"):
    fig, ax = plt.subplots(figsize=(5.7, 5))
    plt.rcParams['font.family'] = 'Arial'
    data = data[data <= np.percentile(data, 98) * 1.5]
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    kdelines = sns.kdeplot(data, ax=ax, color='#253761', linewidth=2)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    formatter = plt.gca().get_yaxis().get_major_formatter()
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.gca().yaxis.get_offset_text().set_fontsize(14)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    formatter = plt.gca().get_xaxis().get_major_formatter()
    plt.gca().xaxis.set_major_formatter(formatter)
    plt.gca().xaxis.get_offset_text().set_fontsize(14)
    plt.xlim(0, np.percentile(data, 95) * 2)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("Contacts", fontsize=24)
    plt.ylabel("Density", fontsize=24)

    x = kdelines.lines[0].get_xdata()
    y = kdelines.lines[0].get_ydata()
    max_idx = np.argmax(y)
    print(max_idx,x[max_idx] * 0.1, x[max_idx] * 1.85, y[max_idx] )
    
    ax.fill_between((x[max_idx] * 0.1, x[max_idx] * 1.85), 
                    0, ax.get_ylim()[1], alpha=0.5 , color='#bcbcbc')
    ax.axvline(x[max_idx] * 0.1, linestyle='--', color='k')
    ax.axvline(x[max_idx] * 1.85, linestyle='--', color='k')

    # plt.plot(x[max_idx], y[max_idx], ms=10, color='r')
    plt.savefig(f'{output}.kde.plot.png', dpi=600, bbox_inches='tight')
    plt.savefig(f'{output}.kde.plot.pdf', dpi=600, bbox_inches='tight')
    

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('cool_file', 
            help='')
    pOpt.add_argument('-o', '--output', type=str,
            default="output", help='output prefix of plot file [default: output]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    cool_file = args.cool_file 

    cool = cooler.Cooler(cool_file)
    bins = cool.bins()[:]
    binsize = cool.binsize
    matrix = cool.matrix(balance=False, sparse=True)[:]

    sum_values = np.array(matrix.sum(axis=1)).T[0]
    
    sum_values = np.array(matrix.sum(axis=1).T[0])
    
    small_bins = bins[bins['end'] - bins['start'] < cool.binsize]
    small_bins_sum_values = sum_values.T[small_bins.index]
    adjust_small_bins_sum_values = small_bins_sum_values.T / \
        ((small_bins['end'] - small_bins['start']) / binsize).values
    sum_values[:, small_bins.index] = adjust_small_bins_sum_values
    sum_values = np.nan_to_num(sum_values)
    # print(sum_values)

    plot(sum_values, args.output)

if __name__ == "__main__":
    main(sys.argv[1:])