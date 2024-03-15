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
    sns.kdeplot(data, ax=ax, color='#253761', linewidth=2)
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
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    cool_file = args.cool_file 

    cool = cooler.Cooler(cool_file)
    bins = cool.bins()[:]
    matrix = cool.matrix(balance=False, sparse=True)[:]
    sum_values = np.array(matrix.sum(axis=1)).T[0]

    plot(sum_values)

if __name__ == "__main__":
    main(sys.argv[1:])