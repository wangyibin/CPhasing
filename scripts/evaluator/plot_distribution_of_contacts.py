#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot the distribution of contacts 
"""

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

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    
    pReq.add_argument('cool', 
            help='')
    pOpt.add_argument('-p', '--percentile', default=95, type=int, 
                      help="percentile [default: %(default)s]")
    pOpt.add_argument('-o', '--output', type=str,
            default='output.hist.plot.png', help='output file [default: output.hist.plot.png]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    cool = cooler.Cooler(args.cool)
    bins = cool.bins()
    matrix = cool.matrix(balance=False, sparse=True)[:]
    sum_array = np.array(matrix.sum(axis=1)).T[0]
    sum_array_nonzero = sum_array[sum_array > 0]
    percentile_value = np.percentile(sum_array_nonzero, args.percentile)
    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax = sns.histplot(sum_array_nonzero, ax=ax)
    ax.axvline(percentile_value, linestyle='--', color='gray')
    ax.tick_params(axis="y", labelsize=16)
    ax.tick_params(axis="x", labelsize=16)
    ax.set_xlabel('Contacts',  fontsize=18)
    ax.set_ylabel('Density', fontsize=18)
    plt.xlim(0, np.percentile(sum_array, 99.99))
    plt.savefig(args.output, dpi=600, bbox_inches='tight')
    plt.savefig(args.output.replace("png", "pdf"), dpi=600, bbox_inches='tight')


    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax = sns.kdeplot(sum_array_nonzero, ax=ax)
    ax.axvline(percentile_value, linestyle='--', color='gray')
    ax.tick_params(axis="y", labelsize=16)
    ax.tick_params(axis="x", labelsize=16)
    ax.set_xlabel('Contacts',  fontsize=18)
    ax.set_ylabel('Density', fontsize=18)
    plt.xlim(0, np.percentile(sum_array_nonzero, 99.99))
    plt.savefig(args.output.replace("hist", "kde"), dpi=600, bbox_inches='tight')
    plt.savefig(args.output.replace("hist", "kde").replace("png", "pdf"), dpi=600, bbox_inches='tight')

if __name__ == "__main__":
    main(sys.argv[1:])
