#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
calculate the cis/trans ratio according by contigs
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

from collections import OrderedDict
from math import sqrt
from joblib import Parallel, delayed

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pOpt.add_argument('-c', '--cool', 
            nargs="+",
            help='Path to cool file',
            required=True)
    pOpt.add_argument('-p', '--contig_pairs',
            help="the pairs of contigs")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    res = OrderedDict()
    
    for cool_file in args.cool:
        prefix = cool_file.rsplit(".", 1)[0]
        cool = cooler.Cooler(cool_file)
        matrix = cool.matrix(balance=False, sparse=True)
        contig_paris = [i.strip().split()[:2] for i in open(args.contig_pairs) if i.strip()]
        
        data = []

        cmd_args = []
        for pair in contig_paris:
            if pair[0] > pair[1]:
                continue
            cmd_args.append(pair)
        
        def func(pair):
            cis1 = matrix.fetch(pair[0]).sum()
            cis2 = matrix.fetch(pair[1]).sum()
            trans = matrix.fetch(*pair).sum()
        
            return trans / sqrt(cis1 * cis2)

        data = Parallel(n_jobs=40)(
            delayed(func)(i) for i in cmd_args
        )

        res[prefix] = list(filter(lambda x: ~np.isinf(x),  data))

    
    res_df = pd.DataFrame(res)
    
    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(5.5, 5))
    boxprops = dict(color='black', linewidth=1.5)
    medianprops=dict(color='black', linewidth=2.5)
    whiskerprops = dict(linestyle='--')
    colors = ['#a83836', "#253761",  '#df8384', '#8dc0ed',]
    # bplot = ax.boxplot(data, 
    #                    showfliers=False, 
    #                    patch_artist=True, 
    #                    notch=True, 
    #                    widths=0.35,
    #                    medianprops=medianprops,
    #                    whiskerprops=whiskerprops,
    #                    boxprops=boxprops)
    ax = sns.violinplot(data=res_df, ax=ax, palette=colors)
    # for patch, color in zip(bplot['boxes'], ['#a83836', '#df8384', '#8dc0ed']):
    #     patch.set_facecolor(color)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Before\n realign", "After\n realign", ], fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=16)
    ax.set_ylabel("Normalized contacts", fontsize=22)
    plt.savefig('boxplot.png', dpi=600, bbox_inches='tight')

        
if __name__ == "__main__":
    main(sys.argv[1:])