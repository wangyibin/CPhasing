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
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import pandas as pd 

import numpy as np 

from collections import OrderedDict
from math import sqrt
from joblib import Parallel, delayed
from pandarallel import pandarallel

logger = logging.getLogger(__name__)

PAF_HADER = ["contig1", "length1", "start1", "end1", "strand",
                 "contig2", "length2", "start2", "end2", "mismatch",
                 "matches", "mapq", "identity"]

def read_paf(paf):
    logger.info(f"Load alignments results `{paf}`")
    df = pd.read_csv(paf, sep='\t', header=None, usecols=range(13),
                        names=PAF_HADER, index_col=None)
    df['identity'] = df['identity'].map(lambda x: x.replace("id:f:", "")).astype('float64')
    df = df.sort_values(['contig2', 'start2'])

    return df 

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')

    pOpt = p.add_argument_group('Optional arguments')
    pOpt.add_argument('-c', '--cool', 
            help='Path to cool file',
            required=True)
    pOpt.add_argument('-p', '--paf', 
            help='Path to wfmash (-m) paf file',
            required=True)
    pOpt.add_argument('--min-identity', type=float, default=0.9,
            help='Minumum identity')
    pOpt.add_argument('--min-length', type=float, default=100000,
            help='Minumum length of alignments')
    pOpt.add_argument('-o', '--output', type=str,
            default="output", help='default `output`')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    cool = cooler.Cooler(args.cool)
    paf = read_paf(args.paf)

    matrix = cool.matrix(balance=False, sparse=True)

    def func(row):
        cis1 = matrix.fetch(f"{row['contig1']}:{row['start1']}-{row['end1']}").sum()
        cis2 = matrix.fetch(f"{row['contig2']}:{row['start2']}-{row['end2']}").sum()
        trans = matrix.fetch(f"{row['contig1']}:{row['start1']}-{row['end1']}", 
                             f"{row['contig2']}:{row['start2']}-{row['end2']}").sum()
        
        if cis1 * cis2 > 0:
            ntrans = trans / sqrt(cis1 * cis2)
        else:
            ntrans = np.nan

        identity = row['identity']

        return pd.Series([ntrans, identity], index=['ntrans', 'identity'])

    paf = paf[(paf['end1'] - paf['start1']) >= args.min_length]
    paf = paf[paf['identity'] >= args.min_identity ]
    
    paf['identity'] = paf['identity'] * 100

    pandarallel.initialize(nb_workers=40, verbose=0, progress_bar=True)
    res = paf.parallel_apply(func, axis=1)
    # .sort_values(['identity'], ascending=False)

    fig, ax = plt.subplots(figsize=(5.7, 5))
    
    sns.regplot(data=res, y='ntrans', x='identity', fit_reg=False,
                scatter_kws={"s": 5, "color": '#bcbcbc', "alpha": 0.5})

    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    plt.tick_params(which='both', width=1.5, length=5)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    sns.despine()
    plt.ylim(0, 2.0)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("Alignment identity (%)", fontsize=20)
    plt.ylabel("Normalized ${trans}$ contacts", fontsize=20)


    plt.savefig(f"{args.output}.png", dpi=600, bbox_inches='tight')
    plt.savefig(f"{args.output}.pdf", dpi=600, bbox_inches='tight')
    
    

if __name__ == "__main__":
    main(sys.argv[1:])