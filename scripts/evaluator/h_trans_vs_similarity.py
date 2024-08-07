#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
evaluate the trans vs similarity by HG002
"""

import argparse
import logging
import os
import os.path as op
import sys

import concurrent.futures
import cooler
import pandas as pd
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from pandarallel import pandarallel

from math import sqrt

PAF_HADER = ["contig1", "length1", "start1", "end1", "strand",
                 "contig2", "length2", "start2", "end2", "mismatch",
                 "matches", "mapq", "identity"]


def chrRangeID(args, axis=1):
    """
    Chrom range transformation.
    Examples:
    --------
    >>> args = ["Chr1", 100, 200]
    >>> chrRangeID(args)
    "Chr1:100-200"
    >>> args = "Chr1:100-200"
    >>> chrRangeID(args, axis=1)
    ("Chr1", "100", "200")
    """
    if axis == 0:
        chrom, start, end = map(str, args)
        return "{}:{}-{}".format(chrom, start, end)
    elif axis == 1:
        chrom, ranges = args.split(':')
        start, end = ranges.split('-')
        return pd.Series((chrom, start, end))
    else:
        return 
    

def chrRangeID_func_1(row):
    return f"{row.chrom1}:{row.start1}-{row.end1}"

def chrRangeID_func_2(row):
    return f"{row.chrom2}:{row.start2}-{row.end2}"


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('paf', 
            help='self align results of window slided of hg002 by wfmash')
    pReq.add_argument('cool', 
            help='cool of hg002, binsize equal to the window of paf')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    df = pd.read_csv(args.paf, sep='\t', header=None, usecols=range(13),
                         names=PAF_HADER, index_col=None)
    df['identity'] = df['identity'].map(lambda x: x.replace("id:f:", "")).astype(float)
    
    df['contig1'] = df['contig1'].str.replace("_sliding", "")
    df['contig2'] = df['contig2'].str.replace("_sliding", "")
    
    res1 = df['contig1'].apply(chrRangeID)
    res2 = df['contig2'].apply(chrRangeID)
    
    df['contig1'] = res1[0]
    df['start1'] = (res1[1]).astype(int) - 1
    df['end1'] = res1[2].astype(int)

    df['contig2'] = res2[0]
    df['start2'] = (res2[1]).astype(int) - 1
    df['end2'] = res2[2].astype(int)
    
    cool = cooler.Cooler(args.cool)
    pixels = cool.pixels(join=True)[:]
    pixels_2 = pixels.copy()
    pixels_2.columns = ['chrom2', 'start2', 'end2', 
                        'chrom1', 'start1', 'end1', 'count']
    
    pixels = pd.concat([pixels, pixels_2], axis=0)
    pandarallel.initialize(nb_workers=10, verbose=0)
    pixels['contig1'] = pixels.parallel_apply(chrRangeID_func_1, axis=1)
    pixels['contig2'] = pixels.parallel_apply(chrRangeID_func_2, axis=1)

    count_db = pixels.set_index(['contig1', 'contig2']).to_dict()['count']


    data = []
    similarity = []
    for i, row in df.iterrows():
        try:

            cis1 = count_db[(f"{row.contig1}:{row.start1}-{row.end1}", 
                                 f"{row.contig1}:{row.start1}-{row.end1}")]
            cis2 = count_db[(f"{row.contig2}:{row.start2}-{row.end2}", 
                                 f"{row.contig2}:{row.start2}-{row.end2}")]
            trans = count_db[(f"{row.contig1}:{row.start1}-{row.end1}", 
                                 f"{row.contig2}:{row.start2}-{row.end2}")]
        except:
            continue 
        
        nh = trans / sqrt(cis1 * cis2)

        if row.identity >= 0.999 and row.identity <= 1.0 and nh < 1.0:
            data.append(nh)

            similarity.append(-np.log10(row.identity))

    
    



    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(5.5, 5))

    scatter_params = dict(color='#209093', s=10)
    line_params = dict(color='#032F49', lw=2)
    legend_elements = [Line2D([0], [0], **line_params)]
    sns.regplot(x=data, y=similarity, scatter_kws=scatter_params, line_kws=line_params)
    
    plt.xlim(0, 1)
    plt.xlabel("Normalized trans contacts", fontsize=16)
    plt.ylabel("Identity", fontsize=16)



    plt.savefig('trans_vs_similarity.replot.png', dpi=600, bbox_inches='tight')


if __name__ == "__main__":
    main(sys.argv[1:])