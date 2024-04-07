#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot the porec concatmers 
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 
import numpy as np 

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.collections import LineCollection, PathCollection
from matplotlib.patches import Rectangle
from itertools import combinations

from cphasing.core import PoreCTable
from cphasing.utilities import humanized2numeric



def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('porec', 
            help='porec table')
    pReq.add_argument('contig', 
                      help="which contig to plot")
    pOpt.add_argument('-b', '--binsize', default="100", 
                      help='binsize')
    pOpt.add_argument('-q', '--min_mapq', default=2, type=int,
                      help='minimum mapping quality of fragments')
    pOpt.add_argument('--no_annotate_chimeric', default=False, action='store_true',
                      help='dont annotate the chimeric position')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    


    
    args = p.parse_args(args)
    
    binsize = humanized2numeric(args.binsize)

    df = pd.read_csv(args.porec, sep='\t', names=PoreCTable.HEADER, 
                        header=None, dtype=PoreCTable.META, 
                        usecols=['read_idx', 'chrom',
                                    'start', 'end', 'mapping_quality'])
    
   
    df = df.query(f'mapping_quality >= {args.min_mapq} ')
    df.drop(['mapping_quality'], axis=1, inplace=True)
    df.set_index('chrom', inplace=True)
    
    df = df.loc[args.contig]
    
    df = df.reset_index()
    df = df.set_index('read_idx')
    
    df_grouped = df.groupby('read_idx')

    df = df.loc[df_grouped['start'].nunique() > 1]
    
    df['position'] = ((df['end'] + df['start']) // 2) // binsize 
    df.drop(['start', 'end'], axis=1, inplace=True)
    
    df = df.loc[df.groupby('read_idx')['position'].nunique() > 1]
    df = df.sort_values('position')
    max_position = max(df['position'])

    idx = 1
    lines = []
    dots = []
    for i, tmp_df in df.groupby('read_idx', sort=False):
        positions = tmp_df.position.values.tolist()
        if len(positions) <= 1:
            continue
        line = [(pos, idx) for pos in positions]
        dots.extend(line)
        line = sorted(line, key=lambda x: x[0]) 
        
        lines.append([line[0], line[-1]])
        idx += 1

    dots_x, dots_y = list(zip(*dots))
    
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)

    line_segments = LineCollection(lines, linewidth=0.3, color="grey")

    ax.add_collection(line_segments)

    ax.scatter(dots_x, dots_y, c='r', s=0.005)
    ax.autoscale()
    ax.set_xlabel(args.contig)
    
    # ax.set_xticks(m)
    plt.yticks([])
    
    sns.despine(trim=True, left=True)
    plt.savefig('output.png', dpi=600)

    
    
    






if __name__ == "__main__":
    main(sys.argv[1:])