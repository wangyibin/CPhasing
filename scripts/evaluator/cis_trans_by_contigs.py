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
    pOpt.add_argument('-c', '--contacts', 
            nargs="+",
            help='Path to contacts file',
            required=True)
    pOpt.add_argument('-p', '--contig_pairs',
            help="the pairs of contigs")
    pOpt.add_argument('-o', '--output', type=str,
            default="cis_trans_contigs.tsv", help='is_trans_contigs.tsv')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    res = []
    
    for contact_file in args.contacts:
        prefix = contact_file.rsplit(".", 1)[0]
        contact_df = pd.read_csv(contact_file, sep='\t', index_col=None, 
                                    header=None)
        contact_df.columns = ['contig1', 'contig2', 'count']
        contact_df = contact_df.set_index(['contig1', 'contig2'])
        contact_db = contact_df.to_dict()['count']
        
        contig_paris = [i.strip().split()[:2] for i in open(args.contig_pairs) if i.strip()]

        data = []

        def func(pair):
            try:
                cis1 = contact_db[(pair[0], pair[0])]
                cis2 = contact_db[(pair[1], pair[1])]
                trans = contact_db[(pair[0], pair[1])]
            except KeyError:
                return 0
            
            if cis1 == 0 or cis2 == 0:
                return 0
            
            return trans / sqrt(cis1 * cis2)

        data = []
        for pair in contig_paris:
            if pair[0] > pair[1]:
                continue
            data.append((tuple(pair), func(pair)))

        
      
        # res[prefix] = list(filter(lambda x: ~np.isinf(x),  data))
        data = dict(data)
        data = pd.DataFrame(data.values(), index=pd.MultiIndex.from_tuples(data.keys()))
 
        data.columns = [prefix]
        data.to_csv(f"{prefix}.trans.tsv", index=True, header=None, sep='\t')
        res.append(data)
    
    res_df = pd.concat(res, axis=1)
    print(res_df)
    if len(res_df.columns) == 1:
        res_df = res_df[res_df > 0]
        res_df = res_df[res_df < .2]
    else:
        res_df = res_df[(res_df < .002).all()]
    res_df.to_csv(args.output, index=True, header=None, sep='\t')
    
    
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