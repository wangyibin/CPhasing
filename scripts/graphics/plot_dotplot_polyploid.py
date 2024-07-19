#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot the dotplot of the simulated polyploids
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

import numpy as np
import pandas as pd

from collections import OrderedDict
from intervaltree import Interval, IntervalTree
from pathlib import Path
from pytools import natsorted
from natsort import natsort_keygen

from cooler.util import binnify
import scipy.stats
from cphasing.agp import import_agp
from cphasing.utilities import read_chrom_sizes


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('agp', 
            help='')
    pReq.add_argument('contigsizes')
    pOpt.add_argument('--title', type=str, default=None)
    pOpt.add_argument('-o', '--output', type=str,
            default=None, help='output file [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    agp_df, _ = import_agp(args.agp)

    binsize = 10000
    contigsizes = read_chrom_sizes(args.contigsizes)
    contigsizes = contigsizes.sort_index(key=natsort_keygen())
    
    chromsizes = contigsizes.copy().reset_index()
    chromsizes['source'] = chromsizes['chrom'].str.split('.', expand=True)[0]

    contig_position = []
    for source, tmp_df in chromsizes.groupby('source'):
        cumsum = np.cumsum(tmp_df['length'])
        start = np.r_[[0], cumsum[:-1]]
        end = cumsum.values
        tmp_df['start'] = start
        tmp_df['end'] = end
        contig_position.append(tmp_df)

    contig_position = pd.concat(contig_position, axis=0)
    contig_position.set_index(['chrom'], inplace=True)
    
    chromsizes = chromsizes.groupby('source').sum()['length']
    chromsizes_df = chromsizes
    chromsizes = chromsizes.to_dict()
    agp_df['source'] = agp_df['id'].str.split(".", expand=True)[0]

    agp_df = agp_df[agp_df.index != agp_df['id']]
    
    res_df = agp_df.groupby(['chrom', 'source'])['tig_end'].sum()
    res_df = res_df[res_df > 0]
    
    res_df = res_df.reset_index().groupby('chrom').apply(lambda x: x.sort_values('tig_end')).reset_index(drop=True)
    res_df = res_df.sort_values('tig_end')

    res_df = res_df.reset_index()[['chrom', 'source']]
  
    res_df.set_index('chrom', inplace=True)
    res_db = res_df.to_dict()['source']

    assign_results = []
    for chrom, tmp_df in agp_df.groupby('chrom'):
        try:
            real_source = res_db[chrom]
        except KeyError:
            continue

        tmp_df['correct'] = tmp_df['source'] == real_source 
        tmp_df['real_source'] = real_source
        
        assign_results.append(tmp_df)
    
    res_df = pd.concat(assign_results, axis=0)


    y_data = []
    chromsizes_df = chromsizes_df.to_frame()
    chromsizes_df['chrom'] = chromsizes_df.index.str[0]

    max_chromsizes = chromsizes_df.groupby('chrom').max()
    for chrom in chromsizes:
        if chrom[0] not in y_data:
            y_data.append(chrom[0])
    
    max_chromsizes_binify = binnify(max_chromsizes['length'], binsize)
    contig_position['start'] = contig_position['start'] // binsize * binsize
    contig_position['end'] = contig_position['end'] // binsize * binsize
    
    ref_edges = max_chromsizes_binify.reset_index().groupby('chrom').apply(lambda x: x.iloc[-1])['index'].values
    ref_position = max_chromsizes_binify.apply(lambda row: row['chrom'] + ":" + str(row['start']) + "-" + str(row['end']), axis=1)

    contig_position['position'] = contig_position.apply(lambda row: row['source'][0] + ":" + str(row['start']) + "-" + str(row['end']), axis=1)
    
    contig_position['query'] = contig_position.apply(lambda row: row['source'][1], axis=1)
    
    contig_position_db = contig_position.to_dict()['position']
    ref_position = ref_position.reset_index().set_index(0)
    x_data = []
    y_data = []
    x_idx = 0
    res_df['index'] = list(range(len(res_df)))
    group_edges = []
    colors = []
    color_db = ["#FBB4AE", "#B3CDE3", ]
    rho_datas = OrderedDict()
  
    for group, tmp_df in res_df.groupby('chrom'):
        tmp_df['position'] = tmp_df['id'].map(contig_position_db.get)
        tmp_x_data = []
        tmp_y_data = []
        for chrom, row in tmp_df.iterrows():
            tmp_data = row['position']
            tmp_chrom, tmp_range = tmp_data.split(":")
            tmp_start, tmp_end = tmp_range.split("-")
            tmp_start, tmp_end = int(tmp_start), int(tmp_end)
            
            tmp_data = np.array(list(range(tmp_start, tmp_end + binsize, binsize)))
            tmp_data = list(map(lambda x: f"{tmp_chrom}:{x[0]}-{x[1]}", list(zip(tmp_data[:-1], tmp_data[1:]))))
            tmp_data = ref_position.reindex(tmp_data)
            tmp_data = tmp_data.dropna()
            tmp_data = tmp_data['index'].values
            
            if row.orientation == "-":
                tmp_data = tmp_data[::-1]

            colors.extend([color_db[row.correct]] * len(tmp_data))
            y_data.extend(list(tmp_data))
            tmp_y_data.extend(list(tmp_data))
            x_data.extend(list(range(x_idx, x_idx + len(tmp_data))))
            tmp_x_data.extend(list(range(x_idx, x_idx + len(tmp_data))))
    

            x_idx += len(tmp_data)
        
        rho_datas[group] = scipy.stats.spearmanr(tmp_x_data, tmp_y_data)[0]

        group_edges.append(x_idx)
        
    plt.rcParams['font.family'] = 'Arial'

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)

    plt.scatter(x_data, y_data, s=0.3, color=color_db[1])
    plt.xticks([])

    ref_half_pos = np.r_[[0], ref_edges]
    ref_half_pos = ref_half_pos[:-1] + (ref_half_pos[1:] - ref_half_pos[:-1])//2

    plt.yticks(ref_half_pos, list(map(lambda x: f"chr{x}", max_chromsizes.index.values.tolist())),
               fontsize=14)
    plt.hlines(ref_edges[:-1], ax.get_xlim()[0], ax.get_xlim()[1], 
                color="#bcbcbc", linewidth=0.5, alpha=0.3)
    plt.vlines(group_edges[:-1], ax.get_ylim()[0], ax.get_ylim()[1], 
                color="#bcbcbc", linewidth=0.5, alpha=0.3)
    plt.xlim(0, max(group_edges))
    plt.ylim(0, max(ref_edges))
    plt.xlabel("Groups", fontsize=16)
    plt.ylabel("Ref", fontsize=16)
    if args.title:
        plt.title(args.title, fontsize=16)
    # ax.spines['bottom'].set_linewidth(1.5)
    # ax.spines['left'].set_linewidth(1.5)
    # ax.spines['top'].set_linewidth(1.5)
    # ax.spines['right'].set_linewidth(1.5)
    ax.tick_params(axis='both', length=5, width=1.5)
    # sns.despine()

    output = Path(args.agp).stem if not args.output else args.output
    plt.savefig(f"{output}.dotplot.png", dpi=600, bbox_inches='tight')
    plt.savefig(f"{output}.dotplot.pdf", dpi=300, bbox_inches='tight')

    rho_df = pd.DataFrame(rho_datas.items(), columns=["group", "rho"]).dropna(axis=0)
    rho_df.to_csv(f"{output}.dotplot.rho.tsv", sep='\t', index=None, header=None)


if __name__ == "__main__":
    main(sys.argv[1:])