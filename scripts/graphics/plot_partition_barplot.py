#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot the evaluation of partition
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

from matplotlib.collections import LineCollection, PathCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

import pandas as pd 
import numpy as np

from collections import OrderedDict
from intervaltree import Interval, IntervalTree
from pathlib import Path
from pytools import natsorted
from natsort import natsort_keygen

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
    pOpt.add_argument('--legend', default=False, action='store_true',
                        help='add legend')
    pOpt.add_argument('-o', '--output', type=str,
            default=None, help='output file [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    agp_df, _ = import_agp(args.agp)
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

    line_pos_db = OrderedDict()
    line_start_db = []
    idx = 0
    for i, chrom in enumerate(chromsizes, 1):
        if chrom[1] not in line_pos_db:
            line_pos_db[chrom[1]] = idx
            idx += 1
        if chrom[0] not in line_start_db:
            line_start_db.append(chromsizes[chrom] + 1000000)

    
    max_idx = max(line_pos_db.values()) + 1
    line_start_db = np.r_[[0], np.cumsum(line_start_db[::max_idx])]

    color_db = ["#FBB4AE", "#B3CDE3", ]
    colors = []
    lines = []
    error_lines = []
    unanchored_lines = []
    stat = res_df.groupby('correct')['source']
    correct_contigs = res_df[res_df['correct'] == True]['id']
    misassignment_contigs = res_df[res_df['correct'] == False]['id']
    
    missing_contigs = list(set(contigsizes.index) - set(res_df['id']))
    

    for group, tmp_df in res_df.groupby('chrom'):
        try:
            chrom = res_db[group]
        except:
            continue
        
        idx = line_pos_db[chrom[1]]
        init_start = line_start_db[int(chrom[0]) - 1]
        existed_error = IntervalTree()
        for i, row in tmp_df.iterrows():
            
            tmp_row = contig_position.loc[row.id]
            contig_start = tmp_row.start
            contig_end = tmp_row.end

            colors.append(color_db[row.correct])
            if row.correct:
                line = [(contig_start+ init_start, idx), 
                        (contig_end + init_start, idx)]
                lines.append(line)
    
            else:
                err_idx = idx + 0.2
                if existed_error.overlap(contig_start+ init_start, contig_end + init_start):
                    err_idx += 0.2

                
                line = [(contig_start+ init_start, err_idx), 
                        (contig_end + init_start, err_idx)]
                
                existed_error.addi(contig_start+ init_start, contig_end + init_start)

                error_lines.append(line)

    back_lines = []
    
    for chrom in chromsizes:
        idx = line_pos_db[chrom[1]]
        init_start = line_start_db[int(chrom[0]) - 1]
        line = [(init_start, idx), (init_start + chromsizes[chrom], idx)]
        back_lines.append(line)
        
    for contig in missing_contigs:
        chrom = contig.split(".")[0]
        
        idx = line_pos_db[chrom[1]]
        init_start = line_start_db[int(chrom[0]) - 1]
        tmp_row = contig_position.loc[contig]
        contig_start = tmp_row.start
        contig_end = tmp_row.end
        line = [(init_start + contig_start, idx), (init_start + contig_end, idx)]
        unanchored_lines.append(line)

    plt.rcParams['font.family'] = 'Arial'

    fig = plt.figure(figsize=(3, 0.12 * (max_idx + 1)))

    ax = fig.add_subplot(111)
    line_segments = LineCollection(back_lines, linewidth=1.5, color="#bcbcbc", alpha=0.15)
    ax.add_collection(line_segments)
    line_segments = LineCollection(unanchored_lines, linewidth=1.5, color="#bcbcbc")
    ax.add_collection(line_segments)
    line_segments = LineCollection(error_lines, linewidth=1.5, color=color_db[0], alpha=0.8)
    ax.add_collection(line_segments)
    line_segments = LineCollection(lines, linewidth=1.5, color=color_db[1])
    ax.add_collection(line_segments)
    
    plt.ylim(max_idx + 1 + 5)
    ax.autoscale()

    ax.set_xticks(line_start_db[:-1] + (line_start_db[1:] - line_start_db[:-1])// 2)
    ax.set_xticklabels(list(map(lambda x: f"chr{x}", [1, 2, 3, 4, 5])), fontsize=6)

    ax.set_yticks(range(max_idx))
    ax.set_yticklabels(line_pos_db.keys(), fontsize=6)
    ax.tick_params(axis='both', length=0)
    # ax.set_ylabel("Haplotypes")
    sns.despine(trim=True, bottom=True, left=True)

    if args.legend:
        custom_legend_lines = [
                    Line2D([0], [0], color=color_db[1], lw=3),
                    Line2D([0], [0], color=color_db[0], lw=3),
                    Line2D([0], [0], color="#bcbcbc", lw=3),]
        ax.legend(custom_legend_lines, ['Correct assignment', 'Misassignment', 'Unanchored sequences'],
                    loc='lower center', bbox_to_anchor=(0.5, 1.02),
            fancybox=False, shadow=False, ncol=3, fontsize='x-small',
            frameon=False)

    output = Path(args.agp).stem if not args.output else args.output
    plt.savefig(f"{output}.png", dpi=600, bbox_inches='tight')
    plt.savefig(f"{output}.pdf", dpi=300, bbox_inches='tight')


    total_length = contigsizes.sum().values[0]
    total_count = len(contigsizes)

    correct_count = len(correct_contigs)
    misassignment_count = len(misassignment_contigs)
    missing_count = len(missing_contigs)

    correct_length = contigsizes.loc[correct_contigs].sum().values[0] if correct_count else 0
    misassignment_length = contigsizes.loc[misassignment_contigs].sum().values[0] if misassignment_count else 0
    missing_length = contigsizes.loc[missing_contigs].sum().values[0] if missing_count else 0

    E_F_correct_contigs = list(filter(lambda x: x.split(".")[0][1] in {"E", "F"}, correct_contigs))
    E_F_misassignments_contigs = list(filter(lambda x: x.split(".")[0][1] in {"E", "F"}, misassignment_contigs))
    E_F_missing_contigs = list(filter(lambda x: x.split(".")[0][1] in {"E", "F"}, missing_contigs))
    E_F_total_contigs = list(filter(lambda x: x.split(".")[0][1] in {"E", "F"}, contigsizes.index.values))

    E_F_correct_count = len(E_F_correct_contigs)
    E_F_misassignments_count = len(E_F_misassignments_contigs)
    E_F_missing_count = len(E_F_missing_contigs)
    E_F_total_count = len(E_F_total_contigs)

    E_F_correct_length = contigsizes.loc[E_F_correct_contigs].sum().values[0] if E_F_correct_count else 0 
    E_F_misassignments_length = contigsizes.loc[E_F_misassignments_contigs].sum().values[0] if E_F_misassignments_count else 0 
    E_F_missing_length = contigsizes.loc[E_F_missing_contigs].sum().values[0] if E_F_missing_count else 0 
    E_F_total_length = contigsizes.loc[E_F_total_contigs].sum().values[0] 
    

    

    print(output, correct_count, misassignment_count, missing_count, total_count, 
            correct_length, misassignment_length, missing_length, total_length, 
            correct_length/ total_length, misassignment_length / total_length, missing_length / total_length,
            E_F_correct_count, E_F_misassignments_count, E_F_missing_count, E_F_total_count,
            E_F_correct_length, E_F_misassignments_length, E_F_missing_length, E_F_total_length, 
            E_F_correct_length / E_F_total_length, E_F_misassignments_length / E_F_total_length, 
            E_F_missing_length / E_F_total_length, sep='\t')
    



if __name__ == "__main__":
    main(sys.argv[1:])