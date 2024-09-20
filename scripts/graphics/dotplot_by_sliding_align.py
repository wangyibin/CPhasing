#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

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
import matplotlib.gridspec as gridspec
import seaborn as sns
import patchworklib as pw
from matplotlib.ticker import ScalarFormatter, MaxNLocator

from collections import OrderedDict
from functools import reduce
from pathlib import Path 
from pytools import natsorted


from cphasing.utilities import (chrom_ticks_convert, 
                                run_cmd, list_flatten,
                                is_empty)




PAF_HADER = ["contig1", "length1", "start1", "end1", "strand",
                 "contig2", "length2", "start2", "end2", "mismatch",
                 "matches", "mapq", "identity"]


def align(ref, query, window=100000, step=50000, threads=10):
    cmd = ["seqkit", "sliding", "-W", str(window), "-s", str(step),
           str(query), "-o", f"slding.{window}.{step}.fasta" ]

    run_cmd(cmd)

    cmd = ["wfmash", "-t", str(threads),  "-m", "-N", str(ref), f"slding.{window}.{step}.fasta", ">",
           "sliding.wfmash.out"]
    os.system(" ".join(cmd))

    os.remove(f"slding.{window}.{step}.fasta")

    return "sliding.wfmash.out"



def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('ref', 
            help='')
    pReq.add_argument('query', 
            help='')
    pOpt.add_argument('-o', '--output', type=str,
            default=None, help='output file [default: %(default)s]')
    pOpt.add_argument('--merge',
                      action='store_false',
                      default=True,
                      help="dont merge plot into one picture")
    pOpt.add_argument('-t', '--threads', type=int, default=10,
            help='number of program threads [default:%(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    if not Path("sliding.wfmash.out").exists() or is_empty("sliding.wfmash.out"):
        paf = align(ref=args.ref, query=args.query)
    else:
        paf = "sliding.wfmash.out"

    df = pd.read_csv(paf, sep='\t', header=None, usecols=range(13),
                         names=PAF_HADER, index_col=None)
    

    df['identity'] = df['identity'].map(lambda x: x.replace("id:f:", "")).astype(float)
    
    _df = df.drop_duplicates(subset=['contig1'], keep='first')
    _df = _df[_df['contig1'].str.startswith('Chr')]
    _df = _df[_df['contig1'].str.startswith('group')]
    total_count = len(_df)
    df = df.drop_duplicates(subset=['contig1'], keep=False)
    total_count = len(df)

    raw_contig_names = df['contig1'].str.split(":", n=2, expand=True)
    raw_contig_name = raw_contig_names[0].str.replace("_sliding", "")
    df['raw_chrom'] = raw_contig_name
    start_end = raw_contig_names[1].str.split("-", expand=True)
    df['raw_start'] = start_end[0].astype(int)
    df['raw_end'] = start_end[1].astype(int)

    ref_chroms = df['contig2'].values.tolist()
    ref_chroms = natsorted(set(ref_chroms))
    ref_idx = OrderedDict()
    idx_db = list_flatten(list(zip(range(23), range(23))))
    for i, chrom in enumerate(ref_chroms):

        ref_idx[chrom] = (i % 2, idx_db[i])
    

    colors = ["k", "#B3CDE3", "grey"]

    axes = OrderedDict()

    fig = plt.figure(figsize=(60, 4.2))
    gs = gridspec.GridSpec(2, 24)

  
    total_swith_error_count = 0
    total_misassembly_count = 0
    ploted_chrom = set()
    for raw_chrom, tmp_df in df.groupby('raw_chrom', sort=False):
        if raw_chrom.startswith('chr'):
            continue
        
        if raw_chrom.startswith('utg'):
            continue

        best_chrom = tmp_df.groupby('contig2')['raw_chrom'].count().idxmax()
        if best_chrom in ploted_chrom:
            continue
        tmp_df = tmp_df[['contig1', 'raw_chrom', 'raw_start', 'raw_end', 'strand', 'contig2', 'start2', 'end2']]
        
        trio = best_chrom.split("_")[1]
        if trio == "PATERNAL":
            tmp_df['color'] = colors[0]
        else:
            tmp_df['color'] = colors[1]

        incorrect_df = tmp_df[tmp_df['contig2'] != best_chrom]
        correct_df = tmp_df[tmp_df['contig2'] == best_chrom]
        if not incorrect_df.empty:
            switch_errors = incorrect_df[incorrect_df['contig2'].str.split("_", expand=True)[0] == best_chrom.split("_")[0]]
            
            misassemblies = incorrect_df[incorrect_df['contig2'].str.split("_", expand=True)[0] != best_chrom.split("_")[0]]
            
            if not switch_errors.empty:
                total_swith_error_count += len(switch_errors)
                if trio == "PATERNAL":
                    tmp_df.loc[switch_errors.index, ['color']] = colors[1]
                else:
                    tmp_df.loc[switch_errors.index, ['color']] = colors[0]
            
            if not misassemblies.empty:
                total_misassembly_count += len(misassemblies)
                tmp_df.loc[misassemblies.index, ['color']] = colors[2]

        if args.merge:
            # ax = pw.Brick(figsize=(3, 3))
            
            ax = fig.add_subplot(gs[ref_idx[best_chrom][0], ref_idx[best_chrom][1]])
            plt.subplots_adjust(wspace=0.8, hspace=0.8)
        else:
            fig, ax = plt.subplots(figsize=(3, 3))

        
        line_list = []
        color_list = []
        for i, row in tmp_df.iterrows():
            
            if row.color == 'grey':
                line_list.append([(row.raw_start, row.raw_end),
                                    (0, 0)])

                color_list.append(row.color)
            else:
                if row.strand == "+":
                    line_list.append([(row.raw_start, row.raw_end),
                                    (row.start2, row.end2)])
                else:
                     line_list.append([(row.raw_start, row.raw_end),
                                    (row.end2, row.start2)])
                color_list.append(row.color)
        
        for i, line in enumerate(line_list):
            plt.plot(line[0], line[1], color=color_list[i])
        
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        _xticks = np.r_[ax.get_xticks()[1:-1]]
        ax.set_xticks(_xticks)
        xticklabels = chrom_ticks_convert(_xticks, add_suffix=False)
        ax.set_xticklabels(xticklabels, fontsize=8, rotation=90)

        _yticks = np.r_[ax.get_yticks()[1:-1]]
        ax.set_yticks(_yticks)
        yticklabels = chrom_ticks_convert(_yticks, add_suffix=False)
        ax.set_yticklabels(yticklabels, fontsize=8)
        plt.xlabel(f"Qry: {raw_chrom} (Mb)", fontsize=10)
        plt.ylabel(f"{best_chrom} (Mb)", fontsize=10)

        if not args.merge:
            output = Path(args.paf).stem if not args.output else args.output
            plt.savefig(f"{output}.{raw_chrom}.dotplot.png", dpi=600, bbox_inches='tight')
            plt.savefig(f"{output}.{raw_chrom}.dotplot.pdf", dpi=300, bbox_inches='tight')

        else:
            axes[ref_idx[best_chrom]] = ax 

        ploted_chrom.add(best_chrom)

    if not args.merge:
        return 
    

    print(f"Misassemblies rate: {total_misassembly_count / total_count}")
    print(f"Switch error rate: {total_swith_error_count/ total_count}")

    fig.savefig("merge.png", dpi=600, bbox_inches='tight')
    fig.savefig("merge.pdf", dpi=600, bbox_inches='tight')



if __name__ == "__main__":
    main(sys.argv[1:])