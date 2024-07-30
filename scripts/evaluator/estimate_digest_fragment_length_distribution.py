#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
estimate the digest fragment length distribution
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
import matplotlib.lines as mlines

from matplotlib.ticker import MaxNLocator
import colormaps as cmaps

bluered_12 = list(map(lambda x: mpl.colors.rgb2hex(x.colors),  list(cmaps.bluered_12)))


import pandas as pd
import numpy as np

from pathlib import Path


from cphasing.utilities import run_cmd

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('fasta', 
            help='')
    pOpt.add_argument('-p', '--pattern', nargs="+",
                    help="Motif of restriction sites, multiple site use comma seperate", 
                    required=True)
    pOpt.add_argument('-o', '--output', type=str,
            default=None, help='output file [default: fasta_prefix.png]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    if args.output:
        prefix = args.output
    else:
        prefix = Path(args.fasta).stem
    

    new_motifs = []
    for motifs in args.pattern:
        tmp_motif = []
        motifs = motifs.split(",")
        for motif in motifs:
            if "N" in motif:
                for s in "AGCT":
                    tmp_motif.append(motif.replace("N", s))
            else:
                tmp_motif.append(motif)

        new_motifs.append(tmp_motif)

    beds = []
    for i, motif in enumerate(new_motifs):
        cmd = ["cphasing-rs", "digest", args.fasta, "-p", ",".join(motif), "-o", f"{i}.motif.bed"]
        run_cmd(cmd)
    
        beds.append(f"{i}.motif.bed")
    
    
    res = []
    stat_res = []
    for i, bed in enumerate(beds):
        df = pd.read_csv(bed, sep='\t', index_col=None, header=None)
        tmp = []
        for chrom, tmp_df in df.groupby(0):
            data = np.r_[[0], tmp_df[1].values]
            data = data[1:] - data[:-1] 
            tmp.extend(data.tolist())

        stat_res.append((args.pattern[i], np.median(tmp), np.mean(tmp)))
        res.append(np.array(tmp))
        os.remove(bed)
    
    stat_res = pd.DataFrame(stat_res, columns=["Site", "Median_Length_(bp)", "Mean_Length_(bp)"])
    stat_res.to_csv(f"{args.output}.RE.stat", sep='\t', index=None, header=True)
    fig, ax = plt.subplots(figsize=(5.7, 5))

    plt.rcParams['font.family'] = 'Arial'
    
    for i in range(len(res)):
        data = res[i]
        data = data[(data < 5000) & (data > 0)]
        # sns.histplot(data.tolist(), color=bluered_12[i], ax=ax, alpha=0.3, stat='density', #linewidth=0,
        #               bins=200)
        sns.kdeplot(data.tolist(), color=bluered_12[i], ax=ax, alpha=0.8, linewidth=2)
    
    legend_items = [mlines.Line2D([], [], color=color, alpha=0.8, marker=None, linewidth=2,
                               markersize=10, label=f'{sample}') for sample, (color, marker) in list(zip(args.pattern, zip(bluered_12, range(len(args.pattern)))))]
    second_legend = ax.legend(handles=legend_items, fontsize=10, loc='best')
    ax.add_artist(second_legend)
    
    ax.set_xlabel("Fragment length (bp)", fontsize=20)
    ax.set_ylabel("Density", fontsize=20)
    
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    plt.xticks(fontsize=18, rotation=45, ha='right')

    plt.yticks(fontsize=18)

    plt.savefig(f'{args.output}.kde.plot.png', dpi=600, bbox_inches='tight')
    plt.savefig(f'{args.output}.kde.plot.pdf', dpi=600, bbox_inches='tight')
    
    plt.cla()
    # for i in range(len(res)):
    #     data = res[i]
    #     data = data[(data < 5000) & (data > 0)]
    #     sns.histplot(data.tolist(), color=bluered_12[i], ax=ax, alpha=0.3, stat='density', #linewidth=0,
    #                   bins=200)
        
    
    # legend_items = [mlines.Line2D([], [], color=color, alpha=0.3, marker=None, linewidth=5,
    #                            markersize=10, label=f'{sample}') for sample, (color, marker) in list(zip(args.pattern, zip(bluered_12, range(len(args.pattern)))))]
    # second_legend = ax.legend(handles=legend_items, fontsize=10, loc='best')
    # ax.add_artist(second_legend)
    
    # ax.set_xlabel("Fragment length (bp)", fontsize=20)
    # ax.set_ylabel("Density", fontsize=20)
    
    # ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    # plt.xticks(fontsize=18, rotation=45, ha='right')

    # plt.yticks(fontsize=18)

    # plt.savefig(f'hist.plot.png', dpi=600, bbox_inches='tight')
    # plt.savefig(f'hist.plot.pdf', dpi=600, bbox_inches='tight')


if __name__ == "__main__":
    main(sys.argv[1:])