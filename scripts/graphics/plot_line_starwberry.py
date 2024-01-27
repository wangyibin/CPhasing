#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot line plot 
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


def plot(data, title, x, y, hue, output):
    colors = {"C-Phasing": "#b02418",
              "HapHiC": "#cb6e7f",
              "HapHiC_remove": "#209093",
              "ALLHiC": "#253761",
              "ALLHiC_pregroup": "#8896ae"}
    maker = {"C-Phasing": "s",
              "HapHiC": "^",
              "ALLHiC": "o",
              "ALLHiC_pregroup": "D"}
    makers = ["s", "^", "o", "D"]
    fig, ax = plt.subplots(figsize=(5.5,5))

    if y == "Wall time (s)":
        data[y] = np.array(data[y])

    plt.rcParams['font.family'] = 'Arial'
    sns.pointplot(data=data, x="N50", y=y, ax=ax, hue=hue, palette=colors)
    plt.setp(ax.collections, alpha=.5) 
    plt.setp(ax.lines, alpha=.3)
    ax.set_xlabel(x, fontsize=24)
    ax.set_ylabel(y, fontsize=24)
    plt.xticks(fontsize=18, rotation=45)

    if y == "Wall time (s)":
        ax.set_yscale("log")

    if np.max(data[y]) > 1000 and y != "Wall time (s)":
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        formatter = plt.gca().get_yaxis().get_major_formatter()
        plt.gca().yaxis.set_major_formatter(formatter)
        plt.gca().yaxis.get_offset_text().set_fontsize(14)

    plt.xlabel("")
    plt.yticks(fontsize=18)
    # plt.title(f"{title}", fontsize=24, fontweight='bold')
    # plt.legend(title="Category", bbox_to_anchor=(1.05, 1), loc='upper left' )
    plt.legend().remove()

    if y == "Group numbers":
        chromosome_number = int(title)*5
        plt.axhline(y=chromosome_number, color='black', linestyle="--")

    plt.savefig(output, dpi=600, bbox_inches='tight')
    plt.savefig(output.replace("png", "pdf"), dpi=600, bbox_inches='tight')


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('tsv', 
            help='evalutation results on simulation')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    df = pd.read_csv(args.tsv, sep='\t', header=0, index_col=None)

    for sample, tmp_df in df.groupby('Sample'):
        for column in tmp_df.columns:
            if column == "Software" or column == "Sample" or \
                column == "N50" or column == "Chromosome numbers":
                continue

            output_middle = column.replace(" ", "_").replace("(", "").replace(")", "")
            sample = sample.replace(" ", "_")
            plot(tmp_df[["Software", "N50", column]], sample, "N50", column, "Software", f"{sample}.{output_middle}.png")

        

if __name__ == "__main__":
    main(sys.argv[1:])