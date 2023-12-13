#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot violin plot 
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


from matplotlib.ticker import MaxNLocator


def plot(data, title, output):
    colors = {"C-Phasing": "#b02418",
              "HapHiC": "#cb6e7f",
              "ALLHiC": "#253761",
              "ALLHiC_pregroup": "#8896ae"}

    softwares = ["C-Phasing", "HapHiC", "ALLHiC", "ALLHiC_pregroup"]
    tmp_res = []
    for software in softwares:
        tmp_res.append(data[data["Software"] == software])
    data = pd.concat(tmp_res, axis=0)
    fig, ax = plt.subplots(figsize=(9, 5))
    plt.rcParams['font.family'] = 'Arial'
    ax = sns.violinplot(data=data, x="N50", y="F1_score", hue="Software",
                        palette=colors, order=["500 kb", "1 mb", "2 mb", "5 mb", "10 mb"])

    x_positions = [0.5, 1.5, 2.5, 3.5] 
    for x_position in x_positions:
        ax.plot([x_position, x_position], [-0.2, 1.2], linestyle='dashed', color='black')
    ax.set_xlabel("N50", fontsize=24)
    ax.set_ylabel("F1 score", fontsize=24)
    plt.xticks(fontsize=18)

    plt.yticks(fontsize=18)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

    plt.title(f"{title}", fontsize=24, fontweight='bold')
    plt.legend(title="Category", bbox_to_anchor=(1.05, 1), loc='upper left' )
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

    # output = f"F_v.violinplot.png"
    # title = "$\it{Fragaria\ virginiana}$"
    output = f"F_c.violinplot.png"
    title = "$\it{Fragaria\ chiloensis}$"
    plot(df[["Software", "N50", "F1_score"]], title, output)

    

if __name__ == "__main__":
    main(sys.argv[1:])

