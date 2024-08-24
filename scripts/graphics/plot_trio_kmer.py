#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot the trio kmer distribution 
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

import pandas as pd 
import numpy as np 


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('txt', 
            help='results from trioeval')
    pOpt.add_argument("--hapinfo", default=None, 
                      help="haplotype information for color and shape contig in picture")
    pOpt.add_argument('-o', '--output', type=str,
            default="trio.plot.png", help='output file [default: "trio.plot.png"]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    df = pd.read_csv(args.txt, sep=r'\s+', header=None, index_col=None, 
                        comment="C", usecols=range(9), engine="python",
                        names=["type", "seq_name", "pat_kmer", "mat_kmer", "pat_pat", "pat_mat", "mat_pat", "mat_mat",  "seq_len"])

    df = df[df["type"] == "S"]
    df = df.drop("type", axis=1)

    # df.columns = ["seq_name", "pat_kmer", "mat_kmer", "pat_pat", "mat_pat"] #, "mat_mat", "seq_len"]
    df['Total k-mers'] = df['pat_kmer'] + df['mat_kmer']
    df["pat_pat"] = df["pat_pat"] / 1000
    df["mat_mat"] = df["mat_mat"] / 1000
    # df['Total k-mers'] = df['pat_kmer'] + df['mat_kmer'] * 1000
    df['Length (kb)'] = df['seq_len'] / 1000
    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(5, 5))

    sns.scatterplot(x="pat_kmer", y="mat_kmer", data=df, ax=ax, 
                    size='seq_len', sizes=(10, 200), color="#B3CDE3")
#     plt.legend([])
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    sns.despine()
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    plt.tick_params(which='both', width=1.5, length=5)
    # plt.xscale("log")
    # plt.yscale("log")
    plt.xlabel("# paternal k-mers (x$10^3$)", fontsize=16)
    plt.ylabel("# maternal k-mers (x$10^3$)", fontsize=16)

    plt.savefig(args.output, dpi=600, bbox_inches="tight")
    plt.savefig(args.output.replace("png", "pdf"), dpi=600, bbox_inches="tight")
    

    

if __name__ == "__main__":
    main(sys.argv[1:])