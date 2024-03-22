#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
estimate the fragment length from in silco digestion genome
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

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bed', 
            help='restriction postiion bed, can generate from `cphasing-rs digest genome.fasta -p GATC -s 0 > GATC.bed')
    pOpt.add_argument('-o', '--output', type=str,
            default="fragment.length.hist.png", help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    bed_df = pd.read_csv(args.bed, sep='\t', index_col=None, header=None, 
                            usecols=range(3), names=["chrom", "start", "end"])
    
    res = []
    for chrom, tmp_df in bed_df.groupby('chrom'):
        tmp_fragment_lengths = tmp_df['start'] - np.r_[0, tmp_df['start']][:-1]

        res.append(tmp_fragment_lengths)
    fragment_lengths = np.concatenate(res)
    fragment_lengths = np.sort(fragment_lengths)
    median_length = np.median(fragment_lengths)
    mean_length = np.mean(fragment_lengths)
    print(f"Median length\t{median_length}\nMean length\t{mean_length}", file=sys.stdout)
    
    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(5, 5))
    sns.histplot(fragment_lengths, ax=ax, color='r', alpha=0.3, stat='density')
    plt.xlim(1, 2000)
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.tick_params(axis='both', labelsize=12)
    plt.xlabel("Fragment length (bp)", fontsize=16)
    plt.ylabel("Density", fontsize=16)
    plt.savefig(args.output, dpi=600, bbox_inches="tight")
    plt.savefig(args.output.replace(".png", ".pdf"), dpi=600, bbox_inches="tight")
    

if __name__ == "__main__":
    main(sys.argv[1:])