#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot the stat of the homo vs nonhomo
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

from itertools import combinations, groupby
from math import sqrt
from statannotations.Annotator import Annotator

from cphasing.utilities import list_flatten

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pOpt.add_argument('-c', '--cool',
            help='Path to cool file',
            required=True)
    pOpt.add_argument("-hl", "--homo_list",
                    
            help="homolog chroms rellationship",
            required=True)
    pOpt.add_argument("--labels", nargs="*", default=None,
                      help="labels of different cool data, the number of labels must equal to cool")
    pOpt.add_argument("--color_index", nargs="*", default=None, 
                      help="color index used for each labels")
    pOpt.add_argument('--wi_test', action="store_true", default=False,
                      help="plot the wi-test, only support two cool file")
    pOpt.add_argument("--rotation", default=0, type=int, 
                      help="rotate the xticks [default: %(default)s]")
    pOpt.add_argument('-o', '--output', type=str,
            default="boxplot.png", help='output picture [default: "boxplot.png"]')
    pOpt.add_argument('-h', '--help', action='help',
        help='show help message and exit.')
    
    args = p.parse_args(args)
    clean_label = []
    for label in args.labels:
        label = label.replace("$", "").replace("\it", "").replace("{", "").replace("}", "")
        clean_label.append(label)

    cool = cooler.Cooler(args.cool)
    matrix = cool.matrix(balance=False, sparse=True)
    homo_chroms = [tuple(sorted(i.strip().split())) for i in open(args.homo_list) if i.strip()]

    chroms = sorted(list_flatten(homo_chroms))
    nonhomo_chroms = []
    homo_contacts = []
    nonhomo_contacts = []

    for chrom_pair in combinations(chroms, 2):
        trans = matrix.fetch(*chrom_pair).sum()
        cis1 = matrix.fetch(chrom_pair[0]).sum()
        cis2 = matrix.fetch(chrom_pair[1]).sum()
        if cis1 == 0 or cis2 == 0:
            continue
        contacts = trans / sqrt(cis1 * cis2)
        if chrom_pair in homo_chroms:
            homo_contacts.append(contacts)
        else:
            nonhomo_chroms.append(chrom_pair)
            nonhomo_contacts.append(contacts)
    

    homo_df = pd.DataFrame(homo_contacts)
    homo_df.columns = [clean_label[0]]
    
    nonhomo_df = pd.DataFrame(nonhomo_contacts)
    nonhomo_df.columns = [clean_label[1]]

    res_df = pd.concat([homo_df, nonhomo_df], axis=0).melt()
    res_df.dropna(axis=0, inplace=True)

    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(3.5, 5))

    colors = ['#df8384', '#8dc0ed', '#a83836',  "#253761",]

    pairs = list(combinations(clean_label, 2))

    if args.color_index:
        new_color = []
        for i in args.color_index:
            new_color.append(colors[int(i)])
    else:
        new_color = colors

    res_df.to_csv(args.output.replace(".png", "") + ".csv", sep='\t', index=True, header=True)
    ax = sns.violinplot(data=res_df, x='variable', y='value', ax=ax, palette=new_color, alpha=1)
    if args.wi_test:
        annotator = Annotator(ax, pairs, data=res_df, x='variable', y='value')
        annotator.configure(test="Mann-Whitney", text_format='full', show_test_name=False, loc='inside')
        annotator.apply_and_annotate()
    
    labels = args.labels.copy() 

    homo_num, nonhomo_num = res_df.groupby('variable')['value'].count().values.tolist()
    labels[0] = f"{labels[0]}\n(n={homo_num})"
    labels[1] = f"{labels[1]}\n(n={nonhomo_num})"
    ax.set_xticks(range(len(args.labels)))
    ax.set_xticklabels(labels, fontsize=8, rotation=args.rotation, ha='center')
    ax.set_xlabel("")
    
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=16)
    ax.set_ylabel("Normalized contacts", fontsize=20)
    plt.savefig(args.output, dpi=600, bbox_inches='tight')
    plt.savefig(args.output.replace("png", "pdf"), dpi=600, bbox_inches='tight')

if __name__ == "__main__":
    main(sys.argv[1:])
