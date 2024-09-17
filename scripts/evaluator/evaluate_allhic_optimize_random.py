#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
evaluate the formulation of allhic optimize 
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
import random
import numpy as np
import pandas as pd 

from collections import defaultdict
from cphasing.core import CountRE
from cphasing.core import PairTable



def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('count_re', 
            help='count RE table of a group')
    pReq.add_argument('pairtable', help="pairtable from allhic extract")
    pOpt.add_argument('-r', '--real', default=None,
                      help='one columns of real contig order')
    pOpt.add_argument('-n', '--number', default=10, type=int, 
                      help='number of random shuffle')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    cr = CountRE(args.count_re, has_header=False)
    length_db = cr.length_db 

    prefix = args.count_re.replace(".txt", "")
    pt = PairTable(args.pairtable)
    count_db = defaultdict( int, pt.data['ObservedLinks'].to_dict(),)
    contigs = cr.contigs 

    if args.real:
        real_contigs = [i.strip().split()[0] for i in open(args.real) if i.strip()]


    random.seed()


    S_values = []
    correct_values = []
    for n in range(args.number):
        random.shuffle(contigs)
    
        midpoint_db = []
        cum_sum = 0 
        for i, contig in enumerate(contigs):
            length = length_db[contig]
            midpoint_db.append(cum_sum + length / 2)
            cum_sum += length 

        S = 0
        for i in range(len(contigs)-1):
            for j in range(i + 1, len(contigs)):
                contig1, contig2 = contigs[i], contigs[j]
                dist = midpoint_db[j] - midpoint_db[i]
                count = count_db[(contig1, contig2)]

                S += count / dist 

        if args.real:
            
            if (contigs == real_contigs) or (contigs[::-1] == real_contigs):
                correct_values.append(S)
            else:
                S_values.append(S)
        else:
            S_values.append(S)


    if args.real:
        midpoint_db = []
        cum_sum = 0 
        for i, contig in enumerate(real_contigs):
            length = length_db[contig]
            midpoint_db.append(cum_sum + length / 2)
            cum_sum += length 

        real_S = 0
        for i in range(len(real_contigs)-1):
            for j in range(i + 1, len(real_contigs)):
                contig1, contig2 = real_contigs[i], real_contigs[j]
                dist = midpoint_db[j] - midpoint_db[i]
                count = count_db[(contig1, contig2)]

                real_S += count / dist 

    else:
        real_S = None 

    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(7, 5))
    # sns.lineplot(range(len(S_values)), S_values, ax=ax, color='#bcbcbc')
    # if args.real:
    #     ax.axhline(real_S, color='r')
    sns.histplot(S_values, alpha=0.7, stat='density')
    if args.real:
        ax.axvline(real_S, color='r')
        if correct_values:
            for v in correct_values:
                ax.axvline(v, color='g')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_ylabel("Density", fontsize=16)
    ax.set_xlabel("Score", fontsize=16)

    plt.savefig(f"{prefix}.evulate_allhic_optimize.png", dpi=600, bbox_inches='tight')
    

if __name__ == "__main__":
    main(sys.argv[1:])