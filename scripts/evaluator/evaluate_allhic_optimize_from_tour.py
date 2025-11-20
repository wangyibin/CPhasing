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

from matplotlib.lines import Line2D

from collections import defaultdict
from pathlib import Path
from cphasing.core import CountRE
from cphasing.core import PairTable, Clm



def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('tour', 
                      help='tour file from allhic optimize')
    pReq.add_argument('count_re', 
                      help='count RE table of a group')
    pReq.add_argument('clm',
                        help='clm file from allhic optimize')
    pReq.add_argument('pairtable', help="pairtable from allhic extract")
    pOpt.add_argument('-r', '--real', default=None,
                      help='one columns of real contig order')
    pOpt.add_argument('-n', '--number', default=10, type=int, 
                      help='number of random shuffle')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    tour_records = []
    idx = 0
    with open(args.tour) as fp:
        for line in fp:
            if line.startswith(">"):
                # if line[1:4] != "GA1" and line[1:4] != "INT":
                #     next_line = next(fp)
                #     continue 
    
                idx += 1 
                next_line = next(fp)
                contigs = list(map(lambda x: x[:-1], next_line.split()))
                tour_records.append(contigs)


    cr = CountRE(args.count_re, has_header=False, minRE=1)
    length_db = cr.length_db 

    prefix = Path(args.count_re).stem.replace(".txt", "")
    pt = PairTable(args.pairtable)
    count_db = defaultdict( int, pt.data['ObservedLinks'].to_dict(),)
    clm = Clm(args.clm)
    count_db = clm.count_db()
    
    
    ## filter count_db to only keep those in length_db
    for key in list(count_db.keys()):
        contig1, contig2 = key
        if (contig1 not in length_db) or (contig2 not in length_db):
            del count_db[key]
    # for key in list(count_db1.keys()):
    #     contig1, contig2 = key
    #     if (contig1 not in length_db) or (contig2 not in length_db):
    #         del count_db1[key]
    

    print(f"Summary of contacts: {sum(count_db.values())} ")
    cr_contigs = cr.contigs 

    if args.real:
        real_contigs = [i.strip().split("\t")[0] for i in open(args.real) if i.strip()]


    random.seed()

    S_values = []
    correct_values = []
    for contigs in tour_records:
        midpoint_db = []
        cum_sum = 0 
        for i, contig in enumerate(contigs):
            length = length_db[contig]
            midpoint_db.append(cum_sum + length/2)
            cum_sum += length 

        S = 0
        for i in range(len(contigs)):
            contig1 = contigs[i]
            for j in range(i + 1, len(contigs)):
                contig2 = contigs[j]
                try:
                    dist = midpoint_db[j] - midpoint_db[i]
                except IndexError:
                    continue
                try:
                    count = count_db[(contig1, contig2)]
                except KeyError:
                    count = 0

                S += count * np.log(dist)
                # Li = float(length_db[contig1])
                # Lj = float(length_db[contig2])
                # s0 = max(1e3, np.median(list(length_db.values())) / 10.0) 
                # d_mid = dist
                # c_len = 0.25 * (Li + Lj)  # ((Li+Lj)/4)
                # d_eff = np.sqrt(d_mid * d_mid + c_len * c_len + s0 * s0)

                # w_len = (min(Li, Lj) / max(1.0, 0.5 * (Li + Lj)))  # [0,1]
                # w_cnt = np.log1p(count)  

                # S += w_cnt * w_len * np.log1p(d_eff / s0)
                
        S_values.append(S)

    ## filter real_contigs 
    real_contigs = list(filter(lambda x: x in contigs, real_contigs))
    if args.real:
        midpoint_db = []
        cum_sum = 0 

        for i, contig in enumerate(real_contigs):
            try:
                length = length_db[contig]
                
            except:
                continue
            midpoint_db.append(cum_sum + length/2)
            cum_sum += length 

        real_S = 0
        for i in range(len(real_contigs)):
            for j in range(i + 1, len(real_contigs)):
                try:
                    dist = midpoint_db[j] - midpoint_db[i]
                except:
                    continue
                contig1, contig2 = real_contigs[i], real_contigs[j]
                try:
                    count = count_db[(contig1, contig2)]
                except KeyError:
                    count = 0

                real_S += count * np.log(dist)
               
                # Li = float(length_db[contig1])
                # Lj = float(length_db[contig2])
                # s0 = max(1e3, np.median(list(length_db.values())) / 10.0) 
                # d_mid = dist
                # c_len = 0.25 * (Li + Lj)  # ((Li+Lj)/4)
                # d_eff = np.sqrt(d_mid * d_mid + c_len * c_len + s0 * s0)

                # w_len = (min(Li, Lj) / max(1.0, 0.5 * (Li + Lj)))  # [0,1]
                # w_cnt = np.log1p(count)  
                # real_S += w_cnt * w_len * np.log1p(d_eff / s0)

    else:
        real_S = None 

    print(len(contigs), len(real_contigs))

    print(S_values)
    print(real_S)
    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(5.5, 5))
    sns.lineplot(data=(S_values), ax=ax, color='grey', alpha=0.3, linewidth=3)
    if args.real:
        ax.axhline(real_S, color='r', linestyle='--')
    
    legend_elements = [Line2D([0], [0], color='grey', lw=3, label="allhic optimize", ),
                       Line2D([0], [0], color='red', linestyle='--', lw=2, label="real" )]


    plt.legend(handles=legend_elements, loc="best")
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    ax.set_ylabel("Score", fontsize=16)
    ax.set_xlabel("Generation-500", fontsize=16)
    
    plt.savefig(f"{prefix}.evulate_allhic_optimize.png", dpi=600, bbox_inches='tight')
    
    


if __name__ == "__main__":
    main(sys.argv[1:])