#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
simulate collapsed contigs

"""

import argparse
import logging
import os
import os.path as op
import sys

import random
import pandas as pd

from collections import defaultdict
from pyfaidx import Fasta
from cphasing.core import AlleleTable

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('alleletable', 
            help='allele table from `cphasing alleles`')
    pReq.add_argument('fasta', 
            help='raw fasta')
    pOpt.add_argument('--ratio', type=float, default=0.05,
            help="the ratio of collapsed contigs")
    pOpt.add_argument('--seed', type=int, default=1212,
            help="random seed of program")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    random.seed(args.seed)
    ratio = args.ratio 

    fasta = Fasta(args.fasta)
    contig_counts = len(fasta.keys())
    
    at = AlleleTable(args.alleletable, sort=False, fmt='allele2')
    high_similarity_data = at.data[at.data['similarity'] >= 0.99]
    high_similarity_data = high_similarity_data[high_similarity_data[1] < high_similarity_data[2]]
    idx1 = abs(high_similarity_data['mz1'] - high_similarity_data['mz2']) / high_similarity_data['mz1'] < 0.1
    idx2 = abs(high_similarity_data['mz1'] - high_similarity_data['mz2']) / high_similarity_data['mz2'] < 0.1
    high_similarity_data = high_similarity_data[pd.concat([idx1, idx2], axis=1).all(axis=1)]


    high_similarity_contig_pairs = high_similarity_data[[1, 2]].values

    collapsed_contig_count = int(ratio * contig_counts)
    if collapsed_contig_count > len(high_similarity_contig_pairs):
        collapsed_contig_count = len(high_similarity_contig_pairs)
    
    collapsed_contigs = set()

    collapsed_contig_db = defaultdict(list)

    while len(collapsed_contigs) < collapsed_contig_count:
        
        collapsed_contig_pair = random.choice(high_similarity_contig_pairs)
        idx = random.choice([0, 1])
        collapsed_contig = collapsed_contig_pair[idx]
        if len(fasta[collapsed_contig]) < 10000:
            continue
        collapsed_contigs.add(collapsed_contig)
        collapsed_contig_db[collapsed_contig_pair[abs(idx - 1)]].append(collapsed_contig)

    for contig in list(fasta.keys()):
        if contig in collapsed_contigs:
            continue 
        if contig in collapsed_contig_db:
            output_id = "_".join(set(collapsed_contig_db[contig]))
            print(f">{contig}|collapse:{output_id}\n{fasta[contig]}")
        else:
            print(f">{contig}\n{fasta[contig]}")

if __name__ == "__main__":
    main(sys.argv[1:])
