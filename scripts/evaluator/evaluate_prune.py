#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
evaluate the prune 
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 

from itertools import product, combinations

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('contigs', 
            help='contig list')
    pReq.add_argument('prunetable',
            help='prune table from kprune')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    contigs = args.contigs
    pt = args.prunetable
    
    contigs_df = pd.read_csv(contigs, sep='\t', usecols=(0,), header=None, index_col=None)

    chrom_contigs = pd.DataFrame(contigs_df[0].str.split(".").values.tolist(), columns=["chrom", "contig"])
    chrom_contigs['contig'] = contigs_df[0]
    chrom_contigs['hap'] = chrom_contigs['chrom'].str[:-1]
    
    res = []
    for i, df in chrom_contigs.groupby('hap'):
        tmp_list = []
        for j, df2 in df.groupby('chrom'):
            tmp_list.append(df2['contig'].tolist())

        for n, m in list(combinations(tmp_list, 2)):
            res.extend(list(set(product(n, m))))
    
    res2 = [(j, i) for i, j in res]
    res.extend(res2)
    
    inter_pairs = set(res)
    
    pt = pd.read_csv(pt, sep='\t', header=None, usecols=(0, 1), index_col=None)
    pt_pairs = set((map(tuple, pt.values.tolist())))

    correct_pairs = inter_pairs & pt_pairs
    incorrect_pairs = pt_pairs - inter_pairs
    loss_pairs = inter_pairs - pt_pairs
    precision = len(correct_pairs) / len(pt_pairs)
    recall = len(correct_pairs) / len(inter_pairs)
    for pair in loss_pairs:
        print("\t".join(pair), file=sys.stdout)
    print(f"Precision: {precision:.4}", file=sys.stderr)
    print(f"Recall: {recall:.4}", file=sys.stderr)

if __name__ == "__main__":
    main(sys.argv[1:])