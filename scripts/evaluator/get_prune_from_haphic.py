#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
get allelic and cross-allelic contig pairs from HapHiC
"""

import argparse
import logging
import os
import os.path as op
import sys

import pickle
import pandas as pd

from itertools import product, combinations


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('full_pkl', 
            help='full pkl in 01.cluster from HapHiC')
    pReq.add_argument('contigs', 
            help='contig list')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    contigs = args.contigs

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
    
    res2 = []
    for i, j in res:
        if i > j:
            res2.append((j, i))
        else:
            res2.append((i, j))
    
    inter_pairs = set(res2)

    data = pickle.load(open(args.full_pkl, 'rb'))
    pt_pairs = set(list(data.keys()))

    correct_pairs = inter_pairs & pt_pairs
    incorrect_pairs = pt_pairs - inter_pairs

    precision = len(correct_pairs) / len(pt_pairs)
    recall = len(correct_pairs) / len(inter_pairs)

    print(f"Precision: {precision:.4}", file=sys.stderr)
    print(f"Recall: {recall:.4}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1:])