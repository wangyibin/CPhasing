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
import gc

import pickle
import pandas as pd

from itertools import product, combinations, permutations


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('raw', 
            help='raw pairs.txt from allhic extract')
    pReq.add_argument('prune',
            help='pruned pairs.txt from ALLHiC_prune and allhic extract')
    pReq.add_argument('contigs', 
            help='contig list')
    pOpt.add_argument('-c', '--contacts',
                    default=None, help="contacts file")
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    contigs = args.contigs
    contacts = args.contacts 
    contacts_dict = {}
    if contacts:
        with open(contacts) as fp:
            for line in fp:
                line_list = line.strip().split()
                contig_pair = (line_list[0], line_list[1])
                # if int(line_list[2]) < 3:
                #     continue
                contacts_dict[contig_pair] = line_list[2]
                
    contigs_df = pd.read_csv(contigs, sep='\t', usecols=(0,), header=None, index_col=None)
    chrom_contigs = pd.DataFrame(contigs_df[0].str.split(".").values.tolist(), columns=["chrom", "contig"])
    chrom_contigs['contig'] = contigs_df[0]
    chrom_contigs['hap'] = chrom_contigs['chrom'].str[:-1]
    chrom_contigs['hap'] = chrom_contigs['hap'].str.split("_").map(lambda x: x[0])
    
    contig_list = chrom_contigs['contig'].values.tolist()
    contig_list.sort()
    total_pairs = set(list(combinations(contig_list, 2)))

    del contig_list 
    gc.collect()
    # print(len(total_pairs))

    inter_pairs = []
    for i, df in chrom_contigs.groupby('hap'):
        tmp_list = []
        for j, df2 in df.groupby('chrom'):
            tmp_list.append(df2['contig'].tolist())

        tmp_list.sort()  
        for n, m in list(combinations(tmp_list, 2)):
            inter_pairs.extend(list(set(product(n, m))))
    
    
    inter_pairs = set(inter_pairs)
    
    intra_pairs = []
    for i, df in chrom_contigs.groupby('chrom'):
        tmp_list = df['contig'].values.tolist()
        tmp_list.sort()
        intra_pairs.extend(list(combinations(tmp_list, 2)))

    intra_pairs = set(intra_pairs)

    nonhomo_inter_pairs = total_pairs - inter_pairs - intra_pairs

    raw_df = pd.read_csv(args.raw, sep='\t', header=None, index_col=None, comment="#")
    raw_df = raw_df[raw_df[8] == "ok"]

    raw_pairs = set(map(tuple, raw_df[[2, 3]].values.tolist()))
    prune_df = pd.read_csv(args.prune, sep='\t', header=None, index_col=None, comment="#")
    prune_df = prune_df[prune_df[8] == "ok"]
    prune_pairs = set(map(tuple, prune_df[[2, 3]].values.tolist()))
    
    pt_pairs = raw_pairs - prune_pairs
    
    if contacts:
        inter_pairs = inter_pairs.intersection(set(contacts_dict.keys()))
        pt_pairs = pt_pairs.intersection(set(contacts_dict.keys()))
        pt_pairs = pt_pairs - nonhomo_inter_pairs
        
    correct_pairs = inter_pairs & pt_pairs
    incorrect_pairs = pt_pairs - inter_pairs

    precision = 1 - len(incorrect_pairs) / len(pt_pairs)
    recall = len(correct_pairs) / len(inter_pairs)
    f1_score = (2 * precision * recall) / (precision + recall)
    for pair in incorrect_pairs:
        print("\t".join(pair), file=sys.stdout)
    print(f"Precision:{precision:.4}")
    print(f"Recall:{recall:.4}")
    print(f"F1 score:{f1_score:.4}")

if __name__ == "__main__":
    main(sys.argv[1:])