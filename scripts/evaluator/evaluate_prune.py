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
import gc
import pandas as pd 

from joblib import Parallel, delayed
from itertools import product, combinations, permutations

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
    pOpt.add_argument('-c', '--contacts',
                    default=None, help="contacts file")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    contigs = args.contigs
    pt = args.prunetable
    contacts = args.contacts 
    contacts_dict = {}
    if contacts:
        with open(contacts) as fp:
            for line in fp:
                line_list = line.strip().split()
                contig_pair = (line_list[0], line_list[1])
                count = float(line_list[2])
                # if count >= 1:
                contacts_dict[contig_pair] = count
                
    
    contigs_df = pd.read_csv(contigs, sep='\t', usecols=(0,), header=None, index_col=None)
   
    chrom_contigs = pd.DataFrame(contigs_df[0].str.split(".").values.tolist(), columns=["chrom", "contig"])
    chrom_contigs['contig'] = contigs_df[0]
    chrom_contigs['hap'] = chrom_contigs['chrom'].str[:-1]
    chrom_contigs['hap'] = chrom_contigs['hap'].str.split("_").map(lambda x: x[0])
    
    contig_list = chrom_contigs['contig'].values.tolist()
    contig_list.sort()
    total_pairs = set(combinations(contig_list, 2))

    del contig_list 
    gc.collect()


    # res = []
    # for i, df in chrom_contigs.groupby('hap'):
    #     tmp_list = []
    #     for j, df2 in df.groupby('chrom'):
    #         tmp_list.append(df2['contig'].tolist())
            
    #     for n, m in list(combinations(tmp_list, 2)):
            
    #         res.extend(list(set(product(n, m))))
    
    # res2 = []
    # for i, j in res:
    #     if i > j:
    #         res2.append((j, i))
    #     else:
    #         res2.append((i, j))
    
    # inter_pairs = set(res2)
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
        tmp_list = sorted(tmp_list)
        intra_pairs.extend(list(combinations(tmp_list, 2)))
   
    intra_pairs = set(intra_pairs)
    
    nonhomo_inter_pairs = total_pairs - inter_pairs - intra_pairs
    
    pt = pd.read_csv(pt, sep='\t', header=None, usecols=(0, 1), index_col=None)
    pt_pairs = set((map(tuple, pt.values.tolist())))
    
    if contacts:
        inter_pairs = inter_pairs.intersection(set(contacts_dict.keys()))
        pt_pairs = pt_pairs.intersection(set(contacts_dict.keys()))
        pt_pairs = pt_pairs - nonhomo_inter_pairs

    correct_pairs = inter_pairs & pt_pairs

    incorrect_pairs = pt_pairs - inter_pairs
    loss_pairs = inter_pairs - pt_pairs
    precision = 1 - len(incorrect_pairs) / len(pt_pairs)
    recall = len(correct_pairs) / len(inter_pairs)
    f1_score = (2 * precision * recall) / (precision + recall)

    for pair in loss_pairs:
        print("\t".join(pair), file=sys.stdout)
    print(f"Precision:{precision:.4}")
    print(f"Recall:{recall:.4}")
    print(f"F1 score:{f1_score:.4}")

if __name__ == "__main__":
    main(sys.argv[1:])