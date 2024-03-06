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

from itertools import product, combinations, permutations


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
                contacts_dict[contig_pair] = line_list[2]
                
    contigs_df = pd.read_csv(contigs, sep='\t', usecols=(0,), header=None, index_col=None)
    chrom_contigs = pd.DataFrame(contigs_df[0].str.split(".").values.tolist(), columns=["chrom", "contig"])
    chrom_contigs['contig'] = contigs_df[0]
    chrom_contigs['hap'] = chrom_contigs['chrom'].str[:-1]
    chrom_contigs['hap'] = chrom_contigs['hap'].str.split("_").map(lambda x: x[-1])

    total_pairs = set(map(tuple, permutations(chrom_contigs['contig'].values.tolist(), 2)))


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
    
    intra_pairs = set()
    for i, df in chrom_contigs.groupby('chrom'):
        for pair in list(permutations(df['contig'].values.tolist(), 2)):
            intra_pairs.add(pair)
    
    nonhomo_inter_pairs = total_pairs - inter_pairs - intra_pairs

    data = pickle.load(open(args.full_pkl, 'rb'))
    full_pairs = set(list(data.keys()))
    all_pairs = set(list(contacts_dict.keys()))
    
    pt_pairs = all_pairs - full_pairs

    if contacts:
        inter_pairs = inter_pairs.intersection(set(contacts_dict.keys()))
        pt_pairs = pt_pairs.intersection(set(contacts_dict.keys()))
        pt_pairs = pt_pairs - nonhomo_inter_pairs
        
    correct_pairs = inter_pairs & pt_pairs

 

    incorrect_pairs = inter_pairs - correct_pairs
    precision = len(correct_pairs) / len(pt_pairs)
    recall = len(correct_pairs) / len(inter_pairs)
    f1_score = (2 * precision * recall) / (precision + recall)
    # for pair in pt_pairs:
    #     print("\t".join(pair), file=sys.stdout)
    print(f"Precision:{precision:.4}")
    print(f"Recall:{recall:.4}")
    print(f"F1 score:{f1_score:.4}")

if __name__ == "__main__":
    main(sys.argv[1:])