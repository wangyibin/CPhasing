#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
simulate chimeric contigs
    1. inner chrom chimeric: 0.2
    2. inter chrom chimeric: 0.8
        a. homo inter chrom chimeric: 0.8
            homo region switch: 0.7
            nonhomo region switch: 0.3
        b. nonhomo inter chrom chimeric: 0.2
"""

import argparse
import logging
import os
import os.path as op
import sys

import random
import pandas as pd
import numpy as np
import scipy.stats

from collections import defaultdict
from pyfaidx import Fasta
from cphasing.core import AlleleTable
from cphasing.utilities import list_flatten


def random_join():
    pass 



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
    pOpt.add_argument('--chimeric_ratio', default=0.05, type=float,
            help="The ratio of chimeric contigs")
    pOpt.add_argument('--inner_chrom_ratio', default=0.2, type=float, 
            help="The ratio of inner chromosome chimeric contigs")
    pOpt.add_argument('--inter_chrom_ratio', default=0.8, type=float,
            help="The ratio of inter chromosome chimeric contigs")
    pOpt.add_argument('--homo_inter', default=0.8, type=float,
            help="The ratio of homo inter chromosome chimeric contigs")
    pOpt.add_argument('--nonhomo_inter', default=0.2, type=float,
            help="The ratio of nonhomo inter chromosome chimeric contigs")
    pOpt.add_argument('--homo_region_switch', default=0.7, type=float,
            help="The ratio of homo inter chromosome chimeric contigs")
    pOpt.add_argument('--nonhomo_region_switch', default=0.3, type=float,
            help="The ratio of homo inter chromosome chimeric contigs")
    pOpt.add_argument('--seed', type=int, default=1212,
            help="random seed of program")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    random.seed(args.seed)

    chimeric_ratio = args.chimeric_ratio
    assert chimeric_ratio <= 0.5, \
        "The ratio of chimeric is set too high, must less than 0.5"

    inner_chrom_ratio, inter_chrom_ratio = args.inner_chrom_ratio, args.inter_chrom_ratio
    assert inner_chrom_ratio + inter_chrom_ratio == 1.0, \
        "The sum of `--inter_chrom_ratio` and `--inner_chrom_ratio` must be equal to 1.0"
    
    homo_inter, nonhomo_inter = args.homo_inter, args.nonhomo_inter 
    assert homo_inter + nonhomo_inter == 1.0, \
        "The sum of `--homo_inter` and `--nonhomo_inter` must be equal to 1.0"

    homo_region_switch, nonhomo_region_switch = args.homo_region_switch, args.nonhomo_region_switch
    assert homo_region_switch + nonhomo_region_switch == 1.0, \
        "The sum of `--homo_region_switch` and `--nonhomo_region_switch` must be equal to 1.0"


    fasta = Fasta(args.fasta)
    contig_list = list(fasta.keys())
    total_contig_count = len(contig_list)
    sim_conitg_pair_count = int(total_contig_count * chimeric_ratio)
    inner_chrom_contig_count = int(inner_chrom_ratio * sim_conitg_pair_count)
    inter_chrom_contig_count = int(inter_chrom_ratio * sim_conitg_pair_count)

    homo_inter_contig_count = int(homo_inter * inter_chrom_contig_count)
    non_homo_inter_contig_count = int(nonhomo_inter * inter_chrom_contig_count)

    homo_region_switch_contig_count = int(homo_region_switch * homo_inter_contig_count)
    non_homo_region_switch_contig_count = int(nonhomo_region_switch * homo_inter_contig_count)

    at = AlleleTable(args.alleletable, sort=False, fmt='allele2')

    high_similarity_data = at.data[at.data['similarity'] >= 0.95]
    high_similarity_data = high_similarity_data[high_similarity_data[1] < high_similarity_data[2]]
    idx1 = abs(high_similarity_data['mz1'] - high_similarity_data['mz2']) / high_similarity_data['mz1'] < 0.2
    idx2 = abs(high_similarity_data['mz1'] - high_similarity_data['mz2']) / high_similarity_data['mz2'] < 0.2
    high_similarity_data = high_similarity_data[pd.concat([idx1, idx2], axis=1).all(axis=1)]

    idx1 = abs(high_similarity_data['mz1'] - high_similarity_data['mzShared']) / high_similarity_data['mz1'] < 0.2
    idx1 = abs(high_similarity_data['mz2'] - high_similarity_data['mzShared']) / high_similarity_data['mz2'] < 0.2
    high_similarity_data = high_similarity_data[pd.concat([idx1, idx2], axis=1).all(axis=1)]

    print(high_similarity_data, file=sys.stderr)
    high_similarity_contig_pairs = set(map(tuple, high_similarity_data[[1, 2]].values.tolist()))


    assert homo_region_switch_contig_count < len(set(list_flatten(high_similarity_contig_pairs))), \
        "not enough homo pair in allele table to simulate, please decrease the ratio of homo switch error"
    
    homo_region_switch_contig_list = set()
    homo_region_switch_contig_pairs = set()

    while len(homo_region_switch_contig_pairs) < homo_region_switch_contig_count//2:

        pair = random.choice(list(high_similarity_contig_pairs))

        if pair in homo_region_switch_contig_pairs:
            continue
        
        if pair[0] in homo_region_switch_contig_list or pair[1] in homo_region_switch_contig_list:
            continue 

        homo_region_switch_contig_pairs.add(pair)

        homo_region_switch_contig_list.add(pair[0])
        homo_region_switch_contig_list.add(pair[1])


    exists_contig_pairs = homo_region_switch_contig_pairs.copy()
    
    exists_contig_list = homo_region_switch_contig_list.copy()
    inner_chrom_contig_pairs = set()
    non_homo_region_switch_contig_pairs = set() 
    non_homo_inter_contig_pairs = set()
    
    while True:
        contig_list = set(contig_list) - exists_contig_list
        if len(contig_list) <= 2:
            break 
        pair = set(random.sample(contig_list, 2))
        contig1, contig2 = pair
        if contig1 > contig2:
            contig2, contig1 = contig1, contig2 

        pair = (contig1, contig2)

        if pair in exists_contig_pairs:
            continue 

        if contig1 in exists_contig_list or contig2 in exists_contig_list:
            continue 

        if pair in high_similarity_contig_pairs:
            continue 

        chrom1 = contig1.rsplit(".", 1)[0]
        chrom2 = contig2.rsplit(".", 1)[0]

        hap1 = chrom1[0]
        hap2 = chrom2[0]

        if hap1 == hap2:
            if chrom1 == chrom2:
                if len(inner_chrom_contig_pairs) < inner_chrom_contig_count:
                    inner_chrom_contig_pairs.add(pair)
                    exists_contig_pairs.add(pair)
                    exists_contig_list.add(contig1)
                    exists_contig_list.add(contig2)

            else:
                if len(non_homo_region_switch_contig_pairs) < non_homo_region_switch_contig_count:
                    non_homo_region_switch_contig_pairs.add(pair)
                    exists_contig_pairs.add(pair)
                    exists_contig_list.add(contig1)
                    exists_contig_list.add(contig2)

                    
        else:
            if len(non_homo_inter_contig_pairs) < non_homo_inter_contig_count:
                non_homo_inter_contig_pairs.add(pair)
                exists_contig_pairs.add(pair)
                exists_contig_list.add(contig1)
                exists_contig_list.add(contig2)

    
        
        

        if (len(non_homo_inter_contig_pairs) >= non_homo_inter_contig_count) and \
            (len(inner_chrom_contig_pairs) >= inner_chrom_contig_count) and \
            (len(non_homo_region_switch_contig_pairs) >= non_homo_region_switch_contig_count):
            break
        
        
    
    total_chimeric_pairs = list_flatten([inner_chrom_contig_pairs, non_homo_inter_contig_pairs, 
        homo_region_switch_contig_pairs, non_homo_region_switch_contig_pairs])
    
    print(f"Inner chrom: {len(inner_chrom_contig_pairs)}\n"
          f"non homoe inter: {len(non_homo_inter_contig_pairs)}\n" 
          f"homo region switch: {len(homo_region_switch_contig_pairs) * 2}\n"
          f"non homo region switch: {len(non_homo_region_switch_contig_pairs)}",
          file=sys.stderr)
    
    type_suffix_list = ["inner_chrom", "inter_nonhomo_chrom", 
                        "inter_homo_switch", "inter_nonhomo_switch"]
    
    exists_contig_list = set()
    for idx, pairs in enumerate([inner_chrom_contig_pairs, 
                       non_homo_inter_contig_pairs, 
                       homo_region_switch_contig_pairs,
                       non_homo_region_switch_contig_pairs]):
        type_suffix = type_suffix_list[idx]
        order_list = scipy.stats.bernoulli.rvs(0.5, size=len(pairs)*5, random_state=args.seed).tolist()

        # (0) ctg_1 left/right, (1) ctg_2 left/right,
        # (2) ctg_1_half revcom/not, (3) ctg_2_half revcom/not
        # (4) ctg_1_half ctg_2_half ordering 
        order_list = np.array(order_list).reshape((len(pairs), 5)).tolist()
        
        for pair, order in zip(pairs, order_list):   

            contig1, contig2 = pair 
            exists_contig_list.add(contig1)
            exists_contig_list.add(contig2)
            contig1_seq = fasta[contig1]
            contig2_seq = fasta[contig2]
            contig1_len, contig2_len = len(contig1_seq), len(contig2_seq)

            contig1_left = contig1_seq[:contig1_len//2]
            contig1_right = contig1_seq[contig1_len//2:]
            contig2_left = contig2_seq[:contig2_len//2]
            contig2_right = contig2_seq[contig2_len//2:]

            if idx == 2:
                new_conitg1 = str(contig1_left) + str(contig2_right)
                new_contig2 = str(contig2_left) + str(contig1_right) 
                new_id1 = f"{contig1}|{contig2};{contig1_len//2};left+|right+;chimeric_{type_suffix}"
                new_id2 = f"{contig2}|{contig1};{contig2_len//2};left+|right+;chimeric_{type_suffix}"

                print(f">{new_id1}\n{new_conitg1}", file=sys.stdout)
                print(f">{new_id2}\n{new_contig2}", file=sys.stdout)

            else:
                strand1 = "+" if order[2] == 1 else "-"
                strand2 = "+" if order[3] == 1 else "-"
                if order[0] == 0 and order[1] == 0:
                    tmp_contig1 = contig1_left if strand1 == "+" else contig1_left[::-1]
                    tmp_contig2 = contig2_left if strand2 == "+" else contig2_left[::-1]
                    new_contig = str(tmp_contig1) + str(tmp_contig2)
                    new_id = f"{contig1}|{contig2};{contig1_len//2};left{strand1}|left{strand2};chimeric_{type_suffix}"
                    contig_frag1 = contig1_right
                    contig_frag_id1 = f"{contig1};right;frag"
                    contig_frag2 = contig2_right 
                    contig_frag_id2 = f"{contig2};right;frag"
                    
                elif order[0] == 0 and order[1] == 1:
                    tmp_contig1 = contig1_left if strand1 == "+" else contig1_left[::-1]
                    tmp_contig2 = contig2_right if strand2 == "+" else contig2_right[::-1]
                    new_contig = str(tmp_contig1) + str(tmp_contig2)
        
                    new_id = f"{contig1}|{contig2};{contig1_len//2};left{strand1}|right{strand2};chimeric_{type_suffix}"
                    contig_frag1 = contig1_right
                    contig_frag2 = contig2_left 
                    contig_frag_id1 = f"{contig1};right;frag"
                    contig_frag2 = contig2_right 
                    contig_frag_id2 = f"{contig2};left;frag"

                elif order[0] == 1 and order[1] == 0:
                    tmp_contig1 = contig1_right if strand1 == "+" else contig1_right[::-1]
                    tmp_contig2 = contig2_left if strand2 == "+" else contig2_left[::-1]
                    new_contig = str(tmp_contig1) + str(tmp_contig2)
               
                    new_id = f"{contig1}|{contig2};{contig1_len//2};right{strand1}|left{strand2};chimeric_{type_suffix}"
                    contig_frag1 = contig1_left
                    contig_frag_id1 = f"{contig1};left;frag"
                    contig_frag2 = contig2_right 
                    contig_frag_id2 = f"{contig2};right;frag"
                else:
                    tmp_contig1 = contig1_right if strand1 == "+" else contig1_right[::-1]
                    tmp_contig2 = contig2_right if strand2 == "+" else contig2_right[::-1]
                    new_contig = str(tmp_contig1) + str(tmp_contig2)
                    
                    new_id = f"{contig1}|{contig2};{contig1_len//2};right{strand1}|right{strand2};chimeric_{type_suffix}"
                    contig_frag1 = contig1_left
                    contig_frag_id1 = f"{contig1};left;frag"
                    contig_frag2 = contig2_left 
                    contig_frag_id2 = f"{contig2};left;frag"
                
                print(f">{new_id}\n{new_contig}", file=sys.stdout)
                print(f">{contig_frag_id1}\n{contig_frag1}", file=sys.stdout)
                print(f">{contig_frag_id2}\n{contig_frag2}", file=sys.stdout)

    else:
        for contig in fasta.keys():
            if contig in exists_contig_list:
                continue

            print(f">{contig}\n{str(fasta[contig])}", file=sys.stdout)



if __name__ == "__main__":
    main(sys.argv[1:])