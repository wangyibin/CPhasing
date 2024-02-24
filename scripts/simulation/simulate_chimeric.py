#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
simulate chimeric contigs
    1. intra chrom chimeric: 0.1
    2. inter chrom chimeric: 0.9
        a. homo inter chrom chimeric: 0.8
            homo region switch: 0.8
            nonhomo region switch: 0.2
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
from pathlib import Path
from cphasing.core import AlleleTable
from cphasing.utilities import list_flatten, run_cmd


logger = logging.getLogger(__name__)


def run_minimap2(fasta, threads=10):
    cmd = f'minimap2 -DP -k 19 -w 19 -m 200 -t {threads} {fasta} {fasta} > genome.paf'
    os.system(cmd)


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    # pReq.add_argument('alleletable', 
    #         help='allele table from `cphasing alleles`')
    pReq.add_argument('fasta', 
            help='raw fasta')
    
    pOpt.add_argument('--chimeric_ratio', default=0.05, type=float,
            help="The ratio of chimeric contigs. [default: %(default)s]")
    pOpt.add_argument('--intra_chrom_ratio', default=0.1, type=float, 
            help="The ratio of intra chromosome chimeric contigs. [default: %(default)s]")
    pOpt.add_argument('--inter_chrom_ratio', default=0.9, type=float,
            help="The ratio of inter chromosome chimeric contigs. [default: %(default)s]")
    pOpt.add_argument('--homo_inter', default=0.8, type=float,
            help="The ratio of homo inter chromosome chimeric contigs. [default: %(default)s]")
    pOpt.add_argument('--nonhomo_inter', default=0.2, type=float,
            help="The ratio of nonhomo inter chromosome chimeric contigs. [default: %(default)s]")
    pOpt.add_argument('--homo_region_switch', default=0.8, type=float,
            help="The ratio of homo region switch chimeric contigs. [default: %(default)s]")
    pOpt.add_argument('--nonhomo_region_switch', default=0.2, type=float,
            help="The ratio of nonhomo region switch chimeric contigs. [default: %(default)s]")
    pOpt.add_argument('--seed', type=int, default=1212,
            help="random seed of program")
    pOpt.add_argument('-t', '--threads', type=int, default=8,
            help='number of program threads [default:%(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    

    random.seed(args.seed)

    chimeric_ratio = args.chimeric_ratio
    assert chimeric_ratio <= 0.5, \
        "The ratio of chimeric is set too high, must less than 0.5"

    intra_chrom_ratio, inter_chrom_ratio = args.intra_chrom_ratio, args.inter_chrom_ratio
    assert intra_chrom_ratio + inter_chrom_ratio == 1.0, \
        "The sum of `--inter_chrom_ratio` and `--intra_chrom_ratio` must be equal to 1.0"
    
    homo_inter, nonhomo_inter = args.homo_inter, args.nonhomo_inter 
    assert homo_inter + nonhomo_inter == 1.0, \
        "The sum of `--homo_inter` and `--nonhomo_inter` must be equal to 1.0"

    homo_region_switch, nonhomo_region_switch = args.homo_region_switch, args.nonhomo_region_switch
    assert homo_region_switch + nonhomo_region_switch == 1.0, \
        "The sum of `--homo_region_switch` and `--nonhomo_region_switch` must be equal to 1.0"


    if not Path("genome.paf").exists():
        run_minimap2(args.fasta, args.threads)
    else:
        logger.warn("Use existed `genome.paf`")
    paf = "genome.paf"
    
    paf = pd.read_csv(paf, usecols=range(11), sep='\t', index_col=None, header=None)

    paf[(paf[0] < paf[5]) & (paf[9] >= 2000) & (paf[9] / paf[10] >= 0.90)]
    all_similar_pairs = paf[[0, 5]]
 
    all_similar_pairs.columns = ['contig1', 'contig2']
    all_similar_pairs = all_similar_pairs.values.tolist()
    all_similar_pairs = set(map(tuple, all_similar_pairs))
    

    paf = paf[(paf[0] < paf[5]) & (paf[9] >= 50000) & (paf[9] / paf[10] >= 0.95) \
                & ((paf[9] / paf[10]) < 1)]
    paf = paf[(paf[2] < 100) | ((paf[1] - paf[3]) < 100)]
    paf = paf[(paf[7] < 100) | ((paf[6] - paf[8]) < 100)]

    def func(row):
        if row[2] < 100:
            row[2] = 0
        if row[1] - row[3] < 100:
            row[3] = row[1]
        
        if row[7] < 100:
            row[7] = 0
        if row[6] - row[8]:
            row[8] = row[6]

        return row 
    
    paf = paf.apply(func, axis=1)

    high_similarity_data = paf[[0, 5, 
                         1, 2, 3, 
                         4, 6, 7, 
                         8]]
    # print(high_similarity_data, file=sys.stderr)
    high_similarity_data.columns = ['contig1', 'contig2',
                             'length1', 'start1', 'end1',
                             'strand', 'length2', 'start2',
                             'end2']
    
    fasta = Fasta(args.fasta)
    contig_list = list(fasta.keys())
    total_contig_count = len(contig_list)
    sim_conitg_pair_count = int(total_contig_count * chimeric_ratio)
    intra_chrom_contig_count = int(intra_chrom_ratio * sim_conitg_pair_count)
    inter_chrom_contig_count = int(inter_chrom_ratio * sim_conitg_pair_count)

    homo_inter_contig_count = int(homo_inter * inter_chrom_contig_count)
    non_homo_inter_contig_count = int(nonhomo_inter * inter_chrom_contig_count)

    homo_region_switch_contig_count = int(homo_region_switch * homo_inter_contig_count)
    non_homo_region_switch_contig_count = int(nonhomo_region_switch * homo_inter_contig_count)

    
    # print(high_similarity_data, file=sys.stderr)
    high_similarity_contig_pairs = set(map(tuple, high_similarity_data[['contig1', 'contig2']].values.tolist()))
    high_similarity_data = high_similarity_data.set_index(['contig1', 'contig2'])

    assert homo_region_switch_contig_count < len(set(list_flatten(high_similarity_contig_pairs))), \
        "not enough homo pair in allele table to simulate, please decrease the ratio of homo switch error"
    
    homo_region_switch_contig_list = set()
    homo_region_switch_contig_pairs = set()

    while len(homo_region_switch_contig_pairs) < homo_region_switch_contig_count//2:

        pair = random.choice(list(high_similarity_contig_pairs))
        contig1, contig2 = pair

        if len(fasta[contig1]) < 25000 or len(fasta[contig2]) < 25000:
            continue
        if pair in homo_region_switch_contig_pairs:
            continue
        
        if pair[0] in homo_region_switch_contig_list or pair[1] in homo_region_switch_contig_list:
            continue 
        
        # tmp_row = high_similarity_data.loc[pair]
        # if not ((tmp_row['start1'] == 0 and tmp_row['start2'] == 0 ) \
        #     or (tmp_row['start1'] !=0 and tmp_row['start2'] != 0)):
        #     continue
        homo_region_switch_contig_pairs.add(pair)

        homo_region_switch_contig_list.add(pair[0])
        homo_region_switch_contig_list.add(pair[1])


    exists_contig_pairs = homo_region_switch_contig_pairs.copy()
    
    exists_contig_list = homo_region_switch_contig_list.copy()
    intra_chrom_contig_pairs = set()
    non_homo_region_switch_contig_pairs = set() 
    non_homo_inter_contig_pairs = set()
    
    while True:
        contig_list = set(contig_list) - exists_contig_list
        if len(contig_list) <= 2:
            break 
        pair = set(random.sample(contig_list, 2))
        contig1, contig2 = pair

        if len(fasta[contig1]) < 25000 or len(fasta[contig2]) < 25000:
            continue

        if contig1 > contig2:
            contig2, contig1 = contig1, contig2 

        pair = (contig1, contig2)

        if pair in exists_contig_pairs:
            continue 

        if contig1 in exists_contig_list or contig2 in exists_contig_list:
            continue 

        if pair in high_similarity_contig_pairs:
            continue 

        
        chrom1, contig1_suffix = contig1.rsplit(".", 1)
        chrom2, contig2_suffix = contig2.rsplit(".", 1)
        
        hap1 = chrom1[0]
        hap2 = chrom2[0]

        if hap1 == hap2:
            if chrom1 == chrom2:
                if abs(int(contig1_suffix[3:]) - int(contig2_suffix[3:])) < 2:
                    continue
                if len(intra_chrom_contig_pairs) < intra_chrom_contig_count:
                    intra_chrom_contig_pairs.add(pair)
                    exists_contig_pairs.add(pair)
                    exists_contig_list.add(contig1)
                    exists_contig_list.add(contig2)

            else:
                if len(non_homo_region_switch_contig_pairs) < non_homo_region_switch_contig_count \
                    and pair not in all_similar_pairs:
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
            (len(intra_chrom_contig_pairs) >= intra_chrom_contig_count) and \
            (len(non_homo_region_switch_contig_pairs) >= non_homo_region_switch_contig_count):
            break
        
        
    
    total_chimeric_pairs = list_flatten([intra_chrom_contig_pairs, 
                                         non_homo_inter_contig_pairs, 
                                            homo_region_switch_contig_pairs, 
                                            non_homo_region_switch_contig_pairs])
    
    logger.info(f"Intra chrom: {len(intra_chrom_contig_pairs)}\n"
          f"Inter chrom:\n"
          f"├──nonhomo inter chrom: {len(non_homo_inter_contig_pairs)}\n" 
          f"└──homo inter chrom:\n"
          f"   ├──homo region switch: {len(homo_region_switch_contig_pairs) * 2}\n"
          f"   └──nonhomo region switch: {len(non_homo_region_switch_contig_pairs)}",
    )
    
    type_suffix_list = ["intra_chrom", "inter_nonhomo_chrom", 
                        "inter_homo_switch", "inter_nonhomo_switch"]
    
    exists_contig_list = set()
    for idx, pairs in enumerate([intra_chrom_contig_pairs, 
                       non_homo_inter_contig_pairs, 
                       homo_region_switch_contig_pairs,
                       non_homo_region_switch_contig_pairs]):
        type_suffix = type_suffix_list[idx]
        order_list = scipy.stats.bernoulli.rvs(0.5, size=len(pairs)*5, 
                                                    random_state=args.seed).tolist()

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

            if idx == 2:
                tmp_row = high_similarity_data.loc[pair]
                contig1_range = (tmp_row['start1'], tmp_row['end1'])
                contig2_range = (tmp_row['start2'], tmp_row['end2'])

                contig1_pos = (contig1_range[1] - contig1_range[0])//2
                contig2_pos = (contig2_range[1] - contig2_range[0])//2
           
                contig1_seq_frag1 = contig1_seq[contig1_range[0]: contig1_range[1]]
                contig1_left = contig1_seq_frag1[:contig1_pos]
                contig1_right = contig1_seq_frag1[contig1_pos:]
                if contig1_range[0] == 0:
                    contig1_seq_frag2 = contig1_seq[contig1_range[1]:]
                    contig1_frag_start = contig1_range[1]
                    contig1_frag_end = contig1_len 
                else:
                    contig1_seq_frag2 = contig1_seq[:contig1_range[0]]
                    contig1_frag_start = 0
                    contig1_frag_end = contig1_range[0]

                contig2_seq_frag1 = contig2_seq[contig2_range[0]: contig2_range[1]]
                contig2_left = contig2_seq_frag1[:contig2_pos]
                contig2_right = contig2_seq_frag1[contig2_pos:]
                if contig2_range[0] == 0:
                    contig2_seq_frag2 = contig2_seq[contig2_range[1]:]
                    contig2_frag_start = contig2_range[1]
                    contig2_frag_end = contig2_len 
                else:
                    contig2_seq_frag2 = contig2_seq[:contig2_range[0]]
                    contig2_frag_start = 0
                    contig2_frag_end = contig2_range[0]

                new_conitg1 = str(contig1_left) + str(contig2_right)
                new_contig2 = str(contig2_left) + str(contig1_right) 
                new_id1 = f"{contig1}|{contig2};{contig1_pos};L+|R+;chimeric_{type_suffix}"
                new_id2 = f"{contig2}|{contig1};{contig2_pos};L+|R+;chimeric_{type_suffix}"
                
                print(f">{new_id1}\n{new_conitg1}", file=sys.stdout)
                print(f">{new_id2}\n{new_contig2}", file=sys.stdout)
                if len(contig1_seq_frag2) > 0:
                    print(f">{contig1};L;frag\n{contig1_seq_frag2};{contig1_frag_start}-{contig1_frag_end}", 
                            file=sys.stdout)
                if len(contig2_seq_frag2) > 0:
                    print(f">{contig2};L;frag\n{contig2_seq_frag2};{contig2_frag_start}-{contig2_frag_end}", 
                          file=sys.stdout)
            else:
                contig1_left = contig1_seq[:contig1_len//2]
                contig1_right = contig1_seq[contig1_len//2:]
                contig2_left = contig2_seq[:contig2_len//2]
                contig2_right = contig2_seq[contig2_len//2:]
                strand1 = "+" if order[2] == 1 else "-"
                strand2 = "+" if order[3] == 1 else "-"
                if strand1 == "-" and strand2 == "-":
                    strand1 = "+"
                if order[0] == 0 and order[1] == 0:
                    tmp_contig1 = contig1_left if strand1 == "+" else contig1_left.complement
                    tmp_contig2 = contig2_left if strand2 == "+" else contig2_left.complement
                    new_contig = str(tmp_contig1) + str(tmp_contig2)
                    new_id = f"{contig1}|{contig2};{contig1_len//2};L{strand1}|L{strand2};chimeric_{type_suffix}"
                    contig_frag1 = str(contig1_right)
                    contig_frag_id1 = f"{contig1};R;frag;{contig1_len//2}-{contig1_len}"
                    contig_frag2 = str(contig2_right)
                    contig_frag_id2 = f"{contig2};R;frag;{contig2_len//2}-{contig2_len}"
                    
                elif order[0] == 0 and order[1] == 1:
                    tmp_contig1 = contig1_left if strand1 == "+" else contig1_left.complement
                    tmp_contig2 = contig2_right if strand2 == "+" else contig2_right.complement
                    new_contig = str(tmp_contig1) + str(tmp_contig2)
        
                    new_id = f"{contig1}|{contig2};{contig1_len//2};L{strand1}|R{strand2};chimeric_{type_suffix}"
                    contig_frag1 = str(contig1_right)
                    contig_frag_id1 = f"{contig1};R;frag;{contig1_len//2}-{contig1_len}"
                    contig_frag2 = str(contig2_left)
                    contig_frag_id2 = f"{contig2};L;frag;0-{contig2_len//2}"

                elif order[0] == 1 and order[1] == 0:
                    tmp_contig1 = contig1_right if strand1 == "+" else contig1_right.complement
                    tmp_contig2 = contig2_left if strand2 == "+" else contig2_left.complement
                    new_contig = str(tmp_contig1) + str(tmp_contig2)
               
                    new_id = f"{contig1}|{contig2};{contig1_len//2};R{strand1}|L{strand2};chimeric_{type_suffix}"
                    contig_frag1 = str(contig1_left)
                    contig_frag_id1 = f"{contig1};L;frag;0-{contig1_len//2}"
                    contig_frag2 = str(contig2_right) 
                    contig_frag_id2 = f"{contig2};R;frag;{contig2_len//2}-{contig2_len}"
                else:
                    tmp_contig1 = contig1_right if strand1 == "+" else contig1_right.complement
                    tmp_contig2 = contig2_right if strand2 == "+" else contig2_right.complement
                    new_contig = str(tmp_contig1) + str(tmp_contig2)
                    
                    new_id = f"{contig1}|{contig2};{contig1_len//2};R{strand1}|R{strand2};chimeric_{type_suffix}"
                    contig_frag1 = str(contig1_left)
                    contig_frag_id1 = f"{contig1};L;frag;0-{contig1_len//2}"
                    contig_frag2 = str(contig2_left)
                    contig_frag_id2 = f"{contig2};L;frag;0-{contig2_len//2}"
                
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