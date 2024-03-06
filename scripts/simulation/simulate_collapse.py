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
from pathlib import Path
from cphasing.core import AlleleTable
from cphasing.utilities import list_flatten, run_cmd


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
    pReq.add_argument('fasta', 
            help='raw fasta')
    pOpt.add_argument('--ratio', type=float, default=0.05,
            help="the ratio of collapsed contigs. [default: %(default)s]")
    pOpt.add_argument('-t', '--threads', type=int, default=8,
            help='number of program threads [default:%(default)s]')
    pOpt.add_argument('--seed', type=int, default=1212,
            help="random seed of program. [default: %(default)s]")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    random.seed(args.seed)
    ratio = args.ratio 

    fasta = Fasta(args.fasta)
    contig_counts = len(fasta.keys())
    
    if not Path("genome.paf").exists():
        run_minimap2(args.fasta, args.threads)

    paf = "genome.paf"
    
    paf = pd.read_csv(paf, usecols=range(11), sep='\t', index_col=None, header=None)

    paf = paf[(paf[0] < paf[5]) & (paf[9] >= 10000) & (paf[9] / paf[10] >= 0.95)]
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
  
    high_similarity_data.columns = ['contig1', 'contig2',
                             'length1', 'start1', 'end1',
                             'strand', 'length2', 'start2',
                             'end2']
    
    high_similarity_contig_pairs = high_similarity_data[["contig1", "contig2"]].values
    high_similarity_data = high_similarity_data.set_index(['contig1', 'contig2'])

    collapsed_contig_count = int(ratio * contig_counts)
    if collapsed_contig_count > len(high_similarity_contig_pairs):
        collapsed_contig_count = len(high_similarity_contig_pairs)
    
    collapsed_contigs = set()

    collapsed_contig_db = defaultdict(list)

    exists_contig_list = set()
    exists_contig_pairs = set()
    while len(collapsed_contigs) < collapsed_contig_count:
        
        collapsed_contig_pair = random.choice(high_similarity_contig_pairs)
        idx = random.choice([0, 1])
        
        if tuple(collapsed_contig_pair) in exists_contig_pairs:
            continue 
        
        collapsed_contig = collapsed_contig_pair[idx]
        contig = collapsed_contig_pair[abs(idx-1)]

        if len(fasta[collapsed_contig]) < 5000:
            continue

        if contig in exists_contig_list or collapsed_contig in exists_contig_list:
            continue
        
        collapsed_contig_seq = fasta[collapsed_contig]
        contig_seq = fasta[contig]
        collapsed_contig_seq_len, contig_len = len(collapsed_contig_seq), len(contig_seq)

        
        exists_contig_pairs.add(tuple(collapsed_contig_pair))
        exists_contig_list.add(collapsed_contig)
        exists_contig_list.add(contig)

        tmp_row = high_similarity_data.loc[tuple(collapsed_contig_pair)]
        if tmp_row.empty:
            continue
        
        try:
             collapsed_contig_range = (tmp_row['start1'].values[0], tmp_row['end1'].values[0])

        except AttributeError:
            collapsed_contig_range = (tmp_row['start1'], tmp_row['end1'])

        try:
            contig_range = (tmp_row['start2'].values[0], tmp_row['end2'].values[0])
        except AttributeError:
            contig_range = (tmp_row['start2'], tmp_row['end2'])

        if idx == 1:
            contig_range, collapsed_contig_range = collapsed_contig_range, contig_range 

      
        collapsed_contig_seq_frag1 = collapsed_contig_seq[collapsed_contig_range[0]: collapsed_contig_range[1]]
       
        if collapsed_contig_range[0] == 0:
            collapsed_contig_seq_frag2 = collapsed_contig_seq[collapsed_contig_range[1]:]
            collapsed_contig_frag_start = collapsed_contig_range[1]
            collapsed_contig_frag_end = collapsed_contig_seq_len 
        else:
            collapsed_contig_seq_frag2 = collapsed_contig_seq[:collapsed_contig_range[0]]
            collapsed_contig_frag_start = 0
            collapsed_contig_frag_end = collapsed_contig_range[0]

        collapsed_contig_frag_id = \
            f"{collapsed_contig}:{collapsed_contig_frag_start}-{collapsed_contig_frag_end};frag"

        contig_seq_frag1 = contig_seq[contig_range[0]: contig_range[1]]
        contig_seq_frag1_id = (f"{contig}:{contig_range[0]}-{contig_range[1]};"
                                f"collapsed|{collapsed_contig}:{collapsed_contig_range[0]}-{collapsed_contig_range[1]}")
        if contig_range[0] == 0:
            contig_seq_frag2 = contig_seq[contig_range[1]:]
            contig_frag_start = contig_range[1]
            contig_frag_end = contig_len 
        else:
            contig_seq_frag2 = contig_seq[:contig_range[0]]
            contig_frag_start = 0
            contig_frag_end = contig_range[0]
        contig_frag_id = f"{contig}:{contig_frag_start}-{contig_frag_end};frag"


        collapsed_contigs.add(collapsed_contig)
        collapsed_contig_db[collapsed_contig_pair[abs(idx - 1)]].append(collapsed_contig)

        print(f">{contig_seq_frag1_id}\n{str(contig_seq_frag1)}", file=sys.stdout)
        if len(collapsed_contig_seq_frag2) > 0:
            print(f">{collapsed_contig_frag_id}\n{str(collapsed_contig_seq_frag2)}", file=sys.stdout)
        if len(contig_seq_frag2) > 0:
            print(f">{contig_frag_id}\n{str(contig_seq_frag2)}", file=sys.stdout)
    
    else:
        for contig in fasta.keys():
            if contig in exists_contig_list:
                continue 
            print(f">{contig}\n{str(fasta[contig])}", file=sys.stdout)

    
    
    
    # for contig in list(fasta.keys()):
    #     if contig in collapsed_contigs:
    #         continue 
    #     if contig in collapsed_contig_db:
    #         output_id = "_".join(set(collapsed_contig_db[contig]))
    #         print(f">{contig}|collapse:{output_id}\n{fasta[contig]}")
    #     else:
    #         print(f">{contig}\n{fasta[contig]}")

if __name__ == "__main__":
    main(sys.argv[1:])
