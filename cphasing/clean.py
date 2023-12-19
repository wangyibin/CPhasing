#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
clean contig by contacts to improve the performance 
and fault tolerance of the algorithm
"""

import argparse
import logging
import os
import os.path as op
import sys

import cooler 

import numpy as np
import pandas as pd 

from collections import defaultdict
from math import sqrt
from joblib import Parallel, delayed


def identify_high_error_h_trans(contacts, contig_pairs, percent=0.98):

    contact_df = pd.read_csv(contacts, sep='\t', index_col=None, 
                                header=None, names=["contig1", "contig2", "count"])

    contact_dict = contact_df.set_index(['contig1', 'contig2']).to_dict()['count']

    contig_pairs = [tuple(i.strip().split()) for i in open(contig_pairs) if i.strip()]
    contig_pairs = list(filter(lambda x: x in contact_dict, contig_pairs))
    # contig_pairs = contact_df[['contig1', 'contig2']]
    # contig_pairs = contig_pairs[contig_pairs['contig1'] != contig_pairs['contig2']]
    # contig_pairs = contig_pairs.values.tolist()
    # contig_pairs = list(map(tuple, contig_pairs))
    def func(pair):
        cis1 = contact_dict.get((pair[0], pair[0]), 0)
        cis2 = contact_dict.get((pair[1], pair[1]), 0)
        trans = contact_dict.get(pair, 0)

        try:
            N = trans / sqrt(cis1 * cis2)
        except:
            N = np.nan

        
        return pair[0], pair[1], N
    
    res = []
    for pair in contig_pairs:
        res.append(func(pair))

    res_df = pd.DataFrame(res)
    res_df = res_df[~res_df[2].isna()]
    up_limit = np.percentile(res_df[2], 98)
    print(res_df[res_df[2] > 0.9])

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('contacts', 
            help='')
    pReq.add_argument('contig_pairs')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    identify_high_error_h_trans(args.contacts, args.contig_pairs)
  

if __name__ == "__main__":
    main(sys.argv[1:])
