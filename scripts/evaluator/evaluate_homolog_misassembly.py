#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
detect the homolog misassemblies
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 

from itertools import combinations

from cphasing.core import AlleleTable 
from cphasing.agp import import_agp

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('agp', 
            help='')
    pReq.add_argument('alleletable')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    agp, _ = import_agp(args.agp)
    at = AlleleTable(args.alleletable, sort=False, fmt="allele2")
    
    df = at.data[at.data['similarity'] >= 0.85].set_index([1, 2])
    
    for chrom, tmp_df in agp.groupby('chrom'):
        if len(tmp_df) <= 1:
            continue
        

   
        res_df = df.reindex(list(combinations(tmp_df['id'].values.tolist(), 2))).dropna()
        if res_df.empty:
            
            continue
        
        print(chrom)
        print(res_df)


if __name__ == "__main__":
    main(sys.argv[1:])