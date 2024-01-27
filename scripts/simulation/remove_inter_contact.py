#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
remove inter contact in tair
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('porec', 
            help='')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    HEADER = ["read_idx", "read_length",
              "read_start", "read_end", 
              "strand", "chrom", "start", 
              "end", "mapping_quality", 
              "identity", "filter_reason"]
    df = pd.read_csv(args.porec, sep='\t', comment='#', header=None, index_col=None, names=HEADER)
    count_df = df.groupby(['read_idx', 'chrom'])['chrom'].count()
    count_df = count_df.rename("count")
    count_df = count_df.reset_index()
    max_idx = count_df.groupby(['read_idx'])['count'].idxmax().reset_index()['count']
    idx_df = count_df.loc[max_idx][['read_idx', 'chrom']]
    
    df = df.set_index(['read_idx', 'chrom'])
    df = df.loc[idx_df.values.tolist()]
    # unique_df = df.groupby('read_idx')['chrom'].unique().map(len) < 2
    # df = df.set_index('read_idx')
    # df = df.loc[unique_df]
    df = df.reset_index()
    df = df[HEADER]
    df.to_csv('remove_iner.porec.gz', sep='\t', index=None, header=None)



    

if __name__ == "__main__":
    main(sys.argv[1:])