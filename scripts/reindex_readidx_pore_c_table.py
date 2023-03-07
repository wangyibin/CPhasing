#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
reindex the readidx of pore-c table 
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
    pReq.add_argument('table', 
            help='pore-c table')
    pReq.add_argument('output',
            help='output pore-c table')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    df = pd.read_parquet(args.table)
    df['read_idx'] = df['read_name'].astype('category').cat.codes
    df['identity'] = 1.0
    df.to_parquet(args.output)

if __name__ == "__main__":
    main(sys.argv[1:])