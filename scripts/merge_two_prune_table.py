#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
merge prune table 1 and prune table 2, which add the region from 2 to 1
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 

from cphasing.core import PruneTable

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('prune1', 
            help='')
    
    pReq.add_argument('prune2')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    pt1 = PruneTable(args.prune1)
    pt2 = PruneTable(args.prune2)

    pt1.data.set_index(['contig1', 'contig2'], inplace=True)
    pt2.data.set_index(['contig1', 'contig2'], inplace=True)

    res = pt1.data.merge(pt2.data, left_index=True, right_index=True, suffixes=("_1", "_2"))
    res = res[["mz1_2", "mz2_2", "mzShared_2", "similarity_2", "type_2"]]
    res = res.reset_index()
    res.columns = PruneTable.Header
    
    res.to_csv(args.output, header=None, index=None, sep='\t')


if __name__ == "__main__":
    main(sys.argv[1:])