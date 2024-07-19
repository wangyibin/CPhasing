#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
rename the name of AT to Chr??g??
"""


import argparse
import logging
import os
import os.path as op
import sys


import pandas as pd 

from collections import OrderedDict
from cphasing.core import ClusterTable



def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('cluster', 
            help='cluster table')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    ct = ClusterTable(args.cluster)
    
    hap_groups = OrderedDict()
    for group in  ct.groups:
        chrom, hap = group 
        
        if chrom not in hap_groups:
            hap_groups[chrom] = []
        
        hap_groups[chrom].append(group)
    
    
    rename_db = OrderedDict()
    for chrom in hap_groups:
        for i, old_chrom in enumerate(hap_groups[chrom], 1):
            rename_db[old_chrom] = f"Chr{chrom}g{i}"
    
    new_data = OrderedDict()
    for chrom in ct.data:
        new_data[rename_db[chrom]] = ct.data[chrom]
    
    ct.data = new_data 

    ct.save(args.output)

if __name__ == "__main__":
    main(sys.argv[1:])