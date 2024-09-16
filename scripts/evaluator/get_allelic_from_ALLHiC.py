#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
get allelic contig pairs from `Allele.ctg.table`
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

from cphasing.core import AlleleTable

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('table', 
            help='Allele.ctg.table')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    at = AlleleTable(args.table)
    
    for contig_pair in at.contig_pairs:
        print("\t".join(contig_pair), file=args.output)



if __name__ == "__main__":
    main(sys.argv[1:])