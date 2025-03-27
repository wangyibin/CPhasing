#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

"""


import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('contacts', 
            help='')
    pReq.add_argument('contig1')
    pReq.add_argument('contig2')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    
    args = p.parse_args(args)

    data = pd.read_csv(args.contacts, sep='\t', header=None, index_col=(0, 1))
    data  = data.to_dict()[2]
    
    try:
        cis1 = data[(args.contig1, args.contig1)]
    except KeyError:
        cis1 = 0
    try:
        cis2 = data[(args.contig2, args.contig2)]
    except KeyError:
        cis2 = 0

    try:
        if args.contig1 > args.contig2:
            trans = data[(args.contig2, args.contig1)]
        else:
            trans = data[(args.contig1, args.contig2)]
    except KeyError:
        trans = 0

    try:
        res = trans / (cis1 + cis2 + trans)
    except ZeroDivisionError:
        res = np.nan


    print(res)




if __name__ == "__main__":
    main(sys.argv[1:])