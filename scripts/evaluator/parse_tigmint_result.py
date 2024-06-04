#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
parser tigmint-long result to break points
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
    pReq.add_argument('bed', 
            help='tigmint bbraktigs.bed')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    df = pd.read_csv(args.bed, sep='\t', header=None, index_col=None)

    df = df[df[3].str.match('.*-[0-9]$')]

    df[4] = df[3].str.rsplit("-", n=1).map(lambda x: x[1]).astype(int)
    
    df = df[df[4] % 2 == 0]
    
    df[1] = df[1] + (df[2] - df[1])//2

    df[[0, 1]].to_csv(args.output, sep='\t', header=None, index=None)

if __name__ == "__main__":
    main(sys.argv[1:])