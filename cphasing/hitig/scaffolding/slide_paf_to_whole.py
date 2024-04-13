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


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('paf', 
            help='')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    df = pd.read_csv(args.paf, sep='\t', header=None, index_col=None)
    df[0] = df[0].str.rsplit("_", n=1).map(lambda x: x[0])
    df.to_csv(args.output, sep='\t', header=None, index=None)
    

if __name__ == "__main__":
    main(sys.argv[1:])