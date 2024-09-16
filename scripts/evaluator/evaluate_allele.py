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


def read_pairs(_file):
    df = pd.read_csv(_file, sep='\t', header=None, index_col=None, usecols=range(2))

    _df = df.copy()
    _df.columns = [1, 0]
    df = pd.concat([df, _df], axis=0)
    df = df[df[0] < df[1]]
    df = df.drop_duplicates(subset=[0, 1], keep='first')
    
    return df 



def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('ground', 
            help='ground truth')
    pReq.add_argument('query',
            help='query allelic contig pairs')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    df1 = read_pairs(args.ground)
    df2 = read_pairs(args.query)


    set1 = set(list(map(lambda x: tuple(x), df1.values.tolist())))
    set2 = set(list(map(lambda x: tuple(x), df2.values.tolist())))

    positive = len(set1.intersection(set2))
    false = len(set2 - set1)

    
    precision = positive / len(set2)
    recall = positive / len(set1) 

    f1 = 2 * recall * precision /  (recall + precision)
    print(precision, recall, f1, sep='\t')
    


if __name__ == "__main__":
    main(sys.argv[1:])