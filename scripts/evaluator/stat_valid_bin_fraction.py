#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
stat the valid bins of a hic matrix by cool format
"""

import argparse
import logging
import os
import os.path as op
import sys


import cooler

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd 
import numpy as np 

from joblib import Parallel, delayed

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('cool', nargs="+",
            help='')
    pOpt.add_argument("-m", default=1, type=float, 
                      help='minimum contact counts')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    params = []
    for cool_file in args.cool:
        params.append((cool_file, args.m))

    def func(cool_file, m):
        cool = cooler.Cooler(cool_file)
        bins = cool.bins()
        matrix = cool.matrix(balance=False, sparse=True)[:]
        sum_array = np.array(matrix.sum(axis=1)).T[0]
        

        filter_sum_array = sum_array[sum_array >= m]

        return filter_sum_array.shape[0] / len(bins[:])
    

    res = Parallel(n_jobs=len(args.cool))(
        delayed(func)(i, j) for i, j in params
    )
    
    print("\t".join(map(str, res)), file=args.output)

if __name__ == "__main__":
    main(sys.argv[1:])