#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot the remove and retain curve
"""


import argparse
import logging
import os
import os.path as op
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd 
import numpy as np 

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('tsv', 
            help='results from `evaluate_retain_and_remove_pipeline.py`')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    df = pd.read_csv(args.tsv, sep='\t', index_col=None, header=None)
    
    sns.scatterplot(data=df, x=1, y=2, hue=3)
    plt.savefig("output.png", dpi=300, bbox_inches='tight')



if __name__ == "__main__":
    main(sys.argv[1:])
