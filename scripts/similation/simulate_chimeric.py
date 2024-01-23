#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
simulate chimeric contigs
    1. inner chrom chimeric: 0.1
    2. inter chrom chimeric:
        a. homo inter chrom chimeric: 0.8
            homo region switch: 0.6
            nonhomo region switch: 0.2
        b. nonhomo inter chrom chimeric: 0.1
"""

import argparse
import logging
import os
import os.path as op
import sys

import random
import pandas as pd

from collections import defaultdict
from pyfaidx import Fasta
from cphasing.core import AlleleTable



def random_join():
    pass 



def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('alleletable', 
            help='allele table from `cphasing alleles`')
    pReq.add_argument('fasta', 
            help='raw fasta')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    at = AlleleTable(args.alleletable, sort=False, fmt='allele2')
    high_similarity_data = at.data[at.data['similarity'] >= 0.99]
    high_similarity_data = high_similarity_data[high_similarity_data[1] < high_similarity_data[2]]
    idx1 = abs(high_similarity_data['mz1'] - high_similarity_data['mz2']) / high_similarity_data['mz1'] < 0.1
    idx2 = abs(high_similarity_data['mz1'] - high_similarity_data['mz2']) / high_similarity_data['mz2'] < 0.1
    high_similarity_data = high_similarity_data[pd.concat([idx1, idx2], axis=1).all(axis=1)]


    high_similarity_contig_pairs = high_similarity_data[[1, 2]].values

    



if __name__ == "__main__":
    main(sys.argv[1:])