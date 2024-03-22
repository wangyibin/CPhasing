#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
evaluate the mapping accuracy of the hitig slide mapping  methods
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 
import pyranges 

from collections import defaultdict



def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('paf', 
            help='mapping record with paf format, which read id must generate from paftools.js pbsim2fq')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    paf = args.paf 
    paf_df = pd.read_csv(paf, sep='\t', index_col=None, header=None, usecols=range(12))
    # paf_df = paf_df[paf_df[16] == "tp:A:P"].iloc[:, :13]

    paf_df[0] = paf_df[0].map(lambda x: x.rsplit("_", 1)[0].split("!"))
    
    mapping_errors = defaultdict(int)
    mapping_corrects = defaultdict(int)
    
    for i, row in paf_df.iterrows():
        read_items = row[0]
        if read_items[1] != row[5]: 
            mapping_errors[read_items[0]] += 1
            continue 
        tstart, tend = row[0][2: 4]
        tstart, tend = int(tstart), int(tend)
        qstart, qend = row[[7, 8]]

        if tstart <= qstart <= tend:
            if qend - tend < 500:
                mapping_corrects[read_items[0]] += 1
            else:
                mapping_errors[read_items[0]] += 1

        else:
            if tstart - qstart < 500:
                mapping_corrects[read_items[0]] += 1
            else:
                mapping_errors[read_items[0]] += 1
    

    count = 0
    error_rates = []
    for read in mapping_corrects:
        error_count = mapping_errors[read]
        correct_count = mapping_corrects[read]
        error_rate = error_count / (correct_count + error_count)
        error_rates.append(error_rate)
        if error_rate > 0:
            print(read, error_rate)
        if error_rate < 0.1:
            count += 1

    print(sum(error_rates)/ len(error_rates), count / len(error_rates), file=sys.stderr)

        


        

if __name__ == "__main__":
    main(sys.argv[1:])
