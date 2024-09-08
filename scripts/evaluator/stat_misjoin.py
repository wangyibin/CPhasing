#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
stat misjoin points from paftools.js misjoin -e 
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

from pathlib import Path
from pyranges import PyRanges
from shutil import which

from cphasing.utilities import cmd_exists

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('tsv', 
            help='tsv from paftools.js misjoin -e ')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    try:
        df = pd.read_csv(args.tsv, sep='\t', header=None, index_col=None,
                        comment="#")
    except pd.errors.EmptyDataError:
        print(f"No. misassemblies: {0}", file=sys.stderr)
        return 
    candicate_res = []
    for read, tmp_df in df.groupby(1, sort=False):
        
        for i in range(0, len(tmp_df), 2):
            tmp_df2 = tmp_df.iloc[i:i+2]

            if (tmp_df2[5] == '-').all():
                tmp_df2 = tmp_df2.sort_values(3, ascending=False)
            

            res1 = (tmp_df2.iloc[0][6], tmp_df2.iloc[0][9] - 1000, tmp_df2.iloc[0][9] )
            res2 = (tmp_df2.iloc[1][6], tmp_df2.iloc[1][8], tmp_df2.iloc[1][8] + 1000)
                
            candicate_res.append(res1)
            candicate_res.append(res2)

    res_df = pd.DataFrame(candicate_res, columns=['Chromosome', 'Start', 'End'])                
    res_df = PyRanges(res_df).merge(count=True).df
    res_df = res_df.sort_values(['Chromosome', 'Start', 'End'])
    output = Path(args.tsv).stem
    res_df.to_csv(f'{output}.misassembly.txt', sep='\t', header=None, index=None)
    print(f"No. misassemblies: {len(res_df)}", file=sys.stderr)

if __name__ == "__main__":
    main(sys.argv[1:])