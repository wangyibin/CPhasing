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
from cphasing.agp import import_agp

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('tsv', 
            help='tsv from paftools.js misjoin -e ')
    pReq.add_argument("agp", help="agp file")
    pOpt.add_argument('-m', '--min_count', default=2, type=int, 
            help='minimum support of read count [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    try:
        df = pd.read_csv(args.tsv, sep='\t', header=None, index_col=None,
                        comment="#")
    except pd.errors.EmptyDataError:
        print(f"No. misassemblies: {0}", file=sys.stderr)
        return 

    df = df[df[0] != "M"]
    
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
    res_df = res_df[res_df['Count'] >= args.min_count]
    output = Path(args.tsv).stem
    res_df.to_csv(f'.{output}.misassembly.txt', sep='\t', header=None, index=None)

    agp_df, _ = import_agp(args.agp)
    chrom_sizes = agp_df.groupby(agp_df.index)['end'].max()
    chrom_sizes.to_csv(f".{output}.sizes", sep='\t', header=None, index=True)

    cmd = f"grep -v -w U {args.agp} | grep '^Chr' | bedtools flank -g .{output}.sizes -i /dev/stdin  -l 2000 -r 0  | bedtools slop -i /dev/stdin -g .{output}.sizes -r 2000 -l 0 | bedtools intersect -a /dev/stdin  -b .{output}.misassembly.txt -wa 2>/dev/null | bedtools merge -i /dev/stdin  > {output}.misassembly.txt"
    
    os.system(cmd)


    if Path(f".{output}.sizes").exists():
        os.remove(f".{output}.sizes")
    if Path(f".{output}.missassembly").exists():
        os.remove(f".{output}.missassembly")


    res_df = pd.read_csv(f"{output}.misassembly.txt", header=None, index_col=None, sep='\t')

    print(f"No. misassemblies: {int(len(res_df))}", file=sys.stderr)

if __name__ == "__main__":
    main(sys.argv[1:])