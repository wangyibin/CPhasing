#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
get allelic contig pairs from the ground truth
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

from natsort import natsort_keygen
from pytools import natsorted
from pyranges import PyRanges



from cphasing.utilities import read_chrom_sizes

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('contigsize', 
            help='the contigsize from simuCTG.py')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    contigsizes = read_chrom_sizes(args.contigsize)
    contigsizes.reset_index(inplace=True)
    contigsizes.sort_values('chrom', key=natsort_keygen(), inplace=True)
    
    raw_contigs = contigsizes['chrom'].str.split(".",n=2, expand=True)
    contigsizes['raw_chrom'] = raw_contigs[0]
    chromsizes = contigsizes.groupby('raw_chrom')['length'].sum().reset_index()
    chromsizes['homo'] = chromsizes['raw_chrom'].str[0]
    
    contigsizes['homo'] = contigsizes['raw_chrom'].str[0]
    
    res = []
    for raw_chrom, tmp_df in contigsizes.groupby('raw_chrom'):
       
        b = np.cumsum(tmp_df['length'])
        a = np.r_[[0], b[:-1]]
        
        tmp_df['start'] = a 
        tmp_df['end'] = b

        res.append(tmp_df)

    res_df = pd.concat(res, axis=0)
    res_df = res_df.rename(columns={"homo": 'Chromosome', 
                           "start": "Start",
                           "end": "End"})
    
    gr = PyRanges(res_df)

    gf = gr.join(gr, slack=0).new_position('intersection').df[['chrom', 'chrom_b']]
    gf = gf[gf['chrom'] < gf['chrom_b']]

    gf.to_csv(args.output, sep='\t', index=None, header=None)

    
    
        
        
 



if __name__ == "__main__":
    main(sys.argv[1:])