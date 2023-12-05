#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Identify the high confidence regions of contigs by Pore-C/Hi-C contacts
"""

import argparse
import logging
import os
import os.path as op
import sys

import cooler 
import numpy as np 
import pandas as pd
import pyranges as pr 

logger = logging.getLogger(__name__)


def hcr_by_contacts(cool_file, output, percent=95, ):
    percent = int(percent)
    cool = cooler.Cooler(cool_file)
    bins = cool.bins()[:]
    matrix = cool.matrix(balance=False, sparse=True)[:]
    sum_values = np.array(matrix.sum(axis=1).T[0])
    sum_values_nonzero = sum_values[sum_values > 0]

    #median = np.median(sum_values)
    max_value = np.percentile(sum_values_nonzero, percent)

    res = np.where(sum_values < max_value)
    
    hcr_regions = bins.loc[res[1]]
    num_hcr_regions = len(hcr_regions)
    logger.info(f"Identified {num_hcr_regions} high-confidence regions")
    
    hcr_regions.columns = ['Chromosome', 'Start', 'End']
    hcr_regions_pr = pr.PyRanges(hcr_regions)
    hcr_regions_pr.merge().df.to_csv(output, sep='\t', index=None, header=None)
    logger.info(f"Successful output HCRs into `{output}`.")

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('cool_file', 
            help='')
    pOpt.add_argument('-p', '--percent', type=int, default=95,
            help='percentile' )
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    hcr_by_contacts(args.cool_file, args.output, args.percent)

if __name__ == "__main__":
    main(sys.argv[1:])