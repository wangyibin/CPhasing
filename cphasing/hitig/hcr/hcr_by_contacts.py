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
    binsize = cool.binsize
    bins = cool.bins()[:]
    
    contigsizes = cool.chromsizes
    matrix = cool.matrix(balance=False, sparse=True)[:]
    sum_values = np.array(matrix.sum(axis=1).T[0])
    logger.debug("Adjusting small bins value ...")
    small_bins = bins[bins['end'] - bins['start'] < cool.binsize]
    small_bins_sum_values = sum_values.T[small_bins.index]
    adjust_small_bins_sum_values = small_bins_sum_values.T / \
        ((small_bins['end'] - small_bins['start']) / binsize).values
    sum_values[:, small_bins.index] = adjust_small_bins_sum_values

    sum_values_nonzero = sum_values[sum_values > 0]

    #median = np.median(sum_values)
    max_value = np.percentile(sum_values_nonzero, percent)
    logger.debug(f"Percent{percent} value is {max_value}")
    res = np.where(sum_values < max_value)
    
    hcr_regions = bins.loc[res[1]]
    total_length = contigsizes.sum()
    num_hcr_regions = len(hcr_regions)
    logger.debug(f"Identify {num_hcr_regions} regions")
    hcr_length = sum(hcr_regions["end"] - hcr_regions["start"])
    logger.info(f"Identified {hcr_length/total_length:.2%} high-confidence regions")
    
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