#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
simulate collapsed contigs from vcf 
"""

import argparse
import logging
import os
import os.path as op
import sys

import random 
import gzip
import numpy as np 
import pandas as pd


from collections import defaultdict
from ncls import NCLS
from pathlib import Path

from cphasing.utilities import run_cmd, read_fasta, xopen


def get_random_region(chrom, chromsizes, min_size=10000000, max_size=15000000):
    chromsizes = chromsizes.copy()
    chromsizes[2] = chromsizes[1].astype(int)
    chromsizes[1] = 0
    chrom_df = chromsizes[chromsizes[0] == chrom]
    if chrom_df.empty:
        return None, None

    start = chrom_df[1].min()
    end = chrom_df[2].max()
    
    size = random.randint(min_size, max_size)
    if end - start < size:
        return None, None

    pos = random.randint(start, end - size)

    return pos, pos + size


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('chromsizes', 
            help='Chromosome sizes file in tab-separated format.')
    pReq.add_argument('vcf', 
            help='')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    chromsizes = args.chromsizes
    chromsizes = pd.read_csv(chromsizes, sep='\t', header=None, index_col=None, usecols=[0, 1])
    vcf = args.vcf

    res = []
    for chrom, _ in chromsizes.groupby(0):
        if random.randint(0, 1) == 1:
            logging.info(f"Skipping {chrom}.")
            continue

        start, end = get_random_region(chrom, chromsizes)
        if start is None:
            logging.warning(f"No valid region found for {chrom}.")
            continue
        
        res.append((chrom, start, end))
    
    res_df = pd.DataFrame(res, columns=['chrom', 'start', 'end'])

    res_df.to_csv(f"{Path(vcf).stem}.collapsed_regions.bed", sep='\t', index=False, header=False)
    collapsed_bed = f"{Path(vcf).stem}.collapsed_regions.bed"
    logging.info(f"Collapsed regions saved to {Path(vcf).stem}.collapsed_regions.tsv")


    # Simulate the collapse of contigs in the VCF file
    output_vcf = f"{Path(vcf).stem}.collapsed.vcf.gz"
    cmd = f"zgrep '#' {vcf} > header"
    os.system(cmd)
    cmd = f"bedtools intersect -a {vcf} -b {collapsed_bed} -v > body"
    os.system(cmd)

    cmd = f"cat header body | bgzip -c > {output_vcf}"
    os.system(cmd)

    cmd = f"tabix -p vcf {output_vcf}"
    os.system(cmd)


if __name__ == "__main__":
    main(sys.argv[1:])

