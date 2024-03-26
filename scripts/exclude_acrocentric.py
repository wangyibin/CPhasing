#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
exclude the acrocentric for Human genome contig cluster
"""


import argparse
import logging
import os
import os.path as op
import sys

import numpy as np
import pandas as pd 


from cphasing.utilities import run_cmd
from shutil import which

logger = logging.getLogger(__name__)

## the acrocentric region were download from https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.0.bed
## and merge the censat region in chr13, chr14, chr15, chr21, chr22
acrocentric_regions = """chr13\t0\t22498291
chr14\t0\t17708240
chr15\t0\t22694129
chr21\t0\t16341849
chr22\t0\t20711065
"""

def output_acrocentric_regions():
    with open("acrocentric.regions.bed", 'w') as out:
        out.write(acrocentric_regions)

    return  "acrocentric.regions.bed"


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('contigs', 
            help='draft assembly')
    pReq.add_argument('chm13', help="the reference of chm13")
    pOpt.add_argument('-t', '--threads', type=int, default=8,
            help='number of program threads [default:%(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    ## run minigraph
    # if which('minigraph') is None:
    #         raise ValueError(f"minigraph: command not found.")
    
    # cmd = ["minigraph", f"-t{args.threads}", "-xasm", "-g10k", "--min-cov-mapq=0", "-j.001", args.chm13, args.contigs, 
            # "|", "awk", "'{print  $6,$8,$9,$1,$3,$4}'", "OFS='\\t'", ">", "map2chm.bed"]
    
    if which('minigraph') is None:
            raise ValueError(f"minigraph: command not found.")
    
    cmd = ["minimap2", "-cxasm5", "--secondary=no", f"-t{args.threads}", args.chm13, args.contigs, 
             "|", "awk", "'{print  $6,$8,$9,$1,$3,$4}'", "OFS='\\t'", ">", "map2chm.bed"]

    os.system(" ".join(cmd))
    acrocentric_bed = output_acrocentric_regions()
    
    cmd = ["bedtools", "intersect", "-a", "map2chm.bed", "-b", acrocentric_bed, "|", 
           "awk", "'$6-$5>5000{print $4,$5,$6}'", "OFS='\\t'", ">", "acrocentric.contigs.bed"]
    
    os.system(" ".join(cmd))
    logger.info(f"Output bed file: `{'acrocentric.contigs.bed'}`")

if __name__ == "__main__":
    main(sys.argv[1:])

