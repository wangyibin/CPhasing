#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
pipeline of hitig
"""

import argparse
import logging
import os
import os.path as op
import sys

from pathlib import Path
from ..utilities import (
    run_cmd, 
    calculate_Nx_from_contigsizes,
    read_chrom_sizes
    )


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def run(
    fasta, fastq,
    min_as=2000, min_mapq=2,
    nhap=4, window=5000, min_windows=10,
    step=1000,
    min_sa=5, edge=2000, min_depth=1,
    cutoff=0.75, hifi=False,
    threads=10, output="output", 
    steps=("1", "2", "3"), skip_steps=[]
):
    from .correct_alignments import workflow 
    from .find_chimeric import ( 
        find_split_alignment, paf2depth, 
        norm_merge_bp, correct_fasta)
    from .cli import hcr 
    
    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)

    if "1" not in skip_steps and "1" in steps:
        logger.info("""#-------------------------------------#
#  Running step 1. correct alignments #
#-------------------------------------#""")
        pafFile = workflow(fasta, fastq, threads, output, window, min_windows, nhap, min_as, min_mapq, hifi)
    else:
        pafFile = output + ".paf"

    if "2" not in skip_steps and "2" in steps:
        logger.info("""#----------------------------------#
#  Running step 2. find chimeric   #
#----------------------------------#""")
        lis = f"{output}.mapq.LIS.gtf"
        corrected_paf = f"{output}.corrected.paf"
        break_point_file = find_split_alignment.workflow(lis, window, min_sa, edge, output)

        contigsizes, depth_file = paf2depth.workflow(corrected_paf, fasta, window, step, output)

        break_pos_file = norm_merge_bp.workflow(break_point_file, depth_file, contigsizes,
                                                    window, min_depth, cutoff, edge, output)

        correct_fasta.workflow(fasta, break_pos_file, output)

    else:
        lis = output + ".mapq.LIS.gtf"
        break_point_file = "output" + ".mergedSplitAlign.txt"

    if "3" not in skip_steps and "3" in steps:
        logger.info("""#----------------------------------#
#  Running step 3. hcr   #
#----------------------------------#""")
        try:
            hcr.main(args=["-f", fasta, 
                           "-l", lis, 
                           "-sa", break_point_file,  
                           "-p", pafFile
                           ])

        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
        pass 