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
    min_as, min_mapq,
    nhap, window, min_windows,
    min_sa, edge, min_depth,
    cutoff, threads, output, 
    steps, skip_steps
):
    from .correct_alignments import workflow 
    from .find_chimeric import ( 
        find_split_alignment, paf2depth, 
        norm_merge_bp, correct_fasta)
    
    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)

    if "1" not in skip_steps and "1" in steps:
        logger.info("""#-------------------------------------#
#  Running step 1. correct alignments #
#-------------------------------------#""")
        workflow(fasta, fastq, threads, output, window, min_windows, nhap, min_as, min_mapq)
    
    if "2" not in skip_steps and "2" in steps:
        logger.info("""#----------------------------------#
#  Running step 2. find chimeric   #
#----------------------------------#""")
        lis = f"{output}.mapq.LIS.gtf"
        corrected_paf = f"{output}.corrected.paf"
        break_point_file = find_split_alignment.workflow(lis, window, min_sa, edge, output)

        depth_file = paf2depth.workflow(corrected_paf, fasta, window, output)

        break_pos_file = norm_merge_bp.workflow(break_point_file, depth_file, 
                                                    window, min_depth, cutoff, output)

        correct_fasta.workflow(fasta, break_pos_file, output)
