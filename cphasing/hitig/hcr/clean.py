#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
clean low coverage contigs of a fasta
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np


from collections import OrderedDict
from pathlib import Path

from ...utilities import read_fasta

logger = logging.getLogger(__name__)


def clean(fasta, low_coverage_df, output, break_pos=None):
    fa_db = read_fasta(fasta)
    low_coverage_contigs = set(low_coverage_df.index.tolist())

    if break_pos:
        break_pos_db = OrderedDict()
        with open(break_pos, 'r') as fp:
            for line in fp:
                if line.strip():
                    line_list = line.strip().split()
                    contig, positions = line_list 
                    positions = positions.split(",")
                    positions = list(map(int, positions))
                    if contig not in break_pos_db:
                        break_pos_db[contig] = []
                    for pos in positions:
                        break_pos_db[contig].append(pos)
    with open(output, 'w') as out:
        for contig in fa_db:
            if contig in low_coverage_contigs:
                if break_pos:
                    if contig.rsplit(":", 1)[0] in break_pos_db:

                        out.write(f">{contig}\n{fa_db[contig]}\n")
                continue

            out.write(f">{contig}\n{fa_db[contig]}\n")

    logger.info(f"Successful output low coverage contigs removed fasta in `{output}`")


class Clean:
    """
    remove false assemblies or false duplication, which false dupliaction will only reatin one copy.
    """
    PAF_HADER = ["contig1", "length1", "start1", "end1", "strand",
                 "contig2", "length2", "start2", "end2", "mismatch",
                 "matches", "mapq", "identity"]
    
    DEPTH_HEADER = ["chrom", "start", "end", "count"]
    def __init__(self, paf, depth):
        self.paf = paf
        self.depth = depth 
    


    def read_paf(self):
        logger.info(f"Load alignments results `{self.paf}`")
        df = pd.read_csv(self.paf, sep='\t', header=None, usecols=range(13),
                         names=self.PAF_HADER, index_col=None)
        df['identity'] = df['identity'].map(lambda x: x.replace("id:f:", "")).astype('float64')
    
        df = df.sort_values(['contig2', 'start2'])

        self.paf_df = df 

        return df 

    def read_depth(self):
        logger.info(f"Load depth file: `{self.depth}`")
        df = pd.read_csv(self.depth, sep='\t', header=None, 
                            names=self.DEPTH_HEADER, index_col=None)
        

        return df 
        
    
    
        
