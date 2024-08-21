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

from collections import OrderedDict

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