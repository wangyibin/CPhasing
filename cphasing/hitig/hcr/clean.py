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

from ...utilities import read_fasta

logger = logging.getLogger(__name__)


def clean(fasta, low_coverage_df, output):
    fa_db = read_fasta(fasta)
    low_coverage_contigs = set(low_coverage_df.index.tolist())
    with open(output, 'w') as out:
        for contig in fa_db:
            if contig in low_coverage_contigs:
                continue

            out.write(f">{contig}\n{fa_db[contig]}\n")

    logger.info(f"Successful output low coverage contigs removed fasta in `{output}`")