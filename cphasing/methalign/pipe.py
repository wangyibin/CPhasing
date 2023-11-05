#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
a pipeline for aligning ont data by four bases and realign data by five bases
"""

import argparse
import logging
import os
import os.path as op
import sys

from ..utilities import run_cmd, cmd_exists

logger = logging.getLogger(__name__)


def split_bam(bam, num, outprefix):
    cmd = f'cphasing-rs splitbam {bam} -n {num} -o {outprefix}'
    pass

def first_mapping():
    cmd = f'dorado '
    pass 

def realign():
    cmd = f'cphasing-rs aligner'
    pass

def run():
    if not cmd_exists("dorado"):
        logger.warning("No such command of `dorado`")
        sys.exit(-1)
    
    if not cmd_exists("cphasing-rs"):
        logger.warning("No such command of `cphasing-rs`")
        sys.exit(-1)
    
    
