#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
mapper of ALLHiC2
"""

import argparse
import logging
import os
import os.path as op
import sys

from allhic2.utilities import run_cmd


def global_mapping(index_name, fastq, threads):
    command = ['hisat2', '-x', index_name, '-U', fastq, 
                '--no-unal',  '-k', '1' '--no-spliced-alignment'
                '--no-softclip', '--threads', str(threads),
                '--reorder' ]

    run_cmd(command)

def trim_fastq():
    pass

def trimmed_mapping():
    command = []



