#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
utility libraries
"""

import logging
import os
import os.path as op
import sys

from pathlib import Path, PosixPath
from subprocess import check_call

logger = logging.getLogger(__name__)

def run_cmd(command):
    """
    run command on shell

    Params:
    --------
    command: `list` or `tuple`, command list

    Returns:
    --------
    None

    Examples:
    --------
    >>> command = ["ls", "-l"]
    >>> run_cmd(command)
    example.txt result.txt
    0
    """
    logger.info('Running command:')
    logger.info('\t' + ' '.join(command))

    return check_call(command, shell=False)


def xopen(infile, mode='r'):
    """
    open file 

    Params:
    --------
    infile: `str`, input file
    mode: `str`, mode of open ["r"]

    Returns:
    --------
    handle: `_io.TextIOWrapper`

    Examples:
    --------
    >>> xopen('input.fastq.gz', 'r')
    
    """
    import gzip

    if not isinstance(infile, PosixPath):
        infile = Path(infile)
    
    if infile.suffix == ".gz":
        handle = gzip.open(infile, mode + 't')
    else:
        handle = open(infile, mode)
    return handle


def list_flatten(list_2d):
    """
    convert 2d list into 1d list

    Params:
    --------
    list_2d: `list` or `array-like`
            2d list [[1, 2, 3], [4, 5, 6]]
    
    Returns:
    --------
    1d list 

    Examples:
    --------
    >>> l = [[1, 2, 3], [4, 5, 6]]
    >>> list_flatten(l)
    [1, 2, 3, 4, 5, 6]
    """

    return [i for item in list_2d for i in item]
    