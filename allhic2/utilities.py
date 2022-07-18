#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
utility libraries
"""

import argparse
import logging
import os
import os.path as op
import sys

from subprocess import check_call

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

    check_call(command, shell=True)


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

    if infile[-3:] == ".gz":
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
    