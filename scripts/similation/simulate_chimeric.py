#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
simulate chimeric contigs
"""

import argparse
import logging
import os
import os.path as op
import sys


def random_join():
    pass 


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('', 
            help='')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

if __name__ == "__main__":
    main(sys.argv[1:])