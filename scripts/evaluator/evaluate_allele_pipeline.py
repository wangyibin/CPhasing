#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
evaluate allele precision and recall
"""


import argparse
import logging
import os
import os.path as op
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import shutil
import pandas as pd 
import tempfile

import pickle 
import pandas as pd


from pathlib import Path
from itertools import combinations, product, permutations
from joblib import Parallel, delayed
from cphasing.utilities import run_cmd

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('contigsize', 
            help='contigsizes')
    pReq.add_argument("alleletable", help="alleletable")

    pReq.add_argument('allhic_alleletable')
    pReq.add_argument('haphic_log', help="log of haphic, `--verbose` should be add")
    pReq.add_argument('n', type=int)
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    cmd = f'~/code/CPhasing/scripts/evaluator/get_allelic_ground_truth.py {args.contigsizes} > ground.truth.allelic'
    os.system(cmd)

    cmd = f'grep -v "#" {args.alleletable} | cut -f 3,4 > cphasing.allelic'
    os.system(cmd)

    cmd = f'~/code/CPhasing/scripts/evaluator/get_allelic_from_ALLHiC {args.allhic_alleletable} > allhic.allelic'
    os.system(cmd)

    cmd = f'~/code/CPhasing/scripts/evaluator/get_allelic_from_haphic_log.py {args.haphic_log} > haphic.allelic'

    os.system(cmd)
    
    res = {}

    for software, allelic in [("C-Phasing", "cphasing.allelic"),
                              ("HapHiC", "haphic.allelic"),
                              ("ALLHiC", "allhic.allelic")]:
        
        cmd = f"~/code/CPhasing/scripts/evaluator/evaluate_allele.py ground.truth.allelic {allelic}"

        res[software] = os.popen(cmd).read()
    
        print(software, res[software], sep='\t')




if __name__ == "__main__":
    main(sys.argv[1:])