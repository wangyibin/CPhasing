#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
stat the chimeric correction by hitig
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 
import numpy as np

from pathlib import Path 

def parse_chimeric_ID(ID):
    pair, pos, orients, _type = ID.split(";")
    contig1, contig2 = pair.split("|")
    pos = int(pos)
    
    return (contig1, contig2, pos, _type)


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('break_pos', 
            help='break_pos.txt from `hitig`')
    pReq.add_argument('contigsizes',
            help='contig sizes')
    pOpt.add_argument('-d', '--dist', default=5000,
            help='minimum distance between break position and real position. [default: %(default)s]')  
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    contigsizes = dict(i.strip().split()[:2] for i in open(args.contigsizes) if i.strip())
    contigsizes = dict(map(lambda x: (x[0], int(x[1])), contigsizes.items()))

    real_chimeric_contigs = [i for i in contigsizes if "chimeric" in i]
    
    break_pos_file = args.break_pos

    correct_breaks = []
    incorrect_breaks = []

    with open(break_pos_file, 'r') as fp:
        for line in fp:
            if not line.strip():
                continue 
            line_list = line.strip().split()
            contig, positions = line_list 
            
            for position in positions.split(","):
                if contig not in real_chimeric_contigs:
                    incorrect_breaks.append((contig, position))
                    continue
                chimeric_info = parse_chimeric_ID(contig)
                if abs(int(position) - chimeric_info[2]) < int(args.dist):
                    correct_breaks.append((contig, position))
                else:
                    incorrect_breaks.append((contig, position))
                
    P = len(correct_breaks) / (len(incorrect_breaks) + len(correct_breaks))
    R = len(correct_breaks) / len(real_chimeric_contigs)
    print(f"Precision:{P}\nRecall{R}", file=sys.stderr)
    




if __name__ == "__main__":
    main(sys.argv[1:])