#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
to test the alignment cross homolog switch chimeric contigs
"""

import argparse
import logging
import os
import os.path as op
import sys

import numpy as np

from Bio import SeqIO 
from pyfaidx import Fasta 


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('fasta', 
            help='mult only two sequences, which represent haplotype 1 and haplotype 2')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    fasta = Fasta(args.fasta)

    length1 = len(fasta[0])
    length2 = len(fasta[1])
    length1_list = np.linspace(0, length1, 6)
    length2_list = np.linspace(0, length2, 6)

    for i in range(1, len(length1_list)):
        r1, l1 = int(length1_list[i-1]), int(length1_list[i])
        r2, l2 = int(length2_list[i-1]), int(length2_list[i])

        tmp_prefix_1 = fasta[0][r1:l1][: (l1-r1)// 2]
        tmp_suffix_1 =  fasta[0][r1:l1][(l1-r1) // 2: ]

        tmp_prefix_2 = fasta[1][r2:l2][: (l2-r2) // 2]
        tmp_suffix_2 =  fasta[1][r2:l2][(l2-r2) // 2: ]

        print(f">hap1!{r1}!{r1 + (l1-r1)//2}!hap2!{r2 + (l2-r2)//2}!{l2}\n{tmp_prefix_1}{tmp_suffix_2}", file=args.output)
        print(f">hap2!{r2}!{r2 + (l2-r2)//2}!hap1!{r1 + (l1-r1)//2}!{l1}\n{tmp_prefix_2}{tmp_suffix_1}", file=args.output)

if __name__ == "__main__":
    main(sys.argv[1:])
