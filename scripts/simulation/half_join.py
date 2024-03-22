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
    tmp_prefix_1 = fasta[0][: length1 // 2]
    tmp_suffix_1 =  fasta[0][length1 // 2: ]

    tmp_prefix_2 = fasta[1][: length2 // 2]
    tmp_suffix_2 =  fasta[1][length2 // 2: ]

    print(f">hap1!0!{length1//2}!hap2!{length2//2}!{length2}\n{tmp_prefix_1}{tmp_suffix_2}", file=args.output)
    print(f">hap2!0!{length2//2}!hap1!{length1//2}!{length1}\n{tmp_prefix_2}{tmp_suffix_1}", file=args.output)

if __name__ == "__main__":
    main(sys.argv[1:])
