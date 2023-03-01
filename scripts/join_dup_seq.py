#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
join duplicated sequences to a sequence
"""

import argparse
import logging
import os
import os.path as op
import sys

from Bio import SeqIO
from collections import defaultdict


def join_dup_seq(fasta, output):
    db = defaultdict(list)
    
    fasta = SeqIO.parse(open(fasta, 'r'), "fasta")
    for record in fasta:
        db[record.name].append(str(record.seq))
    
    for name in db:
        seq = "NNNNN".join(db[name])
        print(f">{name}\n{seq}", file=output)

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('fasta', 
            help='path of fasta')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    join_dup_seq(args.fasta, args.output)

if __name__ == "__main__":
    main(sys.argv[1:])