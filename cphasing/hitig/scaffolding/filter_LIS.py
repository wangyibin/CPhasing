#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
filter LIS
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 



def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('lis', 
            help='LIS gtf from hitig correct-alignments')
    pOpt.add_argument("-s", "--min_lis", default=3, type=int, 
                      help="minimum length of LIS")
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)



    tmp_res = []
    with open(args.lis) as fp:
        for row in fp:
            row = row.strip().split()
            if row[2] == "LIS":
                if tmp_res:
                    if len(tmp_res) >= args.min_lis:
                        for tmp in tmp_res:
                            print("\t".join(map(str, tmp)), file=args.output)
                tmp_res = []
                tmp_res.append(row)
            else:
                if row[7] == "E":
                    continue
                tmp_res.append(row)


        


if __name__ == "__main__":
    main(sys.argv[1:])