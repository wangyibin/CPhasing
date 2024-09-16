#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
get allelic or cross allelic contigs from haphic log, which haphic should add `--verbose` parameters.
"""

import argparse
import os
import os.path as op
import sys
import re 

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('log', 
            help='log file from haphic --verbose')
    pOpt.add_argument('-m', '--min_concordance_ratio', default=0.2, type=float,
                      help="minimum concordance ratio used in haphic [default %(default)s]")
    pOpt.add_argument('--cross', action='store_true', default=False, 
                        help='output cross allelic informations [default: %(default)s]')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    min_concordance_ratio = args.min_concordance_ratio 

    with open(args.log) as fp:
        flag = 'None'
        for line in fp:
            if "remove_allelic_HiC_links" not in line:
                continue 
            
            if "concordance_ratio" in line:
                flag = "allelic"
            elif "non-maximum" in line:
                flag = "cross-allelic"
            else:
                flag = "None"

            line_list = line.strip().split()
            if flag == 'allelic':
                concordance_ratio = float(line_list[-1].split("=")[-1])
                if concordance_ratio > min_concordance_ratio:
                    print("\t".join(line_list[4:6]), file=args.output)
            elif flag == "cross-allelic":
                if args.cross:
                    print("\t".join(line_list[4:6]), file=args.output)

if __name__ == "__main__":
    main(sys.argv[1:])




