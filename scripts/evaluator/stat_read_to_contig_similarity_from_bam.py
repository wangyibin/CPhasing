#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

import pysam 


def get_query_alignment_termini(aln):

    if aln.is_forward:
        # '5' means hard-clipping
        if aln.cigartuples[0][0] == 5:
            hard_clip = aln.cigartuples[0][1]
        else:
            hard_clip = 0
        return aln.query_alignment_start + hard_clip, aln.query_alignment_end + hard_clip
    else:
        if aln.cigartuples[-1][0] == 5:
            hard_clip = aln.cigartuples[-1][1]
        else:
            hard_clip = 0
        read_len = aln.infer_read_length()
        return read_len - aln.query_alignment_end + hard_clip, read_len - aln.query_alignment_start + hard_clip
    

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bam', 
            help='')
    pOpt.add_argument('--no_secondary', action='store_true', default=False,
                      help='do not output secondary')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    format_options = [b'filter=!flag.unmap'] #, ]
    if args.no_secondary:
        format_options.append(b'filter=!flag.secondary')
    
    last_read = ''
    with pysam.AlignmentFile(args.bam, format_options=format_options, threads=8) as fp:
        for alignment in fp:
            if args.no_secondary:
                if alignment.is_secondary:
                    continue
            read_name, flag = alignment.query_name, alignment.flag
            try:
                inconsistent_sites = alignment.get_tag('MA')
            except KeyError:
                inconsistent_sites = 0
            mapq = alignment.mapq
            nm = alignment.get_tag('NM')
            AS = alignment.get_tag('AS')
            try:
                rf = alignment.get_tag('RF')
            except KeyError:
                rf = 'N'
            ref_name = alignment.reference_name
            ref_start, ref_end = alignment.reference_start, alignment.reference_end
            query_start, query_end = get_query_alignment_termini(alignment)

        
            canonical_identity = 1 - (nm / (query_end - query_start))
            methalign_identity = 1 - ((nm + inconsistent_sites) / (query_end - query_start))
            
            print(read_name, flag, mapq, 
                    query_start, query_end, 
                    ref_name, ref_start, ref_end, 
                    nm, AS, inconsistent_sites, rf, 
                    canonical_identity, methalign_identity,
                    sep='\t')



if __name__ == "__main__":
    main(sys.argv[1:])