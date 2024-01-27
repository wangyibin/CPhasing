#!/usr/bin/env python

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import re

from collections import OrderedDict, defaultdict


def create_peusdo_alleles(contig_list, hap_pairs, output):
    regex_chrom = re.compile(r'(.*)\.ctg\d+')
    chrom_contigs_db = defaultdict(list)
    contig2chrom_db = {}
    for contig in contig_list:
        try:
            tmp_chrom = regex_chrom.match(contig).groups()[0]
        except:
            continue

        if tmp_chrom not in hap_pairs:
            continue
        
        chrom_contigs_db[tmp_chrom].append(contig)
        contig2chrom_db[contig] = tmp_chrom

    for contig in contig2chrom_db:
        tmp_chrom = contig2chrom_db[contig]
        hap_chrom = hap_pairs[tmp_chrom]
        other_contig_list = chrom_contigs_db[hap_chrom]
        
        for contig2 in other_contig_list:
            if contig == contig2:
                continue

            print(f"{contig}\t{contig2}", file=output)

    
            



def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('contig_list', 
            help='contig list')
    pReq.add_argument('hap_pairs', 
            help='two columns for haplotype pairs.')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    contig_list = [i.strip().split()[0] for i in open(args.contig_list) if i.strip()]
    hap_pairs_list = [i.strip().split()[:2] for i in open(args.hap_pairs) if i.strip()]

    hap_pairs = dict(hap_pairs_list) 
    hap_pairs.update(dict(list(map(lambda x: x[::-1], hap_pairs_list))))

    create_peusdo_alleles(contig_list, hap_pairs, args.output)

if __name__ == "__main__":
    main(sys.argv[1:])