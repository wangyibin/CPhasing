#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
evaluate retain cis and remove h-trans
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

import pandas as pd
import numpy as np 

from itertools import combinations, product
from scipy.optimize import curve_fit

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('remove_list')
    pReq.add_argument('contigs', 
            help='contig list')
    pReq.add_argument('contacts', 
            help='')
    
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    remove_list = set([tuple(sorted(i.strip().split()[:2])) for i in open(args.remove_list) if i.strip()])
    contigs_df = pd.read_csv(args.contigs, sep='\t', usecols=(0,), header=None, index_col=None)
    chrom_contigs = pd.DataFrame(contigs_df[0].str.split(".").values.tolist(), columns=["chrom", "contig"])
    chrom_contigs['contig'] = contigs_df[0]
    chrom_contigs['hap'] = chrom_contigs['chrom'].str[:-1]
    chrom_contigs['hap'] = chrom_contigs['hap'].str.split("_").map(lambda x: x[0])

    contact_df = pd.read_csv(args.contacts, sep='\t', index_col=None, 
                                    header=None)
    contact_df.columns = ['contig1', 'contig2', 'count']
    contact_df = contact_df.set_index(['contig1', 'contig2'])
    contact_db = contact_df.to_dict()['count']

    inter_pairs = []
    for i, df in chrom_contigs.groupby('hap'):
        tmp_list = []
        for j, df2 in df.groupby('chrom'):
            tmp_list.append(df2['contig'].tolist())

        tmp_list.sort() 
        
        for n, m in list(combinations(tmp_list, 2)):
            inter_pairs.extend(list(set(product(n, m))))
    
    intra_pairs = []
    for i, df in chrom_contigs.groupby('chrom'):
        tmp_list = df['contig'].values.tolist()
        tmp_list = sorted(tmp_list)
        intra_pairs.extend(list(combinations(tmp_list, 2)))
   
    intra_pairs = set(intra_pairs)

    intra_contacts = contact_df.reindex(intra_pairs).dropna()
    inter_contacts = contact_df.reindex(inter_pairs).dropna()
    total_intra_contacts = intra_contacts['count'].sum()
    total_inter_contacts = inter_contacts['count'].sum()
    retain_contacts =  total_intra_contacts -  intra_contacts.reindex(remove_list).dropna()['count'].sum()
    remove_contacts = inter_contacts.reindex(remove_list).dropna()['count'].sum()
    
    # for pair in contact_db:
    #     contact = contact_db[pair]
    #     if pair in intra_pairs:
    #         if pair not in remove_list:
    #             retain_contacts += contact
    #         total_intra_contacts += contact
    #     elif pair in inter_pairs:
    #         if pair in remove_list:
    #             remove_contacts += contact
    #         total_inter_contacts += contact 

    print(f"Removed h-trans error contacts: {remove_contacts / total_inter_contacts:.4f}" )
    print(f"Retain cis contacts: {retain_contacts / total_intra_contacts:.4f}")





if __name__ == "__main__":
    main(sys.argv[1:])