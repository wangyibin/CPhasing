#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
This script is used to break contigs and duplicate them in gff file.
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

import pyranges as pr 
import gffpandas.gffpandas as gffpd
import re
from collections import OrderedDict

from cphasing.agp import import_agp

def get_dup_contig_in_agp(agp_df):
    dup_contigs = agp_df[agp_df['id'].str.match(r'utg.*_d(\d+)')]['id']
    

    return dup_contigs

def get_breaked_contig_in_agp(agp_df):
    breaked_contigs = agp_df[agp_df['id'].str.match(r'utg.*:(\d+)-(\d+)')]['id']
    

    return breaked_contigs

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('agp', 
            help='')
    pReq.add_argument('gff')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    agp_df, _ = import_agp(args.agp)
    dup_contigs = get_dup_contig_in_agp(agp_df)
    breaked_contigs = get_breaked_contig_in_agp(agp_df)

    gff_df = gffpd.read_gff3(args.gff).df
    regex = re.compile(r'ID=(.*?);')
    dup_df_list = []
    for contig in dup_contigs:
        raw_contig, suffix = contig.rsplit("_", 1)

        tmp_df = gff_df[gff_df['seq_id'].isin([raw_contig])]
        if tmp_df.empty:
            continue
        tmp_df = gffpd.Gff3DataFrame(input_df=tmp_df)

        tmp_df.filter_feature_of_type(['gene'])
        for i, row in tmp_df.df.iterrows():
            attribute = row['attributes']
            ID = regex.findall(attribute)
            ID = ID[0].split(".")[0]
            attribute = attribute.replace(ID, ID + "_" + suffix)
            tmp_df.df.loc[i, 'attributes'] = attribute
        
        tmp_df.df['seq_id'] = tmp_df.df['seq_id'].apply(lambda x: x + "_" + suffix)
        
        dup_df_list.append(tmp_df.df)

    if dup_df_list:
        dup_df = pd.concat(dup_df_list)

    regions = []
    for contig in breaked_contigs.values:
        
        raw_contig, region = contig.rsplit(":", 1)
        start, end = region.split("-")
        start = int(start)
        end = int(end)
        
        regions.append((raw_contig, start, end, contig))
    
    region_df = pd.DataFrame(regions, columns=['Chromosome', 'Start', 'End', 'ID'])
    
    ## for gene gff3
    # gff_bed = gff_df.reset_index()[gff_df['type'] == 'gene'][['index', 'seq_id', 'start', 'end']]
    ## for TE gff3
    gff_bed = gff_df.reset_index()[['index', 'seq_id', 'start', 'end']]
    gff_bed = gff_bed[gff_bed['seq_id'].isin(region_df['Chromosome'].values.tolist())]
   
    gff_bed.columns = ['Index', 'Chromosome', 'Start', 'End']
  
    gff_gr = pr.PyRanges(gff_bed)
    region_gr = pr.PyRanges(region_df)
    breaked_df = gff_gr.join(region_gr).new_position('intersection').df

    duplicated_df = breaked_df[breaked_df['Index'].duplicated(keep=False)]

    if duplicated_df.empty:
        print("No breaked gene in gff file.", file=sys.stderr)
    else:
        print("Breaked gene in gff file:", file=sys.stderr)

        need_to_drop_df = gff_df.loc[duplicated_df['Index'].values]
        # need_to_drop_df.to_csv(sys.stderr, sep="\t", index=False, header=False)

        ## for TE gff3
        gff_df = gff_df.drop(duplicated_df['Index'].values)
    
    if dup_df_list:
       
        pd.concat([gff_df, dup_df]).to_csv(sys.stdout, sep="\t", index=False, header=False)
    else:
        gff_df.to_csv(sys.stdout, sep="\t", index=False, header=False)
 


    

if __name__ == "__main__":
    main(sys.argv[1:])
