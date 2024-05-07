#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
stat the porec table
    fragment length
"""




import argparse
import logging
import os
import os.path as op
import sys

import numpy as np
import pandas as pd 


from cphasing.core import PAFTable


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('paf', 
            help='')
    pOpt.add_argument('--low_memory', action='store_true', default=False, 
                        help='use chunk to reduce memory usage for huge paf file.')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    paf = args.paf 
    HEADER = [
        'read_name',
        'read_length',
        'read_start', 
        'read_end', 
        'chrom', 
        'start', 
        'end', 
        'matches', 
        'align_length', 
        'mapping_quality',
        ]
    
    if not args.low_memory:
        pd.set_option('display.float_format',  '{:,.4f}'.format)
        df = pd.read_csv(paf, sep='\t', usecols=[0, 1, 2, 3, 5, 7, 8, 9, 10, 11], names=HEADER, 
                            header=None, dtype=PAFTable.META)
        
        
        df['read_idx'] = df['read_name'].cat.codes
        del df['read_name']
        df_grouped = df.groupby('read_idx', sort=False)
        mapping_reads = len(df_grouped['chrom'])
        mapping_bases = df_grouped['read_length'].head(1).sum()
        
        df_grouped_nunique = df_grouped['read_start'].nunique()
        
        mean_fragment_number = df_grouped_nunique.mean()
        median_fragment_number = df_grouped_nunique.median()
        singleton_fragment_number = (df_grouped_nunique == 1).count()
        df['fragment_length'] = df['read_end'] - df['read_start']
        df['identity'] = df['matches'] / df['align_length']
        mean_fragment_length = df['fragment_length'].mean()
        median_fragment_length = df['fragment_length'].median()
    
    
        total_alignments = len(df)
        df = df.query("mapping_quality >= 1 & identity > 0.8 & fragment_length > 150" )
        df_grouped = df.groupby('read_idx', sort=False)
        unique_alignments = len(df)


        df_grouped_nunique = df_grouped['read_start'].nunique()
        unique_reads = df_grouped_nunique.count()
        unique_mean_fragment_number = df_grouped_nunique.mean()
        unique_median_fragment_number = df_grouped_nunique.median()
        unique_mean_fragment_length = df['fragment_length'].mean()
        unique_median_fragment_length = df['fragment_length'].median()
        unique_singleton_read_number = (df_grouped_nunique == 1).value_counts()[True]
        unique_valid_read_number = (df_grouped_nunique > 1).value_counts()[True]
        unique_high_order_read_number = (df_grouped_nunique > 2).value_counts()[True]
        df_orders = df_grouped_nunique[df_grouped_nunique > 1] 
    
        pairwise_contacts = (df_orders * (df_orders - 1) / 2).sum()
    
    else:
        chunks = pd.read_csv(paf, sep='\t', usecols=[0, 1, 2, 3, 5, 7, 8, 9, 10, 11], names=HEADER, 
                            header=None, dtype=PAFTable.META, chunksize=1e8)
        
        mapping_reads = 0
        mapping_bases = 0
        total_alignments = 0
        fragment_numbers = []
        fragment_lengths = []

        unique_reads = 0
        unique_fragment_numbers = []
        unique_fragment_lengths = []

        unique_singleton_read_number = 0
        unique_valid_read_number = 0
        unique_high_order_read_number = 0
        pairwise_contacts = 0
        
        
        for df in chunks:
        
        
            df['read_idx'] = df['read_name'].cat.codes
            del df['read_name']
            df_grouped = df.groupby('read_idx', sort=False)
            mapping_reads += len(df_grouped['chrom'])
            mapping_bases += df_grouped['read_length'].head(1).sum()
        
            df_grouped_nunique = df_grouped['read_start'].nunique()
            
            fragment_numbers.append(df_grouped_nunique)
         

            df['fragment_length'] = df['read_end'] - df['read_start']
            
            fragment_lengths.append(df['fragment_length'])
    
            df['identity'] = df['matches'] / df['align_length']
            del df['matches'], df['align_length'], df['read_end']
            total_alignments += len(df)

            df = df.query("mapping_quality >= 1 & identity > 0.8 & fragment_length > 150" )
            df_grouped = df.groupby('read_idx', sort=False)
            unique_alignments = len(df)

            df_grouped_nunique = df_grouped['read_start'].nunique()
            unique_reads += df_grouped_nunique.count()
            unique_fragment_numbers.append(df_grouped_nunique)
            unique_fragment_lengths.append(df['fragment_length'])
            
            unique_singleton_read_number += (df_grouped_nunique == 1).value_counts()[True]
            unique_valid_read_number += (df_grouped_nunique > 1).value_counts()[True]
            unique_high_order_read_number += (df_grouped_nunique > 2).value_counts()[True]
            df_orders = df_grouped_nunique[df_grouped_nunique > 1] 
        
            pairwise_contacts += (df_orders * (df_orders - 1) / 2).sum()

        fragment_lengths = pd.concat(fragment_lengths )
        mean_fragment_length = fragment_lengths.mean()
        median_fragment_length = fragment_lengths.median()
        
        fragment_numbers = pd.concat(fragment_numbers)
        mean_fragment_number = fragment_numbers.mean()
        median_fragment_number = fragment_numbers.median()

        unique_fragment_numbers = pd.concat(unique_fragment_numbers)
        unique_mean_fragment_number = unique_fragment_numbers.mean()
        unique_median_fragment_number = unique_fragment_numbers.median()

        unique_fragment_lengths = pd.concat(unique_fragment_lengths)
        unique_mean_fragment_length = unique_fragment_lengths.mean()
        unique_median_fragment_length = unique_fragment_lengths.median()


    res_df = pd.DataFrame()
    res_df.loc["Mapping reads (M)", "counts"] = mapping_reads / 1e6
    res_df.loc["Mapping bases (Gb)", "counts"] = mapping_bases / 1e9
    res_df.loc["Mapping fragments (M)", "counts"] = total_alignments / 1e6
    res_df.loc["Mean fragment number", "counts"] = mean_fragment_number
    res_df.loc["Median fragment number", "counts"] = median_fragment_number
    res_df.loc["Mean fragment length (bp)", "counts"] = mean_fragment_length
    res_df.loc["Median fragment length (bp)", "counts"] = median_fragment_length

    res_df.loc["Unique reads (M)", "counts"] = unique_reads / 1e6
    res_df.loc["Unique fragments (M)", "counts"] = unique_alignments / 1e6
    res_df.loc["Unique mean fragment number", "counts"] = unique_mean_fragment_number
    res_df.loc["Unique median fragment number", "counts"] = unique_median_fragment_number
    res_df.loc["Unique mean fragment length (bp)", "counts"] = unique_mean_fragment_length
    res_df.loc["Unique median fragment length (bp)", "counts"] = unique_median_fragment_length
    res_df.loc["Unique singleton (M)", "counts"] = unique_singleton_read_number / 1e6
    res_df.loc["Unique reads (Order ≥ 2) (M)", "counts"] = unique_valid_read_number / 1e6
    res_df.loc["Unique reads (Order ≥ 3) (M)", "counts"] = unique_high_order_read_number / 1e6 
    res_df.loc["Pairwise contacts (M)", "counts"] = pairwise_contacts / 1e6
    res_df.loc["Valid rate", "counts"] = unique_valid_read_number / mapping_reads

    res_df.to_csv(args.output, sep='\t', index=True, header=None)
   

if __name__ == "__main__":
    main(sys.argv[1:])