#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
evaluate the haplotype phasing by hic based on homologs
"""
import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np
import polars as pl 

from itertools import combinations
from pyranges import PyRanges

class Eval:
    def __init__(self, pairs, bed, output):
        self.pairs = pairs
        self.bed = bed
        self.output = output

    def read_pairs(self):
        dtype = {'chrom1': pl.Categorical, 'chrom2': pl.Categorical, 'mapq': pl.Int8,
                 'pos1': pl.UInt64, 'pos2': pl.UInt64}
        pairs_df = pl.read_csv(self.pairs, separator='\t', has_header=False,
                                comment_prefix="#", columns=[1, 2, 3, 4,],
                                new_columns=['chrom1', 'pos1', 'chrom2', 'pos2'],
                                        dtypes=dtype)
        
        self.pairs_df = pairs_df
    
    def read_bed(self):
        bed_df = pl.read_csv(self.bed, separator='\t', has_header=False,
                                comment_prefix="#", columns=[0, 1, 2, 3],
                                new_columns=['chrom', 'start', 'end', 'gene'])
        bed_df = bed_df.to_pandas()
        bed_df.sort_values(['chrom', 'start'], inplace=True)
        bed_df.reset_index(drop=True, inplace=True)
        self.bed_df = bed_df

    def get_homolog_pairs(self):
        bed_df = self.bed_df
        bed_df['source'] = bed_df['gene'].str.rsplit('.', n=1).str[0]
        def func(tmp_df):
            return list(combinations(tmp_df.index, 2))
        
        res_df = bed_df.groupby('source')[['chrom', 'start', 'end']].apply(lambda x: pd.Series(func(x)))
        
        return res_df[res_df.apply(lambda x: len(x) > 1)]

    def get_region_index_db(self):
        db = self.bed_df.reset_index().set_index(['chrom', 'start', 'end'])['index'].to_dict()
        
        index_to_region = dict(zip(db.values(), db.keys()))

        return index_to_region, db
    

    def run(self):
        self.read_pairs()
        self.read_bed()
      
        homolog_pairs = self.get_homolog_pairs()
        index_to_retions, regions_to_index = self.get_region_index_db()
        
        bed_df = self.bed_df[['chrom', 'start', 'end']].reset_index().rename(columns={"chrom": "Chromosome",
                                                               "start": "Start", 
                                                               "end": "End", "index": "index"})
        
        gr = PyRanges(bed_df.rename(columns={"index": "index1"}))       

        df1 = self.pairs_df[['chrom1', 'pos1']].to_pandas()
        df1['end'] = df1['pos1'] + 1
        df1.columns = ['Chromosome', 'Start', 'End']
        
        df1 = df1.reset_index()
        gr1 = PyRanges(df1)
        # gr.count_overlaps(gr1).df[['index', 'NumberOverlaps']]
        res1 = gr.join(gr1).new_position("intersection").df[['index1', 'index']]
        res1 = res1.drop_duplicates(subset=['index'])  

        gr = PyRanges(bed_df.rename(columns={"index": "index2"}))       
        df1 = self.pairs_df[['chrom2', 'pos2']].to_pandas()
        df1['end'] = df1['pos2'] + 1
        df1.columns = ['Chromosome', 'Start', 'End']
        df1 = df1.reset_index()
        gr1 = PyRanges(df1)
        # gr.count_overlaps(gr1).df[['index', 'NumberOverlaps']]
        res2 = gr.join(gr1).new_position("intersection").df[['index2', 'index']]
        res2 = res2.drop_duplicates(subset=['index'])     
   
        res = pd.concat([res1.set_index('index'), res2.set_index('index')], axis=1).dropna()
        res = res.astype({"index1": 'int64', "index2": "int64"})
        res = res.reset_index().set_index(['index1', 'index2'])
        counts = res.groupby(['index1', 'index2'])['index'].nunique()
        cis_counts = counts.reset_index()
        cis_counts = cis_counts[cis_counts['index1'] == cis_counts['index2']]
        cis_counts_db = cis_counts[['index1', 'index']].set_index('index1').to_dict()['index']
       
        cis_counts['index1'] = cis_counts['index1'].map(index_to_retions.get)
        cis_counts['index2'] = cis_counts['index2'].map(index_to_retions.get)
        cis_counts = cis_counts.sort_values(by=['index1', 'index2'])
        
        htrans_counts = counts.reindex(homolog_pairs.values).dropna().reset_index()
        htrans_counts['cis1'] = htrans_counts['index1'].map(cis_counts_db.get)
        htrans_counts['cis2'] = htrans_counts['index2'].map(cis_counts_db.get)
        htrans_counts = htrans_counts.fillna(0)
        def func(row):
            return row['index'] / (row['index'] + row['cis1'] + row['cis2'])
        
        htrans_counts['index1'] = htrans_counts['index1'].map(index_to_retions.get)
        htrans_counts['index2'] = htrans_counts['index2'].map(index_to_retions.get)
        htrans_counts =  htrans_counts.sort_values(by=['index1', 'index2'])
        htrans_counts['value'] = htrans_counts.apply(func, axis=1)
        

        htrans_counts.to_csv(self.output, sep='\t', index=None, header=None)
 

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('pairs', 
            help='')
    pReq.add_argument('bed',
            help='')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    eval = Eval(args.pairs, args.bed, args.output)
    eval.run()


if __name__ == "__main__":
    main(sys.argv[1:])