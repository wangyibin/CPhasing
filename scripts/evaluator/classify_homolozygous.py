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

import cooler 
from collections import defaultdict
from scipy.sparse import csr_matrix
from cphasing.utilities import read_chrom_sizes
from cphasing.utilities import list_flatten

def read_paf(paf):
    df = pd.read_csv(paf, sep='\t', header=None, index_col=None, usecols=range(16))
    df.columns = ['contig1', 'length1', 'start1', 'end1', 'strand',
                    'contig2', 'length2', 'start2', 'end2', 'matches',
                    'block', 'mapq', 'gi', 'bi', 'md', 'cg']
    df['identity'] = df['matches'] / df['block']
    
    return df 


def read_paf(paf):
    df = pd.read_csv(paf, sep='\t', header=None, index_col=None, usecols=range(19))
    df.columns = ['contig1', 'length1', 'start1', 'end1', 'strand',
                    'contig2', 'length2', 'start2', 'end2', 'matches',
                    'block', 'mapq', 'tp', 'nm', 'cm', 's1', 's2', 'dv', 
                    'cg']
    df['cg'] = df['cg'].str.replace("cg:Z:","")
    df['identity'] = df['matches'] / df['block']
    df = df[df['block'] > 1000]

    return df 


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('paf', 
            help='')
    pReq.add_argument('chromsizes')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    paf = args.paf
    paf_df = read_paf(paf)
    chromsizes = read_chrom_sizes(args.chromsizes)
    binsize = paf_df['length1'].max()
    
    paf_df['contig1'] = paf_df['contig1'].map(lambda x: x.split("_sliding:"))
    paf_df['contig1'] = paf_df['contig1'].map(lambda x: (x[0], int(x[1].split("-")[0]) - 1, int(x[1].split("-")[1])))
    

    paf_df['contig2'] = paf_df['contig2'].map(lambda x: x.split("_sliding:"))
    paf_df['contig2'] = paf_df['contig2'].map(lambda x: (x[0], int(x[1].split("-")[0]) - 1, int(x[1].split("-")[1])))
    

    # identity_db = defaultdict(list)
    # for idx, row in paf_df.iterrows():
    #     identity_db[row.contig1].append(row['identity'])
    #     identity_db[row.contig2].append(row['identity'])
    
    # identity_db2 = defaultdict(int)
    # for i in identity_db:
    #     identity_db2[i] = max(identity_db[i])

    identity_db2 = (
    pd.concat([
        paf_df[['contig1', 'identity']].rename(columns={'contig1': 'contig'}),
        paf_df[['contig2', 'identity']].rename(columns={'contig2': 'contig'})
    ])
    .groupby('contig')['identity']
    .max()
    .to_dict()
    )

    bins = cooler.util.binnify(chromsizes['length'], binsize)
    bins['identity'] = bins.apply(
    lambda row: identity_db2.get((row.chrom, row.start, row.end), 0),
    axis=1
)

    
    # def func(row):
    #     return identity_db2[(row.chrom, row.start, row.end)]
    
    # bins['identity'] = bins.apply(func, axis=1)

    # index_db = bins.reset_index(drop=False).set_index(['chrom', 'start', 'end']).to_dict()['index']
    

    
    # paf_df['contig1'] = paf_df['contig1'].map(index_db.get)
    # paf_df['contig2'] = paf_df['contig2'].map(index_db.get)
    # paf_df = paf_df[['contig1', 'contig2', 'identity']]

    # paf_df1 = paf_df[paf_df['contig1'] < paf_df['contig2']]
    # paf_df2 = paf_df[paf_df['contig1'] > paf_df['contig2']]

    # paf_df2.columns = ['contig2', 'contig1', 'identity']
    # paf_df = pd.concat([paf_df1, paf_df2], axis=0)
    # paf_df = paf_df.drop_duplicates(subset=['contig1', 'contig2'], keep='first')
    
    
    # matrix = csr_matrix((paf_df['identity'].values, (paf_df['contig1'].values, paf_df['contig2'].values)),
    #         shape=(len(bins), len(bins)),
    #             )


    # classify_bins = [0, 0.99, 0.995, 0.999, 1]
    ## translate this to sparse matrix 
    # indices = np.digitize(matrix.data, classify_bins)


    ## 0: no homology
    ## 1: homology with identity < 0.99
    ## 2: homology with identity >= 0.99 and < 0.995
    ## 3: homology with identity >= 0.995 and < 0.999
    ## 4: homology with identity >= 0.999
    ## get the col indices of the matrix, if the col is 4 any row with the col index is homology with identity >= 0.999
    # homologous_bins = {}
    # for i in range(4, -1, -1):
        
    #     tmp_res = np.where(indices == i)[0]
    #     if i < 4: 
    #         ## remove the bins that are already in homologous_bins
    #         tmp_res = np.setdiff1d(tmp_res, list_flatten(list(homologous_bins.values())))
    #     homologous_bins[i] = tmp_res
    
    # res = []
    # existed_idx = []
    # for i in homologous_bins:
    #     res_bins = bins.loc[homologous_bins[i]]
    #     res_bins['type'] = i  
    #     existed_idx.extend(list(homologous_bins[i]))
    #     res.append(res_bins)

    # unexisted_idx = set(bins.index.tolist()) - set(existed_idx)
    # un_res = bins.loc[unexisted_idx]
    # un_res['type'] = 1
    # res.append(un_res)
    # res = pd.concat(res, axis=0)
    # res.sort_values(by=['chrom', 'start'], inplace=True, ascending=True)
    
 
    
    def func(v):
        if v < 0.99:
            return ">1.0%" 
        elif 0.99 <= v < 0.995:
            return "0.5%-1.0%"
        elif 0.995 <= v < 0.999:
            return "0.1%-0.5%"
        elif 0.999 <= v < 0.9999:
            return "0.01%-0.1%"
        else:
            return "<0.01%"
        
    bins['type'] = bins['identity'].map(func)
    res = bins
    res['length'] = res['end'] - res['start']

    stat_res = res.groupby('type').agg({'length': 'sum'})
    res.drop('length', axis=1, inplace=True)
    
    stat_res.reset_index(inplace=True)
    stat_res['type'] = stat_res['type'].astype(str)
    stat_res = stat_res.set_index('type')
    ## reorder the index
    stat_res = stat_res.reindex(['>1.0%', '0.5%-1.0%', '0.1%-0.5%', '0.01%-0.1%', '<0.01%'])
    
    stat_res['length (Mb)'] = stat_res['length'] / 1e6
    stat_res['length (%)'] = stat_res['length'] / chromsizes['length'].sum() * 100

    stat_res.drop('length', axis=1, inplace=True)

    stat_res.to_csv(args.output.name.replace('.tsv', '_stat.tsv'), sep='\t', header=True, index=True)
    res.to_csv(args.output, sep='\t', header=None, index=None)
    
    
   

if __name__ == "__main__":
    main(sys.argv[1:])