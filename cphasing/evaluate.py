#!/usr/bin/env python
# -*- coding:utf-8 -*-


import argparse
import logging
import os
import os.path as op
import sys

import numpy as np
import pandas as pd 

from bisect import bisect_left
from itertools import combinations
from joblib import Parallel, delayed
from scipy.stats import spearmanr

from .core import AlleleTable, ClusterTable
from .utilities import read_chrom_sizes, list_flatten

def import_bed(bed):
    df = pd.read_csv(bed, sep='\t', header=None, usecols=[0, 1, 2, 3])
    df.columns = ['chrom', 'start', 'end', 'gene']
   

    return df 

def allelic_error(cluster_table, alleletable, contigsizes):
    ct = ClusterTable(cluster_table)
    at = AlleleTable(alleletable, sort=False, fmt='allele2')
    df = read_chrom_sizes(contigsizes)

    allelic_pairs = set(map(tuple, at.data.query(f"similarity > 0.85 & mz1 > 10000")[[1, 2]].values.tolist()))
    
    for group in ct.groups:
        contigs = ct.data[group]
        if len(contigs) <= 1:
            continue
        pairs = set(combinations(contigs, 2))
        
        group_length = df.loc[contigs]['length'].sum()

        error_pairs = allelic_pairs & pairs
        error_contigs = list(set(list_flatten(list(error_pairs))))
    
        length = df.loc[error_contigs]['length'].sum()
        print(group, length/group_length)


def load_gff(gff):

    df = pd.read_csv(gff, sep='\t', header=None, comment='#')
    df.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    df = df[df['feature'] == 'gene']
    df['gene'] = df['attribute'].apply(lambda x: x.split(';')[0].split('=')[1])

    return df

def gff2bed(gff):
    df = load_gff(gff)
    df = df[['chrom', 'start', 'end', 'gene']]
    df.reset_index(drop=True, inplace=True)

    return df

def allelic_completeness_index(mono_bed, contig_bed, poly_bed, ploidy):
    """

    """
    mono_bed = import_bed(mono_bed)
    contig_bed = import_bed(contig_bed)

    poly_bed = import_bed(poly_bed)
    
    contig_bed['source'] = contig_bed['gene'].apply(lambda x: x.rsplit('.', 1)[0])
    contig_source_gene_count = contig_bed.groupby(['chrom', 'source'], as_index=False)['gene'].count()
    
    duplicate = contig_source_gene_count[contig_source_gene_count['gene'] > 1]
    tandem_genes = contig_bed.loc[duplicate.index.to_list()]['gene'].values.tolist()

    # bed = import_bed(args.poly_bed)
    bed = poly_bed
    bed = bed[~bed['gene'].isin(tandem_genes)]
    bed['source'] = bed['gene'].apply(lambda x: x.rsplit('.', 1)[0])
    source_counts = bed.groupby('source')['gene'].count()
    source_counts = source_counts[source_counts == ploidy]
    source_counts = source_counts.index.tolist()
    
    mono_bed = mono_bed[mono_bed['gene'].isin(source_counts)]
    mono_gene_chrom_db = {}
    for chrom, tmp_df in mono_bed.groupby('chrom', sort=False):
        # n_hap_count_db[chrom] = len(tmp_df['gene'].drop_duplicates())
        mono_gene_chrom_db[chrom] = tmp_df['gene'].drop_duplicates().tolist()

    bed = bed[bed['source'].isin(source_counts)]
    bed = bed[bed['chrom'].str.contains('Chr')]
   
    bed['hap'] = bed['chrom'].apply(lambda x: x.rsplit('g', 1)[0])
    
    total_N_hap = 0
    total_v = 0
    total_dup_v = 0
    for hap, tmp_df in bed.groupby('hap', sort=False):
        candicate_genes = mono_gene_chrom_db[hap]
        candicate_genes = set(candicate_genes)
        N_hap = len(candicate_genes)    
        total_N_hap += N_hap
        for chrom, tmp_df2 in tmp_df.groupby('chrom', sort=False):
            # print(tmp_df2[tmp_df2['source'].duplicated(keep=False)])
            tmp_df2 = tmp_df2[tmp_df2['source'].isin(candicate_genes)]
            
            v = len(tmp_df2.drop_duplicates(subset='source', keep='first'))
            total_v += v
            dup_v = len(tmp_df2) - v
            total_dup_v += dup_v
            print(chrom, v/N_hap, dup_v/(v + dup_v), sep='\t')

    print("Total", total_v/total_N_hap/ploidy, total_dup_v/(total_v + total_dup_v), sep='\t')

def find_lis(seq1, seq2):
    n = len(seq1)
    m = len(seq2)
    dp = np.zeros((n+1, m+1))
    for i in range(1, n+1):
        for j in range(1, m+1):
            if seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1] + 1
            else:
                dp[i][j] = max(dp[i-1][j], dp[i][j-1])
    return dp[n][m]



def find_lis_fast(seq1, seq2):
    pos_map = {v: i for i, v in enumerate(seq2)}
    mapped = []
    for x in seq1:
        if x in pos_map:
            mapped.append(pos_map[x])

    lis = []
    for val in mapped:
        idx = bisect_left(lis, val)
        if idx == len(lis):
            lis.append(val)
        else:
            lis[idx] = val
    return len(lis)

def synteny_correlation_index(anchor, bed1, bed2, output):
    """

    """

    anchor = pd.read_csv(anchor, sep='\t', header=None, comment='#')
    anchor.columns = ['gene1', 'gene2', 'matches']
    # anchor = anchor[anchor['gene1'].map(lambda x: x.rsplit(".", 1)[0]) == anchor['gene2']]

    bed1 = import_bed(bed1)
    bed2 = import_bed(bed2)
    bed2.set_index('chrom', inplace=True)


    for chrom, df in bed1.groupby('chrom', sort=False):

        if chrom.startswith('utg'):
            continue
        
        df_merged = df.merge(anchor, left_on='gene', right_on='gene1')
     
        tmp_bed2 = bed2[bed2['gene'].isin(df_merged['gene2'].unique())].copy()
        tmp_bed2.sort_values(by='start', inplace=True)
        

        tmp_bed2 = tmp_bed2.reset_index().reset_index().set_index('gene')
        gene2_idx = tmp_bed2['index'].to_dict()


        df_merged['gene2_idx'] = df_merged['gene2'].map(gene2_idx).fillna(-1).astype(int)
        df_merged = df_merged[df_merged['gene2_idx'] != -1].reset_index(drop=True).reset_index()

        value = spearmanr(df_merged['gene2_idx'], sorted(df_merged['gene2_idx'].tolist()))
        lis = find_lis_fast(df_merged['gene2_idx'].tolist(),sorted(df_merged['gene2_idx'].tolist()))
        lis2 = find_lis_fast(df_merged['gene2_idx'].tolist(), sorted(df_merged['gene2_idx'].tolist())[::-1])
        lis = max(lis, lis2)
        lis_score = lis / len(df_merged)
        value = abs(value.correlation)
        print(chrom, value, lis_score, sep='\t', file=output)
        
    

  