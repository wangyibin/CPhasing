#!/usr/bin/env python
# -*- coding:utf-8 -*-


import logging
import gc
import os
import os.path as op
import sys
import warnings

import cooler 
import numpy as np
import pandas as pd


from collections import defaultdict
from scipy.sparse import triu
from pandarallel import pandarallel 

from cphasing.plot import chrRangeID


logger = logging.getLogger(__name__)


def convert_matrix_with_dup_contigs(cool, dup_contig_path, output, threads=4):
    """
    Params:
    --------
    cool: str
        Path to cool file
    dup_contigs: str
        Path of two columns of `raw_contig dup_contig`
    output: str
        Path of output cool
    """
    pandarallel.initialize(nb_workers=threads, verbose=0)

    cool = cooler.Cooler(cool)
    bins = cool.bins()[:]

    bins_chrom_idx_db = bins.reset_index().groupby('chrom')['index'].apply(list).to_dict()
    bins_chrom = bins.reset_index().set_index('chrom')
    max_bin_id = len(bins)

    matrix = cool.matrix(balance=False, sparse=True)[:].tocsr().astype('float64')
 
    dup_contigs_db = defaultdict(list)
    with open(dup_contig_path) as fp:
        for line in fp:
            line_list = line.strip().split()[:2]
            dup_contigs_db[line_list[0]].append(line_list[1])

   
    bin_res = [bins.reset_index()]
    res = []

    for contig in dup_contigs_db:
        dup_contigs = dup_contigs_db[contig]
        cn = len(dup_contigs) + 1
        
        idxes = bins_chrom_idx_db[contig]
        tmp_bins = bins_chrom.loc[contig]
        tmp_raw_matrix = matrix[idxes]
        
       

        for i, dup_contig in enumerate(dup_contigs):
            
            length_of_tmp_bins = len(tmp_bins)
            tmp_bins_2 = tmp_bins.copy()
            tmp_bins_2['index'] = range(max_bin_id, max_bin_id + length_of_tmp_bins)
            
            tmp_bins_2['chrom'] = dup_contig 
            row_bin_idx_db = dict(zip(range(0, length_of_tmp_bins), 
                                      range(max_bin_id, max_bin_id + length_of_tmp_bins)))
            col_bin_idx_db = dict(zip(tmp_bins['index'].values.tolist(), 
                                    tmp_bins_2['index'].values.tolist()))
            bin_res.append(tmp_bins_2)
            max_bin_id += length_of_tmp_bins
            

            tmp_matrix = (tmp_raw_matrix / cn ).tocoo()
            
            tmp_pixels = pd.Series.sparse.from_coo(tmp_matrix)
            tmp_pixels = (
                            tmp_pixels.reset_index()
                                    .rename(columns={'level_0': 'bin1_id', 
                                                    'level_1': 'bin2_id',
                                                    0: 'count'})
                        )

            def get_col_idx(x):
                try:
                    return col_bin_idx_db[x]
                except KeyError:
                    return x
            
            tmp_pixels['bin1_id'] = tmp_pixels['bin1_id'].parallel_apply(lambda x: row_bin_idx_db[x])
            tmp_pixels['bin2_id'] = tmp_pixels['bin2_id'].parallel_apply(get_col_idx)
            res.append(tmp_pixels) 
        
        ## replace 
        rows, cols = triu(matrix).tocsr()[idxes].nonzero()
        matrix[rows, cols] /= float(cn)
        rows, cols = triu(matrix).tocsr()[:, idxes].nonzero()
        matrix[rows, cols] /= float(cn)

 

    pixels = pd.Series.sparse.from_coo(triu(matrix).tocoo())
    
    pixels = pixels.reset_index()
    pixels.columns = ['bin1_id', 'bin2_id', 'count']
    
    del matrix
    gc.collect()
    new_pixels = pd.concat(res, axis=0)
    
    new_pixels_up = new_pixels[new_pixels['bin1_id'] <= new_pixels['bin2_id']]
    new_pixels_lower = new_pixels[new_pixels['bin1_id'] > new_pixels['bin2_id']]
    new_pixels_lower.columns = ['bin2_id', 'bin1_id', 'count']
    
    new_pixels = pd.concat([pixels, new_pixels_up, new_pixels_lower], axis=0)

    new_pixels = new_pixels.drop_duplicates(subset=['bin1_id', 'bin2_id'])

    bins = pd.concat(bin_res, axis=0)
    bins = bins.drop(['index'], axis=1)

    cooler.create_cooler(output, bins, new_pixels, dtypes={'count': 'float64'})

    
def get_dup_prune_pairs():
    pass
    

    
if __name__ == "__main__":  
    cool = sys.argv[1]
    dup_contigs_path = sys.argv[2]
    output = sys.argv[3]
    convert_matrix_with_dup_contigs(cool, dup_contigs_path, output)
        
        
        