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

try:
    from .algorithms.hypergraph import extract_incidence_matrix2 
    from .utilities import run_cmd, list_flatten
except ImportError:
    from cphasing.utilities import run_cmd, list_flatten

logger = logging.getLogger(__name__)


class CollapseContigs:
    def __init__(self, cool_path):
        self.cool_path = cool_path 
        self.cool = cooler.Cooler(self.cool_path)
        
        pass 

class CollapsedRescue:
    """
    Rescue the collapsed contigs into a group
    """
    def __init__(self, HG, clustertable, alleletable,
                    collapsed_contigs, allelic_similarity: float=.85):
        
        self.HG = HG 
        self.clustertable = clustertable
        self.alleletable = alleletable 
        self.alleletable.data = self.alleletable.data[
            self.alleletable.data['similarity'] >= allelic_similarity]
    
        self.collapsed_contigs = collapsed_contigs
        self.allelic_similarity = allelic_similarity
        
        self.hap_groups = self.clustertable.hap_groups
        self.H = HG.incidence_matrix()
        self.vertices = self.HG.nodes 

    @property
    def vertices_idx(self):
        return dict(zip(self.vertices, 
                        range(len(self.vertices))))
    
    @property
    def idx_to_vertices(self):
        idx_to_vertices = dict(zip(range(len(self.vertices)), self.vertices))
        
        return idx_to_vertices

    @staticmethod
    def _rescue(A, ):
        pass

    def rescue(self):
        vertices_idx = self.vertices_idx
        idx_to_vertices = self.idx_to_vertices
        new_cluster_data = {}
        for k, hap_group in enumerate(self.hap_groups):
            groups = self.hap_groups[hap_group]
            if len(groups) < 2:
                continue 
            contigs = list_flatten(groups)
            contigs_idx = list(map(vertices_idx.get, contigs))
            sub_old2new_idx = dict(zip(contigs_idx, range(len(contigs_idx))))
            sub_alleletable = self.alleletable.data[
                self.alleletable.data[1].isin(contigs) & self.alleletable.data[2].isin(contigs)]
            sub_alleletable[1] = sub_alleletable[1].map(vertices_idx.get).map(sub_old2new_idx.get)
            sub_alleletable[2] = sub_alleletable[2].map(vertices_idx.get).map(sub_old2new_idx.get)
            P_allelic_idx = [sub_alleletable[1], sub_alleletable[2]]
            sub_alleletable.set_index([1], inplace=True)

            sub_H, _ = extract_incidence_matrix2(self.H, contigs_idx)
            sub_A = self.HG.clique_expansion_init(sub_H, P_allelic_idx=P_allelic_idx, allelic_factor=0)
           
            # sub_allelic = set(map(tuple, sub_alleletable[[1, 2]].values.tolist()))

            groups_idx = list(map(lambda x: list(map(vertices_idx.get, x)), groups))
            groups_new_idx = list(map(lambda x: list(map(sub_old2new_idx.get, x)), groups_idx))
            groups_new_idx_db = {}
            for i in range(len(groups_new_idx)):
                groups_new_idx_db.update(dict(zip(groups_new_idx[i], [i] * len(groups_new_idx[i]))))
            
            sub_collapsed_contigs = self.collapsed_contigs.reindex(contigs).dropna(axis=0)
            sub_collapsed_contigs_idx = list(map(vertices_idx.get, sub_collapsed_contigs.index.tolist()))
            sub_collapsed_contigs_idx_new = list(map(sub_old2new_idx.get, sub_collapsed_contigs_idx))

            res = []
           

            for j in sub_collapsed_contigs_idx_new:
                try:
                    tmp_allelic_table = sub_alleletable.loc[j]
                    try:
                        tmp_allelic_table.set_index([2], inplace=True)
                    except:
                        tmp_allelic_table.to_frame().set_index([2], inplace=True)

                except KeyError:
                    tmp_allelic_table = []
                shared_similarity = np.zeros(len(groups))
                tmp_res = []
                for i in range(len(groups)):
                    
                    if i == groups_new_idx_db[j]:
                        tmp_res.append(0)
                    else:
                        if len(tmp_allelic_table) > 1:
                            tmp = tmp_allelic_table.reindex(groups_new_idx[i]).dropna()
                            if len(tmp) > 1:
                                shared_similarity[i] = tmp['mzShared'].sum()

                        tmp_res.append(sub_A[j, groups_new_idx[i]].mean())

                groups_new_idx[np.argmax(tmp_res)].append(j)
                
                res.append(tmp_res)
            
            
            new_groups = []
            for group_idx in groups_new_idx:
                tmp = list(map(lambda x: contigs_idx[x], group_idx))
                tmp = list(map(idx_to_vertices.get, tmp))
                new_groups.append(tmp)

            for k, group in enumerate(new_groups):
                new_cluster_data[f'{hap_group}g{k+1}'] = group 
            
        
        self.clustertable.data = new_cluster_data
        self.clustertable.save("test.clusters.txt")

    

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
        
        
        