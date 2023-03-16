#!/usr/bin/env python
# -*- coding:utf-8 -*-

import logging
import gc
import os 
import sys

import numpy as np
import pandas as pd


from collections import OrderedDict 
from itertools import permutations
from joblib import Parallel, delayed
from scipy.sparse import hstack
from pathlib import Path
from pprint import pformat
from pyfaidx import Fasta

from .algorithms.hypergraph import (
    HyperGraph,
    IRMM,
    extract_incidence_matrix,
    extract_incidence_matrix2, 
    remove_incidence_matrix
    )
from .utilities import list_flatten, get_contigs
from ._config import *

logger = logging.getLogger(__name__)

class HyperPartition:
    """
    Method of contigs partition based on hypergraph partition.

    Params:
    --------
    edges: list
        Hyperedges
    k: str
        Number of groups. Set to k1:k2 to perform multipartition.

    """
    def __init__(self, edges, 
                    contigsizes,
                    prune=None,
                    min_contacts=3,
                    min_length=10000, 
                    resolution1=0.8,
                    resolution2=0.6,
                    min_scaffold_length=10000,
                    threshold=0.01,
                    max_round=1, 
                    threads=4,
                    chunksize=None
                    ):
        
        self.edges = edges
        self.contigsizes = contigsizes ## dataframe
        self.contig_sizes = self.contigsizes.to_dict()['length'] ## dictionary
        self.prune = prune 
        self.min_contacts = min_contacts
        self.min_length = min_length
        self.resolution1 = resolution1
        self.resolution2 = resolution2
        self.min_scaffold_length = int(min_scaffold_length)
        self.threshold = threshold
        self.max_round = max_round
        self.threads = threads
        self.chunksize = int(chunksize) if chunksize else None

        self.contigs = self.contigsizes.index.values.tolist()
        self.H, self.vertices = self.get_hypergraph()
        
        ## remove edges
        del self.edges, edges 
        gc.collect()

        self.filter_hypergraph()
        # self.NW = self.get_normalize_weight()

        if prune:
            self.P_idx, self.prune_pair_df = self.get_prune_pairs()
        else:
            self.P_idx, self.prune_pair_df = None, None

    @property
    def vertices_idx(self):
        return dict(zip(self.vertices, 
                        range(len(self.vertices))))
    @property
    def idx_to_vertices(self):
        idx_to_vertices = dict(zip(range(len(self.vertices)), self.vertices))
        
        return idx_to_vertices
    
    @property
    def vertices_idx_sizes(self):
        """
        through vertices idx to get contig size
        """
        vertices_idx_sizes = {}
        for i, j in self.vertices_idx.items():
            vertices_idx_sizes[j] = self.contig_sizes[i]
        
        return vertices_idx_sizes

    def get_prune_pairs(self):
        
        vertices_idx = self.vertices_idx

        pair_df = pd.read_csv(self.prune, sep='\t', 
                                header=None, index_col=None)
        pair_df[0] = pair_df[0].map(lambda x: vertices_idx.get(x, np.nan))
        pair_df[1] = pair_df[1].map(lambda x: vertices_idx.get(x, np.nan))
        pair_df = pair_df.dropna(axis=0)
        pair_df2 = pair_df.reindex(columns=[1, 0])
        pair_df2.columns = [0, 1]
        pair_df = pd.concat([pair_df, pair_df2], axis=0)
        pair_df = pair_df.drop_duplicates(subset=[0, 1])
  
        P_idx = [pair_df[0], pair_df[1]]
       
        return P_idx, pair_df
    
    def get_hypergraph(self):
        if self.chunksize:
            from .algorithms.hypergraph import generate_hypergraph
            args = [] 
            pesudo_idx = []
            idx = 0
            for i in range(0, len(self.edges), self.chunksize):
                pesudo_idx.append(i + idx)
                chunk_edges = self.edges[i: i + self.chunksize]
                args.append([self.contigs] + chunk_edges)
                idx += 1
            
            
            res = Parallel(n_jobs=min(self.threads, len(args)))(
                        delayed(generate_hypergraph)(i) for i in args)
            
            _H_list = []
            for _H, _vertices in res:
                order = list(map(_vertices.index, self.contigs))
        
                _H = _H[order] 
                _H_list.append(_H)
            
            H = hstack(_H_list)
            real_idx = np.arange(H.shape[1])
            real_idx = real_idx[np.where(~np.isin(real_idx, pesudo_idx))]
        
            H = H[:, real_idx]
            non_zero_contig_idx = np.where(H.sum(axis=1).T != 0)[-1]
            H = H[non_zero_contig_idx]
            vertices = np.array(self.contigs)[non_zero_contig_idx]

        else:
            HG = HyperGraph(self.edges)
            H = HG.incidence_matrix(self.min_contacts)
            vertices = HG.nodes

            del HG 
            gc.collect()

        logger.info(f"Generated hypergraph that containing {H.shape[0]} vertices"    
                    f" and {H.shape[1]} hyperedges.")

        return H, vertices

    def filter_hypergraph(self):
        ## remove too short contigs
        vertices_idx = self.vertices_idx
        
        short_contigs = self.contigsizes[self.contigsizes < self.min_length]

        short_contig_idx = []
        for i in short_contigs:
            try:
                short_contig_idx.append(vertices_idx[i])
            except KeyError:
                continue
        
        if len(short_contig_idx) == 0:
            return

        logger.info(f"Total {len(short_contig_idx)} contigs were removed, "
                        "because it's length too short.")
        logger.info("start to remove ...")
        self.H, _ = remove_incidence_matrix(self.H, short_contig_idx)
        logger.info("done")
        self.vertices = np.delete(self.vertices, short_contig_idx)

    def get_normalize_weight(self):
        contig_sizes = self.contig_sizes
        contig_sizes_df = pd.DataFrame(contig_sizes, index=['length']).T
        vertices_length = contig_sizes_df.loc[self.vertices]

        a = vertices_length['length'].astype('float32')
        
        NW = (a.max() ** 2) / np.outer(a, a)
 
        return NW

    def single_partition(self):
        logger.info("Starting cluster ...")
        self.K = IRMM(self.H, #self.NW, 
                        self.P_idx,
                        self.resolution1, 
                        self.threshold, self.max_round,
                        threads=self.threads)
        logger.info("Cluster done.")


        # self.K = list(filter(lambda x: len(x) > 1, self.K))
        self.K = list(map(list, self.K))
        # self.K = sorted(self.K, key=lambda x: len(x), reverse=True)
        self.K = self.filter_cluster()
        
        self.K = sorted(self.K, 
                        key=lambda x: self.contigsizes.iloc[x]['length'].sum(), 
                        reverse=True)
        length_contents = pformat(list(map(
            lambda  x: "{:,}".format(self.contigsizes.iloc[x]['length'].sum()), self.K)))
        logger.info(f"Hyperpartition result {len(self.K)} groups:\n"
                    f"{length_contents}")
        
        return self.K  
        
    @staticmethod
    def _multi_partition(k, prune_pair_df, H, contigsizes, #NW, 
                         resolution, threshold, max_round, num):
        """
        single function for multi_partition.
        """
        k = np.array(list(k))
        sub_H, _ = extract_incidence_matrix2(H, k)
        
        del H 
        gc.collect() 
        
        #sub_NW = NW[list(k)][:, list(k)]

        sub_old2new_idx = dict(zip(k, range(len(k))))
        sub_new2old_idx = dict(zip(range(len(k)), k))

        # sub_vertices = list(np.array(self.vertices)[k])
        # sub_vertives_idx = dict(zip(sub_vertices, range(len(sub_vertices))))

        sub_prune_pair_df = prune_pair_df.reindex(list(permutations(k, 2)))
        sub_prune_pair_df = sub_prune_pair_df.dropna().reset_index()
    
        sub_prune_pair_df[0] = sub_prune_pair_df[0].map(lambda x: sub_old2new_idx[x])
        sub_prune_pair_df[1] = sub_prune_pair_df[1].map(lambda x: sub_old2new_idx[x])

        sub_P_idx = [sub_prune_pair_df[0], sub_prune_pair_df[1]]
        
        new_K = IRMM(sub_H, #sub_NW, 
                        sub_P_idx, 
                        resolution, threshold, 
                        max_round, threads=1, 
                        outprefix=num)
        
        ## remove single group
        # new_K = list(filter(lambda x: len(x) > 1, new_K))
        logger.info(f"Cluster Statistics: {list(map(len, new_K))}")
        new_K = list(map(lambda x: list(map(lambda y: sub_new2old_idx[y], x)), new_K))
        new_K = sorted(new_K, key=lambda x: contigsizes.iloc[x]['length'].sum(), reverse=True)

        return new_K

    def multi_partition(self):
        """
        multiple partition for autopolyploid.
        """
        logger.info("Starting first partition ...")
        prune_pair_df = self.prune_pair_df.reset_index().set_index([0, 1])

        self.K = IRMM(self.H, #self.NW, 
                      None, self.resolution1, self.threshold, 
                        self.max_round, threads=self.threads)
        
        self.K = list(map(list, self.K))
        self.K = self.filter_cluster()
        self.to_cluster(f'first.clusters.txt')

        length_contents = pformat(list(map(
            lambda  x: "{:,}".format(self.contigsizes.iloc[x]['length'].sum()), self.K)))
        logger.info(f"First hyperpartition resulted {len(self.K)} groups:\n"
                    f"{length_contents}")
        

        logger.info("Starting second hyperpartition ...")
    
        args = []
        for num, k in enumerate(self.K, 1):
            args.append((k, prune_pair_df, self.H, self.contigsizes,#self.NW, 
                        self.resolution2, self.threshold, self.max_round, num))
            # results.append(HyperPartition._multi_partition(k, prune_pair_df, self.H, #self.NW, 
            #             self.resolution2, self.threshold, self.max_round, num))
            
        results = Parallel(n_jobs=min(self.threads, len(args) + 1))(
                        delayed(HyperPartition._multi_partition)
                                (i, j, k, l, m, n, o, p) 
                                    for i, j, k, l, m, n, o, p in args)

        self.K = list_flatten(results)
        self.K = self.filter_cluster()
        length_contents = pformat(list(map(
            lambda  x: "{:,}".format(self.contigsizes.iloc[x]['length'].sum()), self.K)))
        logger.info(f"Second hyperpartition resulted {len(self.K)} groups: \n"
                    f"{length_contents}")

    def post_check(self):
        pass

    @staticmethod
    def get_k_size(k, contigsizes, idx_to_vertices):
        k = list(map(idx_to_vertices.get, k))
        
        return contigsizes.loc[k].sum().values[0]
    
    def filter_cluster(self):
        logger.info(f"Removed scaffolding less than {self.min_scaffold_length} in length.")
        return list(
            filter(lambda x: HyperPartition.get_k_size(
                    x, self.contigsizes, self.idx_to_vertices) \
                        >= self.min_scaffold_length, 
                    self.K))
        
    def to_cluster(self, output):
        idx_to_vertices = self.idx_to_vertices

        clusters = list(map(lambda y: list(
                        map(lambda x: idx_to_vertices[x], y)), 
                        self.K))
        with open(output, 'w') as out:
            for i, group in enumerate(clusters, 1):
                print(f'group{i}\t{len(group)}\t{" ".join(group)}', 
                        file=out)

        logger.info(f"Successful output hyperpartition results in `{output}`.")
