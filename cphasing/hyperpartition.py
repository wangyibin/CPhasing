#!/usr/bin/env python
# -*- coding:utf-8 -*-

import logging
import os 
import sys

import numpy as np
import pandas as pd


from collections import OrderedDict 
from itertools import permutations
from joblib import Parallel, delayed
from scipy.sparse import hstack
from pathlib import Path
from pyfaidx import Fasta

from .algorithms.hypergraph import (
    IRMM,
    extract_incidence_matrix, 
    remove_incidence_matrix
    )
from .utilities import listify, list_flatten
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
                    fasta,
                    prune=None,
                    min_length=10000, 
                    resolution1=0.8,
                    resolution2=0.6,
                    threshold=0.01,
                    max_round=10, 
                    threads=4,
                    chunksize=400000
                    ):
        
        self.edges = edges

        self.fasta = fasta
        self.prune = prune 

        self.min_length = min_length
        self.resolution1 = resolution1
        self.resolution2 = resolution2
        self.threshold = threshold
        self.max_round = max_round
        self.threads = threads
        self.chunksize = int(chunksize) 

        self.contigs = self.get_contigs()
        self.H, self.vertices = self.get_hypergraph()
        
        self.filter_hypergraph()
        self.NW = self.get_normalize_weight()

        if prune:
            self.P_idx, self.prune_pair_df = self.get_prune_pairs()
        else:
            self.P_idx, self.prune_pair_df = None, None

    @property
    def vertices_idx(self):
        return dict(zip(self.vertices, 
                        range(len(self.vertices))))
    @property
    def contig_sizes(self):
        fasta = Fasta(self.fasta)
        contig_size_db = OrderedDict(list(map(lambda x: (x.name, len(x)), list(fasta))))
        
        return contig_size_db 

    def get_contigs(self):
        fasta = Fasta(self.fasta)
        contigs = list(map(lambda x: x.name, fasta))
        lengths = list(map(len, list(fasta)))
        length_db = dict(zip(contigs, lengths))
        contigs = sorted(contigs, key=lambda x: length_db[x], reverse=True)

        return contigs

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
        from .algorithms.hypergraph import generate_hypergraph
        if self.chunksize:
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
            H, vertices = generate_hypergraph(self.edges)

        logger.info(f"Generated hypergraph that containing {H.shape[0]} vertices"    
                    f" and {H.shape[1]} hyperedges.")

        return H, vertices

    def filter_hypergraph(self):
        ## remove too short contigs
        contig_sizes = self.contig_sizes
        vertices_idx = self.vertices_idx
        short_contigs = [i for i in contig_sizes
                                if contig_sizes[i] < self.min_length]

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
        self.H, _ = remove_incidence_matrix(self.H, short_contig_idx)
        self.vertices = np.delete(self.vertices, short_contig_idx)

    def get_normalize_weight(self):
        contig_sizes = self.contig_sizes
        contig_sizes_df = pd.DataFrame(contig_sizes, index=['length']).T
        vertices_length = contig_sizes_df.loc[self.vertices]

        a = vertices_length['length'].astype('float32')
        
        NW = (a.max() ** 2) / np.outer(a, a)
 
        return NW

    
    @classmethod
    def _multi_partition(self, k, prune_pair_df, H, NW, resolution, threshold, max_round, num):
        """
        single function for multi_partition.
        """
        k = np.array(list(k))
        sub_H, _ = extract_incidence_matrix(H, k)
        sub_NW = NW[list(k)][:, list(k)]

        sub_old2new_idx = dict(zip(k, range(len(k))))
        sub_new2old_idx = dict(zip(range(len(k)), k))

        # sub_vertices = list(np.array(self.vertices)[k])
        # sub_vertives_idx = dict(zip(sub_vertices, range(len(sub_vertices))))

        sub_prune_pair_df = prune_pair_df.reindex(list(permutations(k, 2)))
        sub_prune_pair_df = sub_prune_pair_df.dropna().reset_index()
    
        sub_prune_pair_df[0] = sub_prune_pair_df[0].map(lambda x: sub_old2new_idx[x])
        sub_prune_pair_df[1] = sub_prune_pair_df[1].map(lambda x: sub_old2new_idx[x])

        sub_P_idx = [sub_prune_pair_df[0], sub_prune_pair_df[1]]
        
        new_K = IRMM(sub_H, sub_NW, 
                        sub_P_idx, 
                        resolution, threshold, 
                        max_round, threads=1, 
                        outprefix=num)
        
        ## remove single group
        new_K = filter(lambda x: len(x) > 1, new_K)
        new_K = list(map(lambda x: list(map(lambda y: sub_new2old_idx[y], x)), new_K))
        
        return new_K

    def multi_partition(self):
        """
        multiple partition for autopolyploid.
        """
        prune_pair_df = self.prune_pair_df.reset_index().set_index([0, 1])

        self.K = IRMM(self.H, self.NW, None, 
                        self.resolution1, self.threshold, 
                        self.max_round, threads=self.threads)
        self.K = filter(lambda x: len(x) > 1, self.K)

        logger.info("multi-partition ...")
        # results = []
        args = []
        for num, k in enumerate(self.K, 1):
            args.append((k, prune_pair_df, self.H, self.NW, 
                        self.resolution2, self.threshold, self.max_round, num))
           
        results = Parallel(n_jobs=min(self.threads, len(args)))(
                        delayed(self._multi_partition)
                                (i, j, k, l, m, n, o, p) 
                                    for i, j, k, l, m, n, o, p in args)

        self.K = list_flatten(results)
        logger.info(f"Multi partition results: {list(map(len, self.K))}")

    def single_partition(self):
        
        logger.info("Starting to cluster ...")
        self.K = IRMM(self.H, self.NW, self.P_idx,
                        self.resolution1, 
                        self.threshold, self.max_round,
                        threads=self.threads)
        logger.info("Cluster done.")


        self.K = filter(lambda x: len(x) > 1, self.K)
        self.K = sorted(self.K, key=lambda x: len(x), reverse=True)
        
        return self.K

    def to_cluster(self, output):
        idx_to_vertices = dict(zip(range(len(self.vertices)), self.vertices))

        with open('idx_to_vertices.list', 'w') as out:
            for idx, contig in idx_to_vertices.items():
                print(idx, contig, sep='\t', file=out)
                
        clusters = list(map(lambda y: list(
                        map(lambda x: idx_to_vertices[x], y)), 
                        self.K))
        with open(output, 'w') as out:
            for i, group in enumerate(clusters, 1):
                print(f'group{i}\t{len(group)}\t{" ".join(group)}', 
                        file=out)

        logger.info(f"Successful output cluster results in `{output}`.")
