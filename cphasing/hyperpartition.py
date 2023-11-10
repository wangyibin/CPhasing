#!/usr/bin/env python
# -*- coding:utf-8 -*-

import logging
import gc
import os 
import sys

import numpy as np
import igraph as ig
import pandas as pd


from collections import defaultdict, OrderedDict 
from itertools import (
    permutations, 
    combinations, 
    product, 
    groupby
    )
from joblib import Parallel, delayed
from scipy.sparse import hstack, csr_matrix
from pathlib import Path
from pprint import pformat

from .algorithms.hypergraph import (
    HyperGraph,
    IRMM,
    extract_incidence_matrix2, 
    )
from .core import (
    AlleleTable, 
    ClusterTable,
    PruneTable
)
from .kprune import KPruneHyperGraph
from .utilities import list_flatten
from ._config import *

logger = logging.getLogger(__name__)

class HyperPartition:
    """
    Method of contigs partition based on hypergraph partition.

    Params:
    --------
    edges: HyperEdges
        Hyperedges
    contigsizes: pd.DataFrame
        contig's sizes 
    prune: None (default) or list
        a list contain several contig pairs that should be prune
    allelic_factor: int, default is -1
        factor of the weight of allelic
    allelic_similarity: float, default is 0.8
        the minimum value of similarity
    whitelist: None (default) or list
        a list contain several contigs that only use it to hyperpartition
    blacklist: None (default) or list
        a list contain several contigs that don't use it to hyperpartition
    min_contacts: int
        minimum contacts of contig pairs 
    min_length: int
        minimum length of contigs
    resolution1: float
        the resolution of first partition
    resolution2: float
        the resolution of second partition
    min_scaffold_length: int 
        minimum length of scaffoldings
    threshold: int
        the threshold for IRMM algorithms, default not adapt IRMM
    max_round: int
        max round of IRMM, if set 1, we only partition through louvain algorithms
    threads: int
        number of threads
    chunksize: None (default) or int
        not use

    """
    def __init__(self, edges, 
                    contigsizes,
                    ultra_long=None,
                    ul_weight=1.0,
                    k=None,
                    alleletable=None,
                    prunetable=None,
                    allelic_factor=-1,
                    allelic_similarity=0.8,
                    min_allelic_overlap=0.1,
                    ultra_complex=5.0,
                    whitelist=None,
                    blacklist=None,
                    min_contacts=3,
                    min_length=10000, 
                    resolution1=1.0,
                    resolution2=1.0,
                    min_weight=1.0,
                    min_scaffold_length=10000,
                    threshold=0.01,
                    max_round=1, 
                    threads=4,
                    chunksize=None
                    ):
        
        self.edges = edges
        self.contigsizes = contigsizes ## dataframe
        self.ultra_long = ultra_long
        self.ul_weight = ul_weight
        self.k = k

        self.prunetable = prunetable
        self.alleletable = AlleleTable(alleletable, sort=False, fmt='allele2') if alleletable else None
        self.allelic_factor = allelic_factor
    
        self.allelic_similarity = allelic_similarity
        self.min_allelic_overlap = min_allelic_overlap
        self.ultra_complex = ultra_complex

        self.whitelist = whitelist
        self.blacklist = blacklist
        self.min_contacts = min_contacts
        self.min_length = min_length

        self.resolution1 = resolution1
        self.resolution2 = resolution2
        self.min_weight = min_weight
        self.min_scaffold_length = int(min_scaffold_length)

        self.threshold = threshold
        self.max_round = max_round
        self.threads = threads
        self.chunksize = int(chunksize) if chunksize else None

        self.contig_sizes = self.contigsizes.to_dict()['length'] ## dictionary
        self.contigs = self.contigsizes.index.values.tolist()
        self.HG = HyperGraph(self.edges)
        
        self.K = []
        ## remove edges
        del self.edges, edges 
        gc.collect()

        self.filter_hypergraph()
        
        self.H, self.vertices = self.get_hypergraph()
        # self.NW = self.get_normalize_weight()
        if self.prunetable:
            self.P_allelic_idx, self.P_weak_idx, self.prune_pair_df = self.get_prune_pairs()
        elif self.alleletable:
            self.P_allelic_idx, self.P_weak_idx, self.prune_pair_df = self.get_prune_pairs_kprune()
        else:
            self.P_allelic_idx, self.P_weak_idx, self.prune_pair_df = None, None, None

        if self.P_allelic_idx:
            self.allelic_idx_set = set(map(tuple, 
                            self.prune_pair_df[(self.prune_pair_df['type'] == 0) & 
                                    (self.prune_pair_df['similarity'] >= self.allelic_similarity)]
                            [['contig1', 'contig2']].values)
                            )
        # if allelic_factor is not None:
        #     self.allelic_factor_df = self.prune_pair_df[
        #         (self.prune_pair_df['type'] == 0)][['contig1', 'contig2', 'similarity']]

        #     self.allelic_factor_df['factor'] = np.log10(1 - self.allelic_factor_df['similarity'])
        #     self.allelic_factor_df = self.allelic_factor_df.fillna(0)
            
        #     _allelic_factor = np.ones(shape=(self.H.shape[0], self.H.shape[0]))
        #     _allelic_factor[self.allelic_factor_df['contig1'], self.allelic_factor_df['contig2']] = self.allelic_factor_df['factor']

        #     self.allelic_factor = csr_matrix(_allelic_factor, 
        #                                     shape=(self.H.shape[0], self.H.shape[0]))
        # else:
        #     _allelic_factor = np.ones(shape=(self.H.shape[0], self.H.shape[0]))
        #     self.allelic_factor = csr_matrix(_allelic_factor, dtype=np.int8)

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
            
            H = self.HG.incidence_matrix(self.min_contacts)
            vertices = self.HG.nodes

        logger.info(f"Generated hypergraph that containing {H.shape[0]} vertices"    
                    f" and {H.shape[1]} hyperedges.")

        return H, vertices

    def filter_hypergraph(self):
        """
        remove too short contigs and according the whitelist or blacklist 
        to remove contigs
        """
        short_contigs = self.contigsizes[self.contigsizes < self.min_length]
        short_contigs = short_contigs.dropna().index.tolist()
      
        remove_contigs = set()
        remove_contigs.update(set(short_contigs))

        if self.whitelist:
            remove_contigs.update(set(self.contigs) - set(self.whitelist))

        if self.blacklist:
            remove_contigs.update(set(self.blacklist))
        
        if len(remove_contigs) == 0:
            return

        logger.info(f"Total {len(remove_contigs)} contigs were removed,")
        logger.info(f"\tbecause it's length too short (<{self.min_length}) or your specified.")

        self.HG.remove_rows(remove_contigs)

    def get_prune_pairs_kprune(self):
        """
        get prune information by execute kprune on hypergraph
        """
        vertices_idx = self.vertices_idx
        A = HyperGraph.clique_expansion_init(self.H)
        kph = KPruneHyperGraph(self.alleletable, A, vertices_idx)
        kph.run() 
        kph.prunetable = kph.get_prune_table() 
        
        P_allelic_idx_df = pd.DataFrame(kph.allelic_prune_list)
        P_allelic_idx_df.columns = ['contig1', 'contig2']
        P_weak_idx_df = pd.DataFrame(kph.weak_prune_list)
        P_weak_idx_df.columns = ['contig1', 'contig2']
        
        return [P_allelic_idx_df['contig1'], P_allelic_idx_df['contig2']], \
                [P_weak_idx_df['contig1'], P_weak_idx_df['contig2']], kph.prunetable

    def get_prune_pairs(self):
        
        vertices_idx = self.vertices_idx
        
        prunetable = PruneTable(self.prunetable)
        pair_df = prunetable.data
        pair_df['contig1'] = pair_df['contig1'].map(lambda x: vertices_idx.get(x, np.nan))
        pair_df['contig2'] = pair_df['contig2'].map(lambda x: vertices_idx.get(x, np.nan))
        pair_df = pair_df.dropna(axis=0).astype({'contig1': 'int', 'contig2': 'int', 'type': 'int'})

        prunetable.data = pair_df
        prunetable.symmetric_table()

        pair_df = prunetable.data
        # pair_df = pd.concat([pair_df, pair_df2], axis=0)
        pair_df = pair_df.drop_duplicates(subset=['contig1', 'contig2'])
        pair_df = pair_df.reset_index(drop=True)
        
        # P_idx = [pair_df[0], pair_df[1]]
        # tmp_df = pair_df[(pair_df['type'] == 0)  & (pair_df['similarity'] >= self.allelic_similarity)]
        tmp_df = pair_df[pair_df['type'] == 0]
        P_allelic_idx = [tmp_df['contig1'], tmp_df['contig2']]

        # tmp_df = pair_df[(pair_df['type'] == 1) | \
        #                   ((pair_df['type'] == 0) & (pair_df['similarity'] < self.allelic_similarity))]
        tmp_df = pair_df[pair_df['type'] == 1]
        P_weak_idx = [tmp_df['contig1'], tmp_df['contig2']]

        return P_allelic_idx, P_weak_idx, pair_df

    def get_normalize_weight(self):
        contig_sizes = self.contigsizes
        # contig_sizes_df = pd.DataFrame(contig_sizes, index=['length']).T
        vertices_length = contig_sizes.loc[self.vertices]

        a = vertices_length['length'].astype('float32')
        
        NW = np.log10((a.max() ** 2 ) / np.outer(a, a))
        NW = np.ones(NW.shape)
        # print(NW)
        return NW

    @staticmethod
    def is_error(group, allelic_idx_set, vertices_idx_sizes, max_error_rate=0.2):
        combined_contig_pair = set(combinations(group, 2))
        allelic = combined_contig_pair & allelic_idx_set
        if not allelic:
            return False 
        
        allelic_contigs = list(zip(*list(allelic)))[0]
      
        group_length = vertices_idx_sizes.loc[group].sum().values[0]
        allelic_length = vertices_idx_sizes.loc[list(allelic_contigs)].sum().values[0]

        error_rate = allelic_length / group_length

        if error_rate >= max_error_rate:
            return True 
        else: 
            return False

    def single_partition(self, k=None):
        """
        Single partition

        Params:
        -------
        H: scipy.sparse.csr_matrix
            Hypergraph incidence matrix
        P_allelic_idx: list
            Allelic pairs index
        P_weak_idx: list
            Weak pairs index
        resolution1: float
            Resolution of first partition
        threshold: float
            Threshold of first partition
        max_round: int  
            Max round of first partition
        threads: int    
            Threads of first partition

        Returns:
        --------
        cluster_assignments: list
            Cluster assignments
        K: list
            Cluster index
        
        Examples:
        --------
        >>> cluster_assignments, K = single_partition(H, P_allelic_idx, P_weak_idx,
        """
        logger.info("Start hyperpartition ...")
        A, self.cluster_assignments, self.K = IRMM(self.H, #self.NW, 
                                                self.P_allelic_idx,
                                                self.P_weak_idx,
                                                self.allelic_factor,
                                                self.resolution1, 
                                                self.min_weight,
                                                self.threshold, 
                                                self.max_round,
                                                threads=self.threads)

        # self.K = list(filter(lambda x: len(x) > 1, self.K))
        self.K = list(map(list, self.K))
        # self.K = sorted(self.K, key=lambda x: len(x), reverse=True)
        self.K = self.filter_cluster()
        
        vertices_idx_sizes = self.vertices_idx_sizes
        vertices_idx_sizes = pd.DataFrame(vertices_idx_sizes, index=['length']).T
        
        args = []
        if self.ultra_complex:
            _results = []
            for num, group in enumerate(self.K):
                if len(group) > 1 and HyperPartition.is_error(group, self.allelic_idx_set, vertices_idx_sizes):
                    args.append((group, self.k[0], self.prune_pair_df, self.H, vertices_idx_sizes, self.ultra_complex, 
                                self.min_weight, self.allelic_similarity, self.allelic_factor, self.min_scaffold_length,
                                self.threshold, self.max_round, num))
                else:
                    _results.append(group)

            _results2 = Parallel(n_jobs=min(self.threads, len(args) + 1))(
                            delayed(HyperPartition._incremental_partition) 
                             (i, j, _k, l, m, n, o, p, q, r, s, t, u) 
                                for i, j, _k, l, m, n, o, p, q, r, s, t, u in args) 
            _, _results2 = zip(*_results2)
            _results2 = list_flatten(_results2)
            _results.extend(_results2)

            self.K = _results 

        if k:
            logger.info(f"Merging {len(self.K)} groups into {k} groups ...")
            self.K = HyperPartition._merge(A, self.K, vertices_idx_sizes, k, 
                                            self.prune_pair_df, self.allelic_similarity,
                                            self.min_allelic_overlap)
            logger.info("Merging done.")

        self.K = sorted(self.K, 
                        key=lambda x: vertices_idx_sizes.loc[x]['length'].sum(), 
                        reverse=True)
        
        length_contents = pformat(list(map(
            lambda  x: "{:,}".format(vertices_idx_sizes.loc[x]['length'].sum()), self.K)))
        logger.info(f"Hyperpartition result {len(self.K)} groups:\n"
                    f"{length_contents}")
        
        return self.K  

    @staticmethod
    def _incremental_partition(K, k, prune_pair_df, H, vertices_idx_sizes, #NW, 
                            resolution, min_weight=1, allelic_similarity=0.8, 
                            min_allelic_overlap=0.1, allelic_factor=-1,
                           min_scaffold_length=10000, threshold=0.01, max_round=1, num=None):
        """
        single function for incremental_partition.
        """
        K = np.array(list(K))
        sub_H, _ = extract_incidence_matrix2(H, K)
        
        del H 
        gc.collect() 
        
        #sub_NW = NW[list(k)][:, list(k)]

        sub_old2new_idx = dict(zip(K, range(len(K))))
        sub_new2old_idx = dict(zip(range(len(K)), K))
        
        sub_vertices_idx_sizes = vertices_idx_sizes.loc[K]
        sub_vertices_new_idx_sizes = sub_vertices_idx_sizes
        sub_vertices_new_idx_sizes.index = sub_vertices_idx_sizes.index.map(lambda x: sub_old2new_idx[x])

        # sub_vertices = list(np.array(self.vertices)[k])
        # sub_vertives_idx = dict(zip(sub_vertices, range(len(sub_vertices))))
        if prune_pair_df is not None:
            sub_prune_pair_df = prune_pair_df.reindex(list(permutations(K, 2)))
            sub_prune_pair_df = sub_prune_pair_df.dropna().reset_index()
        
            sub_prune_pair_df['contig1'] = sub_prune_pair_df['contig1'].map(lambda x: sub_old2new_idx[x])
            sub_prune_pair_df['contig2'] = sub_prune_pair_df['contig2'].map(lambda x: sub_old2new_idx[x])
            
            tmp_df = sub_prune_pair_df[sub_prune_pair_df['type'] == 0]
            sub_P_allelic_idx = [tmp_df['contig1'], tmp_df['contig2']]
            tmp_df = sub_prune_pair_df[sub_prune_pair_df['type'] == 1]
            sub_P_weak_idx = [tmp_df['contig1'], tmp_df['contig2']]
            
            sub_prune_pair_df = sub_prune_pair_df.reset_index()
        else:
            sub_prune_pair_df = None
            sub_P_allelic_idx = None
            sub_P_weak_idx = None
        
        # sub_allelic_factor_df = sub_prune_pair_df[
        #     (sub_prune_pair_df['type'] == 0)][['contig1', 'contig2', 'similarity']]
        
        # sub_allelic_factor_df['factor'] = np.log10(1 - sub_allelic_factor_df['similarity'])
        # sub_allelic_factor_df = sub_allelic_factor_df.fillna(0)
        
        # _sub_allelic_factor = np.ones(shape=(sub_H.shape[0], sub_H.shape[0]))
        # _sub_allelic_factor[
        #     sub_allelic_factor_df['contig1'], sub_allelic_factor_df['contig2']
        #             ] = sub_allelic_factor_df['factor']

        # sub_allelic_factor = csr_matrix(_sub_allelic_factor)

        A, cluster_assignments, new_K = IRMM(sub_H, #sub_NW, 
                                            sub_P_allelic_idx, 
                                            sub_P_weak_idx,
                                            allelic_factor,
                                            resolution, 
                                            min_weight,
                                            threshold, 
                                            max_round, 
                                            threads=1, 
                                            outprefix=num)

        ## remove the scaffold that is too short
        new_K = list(filter(
                        lambda x: sub_vertices_new_idx_sizes.loc[x].sum().values[0] \
                            >= min_scaffold_length, new_K)
        )
        new_K = list(map(list, new_K))
        if k:
            new_K = HyperPartition._merge(A, new_K, sub_vertices_new_idx_sizes, k, 
                                            sub_prune_pair_df, allelic_similarity, 
                                            min_allelic_overlap)
        
        new_K = list(map(lambda x: list(map(lambda y: sub_new2old_idx[y], x)), new_K))
        new_K = sorted(new_K, key=lambda x: vertices_idx_sizes.loc[x]['length'].sum(), reverse=True)

        return cluster_assignments, new_K

    def incremental_partition(self, k, first_cluster=None):
        """
        incremental partition for autopolyploid.
        """
        logger.info("Starting first partition ...")
        if self.prunetable:
            prune_pair_df = self.prune_pair_df.reset_index().set_index(['contig1', 'contig2'])
        else:
            prune_pair_df = None

        vertices_idx_sizes = self.vertices_idx_sizes
        vertices_idx_sizes = pd.DataFrame(vertices_idx_sizes, index=['length']).T

        if not first_cluster:
            A, _, self.K = IRMM(self.H, #self.NW, 
                        None, None, self.allelic_factor, self.resolution1, 
                            self.min_weight, self.threshold, 
                            self.max_round, threads=self.threads)

            self.K = list(map(list, self.K))
            self.K = self.filter_cluster()
            
            if k[0]:
                self.K = HyperPartition._merge(A, self.K, vertices_idx_sizes, k[0])

            self.K = sorted(self.K, key=lambda x: vertices_idx_sizes.loc[x]['length'].sum(), reverse=True)
            self.to_cluster(f'first.clusters.txt')

            length_contents = pformat(list(map(
                lambda  x: "{:,}".format(vertices_idx_sizes.loc[x]['length'].sum()), self.K)))
            logger.info(f"First hyperpartition resulted {len(self.K)} groups:\n"
                        f"{length_contents}")

        else:
            vertices_idx = self.vertices_idx
            first_cluster = list(map(lambda y: list(filter(lambda x: x in vertices_idx, y)), first_cluster))
            for sub_group in first_cluster:
                self.K.append(list(map(lambda x: vertices_idx[x], sub_group)))

                
            logger.info(f"Load the first results from exists file.")

        logger.info("Starting second hyperpartition ...")

        # prune_pair_df = prune_pair_df.reset_index()
        args = []
        for num, sub_k in enumerate(self.K, 1):
            args.append((sub_k, k[1], prune_pair_df, self.H, vertices_idx_sizes,#self.NW, 
                        self.resolution2, self.min_weight, self.allelic_similarity, 
                        self.min_allelic_overlap, self.allelic_factor, 
                        self.min_scaffold_length, self.threshold, self.max_round, num))
            # results.append(HyperPartition._incremental_partition(k, prune_pair_df, self.H, #self.NW, 
            #             self.resolution2, self.threshold, self.max_round, num))
            
        results = Parallel(n_jobs=min(self.threads, len(args) + 1))(
                        delayed(HyperPartition._incremental_partition)
                                (i, j, _k, l, m, n, o, p, q, r, s, t, u, v) 
                                    for i, j, _k, l, m, n, o, p, q, r, s, t, u, v in args)
        
        self.cluster_assignments, results = zip(*results)
        self.inc_chr_idx = []
        for i, j in enumerate(results):
            for _ in j:
                self.inc_chr_idx.append(i)

        self.K = list_flatten(results)
        self.K = self.filter_cluster()

        length_contents = pformat(list(map(
            lambda  x: "{:,}".format(vertices_idx_sizes.loc[x]['length'].sum()), self.K)))
        logger.info(f"Second hyperpartition resulted {len(self.K)} groups: \n"
                    f"{length_contents}")


    @staticmethod
    def _merge(A, K, vertices_idx_sizes, k=None, 
                prune_pair_df=None, allelic_similarity=0.85,
                min_allelic_overlap=0.1):
        if not k:
            return K 
        
        if prune_pair_df is not None:
            allelic_idx_set = set(map(tuple, 
                            prune_pair_df[(prune_pair_df['type'] == 0) & 
                                                (prune_pair_df['similarity'] >= allelic_similarity)]
                            [['contig1', 'contig2']].values)
                            )

        iter_round = 0 
        flag = 1
      
        while k and len(K) > k:
            if iter_round > (len(K) - k + 10):
                break
            
            current_group_number = len(K)
            value_matrix = np.zeros(shape=(current_group_number, current_group_number))
            flag_matrix = np.ones(shape=(current_group_number, current_group_number))
            res = {}
            for i, group1 in enumerate(K):
                group1 = list(group1)
                for j in range(i + 1, current_group_number):
                    group2 = list(K[j])

                    group1_length = vertices_idx_sizes.loc[list(group1)].sum().values[0]
                    group2_length = vertices_idx_sizes.loc[list(group2)].sum().values[0]
                    
                    if prune_pair_df is not None:
                        product_contig_pair = set(product(group1, group2))
                        allelic = product_contig_pair & allelic_idx_set
                        if allelic:
                            tmp1, tmp2 = list(map(set, zip(*allelic)))
                            tmp1_length = vertices_idx_sizes.loc[list(tmp1)].sum().values[0]
                            tmp2_length = vertices_idx_sizes.loc[list(tmp2)].sum().values[0]
                            overlap1 = tmp1_length / group1_length
                            overlap2 = tmp2_length / group2_length

                            if overlap1 > min_allelic_overlap or overlap2 > min_allelic_overlap:
                                flag = 0
                        
                    value = A[group1, ][:,group2 ].sum() 
                    value_matrix[i, j] = value
                    flag_matrix[i, j] = flag
                    flag = 1

            total_value = value_matrix.sum()
        
            value_matrix = value_matrix + value_matrix.T - np.diag(value_matrix.diagonal())
         
            for i in range(len(K)):
                i_value = value_matrix[i].sum()
                for j in range(i+1, len(K)):
                    j_value = value_matrix[j].sum()
                    value =  value_matrix[i, j]
                    flag = flag_matrix[i, j]
                    if flag:
                        Q = value - (j_value * i_value) / total_value
                    else:
                        Q = - 2**64

                    res[(i, j)] = Q
            
            i1, i2 = max(res,  key=lambda x: res[x])

            if max(res) == 0 or max(res) == 2**64:
                continue
            
            group1 = K[i1]
            group2 = K[i2]
    
            group1.extend(group2)
            
            K[i1] = group1 
            K.pop(i2)

            iter_round += 1
        
        return K

    def merge_cluster(self, merge_cluster):
        """
        merge several clusters into k groups
        """
        vertices_idx = self.vertices_idx
        vertices_idx_sizes = self.vertices_idx_sizes
        vertices_idx_sizes = pd.DataFrame(vertices_idx_sizes, index=['length']).T
        A = self.HG.clique_expansion_init(self.H, self.P_allelic_idx, self.P_weak_idx, 
                                          self.allelic_factor, self.min_weight)
        _K = list(merge_cluster.data.values())
        self.K = list(list(map(lambda x: vertices_idx[x], i)) for i in _K)
        self.K = sorted(self.K, 
                        key=lambda x: vertices_idx_sizes.loc[x]['length'].sum(), 
                        reverse=True)
        
        self.K = HyperPartition._merge(A, self.K, vertices_idx_sizes, self.k[0],
                                        self.prune_pair_df, self.allelic_similarity,
                                        self.min_allelic_overlap)

        self.K = sorted(self.K, 
                        key=lambda x: vertices_idx_sizes.loc[x]['length'].sum(), 
                        reverse=True)
        
        length_contents = pformat(list(map(
            lambda  x: "{:,}".format(vertices_idx_sizes.loc[x]['length'].sum()), self.K)))
        logger.info(f"Merge cluster result {len(self.K)} groups:\n"
                    f"{length_contents}")
        

        return self.K  
        

    def post_check(self):
        """
        check the partition results remove the allelic 
            mis-assembly and add it to correctly groups
        """
        cluster_assignment = self.cluster_assignments
        allelic_idx_set = set(map(tuple, 
                        self.prune_pair_df[(self.prune_pair_df['type'] == 0) & 
                                            (self.prune_pair_df['similarity'] >= 0.85)]
                        [['contig1', 'contig2']].values)
                        )
        idx_to_vertices = self.idx_to_vertices
        graph = cluster_assignment.graph
      
        for i, group in enumerate(self.K):
            res = set(combinations(group, 2)).intersection(allelic_idx_set)
            res_flatten = set(list_flatten(res))

            print(len(res_flatten))
            _, new_K = HyperPartition._incremental_partition(res_flatten, 
                                                          self.prune_pair_df, 
                                                          self.H, self.contigsizes,
                                                          self.resolution1, num=None)
            group = np.array(group)
            
            modularity_results = []
            for new_group in new_K:
                tmp_modularity = []
                for j in range(len(self.K)):
                    new_membership = np.array(cluster_assignment.membership)
                    new_membership[group[np.isin(group, new_group)]] = j
                    tmp_modularity.append(graph.modularity(new_membership))

                modularity_results.append(np.array(tmp_modularity))

            modularity_results = np.array(modularity_results)
            drop_groups = []
            if all(modularity_results[:, i] == np.max(modularity_results[:, i]) ):
                print("equal")
                continue
            
            retain_idx = np.argmax(modularity_results[:, i])
            
            for n, new_group in enumerate(new_K):
                if n != retain_idx:
                    other_indexes = list(range(len(self.K)))
                    other_indexes.pop(i)
                    modularity_results[n][i] = 0 
                    rescue_to_idx = np.argmax(modularity_results[n])
                    print(modularity_results[n])
                    print(rescue_to_idx)
                    rescued_group = self.K[rescue_to_idx]
                    rescued_group.extend(new_group)
                    self.K[rescue_to_idx] =  rescued_group

                    retain_group = self.K[i]
                    print(len(retain_group), len(new_group))
                    retain_group = [contig for contig in retain_group if contig not in new_group]
                    print(len(retain_group)), len(new_group)
                    self.K[i] = retain_group
                # self.K[i] = list(filter(lambda x: x not in new_group, self.K[i])) 
                # print(list(idx_to_vertices[j] for j in new_group))
            # clusters = list(map(lambda y: list(
            #             map(lambda x: idx_to_vertices[x], y)), new_K))

        print(list(map(len, self.K)))
            

    @staticmethod
    def get_k_size(k, contigsizes, idx_to_vertices):
        k = list(map(idx_to_vertices.get, k))
        
        return contigsizes.loc[k].sum().values[0]
    
    def filter_cluster(self):
        logger.info(f"Removed scaffolding less than {self.min_scaffold_length} in length.")

        try: 
            self.inc_chr_idx 
            pop_idx = []
            _K = []
            for i, group in enumerate(self.K):
                _size = HyperPartition.get_k_size(
                            group, self.contigsizes, self.idx_to_vertices)
                if _size >= self.min_scaffold_length:
                    _K.append(group)
                else:
                    pop_idx.append(i)
            
            _new_inc_chr_idx = []
            for i, idx in enumerate(self.inc_chr_idx):
                if i not in pop_idx:
                    _new_inc_chr_idx.append(idx)
            self.inc_chr_idx = _new_inc_chr_idx

        except AttributeError:
            _K = list(
                filter(lambda x: HyperPartition.get_k_size(
                        x, self.contigsizes, self.idx_to_vertices) \
                            >= self.min_scaffold_length, 
                        self.K))
        
        return _K 
        
    def to_cluster(self, output):
        idx_to_vertices = self.idx_to_vertices

        clusters = list(map(lambda y: list(
                        map(lambda x: idx_to_vertices[x], y)), 
                        self.K))
        with open(output, 'w') as out:
            try:
                db = defaultdict(lambda :1)
                for i, chr_idx in enumerate(self.inc_chr_idx):
                    group = clusters[i]
                    
                    hap_idx = db[chr_idx]
                    print(f'Chr{chr_idx + 1:0>2}g{hap_idx}\t{len(group)}\t{" ".join(group)}', 
                        file=out)
                    db[chr_idx] += 1

            except AttributeError:
                for i, group in enumerate(clusters, 1):
                    print(f'Chr{i:0>2}\t{len(group)}\t{" ".join(group)}', 
                        file=out)

        logger.info(f"Successful output hyperpartition results in `{output}`.")
