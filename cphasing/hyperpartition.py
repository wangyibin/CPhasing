#!/usr/bin/env python
# -*- coding:utf-8 -*-

import logging
import gc
import os 
import sys

import math
import numpy as np
import igraph as ig
import pandas as pd


from collections import defaultdict, OrderedDict 
from itertools import (
    permutations, 
    combinations, 
    product, 
    groupby,
    cycle
    )
from joblib import Parallel, delayed, parallel_backend
from scipy.sparse import hstack, csr_matrix
from pathlib import Path
from pprint import pformat
from rich.console import Console
from rich.table import Table

from . import console_html 
from .algorithms.hypergraph import (
    HyperGraph,
    IRMM,
    extract_incidence_matrix2, 
    )
from .algorithms.scaffolding import (
    raw_sort
)
from .core import (
    AlleleTable, 
    ClusterTable,
    PruneTable
)
from .kprune import KPruneHyperGraph
from .utilities import list_flatten, run_cmd, to_humanized2
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
                    fasta=None,
                    alleles_kmer_size=19,
                    alleles_window_size=19,
                    alleles_minimum_similarity=0.5,
                    alleles_diff_thres=0.2,
                    alleles_trim_length=25000,
                    alleletable=None,
                    prunetable=None,
                    normalize=False,
                    contacts=None,
                    kprune_norm_method="auto",
                    allelic_factor=-1,
                    cross_allelic_factor=0.3,
                    allelic_similarity=0.8,
                    min_allelic_overlap=0.1,
                    ultra_complex=5.0,
                    disable_merge_in_first=False,
                    exclude_group_to_second=None,
                    exclude_group_from_first=None,
                    whitelist=None,
                    blacklist=None,
                    min_contacts=3,
                    min_length=10000, 
                    resolution1=1.0,
                    resolution2=1.0,
                    init_resolution1=0.8,
                    init_resolution2=0.8,
                    min_weight=1.0,
                    min_cis_weight=1.0,
                    min_quality1=1,
                    min_quality2=2,
                    mapq_filter_1_to_2=False,
                    min_scaffold_length=10000,
                    is_remove_misassembly=False,
                    is_recluster_contigs=False,
                    refine=False,
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
        
        self.fasta = fasta
        self.alleles_kmer_size = alleles_kmer_size
        self.alleles_window_size = alleles_window_size
        self.alleles_minimum_similarity = alleles_minimum_similarity
        self.alleles_diff_thres = alleles_diff_thres
        self.alleles_trim_length = alleles_trim_length

        self.prunetable = prunetable
        # self.alleletable = AlleleTable(alleletable, sort=False, fmt='allele2') if alleletable else None
        if alleletable:
            self.alleletable = str(Path(alleletable).absolute())
        else:
            self.alleletable = alleletable

        self.normalize = normalize
        self.contacts = contacts
        self.kprune_norm_method = kprune_norm_method
            
        self.allelic_factor = allelic_factor
        self.cross_allelic_factor = cross_allelic_factor
    
        self.allelic_similarity = allelic_similarity
        self.min_allelic_overlap = min_allelic_overlap
        self.ultra_complex = ultra_complex
        self.disable_merge_in_first = disable_merge_in_first
        self.exclude_group_to_second = exclude_group_to_second
        self.exclude_group_from_first = exclude_group_from_first
        self.whitelist = whitelist
        self.blacklist = blacklist
        self.min_contacts = min_contacts
        self.min_length = min_length

        self.resolution1 = resolution1
        self.resolution2 = resolution2
        self.init_resolution1 = init_resolution1
        self.init_resolution2 = init_resolution2
        self.min_weight = min_weight
        self.min_cis_weight = min_cis_weight
        self.min_quality1 = min_quality1
        self.min_quality2 = min_quality2
        self.min_scaffold_length = int(min_scaffold_length)
        self.is_remove_misassembly = is_remove_misassembly
        self.is_recluster_contigs = is_recluster_contigs
        self.refine = refine
        logger.debug(f"is remove misassembly: {self.is_remove_misassembly}")
        self.threshold = threshold
        self.max_round = max_round
        self.threads = threads
        self.chunksize = int(chunksize) if chunksize else None
        
        self.log_dir = "logs"
        Path(self.log_dir).mkdir(exist_ok=True)

        self.contig_sizes = self.contigsizes.to_dict()['length'] ## dictionary
        self.contigs = self.contigsizes.index.values.tolist()
        self.K = []

        if not mapq_filter_1_to_2:
            logger.debug(f"min_quality1: {self.min_quality1}")
            self.HG = HyperGraph(self.edges, min_quality=self.min_quality1)
            # if self.min_quality2 > self.min_quality1:
            #     logger.setLevel(logging.WARNING)
            #     tmp_HG = HyperGraph(self.edges, min_quality=self.min_quality2)
            #     tmp_HG.incidence_matrix(min_contacts=self.min_contacts)
            #     remove_contigs = set(tmp_HG.remove_contigs)
            #     logger.setLevel(logging.INFO)
            #     del tmp_HG
            # else:
            #     remove_contigs = set()
        else:
            logger.debug(f"min_quality2: {self.min_quality2}")
            self.HG = HyperGraph(self.edges, min_quality=self.min_quality2)
            # remove_contigs = set()

        ## remove edges
        del edges 
        gc.collect()

        self.filter_hypergraph()
        self.H, self.vertices = self.get_hypergraph()
        if self.normalize:
            self.NW = self.get_normalize_weight()
        else:
            self.NW = None

        # if self.HG.mapq.size == 0:
        #     del self.HG.edges 
        #     gc.collect()

        if self.prunetable:
            self.P_allelic_idx, self.P_weak_idx, self.prune_pair_df = self.get_prune_pairs()

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
        vertices_idx_sizes = {j: self.contig_sizes[i] 
                                for i, j in self.vertices_idx.items()}
        
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

        logger.info(f"Generated filtered hypergraph that containing {H.shape[0]:,} vertices"    
                    f" and {H.shape[1]:,} hyperedges.")

        return H, vertices

    def filter_hypergraph(self, remove_contigs=set()):
        """
        remove too short contigs and according the whitelist or blacklist 
        to remove contigs
        """
        short_contigs = self.contigsizes[self.contigsizes < self.min_length]
        short_contigs = short_contigs.dropna().index.tolist()

        remove_contigs.update(set(short_contigs))

        if self.whitelist:
            remove_contigs.update(set(self.contigs) - set(self.whitelist))

        if self.blacklist:
            remove_contigs.update(set(self.blacklist))
        
        if len(remove_contigs) == 0:
            return
    
        logger.info(f"Total {len(remove_contigs):,} contigs were removed,")
        logger.info(f"\tbecause it's length too short (<{self.min_length}) "
                        "or your specified in blacklist or not in whitelist.")

        self.HG.remove_rows(remove_contigs)

    def get_prune_pairs_kprune(self):
        """
        get prune information by execute kprune on hypergraph
        """
        vertices_idx = self.vertices_idx
        if self.normalize:
            NW = self.get_normalize_weight()
        else:
            NW = None

        A = HyperGraph.clique_expansion_init(self.H, NW=NW)
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
        # tmp_df = pair_df[pair_df['type'] == 0]
        # P_allelic_idx = [tmp_df['contig1'], tmp_df['contig2']]

        P_allelic_idx = [pair_df.loc[pair_df['type'] == 0, 'contig1'], pair_df.loc[pair_df['type'] == 0, 'contig2']]
        
        # tmp_df = pair_df[(pair_df['type'] == 1) | \
        #                   ((pair_df['type'] == 0) & (pair_df['similarity'] < self.allelic_similarity))]
        # tmp_df = pair_df[pair_df['type'] == 1]
        # P_weak_idx = [tmp_df['contig1'], tmp_df['contig2']]

        P_weak_idx = [pair_df.loc[pair_df['type'] == 1, 'contig1'], pair_df.loc[pair_df['type'] == 1, 'contig2']]

        return P_allelic_idx, P_weak_idx, pair_df

    @staticmethod
    def get_prune_pairs_func(prunetable, vertices_idx):
        
        prunetable = PruneTable(prunetable)
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
        # tmp_df = pair_df[pair_df['type'] == 0]
        # P_allelic_idx = [tmp_df['contig1'], tmp_df['contig2']]

        P_allelic_idx = [pair_df.loc[pair_df['type'] == 0, 'contig1'], pair_df.loc[pair_df['type'] == 0, 'contig2']]
        
        # tmp_df = pair_df[(pair_df['type'] == 1) | \
        #                   ((pair_df['type'] == 0) & (pair_df['similarity'] < self.allelic_similarity))]
        # tmp_df = pair_df[pair_df['type'] == 1]
        # P_weak_idx = [tmp_df['contig1'], tmp_df['contig2']]

        P_weak_idx = [pair_df.loc[pair_df['type'] == 1, 'contig1'], pair_df.loc[pair_df['type'] == 1, 'contig2']]

        return P_allelic_idx, P_weak_idx, pair_df

    def get_normalize_weight(self):
        contig_sizes = self.contigsizes
       
        vertices_length = contig_sizes.loc[self.vertices]

        a = vertices_length['length'].astype('float32')
        
        # NW = np.log10((a.max() ** 2 ) / np.outer(a, a))
        NW = 1 / np.outer(a, a)
        # NW = np.ones(NW.shape)
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

    def single_partition(self, k=None, sort_by_length=True, sort_group=False):
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
        if self.alleletable and not self.prunetable:
            self.kprune(self.alleletable, contacts=self.contacts)


        A = self.HG.clique_expansion_init(self.H, self.P_allelic_idx, self.P_weak_idx, 
                                          NW=self.NW,
                                          allelic_factor=self.allelic_factor, 
                                          min_weight=self.min_weight)
        
        dia = A.diagonal()
        raw_contig_counts = A.shape[0]
        retain_idx = np.where(dia > self.min_cis_weight )[0]
        contig_counts = len(retain_idx)

        if len(retain_idx) < raw_contig_counts:
            A = A[retain_idx, :][:, retain_idx]
            self.H, _, _ = extract_incidence_matrix2(self.H, retain_idx)
            self.vertices = self.vertices[retain_idx]

            logger.info(f"Removed {raw_contig_counts - contig_counts:,} contigs that self edge weight < {self.min_cis_weight} (--min-cis-weight).")

        idx_to_vertices = self.idx_to_vertices

        if self.resolution1 == -1:
            result_K_length = 0
            tmp_resolution = self.init_resolution1
            if not k:
                logger.warning("To automatic search best resolution, the `-n` must be specified.")
            while result_K_length < k:
                # logger.info(f"Automatic search the best resolution ... {tmp_resolution:.1f}")
                raw_A, A, self.cluster_assignments, self.K = IRMM(self.H, 
                                                                  A, self.NW, 
                                                        self.P_allelic_idx,
                                                        self.P_weak_idx,
                                                        self.allelic_factor,
                                                        self.cross_allelic_factor,
                                                        tmp_resolution, 
                                                        self.min_weight,
                                                        self.threshold, 
                                                        self.max_round,
                                                        threads=self.threads)
                self.K = list(map(list, self.K))
                self.K = self.filter_cluster(verbose=0)
                result_K_length = len(self.K)
                logger.info(f"Generated `{result_K_length}` groups at resolution `{tmp_resolution}`.")
                tmp_resolution += 0.2

                
        else:
            raw_A, A, self.cluster_assignments, self.K = IRMM(self.H, A, self.NW, 
                                                        self.P_allelic_idx,
                                                        self.P_weak_idx,
                                                        self.allelic_factor,
                                                        self.cross_allelic_factor,
                                                        self.resolution1, 
                                                        self.min_weight,
                                                        self.threshold, 
                                                        self.max_round,
                                                        threads=self.threads)


        # self.K = list(filter(lambda x: len(x) > 1, self.K))
        self.K = list(map(list, self.K))
        self.K = self.filter_cluster()
        # self.K = sorted(self.K, key=lambda x: len(x), reverse=True)
        

        vertices_idx_sizes = self.vertices_idx_sizes
        vertices_idx_sizes = pd.DataFrame(vertices_idx_sizes, index=['length']).T
        
        args = []
        if self.ultra_complex:
            _results = []
            for num, group in enumerate(self.K):
                if len(group) > 1 and HyperPartition.is_error(group, self.allelic_idx_set, vertices_idx_sizes):
                    args.append((group, self.k[0], self.prune_pair_df, self.H, 
                                 vertices_idx_sizes, self.ultra_complex, self.min_weight, 
                                 self.allelic_similarity, self.allelic_factor, 
                                self.cross_allelic_factor, self.min_scaffold_length,
                                self.threshold, self.max_round, num))
                    HyperPartition._incremental_partition(*args[-1])
                else:
                    _results.append(group)
            # _results2 = Parallel(n_jobs=min(self.threads, len(args) + 1))(
            #                 delayed(HyperPartition._incremental_partition) 
            #                  (i, j, _k, l, m, n, o, p, q, r, s, t, u, v) 
            #                     for i, j, _k, l, m, n, o, p, q, r, s, t, u, v in args) 
            _, _results2 = zip(*_results2)
            _results2 = list_flatten(_results2)
            _results.extend(_results2)

            self.K = _results 
        
        if k and len(self.K) > k:  
            logger.info(f"Merging {len(self.K)} groups into {k} groups ...")
            self.K = HyperPartition._merge(A, self.K, vertices_idx_sizes, k, 
                                            self.prune_pair_df, self.allelic_similarity,
                                            self.min_allelic_overlap, method='sum')

        # self.remove_misassembly()

        if self.is_recluster_contigs:
            self.K = HyperPartition.recluster_contigs(self.K, k, A, raw_A,
                                            vertices_idx_sizes, 
                                            None, 
                                            self.prune_pair_df, self.allelic_similarity)
        if sort_group:
            vertices_length = self.contigsizes.loc[self.vertices]['length'].astype('float32')
            self.K = raw_sort(self.K, A, vertices_length, threads=self.threads)

        if sort_by_length:
            self.K = sorted(self.K, 
                            key=lambda x: vertices_idx_sizes.loc[x]['length'].sum(), 
                            reverse=True)
    
        group_info = [i for i in range(1, len(self.K) + 1)]
        length_contents = list(map(
            lambda  x: "{:,}".format(vertices_idx_sizes.loc[x]['length'].sum()), self.K))
        # length_contents = pformat(list(zip(group_info, length_contents)))
        logger.info(f"Hyperpartition result {len(self.K)} groups:")
                    # f"{length_contents}")

        table = Table(header_style=None)
        table.add_column("GroupIdx", justify="center")
        table.add_column("Length", justify="right")
        for i, j in zip(group_info, length_contents):
            table.add_row(str(i), j)

        with console_html.capture() as capture:
            console_html.print(table)

        logger.info(capture.get().strip())
        
        
        return self.K  

    @staticmethod
    def _incremental_partition(
                                # raw_K, raw_A, raw_idx_to_vertices, 
                                K, 
                                idx_to_vertices, 
                                k, alleletable, prune_pair_df,
                                H, vertices_idx_sizes, NW, resolution, 
                                output_dir="./", init_resolution=0.8, min_weight=1, 
                                min_cis_weight = 1.0, allelic_similarity=0.8, 
                                min_allelic_overlap=0.1, allelic_factor=-1, cross_allelic_factor=0.0,
                                is_remove_misassembly=False, is_recluster_contigs=False, refine=False,
                                min_scaffold_length=10000, threshold=0.01, max_round=1, num=None,
                                kprune_norm_method='cis', threads=1):
        """
        single function for incremental_partition.
        """
        if k == 1:
            return None, None, [K]
        
        if NW is not None:
            merge_method = "sum"
        else:
            merge_method = "mean"

        K = np.array(list(K))
        if len(K) == 1:
            return None, None, [K]

        sub_H, _, sub_edge_idx = extract_incidence_matrix2(H, K)
        
        ## remove low weigth contigs
        sub_A = HyperGraph.clique_expansion_init(sub_H, min_weight=min_weight)
        dia = sub_A.diagonal()
        raw_contig_counts = len(K)
        K = K[np.where(dia > min_cis_weight)[0]]
        contig_counts = len(K)
        if contig_counts == 0:
            return None, None, None
        if (raw_contig_counts - contig_counts) > 0:
            logger.info(f"Removed {raw_contig_counts - contig_counts:,} contigs that self edge weight < {min_cis_weight} (--min-cis-weight).")
            sub_H, _, sub_edge_idx = extract_incidence_matrix2(H, K)


        del H 
        gc.collect() 

        if NW is not None:
            sub_NW = NW[list(K)][:, list(K)]
        else:
            sub_NW = None

        sub_A = HyperGraph.clique_expansion_init(sub_H, NW=sub_NW, min_weight=min_weight)

        sub_old2new_idx = dict(zip(K, range(len(K))))
        
        sub_vertices_idx_sizes = vertices_idx_sizes.reindex(K)
        sub_vertices_new_idx_sizes = sub_vertices_idx_sizes
        sub_vertices_new_idx_sizes.index = sub_vertices_idx_sizes.index.map(sub_old2new_idx.get)

        sub_vertices = np.array(list(idx_to_vertices[i] for i in K))
        sub_vertives_idx = dict(zip(sub_vertices, range(len(sub_vertices))))
        
        if alleletable:
            sub_output_contacts = f"{output_dir}/kprune_workdir/{num}.hypergraph.expansion.contacts"
            HyperGraph.to_contacts(sub_A, sub_vertices, NW=sub_NW, min_weight=min_weight,
                                output=sub_output_contacts)
            
            sub_P_allelic_idx, sub_P_weak_idx, sub_prune_pair_df = HyperPartition.kprune_func(
                                    alleletable, sub_output_contacts,
                                    sub_vertives_idx, 
                                    output_dir=f"{output_dir}/kprune_workdir",
                                    kprune_norm_method=kprune_norm_method,
                                    threads=threads*2)


        elif  prune_pair_df is not None:      
            sub_prune_pair_df = prune_pair_df.reindex(list(permutations(K, 2))).dropna().reset_index()

            sub_prune_pair_df['contig1'] = sub_prune_pair_df['contig1'].map(sub_old2new_idx.get)
            sub_prune_pair_df['contig2'] = sub_prune_pair_df['contig2'].map(sub_old2new_idx.get)
            
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
        # _A = HyperGraph.clique_expansion_init(sub_H, None, None, 
        #                                   sub_NW,
        #                                   allelic_factor, min_weight)
        if resolution < 0.0 and k != 0:
            tmp_resolution = init_resolution
            result_K_length = 0
            auto_round = 1

            while result_K_length < k:
                if tmp_resolution > 10.0 or auto_round > 50:
                    break
                # logger.info(f"Automaticly to search best resolution ... {tmp_resolution}")
                raw_A, A, cluster_assignments, new_K = IRMM(sub_H, sub_A, sub_NW, 
                                                    sub_P_allelic_idx, 
                                                    sub_P_weak_idx,
                                                    allelic_factor,
                                                    cross_allelic_factor,
                                                    tmp_resolution,
                                                    min_weight,
                                                    threshold, 
                                                    max_round, 
                                                    threads=1, 
                                                    outprefix=num)
                new_K = list(filter(
                                lambda x: sub_vertices_new_idx_sizes.reindex(x).sum().values[0] \
                                    >= min_scaffold_length, new_K)
                )
                # filtered_K = list(filter(
                #                 lambda x: sub_vertices_new_idx_sizes.reindex(x).sum().values[0] \
                #                     < min_scaffold_length, new_K)
                # )
              
                new_K = list(map(list, new_K))
                tmp_resolution += 0.2  
                result_K_length = len(new_K)
                auto_round += 1
             
        else:
            raw_A, A, cluster_assignments, new_K = IRMM(sub_H, sub_A, sub_NW, 
                                                    sub_P_allelic_idx, 
                                                    sub_P_weak_idx,
                                                    allelic_factor,
                                                    cross_allelic_factor,
                                                    resolution,
                                                    min_weight,
                                                    threshold, 
                                                    max_round, 
                                                    threads=1, 
                                                    outprefix=num)

        
    
        
        ## remove the scaffold that is too short
        # before_filter_K = new_K.copy()
        # new_K = list(map(list, new_K))
        new_K = list(map(list, filter(
                        lambda x: sub_vertices_new_idx_sizes.loc[list(x)].sum().values[0] \
                            >= min_scaffold_length, new_K)))

        filtered_K = list(filter(
                                lambda x: sub_vertices_new_idx_sizes.reindex(x).sum().values[0] \
                                    < min_scaffold_length, new_K)
                )
        
        # removed_idxes = set(list_flatten(before_filter_K)) - set(list_flatten(new_K))
        # if removed_idxes:
        #     raw_A = np.delete(raw_A, removed_idxes, axis=0)
        #     raw_A = np.delete(raw_A, removed_idxes, axis=1)

        #     A = np.delete(A, removed_idxes, axis=0)
        #     A = np.delete(A, removed_idxes, axis=1)
        #     new_K = [k for i, k in enumerate(new_K) if i not in removed_idxes]
            # for idx in removed_idxes:
            #     if idx in sub_vertices_new_idx_sizes:
            #         del sub_vertices_new_idx_sizes[idx]

            ## update A and raw_A and new_K and sub_vertices_new_idx_sizes and 
        
        
        if k and len(new_K) > k:
            new_K = sorted(new_K, key=lambda x: sub_vertices_new_idx_sizes.loc[x]['length'].sum())
            logger.info(f"FirstGroup{num}: Merging {len(new_K)} groups into {k} groups ...")
            new_K = HyperPartition._merge(A, new_K, sub_vertices_new_idx_sizes, k, 
                                            sub_prune_pair_df, allelic_similarity, 
                                            min_allelic_overlap, method=merge_method)
        

        if sub_prune_pair_df is not None and is_remove_misassembly:
            new_K = HyperPartition._remove_misassembly(new_K, sub_H, sub_prune_pair_df, 
                                                        allelic_similarity=allelic_similarity)
        
        if is_recluster_contigs:
            new_K = HyperPartition.recluster_contigs(new_K, k, A, raw_A,
                                         sub_vertices_new_idx_sizes, 
                                          None, 
                                         sub_prune_pair_df, allelic_similarity)
        
        if refine:
            logger.info("Refining misassembies of allelic contig pairs, which partition into the same group")
            new_K = HyperPartition.refine_allelic_errors(new_K, k, A, 
                                                        sub_vertices_new_idx_sizes, 
                                                        sub_prune_pair_df,
                                                        min_weight=min_weight)
        # vertices_length = sub_vertices_new_idx_sizes['length']
        # new_K = raw_sort(new_K, A, vertices_length, threads=1)
        sub_new2old_idx = dict(zip(range(len(K)), K))
        new_K = list(map(lambda x: list(map(lambda y: sub_new2old_idx[y], x)), new_K))

        # new_flatten_K = list_flatten(new_K)
        # if len(new_flatten_K) < len(raw_K):

        #     new_vertices = list(map(lambda x: list(map(idx_to_vertices.get, x)), new_K))
        #     new_flatten_vertices = list_flatten(new_vertices)
            
        #     raw_vertices = list(map(raw_idx_to_vertices.get, raw_K))
        #     raw_vertices_to_idx = dict(items[::-1] for items in raw_idx_to_vertices.items())
            
        #     new_with_raw_idx = list(map(lambda x: list(map(raw_vertices_to_idx.get, x)), new_vertices))

        #     unanchor_vertices = list(filter(lambda x: x not in set(new_flatten_vertices), set(raw_vertices)))
            
        #     sub_raw_new2old_idx = dict(zip(raw_K, range(len(raw_K))))
        #     # sub_raw_A = raw_A[raw_K, :][:, raw_K]
        #     # unanchor_list_new_idx = list(map(sub_raw_new2old_idx.get, unanchor_list))
        #     unanchor_idx = list(map(raw_vertices_to_idx.get, unanchor_vertices))
        #     # unanchor_idx = list(map(sub_raw_new2old_idx.get, unanchor_idx))

        #     # new_with_raw_idx = list(map(lambda x: list(map(sub_raw_new2old_idx.get, x)), new_with_raw_idx))
        #     anchor_res_db = {}
        #     for idx in unanchor_idx:
        #         res_matrix = np.zeros(len(new_with_raw_idx))
        #         for j, g in enumerate(new_with_raw_idx):
        #             # print(raw_idx_to_vertices[idx])
        #             # print(new_vertices[j])
        #             res_matrix[j] = raw_A[idx, :][:, g].sum()
                
        #         anchor_g = np.argmax(res_matrix)
        #         anchor_res_db[idx] = anchor_g
            
        #     for idx in anchor_res_db:
        #         new_with_raw_idx[anchor_res_db[idx]].append(idx)
            
        #     new_K = new_with_raw_idx
            

        new_K = sorted(new_K, key=lambda x: vertices_idx_sizes.reindex(x).dropna()['length'].sum(), reverse=True)

        return A, cluster_assignments, new_K
    

    def kprune(self, alleletable, first_cluster_file=None, contacts=None, is_run=True):
        if is_run:
            if not contacts or not Path(contacts).exists():
                A = HyperGraph.clique_expansion_init(self.H)
                contacts = "hypergraph.expansion.contacts"
                HyperGraph.to_contacts(A, self.vertices , NW=self.NW, min_weight=self.min_weight, output=contacts)
       
        kprune_output_file = "hypergraph.prune.table"
        if is_run:
            if not first_cluster_file:
                cmd = ["cphasing-rs", "kprune", alleletable, contacts, 
                    kprune_output_file, "-n", self.kprune_norm_method, 
                     "-t", str(self.threads)]
            else:
                cmd = ["cphasing-rs", "kprune", alleletable, contacts, kprune_output_file,
                        "-n", self.kprune_norm_method, "-t", str(self.threads), 
                        "-f", first_cluster_file]
        
            logger.info("Generating the prune table ...")
            flag = run_cmd(cmd, log=f"{self.log_dir}/hyperpartition_kprune.log")
            assert flag == 0, "Failed to execute command, please check log."

            logger.info(f"Prune table stored in `{kprune_output_file}`")
        else:
            logger.warning(f"Load exists prune table `{kprune_output_file}`")
            
        self.prunetable = kprune_output_file
        

        self.P_allelic_idx, self.P_weak_idx, self.prune_pair_df = self.get_prune_pairs()

    @staticmethod
    def alleles(fasta, first_cluster, 
                k=19, w=19, m=0.5, d=0.2, c=60, 
                tl=25000, threads=4):
        from .cli import alleles
        try:
            alleles.main(
                args=[
                    "-f",
                    str(fasta),
                    "-fc",
                    str(first_cluster),
                    "-t",
                    str(threads),
                    "-c", str(c),
                    "-k", str(k),
                    "-w", str(w),
                    "-m", str(m),
                    "-d", str(d),
                    "-tl", str(tl),
                ],
                prog_name='alleles'
            )

        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
        

        fasta_prefix = Path(Path(fasta).name).with_suffix("")
        while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
            fasta_prefix = fasta_prefix.with_suffix("")

        return f"{fasta_prefix}.allele.table"

    @staticmethod
    def kprune_func(alleletable, 
                    contacts,
                    vertices_idx,
                    output_dir="./",
                    first_cluster_file=None, 
                    kprune_norm_method='cis',
                    threads=2,
                    ):

        prefix = Path(Path(contacts).stem).stem 

        prunetable =f"{output_dir}/{prefix}.prune.table"
        
        if not first_cluster_file:
            cmd = ["cphasing-rs", "kprune", alleletable, contacts, 
                prunetable, "-n", kprune_norm_method, 
                    "-t", str(threads)]
        else:
            cmd = ["cphasing-rs", "kprune", alleletable, contacts, prunetable,
                    "-n", kprune_norm_method, "-t", str(threads), 
                    "-f", first_cluster_file]
        logger.info("Generating the prune table ...")
        flag = run_cmd(cmd, log=f"{output_dir}/{prefix}.hyperpartition_kprune.log")
        assert flag == 0, "Failed to execute command, please check log."

        logger.info(f"Prune table stored in `{prunetable}`")   

        P_allelic_idx, P_weak_idx, prune_pair_df = HyperPartition.get_prune_pairs_func(prunetable, vertices_idx)

        return P_allelic_idx, P_weak_idx, prune_pair_df

    def incremental_partition(self, k, first_cluster=None):
        """
        incremental partition for autopolyploid.
        """
        if not first_cluster:
            logger.info("Starting first round partition ...")
        
        self.P_allelic_idx, self.P_weak_idx = None, None
        A = HyperGraph.clique_expansion_init(self.H, self.P_allelic_idx, self.P_weak_idx, 
                                          NW=self.NW,
                                          allelic_factor=self.allelic_factor, 
                                          min_weight=self.min_weight)
        
     
        dia = A.diagonal()
        
        raw_contig_counts = A.shape[0]
        retain_idx = np.where(dia > self.min_cis_weight )[0]
        contig_counts = len(retain_idx)
        
        if len(retain_idx) < raw_contig_counts:
            A = A[retain_idx, :][:, retain_idx]
            self.H, _, _ = extract_incidence_matrix2(self.H, retain_idx)
            self.vertices = self.vertices[retain_idx]
            logger.info(f"Removed {raw_contig_counts - contig_counts:,} contigs that self edge weight < {self.min_cis_weight} (--min-cis-weight).")

        vertices_idx_sizes = self.vertices_idx_sizes
        vertices_idx_sizes = pd.DataFrame(vertices_idx_sizes, index=['length']).T
        
        if not first_cluster:
            
            if self.resolution1 < 0 :
                result_K_length = 0
                tmp_resolution = self.init_resolution1
                auto_round = 1
                while result_K_length < k[0] and k[0] != 0:
                    if tmp_resolution > 10.0 or auto_round > 50:
                        break
                    logger.info(f"Automatic search for best resolution ... {tmp_resolution:.1f}")
                    _, A, _, self.K = IRMM(self.H, A, self.NW, 
                            None, None, self.allelic_factor, 
                                self.cross_allelic_factor, tmp_resolution, 
                                self.min_weight, self.threshold, 
                                self.max_round, threads=self.threads)

                    self.K = list(map(list, self.K))
                    self.K = self.filter_cluster(verbose=0)
                    result_K_length = len(self.K)
                    logger.info(f"Generated `{result_K_length}` groups at resolution `{tmp_resolution:.1f}`.")
                    tmp_resolution += 0.2
                    auto_round += 1
            else:
                _, A, _, self.K = IRMM(self.H, A, self.NW, 
                            None, None, self.allelic_factor, 
                                self.cross_allelic_factor, self.resolution1, 
                                self.min_weight, self.threshold, 
                                self.max_round, threads=self.threads)

            self.K = list(map(list, self.K))
            self.K = self.filter_cluster()
            
            if k[0] and len(self.K) > k[0]:
                if self.disable_merge_in_first:
                    self.K = sorted(self.K, key=lambda x: vertices_idx_sizes.loc[x]['length'].sum(), reverse=True)
                    logger.info(f"Discarded `{len(self.K) - k[0]}` smallest groups (--disable-merge-in-first).")
                    self.K = self.K[:k[0]]
                else:
                    self.K = HyperPartition._merge(A, self.K, vertices_idx_sizes, k[0], method='sum')

            self.K = sorted(self.K, key=lambda x: vertices_idx_sizes.loc[x]['length'].sum(), reverse=True)

            if self.exclude_group_from_first:
                for i in self.exclude_group_from_first:
                    try:
                        self.K[i - 1] = []  
                    except IndexError:
                        logger.warning(f"Group `{i}` exceed the numbers of first cluster. skipped")
                    else:
                        logger.info(f"Exclude group `{i}` from first cluster, and this group will not retain in any group.")

                self.K = list(filter(lambda x: len(x) > 0, self.K))

            assert len(self.K) > 0, "Couldn't run first cluster."
            self.to_cluster(f'first.clusters.txt')
            first_cluster_file = f'first.clusters.txt'

            length_contents = list(map(
                lambda  x: "{:,}".format(vertices_idx_sizes.loc[x]['length'].sum()), self.K))
            
            first_group_info = [i for i in range(1, len(self.K) + 1)]
            table = Table(header_style=None)
            table.add_column("GroupIdx", justify="left")
            table.add_column("Length", justify="right")
            total_length = 0
            for i, j in zip(first_group_info, length_contents):
                table.add_row(*map(str, [i, j]))
                total_length += int(j.replace(",", ""))
            
            # table.add_row("Total", f"{total_length:,}")
            
            # first_length_contents = pformat(list(zip(first_group_info, length_contents)))

            logger.info(f"First hyperpartition resulted in {len(self.K)} groups:")
                        # f"{first_length_contents}")
            with console_html.capture() as capture:
                console_html.print(table)
            logger.info(capture.get())   
            mapq_filter = True
        else:
            vertices_idx = self.vertices_idx
            first_cluster_list = list(first_cluster.data.values())
            first_cluster_list = list(map(lambda y: list(
                                    filter(lambda x: x in vertices_idx, y)), first_cluster_list))
            for sub_group in first_cluster_list:
                self.K.append(list(map(lambda x: vertices_idx[x], sub_group)))
            first_cluster_file = first_cluster.filename
            mapq_filter = False
            logger.debug(f"Load the first cluster results from exists file `{first_cluster_file}`.")

        if self.exclude_group_to_second:
            for i in self.exclude_group_to_second:
                if i > len(self.K):
                    logger.warning(f"Exclude group `{i}` exceed the numbers of first cluster. skipped")

        logger.info("Starting second round hyperpartition ...")

        # raw_A = HyperGraph.clique_expansion_init(self.H, P_allelic_idx=self.P_allelic_idx, allelic_factor=0)
        # raw_K = self.K.copy()
        # raw_vertices = self.vertices.copy()
        # raw_idx_to_vertices = self.idx_to_vertices
        if self.HG.edges.mapq.size and (self.min_quality2 > self.min_quality1) and mapq_filter :
            idx_to_vertices = self.idx_to_vertices

            tmp_K = list(map(lambda y: list(
                            map(lambda x: idx_to_vertices[x], y)), 
                            self.K))
            self.HG = HyperGraph(self.edges, min_quality=self.min_quality2)

            del self.edges
            gc.collect()
            self.filter_hypergraph()
            self.H, self.vertices = self.get_hypergraph()
            
            if self.normalize:
                self.NW = self.get_normalize_weight()
            else:
                self.NW = None

            vertices_idx_sizes = pd.DataFrame(self.vertices_idx_sizes, index=['length']).T
            vertices_idx = self.vertices_idx
        
            tmp_K = list(map(lambda x: list(
                            filter(lambda y: y not in self.HG.remove_contigs, x)), tmp_K))
           
            tmp_K = list(map(lambda x: list(
                            filter(lambda y: y in vertices_idx, x )), tmp_K))
            self.K = list(map(lambda x: list(map(lambda y: vertices_idx[y], x)), tmp_K))

        # W = ((self.HG.mapq + 1) / 60).astype(np.float32)
        if self.prunetable:
            self.P_allelic_idx, self.P_weak_idx, self.prune_pair_df = self.get_prune_pairs()
            prune_pair_df = self.prune_pair_df.reset_index().set_index(['contig1', 'contig2'])

        # elif self.alleletable:
        #     # is_run = False if isinstance(first_cluster, ClusterTable) else True
        #     is_run = True
        #     self.kprune(self.alleletable, first_cluster_file, contacts=self.contacts, is_run=is_run)
        #     prune_pair_df = self.prune_pair_df.reset_index().set_index(['contig1', 'contig2'])

        else:
            prune_pair_df = None

        ## Second round cluster
        sub_threads = 1 
        if self.threads // (len(self.K) + 1) < 2:
            sub_threads = 1
        else:
            sub_threads = self.threads // (len(self.K) + 1)

        args = []
        results = []
        self.exclude_groups = []
        idx_to_vertices = self.idx_to_vertices
  
        if self.fasta is not None and self.alleletable is None and self.prunetable is None:
            # if k[1] is not None or k[1] != 0:
            #     c = 5 * k[1]
            # else:
            #     c = 60
            self.alleletable = str(Path(self.alleles(self.fasta, first_cluster_file, 
                                                    k=self.alleles_kmer_size,
                                                    w=self.alleles_window_size,
                                                    d=self.alleles_diff_thres,
                                                    m=self.alleles_minimum_similarity,
                                                    tl=self.alleles_trim_length,
                                                    # c=c,
                                                    threads=self.threads)).absolute())


            
        Path("kprune_workdir").mkdir(exist_ok=True)
        current_dir = Path.cwd()

        for num, sub_k in enumerate(self.K, 1):
            if isinstance(k[1], dict):
                sub_group_number = int(k[1][num - 1])
            else:
                sub_group_number = int(k[1])

            if self.exclude_group_to_second:
                if num in self.exclude_group_to_second:
                    self.exclude_groups.append(sub_k)
                    continue
    
            # sub_raw_k = raw_K[num - 1]
            args.append((
                        # sub_raw_k, raw_A, raw_idx_to_vertices, 
                        sub_k, 
                        idx_to_vertices, 
                        sub_group_number, self.alleletable, prune_pair_df,
                        self.H, vertices_idx_sizes, self.NW,
                        self.resolution2, current_dir, self.init_resolution2, 
                        self.min_weight, self.min_cis_weight,
                        self.allelic_similarity,  self.min_allelic_overlap, 
                        self.allelic_factor, self.cross_allelic_factor, self.is_remove_misassembly,
                        self.is_recluster_contigs, self.refine,
                        self.min_scaffold_length, self.threshold, self.max_round, num,
                        self.kprune_norm_method, sub_threads))
            
            # results.append(HyperPartition._incremental_partition(args[-1])
 
        with parallel_backend('loky', n_jobs=min(self.threads, len(args))):
            try:
                results = Parallel(n_jobs=min(self.threads, len(args)), return_as="generator")(
                            delayed(HyperPartition._incremental_partition)(*a) for a in args)
            except TypeError:
                results = Parallel(n_jobs=min(self.threads, len(args)))(
                            delayed(HyperPartition._incremental_partition)(*a) for a in args)
            results = list(filter(lambda x: x[2] is not None, results))
 

        self.sub_A_list, self.cluster_assignments, results = zip(*results)
        self.inc_chr_idx = []
        for i, j in enumerate(results):
            for _ in j:
                self.inc_chr_idx.append(i)
        
 
        second_group_info = [f"{i}g{j}" for i, res in enumerate(results, 1) 
                                            for j, _ in enumerate(res, 1) ]

        if self.exclude_group_to_second:
            exclude_group_init_idx = len(results) + 1

        self.K = results 

        # self.vertices = raw_vertices
        self.K = list_flatten(results)

        self.K = self.filter_cluster()

        if self.exclude_group_to_second:
            self.K.extend(self.exclude_groups)

        length_contents = list(map(
            lambda  x: "{:,}".format(int(vertices_idx_sizes.reindex(x).dropna()['length'].sum())), self.K))
        second_length_contents = list(zip(second_group_info, length_contents))
        
        if self.exclude_group_to_second:
            for length_content in length_contents[len(second_length_contents):]:
                second_length_contents.append((str(exclude_group_init_idx), length_content))
                exclude_group_init_idx += 1

        table = Table(header_style=None)
        table.add_column("Group", justify="left")
        table.add_column("Length", justify="right")
        total_length = 0
        for i, j in second_length_contents:
            table.add_row(*map(str, [i, j]))
            total_length += int(j.replace(",", ""))
        table.add_row("Total", f"{total_length:,}")
        second_length_contents = pformat(second_length_contents)
        logger.info(f"Second hyperpartition resulted in {len(self.K)} groups:")
                    # f"{second_length_contents}")
        with console_html.capture() as capture:
            console_html.print(table)

        logger.info(capture.get())
    

    @staticmethod
    def _merge(A, K, vertices_idx_sizes, k=None, 
                prune_pair_df=None, allelic_similarity=0.85,
                min_allelic_overlap=0.1, method='mean'):
        if not k:
            return K 
        
        # if prune_pair_df is None:
        #     return K
        # tmp_df = prune_pair_df[(prune_pair_df['type'] == 0) & 
        #                         (prune_pair_df['similarity'] >= allelic_similarity)][['contig1', 'contig2']]
        # P_allelic_idx = [tmp_df['contig1'], tmp_df['contig2']]
        
        # if P_allelic_idx:
        #     A = A.tolil()
        #     A[P_allelic_idx[0], P_allelic_idx[1]] = -1
        #     A = A.tocsr()
        
        if prune_pair_df is not None:
            allelic_pair_df = prune_pair_df[prune_pair_df['type'] == 0]
            allelic_idx_set = set(map(tuple, 
                            prune_pair_df[(prune_pair_df['type'] == 0) & 
                                                (prune_pair_df['similarity'] >= allelic_similarity)]
                            [['contig1', 'contig2']].values)
                            )
            mz_df = allelic_pair_df[['contig1', 'mz1']]
            mz_df = mz_df.drop_duplicates(subset=['contig1']).set_index('contig1')
            allelic_pair_df = allelic_pair_df.set_index(['contig1', 'contig2'])
            


        iter_round = 0 
        flag = 1
      
        while k and len(K) > k:
            current_group_number = len(K)
            if iter_round > (current_group_number - k + 50):
                break
            value_matrix = np.zeros(shape=(current_group_number, current_group_number))
            flag_matrix = np.ones(shape=(current_group_number, current_group_number))
            res = {}
            for i, group1 in enumerate(K):
                group1 = list(group1)
                for j in range(i + 1, current_group_number):
                    group2 = list(K[j])
                    
                    group1_length = vertices_idx_sizes.reindex(group1).sum().values[0]
                    group2_length = vertices_idx_sizes.reindex(group2).sum().values[0]
                    if prune_pair_df is not None:
                        product_contig_pair = set(product(group1, group2))
                        allelic = product_contig_pair & allelic_idx_set
                        if allelic:
                            tmp1, tmp2 = list(map(list, map(set, zip(*allelic))))
                            tmp1_length = vertices_idx_sizes.loc[tmp1].sum().values[0]
                            tmp2_length = vertices_idx_sizes.loc[tmp2].sum().values[0]
                            overlap1 = tmp1_length / group1_length
                            overlap2 = tmp2_length / group2_length
                            
                            # tmp_mzShared = allelic_pair_df.reindex(product_contig_pair)['mzShared'].sum()
                            # tmp1_mz = mz_df.reindex(group1)['mz1'].sum()
                            # tmp2_mz = mz_df.reindex(group2)['mz1'].sum()
                           
                            # overlap1 = tmp_mzShared  / tmp1_mz 
                            # overlap2 = tmp_mzShared  / tmp2_mz
                            # print(i, j, group1_length, group2_length, overlap1, overlap2)

                            if overlap1 > min_allelic_overlap or overlap2 > min_allelic_overlap:
                                flag = 0
                    
                    if method == "mean":
                        value = A[group1, ][:,group2 ].mean() 
                    else:
                        value = A[group1, ][:,group2 ].sum() 

                    value_matrix[i, j] = value
                    flag_matrix[i, j] = flag
                    flag = 1
            
            value_matrix = value_matrix + value_matrix.T - np.diag(value_matrix.diagonal())
            total_value = np.triu(value_matrix).sum()


            for i in range(current_group_number):
                i_value = value_matrix[i].sum()
                for j in range(i+1, current_group_number):
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

    def merge_cluster(self, groups):
        """
        merge several clusters into k groups
        """
        vertices_idx = self.vertices_idx
        vertices_idx_sizes = self.vertices_idx_sizes
        vertices_idx_sizes = pd.DataFrame(vertices_idx_sizes, index=['length']).T

        A = self.HG.clique_expansion_init(self.H, self.P_allelic_idx, self.P_weak_idx, 
                                          NW=self.NW,
                                          allelic_factor=self.allelic_factor, 
                                          min_weight=self.min_weight)
        
        _K = groups
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
    
    def phased_merge_cluster(self, merge_cluster):
        """
        merge the phased cluster
        """

        ## find hap 
        hap_group = OrderedDict()
        for group in merge_cluster.groups:
            hap = group.rsplit("g")[0]
            if hap not in hap_group:
                hap_group[hap] = []
            hap_group[hap].append(group)
    
    def remove_misassembly(self):
        idx_to_vertices = self.idx_to_vertices
        tmp_df = self.prune_pair_df[(self.prune_pair_df['type'] == 0) & 
                                    (self.prune_pair_df['similarity'] >= self.allelic_similarity)][['contig1', 'contig2']]
        P_allelic_idx = [tmp_df['contig1'], tmp_df['contig2']]
        A = HyperGraph.clique_expansion_init(self.H, P_allelic_idx, allelic_factor=0)

        allelic_idx_set = set(map(tuple, tmp_df.values))
        cadinate_missassembly_groups = []
        removed_missassembly_groups = []
        raw_group_idx_db = {}
        for i, group in enumerate(self.K):
            group.sort()
            
            raw_group_idx_db.update(dict(zip(group, cycle([i]))))
            res = set(combinations(group, 2)).intersection(allelic_idx_set)

            cadinate_missassembly_groups.append(list(map(np.array, res)))
            removed_missassembly_groups.append(list(set(group) - set(list_flatten(list(res)))))
       

        logger.debug(f"Find {sum(map(len, cadinate_missassembly_groups))} pairs cadinate homolog misassemblies..")
        logger.debug(len(removed_missassembly_groups))
        corrected_idx = {}
        for j, mis_groups in enumerate(cadinate_missassembly_groups):
            for mis in mis_groups:
                if len(mis) == 0:
                    continue
                values = np.zeros(shape=(2, len(self.K)))
                
                for k in range(values.shape[1]):
                    values[0, k] = A[mis[0], :][:, removed_missassembly_groups[k]].sum(axis=1)
                    values[1, k] = A[mis[1], :][:, removed_missassembly_groups[k]].sum(axis=1)
    
                
                values_argmax = values.argmax(axis=1)
                res = values_argmax == j
                logger.debug(np.array(list(map(idx_to_vertices.get, mis))))
                logger.debug(values)

                if res.all():
                    continue
                
                if not res.any():
                    corrected_idx[mis[0]] = values_argmax[~res][0]
                    corrected_idx[mis[1]] = values_argmax[~res][0]
                    # removed_missassembly_groups[corrected_idx[mis[0]]].add(mis[0])
                    # removed_missassembly_groups[corrected_idx[mis[1]]].add(mis[1])
                    continue
               
                corrected_idx[mis[~res][0]] = values_argmax[~res][0]
                # removed_missassembly_groups[corrected_idx[mis[~res][0]]].add(mis[~res][0])
                
                
        logger.debug(f"Corrected {len(corrected_idx)} misassemblies.")
        for idx in corrected_idx:
            # logger.debug(" ".join(map(str, [idx, corrected_idx[idx], raw_group_idx_db[idx]])))
            self.K[raw_group_idx_db[idx]].remove(idx)
            self.K[corrected_idx[idx]].append(idx)

    @staticmethod
    def _remove_misassembly(K, sub_H, prune_pair_df=None, allelic_similarity=0.85):
        
        if prune_pair_df is None:
            return K
        tmp_df = prune_pair_df[(prune_pair_df['type'] == 0) & 
                                (prune_pair_df['similarity'] >= allelic_similarity)][['contig1', 'contig2']]
        P_allelic_idx = [tmp_df['contig1'], tmp_df['contig2']]
        allelic_idx_set = set(map(tuple, tmp_df.values))
        # A = sub_raw_A.tolil()
        # A[P_allelic_idx[0], P_allelic_idx[1]] = 0
        # A = A.tocsr()
        A = HyperGraph.clique_expansion_init(sub_H, P_allelic_idx, allelic_factor=0)
        cadinate_missassembly_groups = []
        removed_missassembly_groups = []
        raw_group_idx_db = {}
        for i, group in enumerate(K):
            group.sort()
            
            raw_group_idx_db.update(dict(zip(group, cycle([i]))))
            res = set(combinations(group, 2)).intersection(allelic_idx_set)

            cadinate_missassembly_groups.append(list(map(np.array, res)))
            removed_missassembly_groups.append(set(group) - set(list_flatten(list(res))))
    
        logger.debug(f"Find {sum(map(len, cadinate_missassembly_groups))} pairs cadinate homolog misassemblies..")
        corrected_idx = {}
        for j, mis_groups in enumerate(cadinate_missassembly_groups):
            for mis in mis_groups:
                if len(mis) == 0:
                    continue
                values = np.zeros(shape=(2, len(K)))
                A_mis = A[mis, :]
                for k in range(values.shape[1]):
                    values[:, k] = A_mis[:, list(removed_missassembly_groups[k])].sum()

                values_argmax = values.argmax(axis=1)
                res = values_argmax == j
                # logger.debug(np.array(list(map(idx_to_vertices.get, mis)))[~res])
                if res.all():
                    continue
                
                if not res.any():
                    corrected_idx[mis[0]] = values_argmax[~res][0]
                    corrected_idx[mis[1]] = values_argmax[~res][0]
                    continue
               
                corrected_idx[mis[~res][0]] = values_argmax[~res][0]
                
                
        logger.debug(f"Corrected {len(corrected_idx)} misassemblies.")
        for idx in corrected_idx:
            logger.debug(" ".join(map(str, [idx, corrected_idx[idx], raw_group_idx_db[idx]])))
            K[raw_group_idx_db[idx]].remove(idx)
            K[corrected_idx[idx]].append(idx)
        
        return K
    
    @staticmethod
    def refine_allelic_errors(K, k, A, idx_size, prune_pair_df, 
                              min_weight, allelic_similarity=0.99, 
                              rel_tol=0.2):
        if prune_pair_df is None:
            return K 
        if len(K) == 1:
            return K
    
        
        tmp_df = prune_pair_df[(prune_pair_df['type'] == 0) &
                                (prune_pair_df['similarity'] >= allelic_similarity)]
        
        idx_size_db = idx_size.to_dict()['length']

        tmp_df['length1'] = tmp_df['contig1'].map(idx_size_db)
        tmp_df['length2'] = tmp_df['contig2'].map(idx_size_db)
    
        high_similar_pairs = set()
        for idx, row in tmp_df.iterrows():
            contig1, contig2 = row['contig1'], row['contig2']
            mz1, mz2 = row['mz1'], row['mz2']
            length1, length2 = row['length1'], row['length2']
            if math.isclose(mz1, mz2, rel_tol=rel_tol) or \
                math.isclose(length1, length2, rel_tol=rel_tol):
                if contig1 > contig2:
                    contig1, contig2 = contig2, contig1
                high_similar_pairs.add((contig1, contig2))


        if not high_similar_pairs:
            return K
        
        retain_res = OrderedDict()
        candicate_error_contig_pairs = OrderedDict()
        all_candicate_error_contigs = set()
        all_candicate_error_contig_db = defaultdict(set)
        for idx, group in enumerate(K):
            group_set = set(combinations(sorted(group), 2))
            candicate_allelic_error_pairs = group_set.intersection(high_similar_pairs)
            candicate_error_contigs = set(list_flatten(candicate_allelic_error_pairs))
            all_candicate_error_contigs.update(candicate_error_contigs)
            for pair in candicate_allelic_error_pairs:
                all_candicate_error_contig_db[pair[0]].add(pair[0])
                all_candicate_error_contig_db[pair[0]].add(pair[1]
                                                           )

            retain_contigs = set(group).difference(candicate_error_contigs)
            retain_res[idx] = retain_contigs
    
            for pair in candicate_allelic_error_pairs:
                if pair not in candicate_error_contig_pairs:
                    candicate_error_contig_pairs[pair] = idx
        
        for pair in candicate_error_contig_pairs:
            idx = candicate_error_contig_pairs[pair]
            retain_contigs = retain_res[idx]
            other_idxes = sorted(set(list(retain_res)).difference([idx]))

            counts_list = []
            for contig in pair:
                
                interaction_pairs = list(product([int(contig)], retain_contigs))
                
                interaction_pairs = [(c1, c2) if c1 < c2 else (c2, c1) for c1, c2 in interaction_pairs] 
                contig1_idxes, contig2_idxes = list(zip(*interaction_pairs))
                
                sub_A = A[contig1_idxes, :][:, contig2_idxes]

                counts_list.append(sub_A.sum())
            
            max_idx = pair[np.argmax(counts_list)]
            min_idx = pair[np.argmin(counts_list)]

            retain_res[idx].add(max_idx)

            other_counts_list = []
            for other_idx in other_idxes:
                other_retain_contigs = retain_res[other_idx]
                if len(other_retain_contigs) == 0:
                    continue
                interaction_pairs = list(product([min_idx], other_retain_contigs))
                interaction_pairs = [(c1, c2) if c1 < c2 else (c2, c1) for c1, c2 in interaction_pairs] 
                contig1_idxes, contig2_idxes = list(zip(*interaction_pairs))
                if len(contig1_idxes) == 0 or len(contig2_idxes) == 0:
                    other_counts_list.append(0)
                else:
                    sub_A = A[contig1_idxes,  contig2_idxes]
                    other_counts_list.append(sub_A.sum())
            
            other_counts_db = dict(zip(other_idxes, other_counts_list))

            other_counts_list = sorted(other_counts_db, 
                                       key=lambda x: other_counts_db[x],
                                       reverse=True)

            for i in range(len(other_counts_list)):
                max_chrom = other_counts_list[i]
                try:
                    second_chrom = other_counts_list[i+1]
                except:
                    second_chrom = None
                
                tmp_error_contigs = all_candicate_error_contig_db[min_idx]
                tmp_counts = len(retain_res[max_chrom].intersection(tmp_error_contigs))
                if tmp_counts > 0:
                    continue 
                if second_chrom:
                    if (other_counts_db[max_chrom] > min_weight) and (other_counts_db[max_chrom] != other_counts_db[second_chrom]):
                        retain_res[max_chrom].add(min_idx)
                        break
                else:
                    if other_counts_db[max_chrom] > min_weight:
                        retain_res[max_chrom].add(min_idx)
                        break

        new_K = []
        for idx in sorted(retain_res):
            new_K.append(sorted(retain_res[idx]))
        

        return new_K 
         

    def rescue(self):
        """
        rescue unanchor contigs 
        """
        pass 

    @staticmethod
    def recluster_contigs(K, k, A, raw_A, idx_size, idx2vertices=None,
                          prune_pair_df=None, allelic_similarity=0.85,
                          delta_modularity=0.1, is_recluster_high_identity=False):
        """
        recluster contigs 
        """
        if prune_pair_df is None:
            return K
        
        if len(K) == 1:
            return K

        if len(K) < k:
            return K
        
        all_idx = list_flatten(K)

        tmp_df = prune_pair_df[(prune_pair_df['type'] == 0) & 
                                (prune_pair_df['similarity'] >= 0.99)] # [['contig1', 'contig2']]
        idx_size_db = idx_size.to_dict()['length']

        tmp_df.set_index(['contig1', 'contig2'], inplace=True)
        candicate_contig_assignment = {}
        all_candicate_contigs = set()
        new_K = K.copy()
        for i, g in enumerate(K):
            g_len = sum(idx_size_db[c] for c in g)
            candicate_allele_df = tmp_df.reindex(list(combinations(g, 2))).dropna()
            
            candicate_contigs = set(list_flatten(candicate_allele_df.index))
            for _c in candicate_allele_df.index.tolist():
                candicate_contig_assignment[_c] = i
            
            candicate_len = sum(idx_size_db[c] for c in candicate_contigs)
            candicate_error_rate = candicate_len / g_len

            if candicate_error_rate > 0.3:
                pass 
            
            new_K[i] = list(set(g) - candicate_contigs)
            all_candicate_contigs.update(candicate_contigs)
        
        # for contig in sorted(all_candicate_contigs, key=idx_size_db.get, reverse=True):
        #     m = np.zeros(len(new_K))
        #     for i, g in enumerate(new_K):    
        #         conflict_df = tmp_df.reindex(list(product([contig], g))).dropna()
        #         if len(conflict_df) >= 1:
        #             m[i] = 0
        #         else:
        #             m[i] = A[contig, g].mean()

        #     new_K[m.argmax()].append(contig)
            
        
        if is_recluster_high_identity:

            reclustered_contigs = set()
            for contig_pair in candicate_contig_assignment:
                contig1, contig2 = contig_pair 
                ori_g = candicate_contig_assignment[contig_pair]
                m = np.zeros(shape=(2, len(K)))
                edges = []
                weights = []

                for i, g in enumerate(new_K):
                    

                    conflict_df1 = tmp_df.reindex(list(product([contig1], g))).dropna()

                    if len(conflict_df1) >= 1:
                        m[0, i] = 0
                    else:
                        m[0, i] = A[contig1, g].sum()
                        
                    conflict_df2 = tmp_df.reindex(list(product([contig2], g))).dropna()

                    if len(conflict_df2) >= 1:
                        m[1, i] = 0
                    else:
                        m[1, i] = A[contig2, g].sum()
                    

                    edges.append((i, len(K)))
                    weights.append(m[0, i])
                    edges.append((i, len(K) + 1))
                    weights.append(m[1, i])
                
                types = [0] * len(K) + [1] * 2

                g = ig.Graph.Bipartite(types, edges=edges)
                g.es['weight'] = weights
                matching = g.maximum_bipartite_matching(weights='weight')
                
                matching_results = []
                for j in range(len(K)):
                    matching_results.append((j, matching.match_of(j)))

                matching_results = list(filter(lambda x: x[0] is not None and x[1] is not None, matching_results))

                if len(matching_results) <= 1:
                    # if contig1 not in reclustered_contigs:
                    #     new_K[ori_g].append(contig1)
                    #     reclustered_contigs.add(contig1)
                    # if contig2 not in reclustered_contigs:
                    #     new_K[ori_g].append(contig2)
                    #     reclustered_contigs.add(contig2)
                    continue
                
                if contig1 not in new_K[matching_results[0][1] - len(K)] and \
                    contig1 not in new_K[matching_results[1][1] - len(K)]:
                    if contig1 in reclustered_contigs:
                        continue
                    reclustered_contigs.add(contig1)
                    new_K[matching_results[0][1] - len(K)].append(contig1)
                if contig2 not in new_K[matching_results[0][1] - len(K)] and \
                    contig2 not in new_K[matching_results[1][1] - len(K)]:
                    if contig2 in reclustered_contigs:
                        continue
                    reclustered_contigs.add(contig2)
                    new_K[matching_results[1][1] - len(K)].append(contig2)

        K = new_K.copy()
        
        
        # if A.ndim > 1 and (A.shape[0] != len(all_idx)):
        #     loss_idx = list(set(range(A.shape[0])) - set(all_idx))
        #     if loss_idx:
        #         mask = np.ones(A.shape[0], dtype=bool)
        #         mask[loss_idx] = False
        #         A = A[mask][:, mask]
        
        # if raw_A.ndim > 1 and (raw_A.shape[0] != len(all_idx)):
        #     loss_idx = list(set(range(raw_A.shape[0])) - set(all_idx))
        #     if loss_idx:
        #         mask = np.ones(raw_A.shape[0], dtype=bool)
        #         mask[loss_idx] = False
        #         raw_A = raw_A[mask][:, mask]

        # # raw_A = raw_A.tolil()
        # # if P_allelic_idx:
        # #     raw_A = raw_A.tolil()
        # #     raw_A[P_allelic_idx[0], P_allelic_idx[1]] = 0
        # #     raw_A = raw_A.tocsr()

        # old2new = dict(zip(list_flatten(K), range(len(all_idx))))
        # new2old = dict(zip(range(len(all_idx)), all_idx))
        
      
        # membership = np.zeros(len(list_flatten(K)), dtype=np.int32)
        # graph = ig.Graph.Weighted_Adjacency(A, mode='undirected', loops=False)

        # K = list(map(lambda x: list(map(old2new.get, x)), K))
        # i = 0 
        # for gn, l in enumerate(K):
        #     for j in l:
        #         membership[i] = gn
        #         i += 1
   
        # init_membership = membership.copy()
        # init_modularity = graph.modularity(init_membership)

        # logger.debug(f"Before recluster the modularity is :`{init_modularity:.4f}`")
        # for j, ori_i in enumerate(init_membership):
        #     tmp_memberhip = init_membership.copy()
            
        #     for i in range(len(K)):
        #         if i == ori_i:
        #             continue

        #         tmp_memberhip[j] = i  
        #         modularity = graph.modularity(tmp_memberhip)
        #         is_recluster = True if (modularity - init_modularity) > delta_modularity else False 
    
        #         if is_recluster:
        #             if idx2vertices:
        #                 logger.debug(f"Change `{idx2vertices[j]}` from {ori_i} to {i} group")
        #             else:
        #                 logger.debug(f"Change `{j}` from {ori_i} to {i} group")
                    
        #             membership[j] = i

        # final_modularity = graph.modularity(membership)

        # logger.debug(f"Final modularity is {final_modularity:.4f}")
        
        # new_K = [[] for _ in range(max(membership) + 1)]
        # for i, j in enumerate(membership):
        #     new_K[j].append(new2old[i])
        
        new_K = list(filter(lambda x: len(x) > 0, new_K))

        return new_K

    def rescue_collapsed(self, is_duplicate=True):
        """
        rescue collapsed contigs

        
        """
        pass 
    # def remove_misassembly_hap(self):
    #     if not self.prune_pair_df:
    #         return
        
    #     idx_to_vertices = self.idx_to_vertices
    #     A = HyperGraph.clique_expansion_init(self.H)
    #     logger.debug(A.shape)
    #     allelic_idx_set = set(map(tuple, 
    #                     self.prune_pair_df[(self.prune_pair_df['type'] == 0) & 
    #                                         (self.prune_pair_df['similarity'] >= self.allelic_similarity)]
    #                     [['contig1', 'contig2']].values)
    #                     )
    #     for i, hap_groups in enumerate(self.K):
    #         if len(hap_groups) < 1:
    #             continue 
            
    #         candicate_misassembly_groups = []
    #         removed_misassembly_hap_groups = []
    #         for group in hap_groups:
    #             group.sort()
                
    #             res = set(combinations(group, 2)).intersection(allelic_idx_set)
    #             candicate_misassembly_groups.append(res)
    #             removed_misassembly_hap_groups.append(set(group) - res)
            
    #         for j, candicate_misassemblies in enumerate(candicate_misassembly_groups):
    #             if not candicate_misassemblies:
    #                 continue
                
    #             for mis in candicate_misassemblies:
    #                 logger.debug(list(map(idx_to_vertices.get, mis)))
                

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
    
    def filter_cluster(self, verbose=1):
        if verbose == 1: 
            logger.info(f"Removed groups less than {to_humanized2(self.min_scaffold_length)} in length. (--min-scaffold-length)")

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
                last_chr_idx = 0
                output_group_number = 0 
                for i, chr_idx in enumerate(self.inc_chr_idx):
                    group = clusters[i]
                    hap_idx = db[chr_idx]

                    print(f'Chr{chr_idx + 1:0>2}g{hap_idx}\t{len(group)}\t{" ".join(group)}', 
                        file=out)
                    db[chr_idx] += 1
                    output_group_number += 1
                    last_chr_idx = chr_idx + 1
                else:
                    for i, group in enumerate(clusters[output_group_number:]):
                        print(f'Chr{i + 1 + last_chr_idx:0>2}\t{len(group)}\t{" ".join(group)}', 
                            file=out)

            except AttributeError:
                for i, group in enumerate(clusters, 1):
                    print(f'Chr{i:0>2}\t{len(group)}\t{" ".join(group)}', 
                        file=out)

        logger.info(f"Successful output hyperpartition results in `{output}`.")
