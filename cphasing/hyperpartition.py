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
import polars as pl
import threading

from contextlib import contextmanager
from collections import defaultdict, OrderedDict 
from functools import cached_property
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
from .utilities import (
    list_flatten, 
    run_cmd, 
    to_humanized2, 
    parse_split_contigs,
    is_file_changed,
    get_fasta_from_split_contig,
)
from ._config import *

logger = logging.getLogger(__name__)

# from line_profiler import profile

@contextmanager
def silence_logs(logger_names, level=logging.WARNING):
    loggers = [logging.getLogger(name) for name in logger_names]
    original_levels = [log.level for log in loggers]
    try:
        for log in loggers:
            log.setLevel(level)
        yield
    finally:
        for log, orig_level in zip(loggers, original_levels):
            log.setLevel(orig_level)

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
                    split=False,
                    contacts=None,
                    kprune_norm_method="auto",
                    count_re=None,
                    allelic_factor=-1,
                    cross_allelic_factor=0.3,
                    allelic_similarity=0.8,
                    min_allelic_overlap=0.1,
                    ultra_complex=5.0,
                    disable_merge_in_first=False,
                    merge_use_allele=False,
                    allelic_positive_factor=3.0,
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
                    cluster_method='louvain',
                    min_quality1=1,
                    min_quality2=2,
                    mapq_filter_1_to_2=False,
                    min_scaffold_length=10000,
                    is_remove_misassembly=False,
                    is_recluster_contigs=False,
                    disable_recluster_by_linkage=False,
                    refine=False,
                    collapsed_contigs_df=None,
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
        self.split = split
        self.contacts = contacts
        self.kprune_norm_method = kprune_norm_method
        self.count_re = count_re
            
        self.allelic_factor = allelic_factor
        self.cross_allelic_factor = cross_allelic_factor
    
        self.allelic_similarity = allelic_similarity
        self.min_allelic_overlap = min_allelic_overlap
        self.ultra_complex = ultra_complex
        self.disable_merge_in_first = disable_merge_in_first
        self.merge_use_allele = merge_use_allele
        self.allelic_positive_factor = allelic_positive_factor if min_quality1 != 0 else 1.0
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
        self.cluster_method = cluster_method
        self.min_quality1 = min_quality1
        self.min_quality2 = min_quality2
        self.min_scaffold_length = int(min_scaffold_length)
        self.is_remove_misassembly = is_remove_misassembly
        self.is_recluster_contigs = is_recluster_contigs
        self.is_recluster_by_linkage = False if disable_recluster_by_linkage else True
        
        self.refine = refine
        self.collapsed_contigs_df = collapsed_contigs_df
        logger.debug(f"is remove misassembly: {self.is_remove_misassembly}")
        self.threshold = threshold
        self.max_round = max_round
        self.threads = threads
        self.chunksize = int(chunksize) if chunksize else None
        
        self.log_dir = "logs"
        Path(self.log_dir).mkdir(exist_ok=True)

        self.contig_sizes = self.contigsizes.to_dict()['length'] ## dictionary

        if self.split:
            self.parse_trimmed_contig()
        self.contigs = self.contigsizes.index.values.tolist()
        self.K = []

        self._cached_raw_contigsizes = None
        if self.split:
            self.get_raw_contigsizes()

        if not mapq_filter_1_to_2:
            logger.debug(f"min_quality1: {self.min_quality1}")
            self.HG = HyperGraph(self.edges, min_quality=self.min_quality1, mapq_filter=False)
        else:
            logger.debug(f"min_quality2: {self.min_quality2}")
            if self.min_quality2 > self.min_quality1:
                mapq_filter = True
            else:
                mapq_filter = False
            self.HG = HyperGraph(self.edges, min_quality=self.min_quality2, mapq_filter=mapq_filter)

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
            self.cis_counts = None
        else:
            
            H = self.HG.incidence_matrix(self.min_contacts)
            # cis_df = pl.DataFrame({"col": self.HG.col, "count": self.HG.count})

            # grouped = cis_df.group_by("col").agg(pl.col("count").first())

            # grouped = grouped.sort("col")
            # # filter row not in self.HG.remove_contig_idx
            # grouped = grouped.filter(~pl.col("col").is_in(pl.Series(self.HG.remove_contig_idx)))
            # grouped = grouped.to_pandas()
            # grouped.reset_index(drop=True, inplace=True)
            # grouped.drop('col', axis=1, inplace=True)
            # self.cis_counts = grouped['count'].values
            self.cis_counts = None
            self.HG.clear()

            vertices = self.HG.nodes

        logger.info(f"Generated filtered hypergraph that containing {H.shape[0]:,} vertices"    
                    f" and {H.shape[1]:,} hyperedges.")

        return H, vertices

    def filter_hypergraph(self, remove_contigs=set()):
        """
        remove too short contigs and according the whitelist or blacklist 
        to remove contigs
        """
        logger.debug("Start to filter hypergraph ...")
        short_contigs = self.contigsizes[self.contigsizes < self.min_length]
        short_contigs = short_contigs.dropna().index.tolist()

        remove_contigs.update(set(short_contigs))

        if self.whitelist:
            if not self.split:
                remove_contigs.update(set(self.contigs) - set(self.whitelist))
            else:
                _remov_contigs = set()
                for x in self.contigs:
                    if parse_split_contigs(x)[0] not in self.whitelist:
                        _remov_contigs.add(x)
                remove_contigs.update(_remov_contigs)

        if self.blacklist:
            if not self.split:
                remove_contigs.update(set(self.blacklist))
            else:
                _remov_contigs = set()
                for x in self.contigs:
                    if parse_split_contigs(x)[0] in self.blacklist:
                        _remov_contigs.add(x)
                remove_contigs.update(_remov_contigs)
        
        if len(remove_contigs) == 0:
            return
        
        self.HG.remove_rows(remove_contigs)
        logger.info(f"Total {len(remove_contigs):,} contigs were removed,")
        logger.info(f"\tbecause it's length too short (<{self.min_length}) "
                        "or your specified in blacklist or not in whitelist.")


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
        pair_df = pair_df[pair_df['contig1'].isin(self.vertices) & \
                            pair_df['contig2'].isin(self.vertices)]
        if pair_df.empty:
            logger.warning("No valid prune pairs found in the provided prunetable.")
            return None, None, None
        
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
        """
        Get normalize weight matrix

        Supported methods:
        -----------------
        none: no normalization
        geom/sqrt_len: 1/sqrt(Li*Lj)
        density/len: 1/(Li*Lj)
        log/loglen: log compression
        len_deg/ld: length and degree combined normalization
        len_deg: 1/((Li^alpha)*(Di^beta))

        Params:
        -------
        method: str or dict
            normalization method

        Returns:
        --------
        NW: np.ndarray or None
            normalize weight matrix
        """
        method = "svc"
        if isinstance(method, bool):
            method = "geom" if method else "none"

        alpha = 1.0
        beta = 1.0
        if isinstance(method, dict):
            alpha = float(method.get("alpha", 1.0))
            beta = float(method.get("beta", 1.0))
            method = method.get("method", "len_deg")

        a = self.contigsizes.loc[self.vertices]['length'].astype(np.float64).values
        a[a <= 0] = 1.0 
        deg = np.asarray(self.H.sum(axis=1)).ravel().astype(np.float64)
        deg[deg <= 0] = 1.0

        eps = 1e-9
        if method in ("none", None):
            return None

        if method in ("geom", "sqrt_len"):
            w = 1.0 / np.sqrt(a)
        elif method in ("density", "len"):
            w = 1.0 / a
        elif method in ("log", "loglen"):
            med = float(np.median(a))
            scale = np.log1p(a / max(med, 1.0))
            w = 1.0 / np.sqrt(np.clip(scale, eps, None))
        elif method in ("len_deg", "ld"):
            deg = np.asarray(self.H.sum(axis=1)).ravel().astype(np.float64)
            deg[deg <= 0] = 1.0
            w = 1.0 / (np.power(a, alpha/2.0) * np.power(deg, beta/2.0))
        elif method in ("vc",):
            w = 1.0 / np.sqrt(deg)

        elif method in ("svc", "smoothed_vc"):
            global_density = np.sum(deg) / np.sum(a)
            expected_deg = a * global_density
            smoothed_deg = deg + 0.5 * expected_deg
            w = 1.0 / np.sqrt(smoothed_deg)

        elif method in ("hybrid", ):
            a_n = a / max(float(np.median(a)), 1.0)
            d_n = deg / max(float(np.median(deg)), 1.0)
            w = 1.0 / np.sqrt(np.power(np.clip(a_n, eps, None), alpha) *
                              np.power(np.clip(d_n, eps, None), beta))
        else:
            w = 1.0 / np.sqrt(a)

        NW = np.outer(w, w).astype(np.float32)

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

    def single_partition_star(self, k=None, sort_by_length=True, sort_group=False, is_kprune=False):
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


        if is_kprune:
            if self.fasta is not None and self.alleletable is None and self.prunetable is None:
                c = 100
                fasta_prefix = Path(Path(self.fasta).name)
                while fasta_prefix.suffix in {".fasta", "gz", "fa", '.fa', '.gz', ".split"}:
                    fasta_prefix = fasta_prefix.with_suffix("")
                    
                if self.split:
                    
                    _contigsizes = self.contigsizes[self.contigsizes['length'] > self.min_length]
                    from .utilities import get_fasta_from_split_contig
                    used_fasta = get_fasta_from_split_contig(
                        self.fasta,
                        _contigsizes,
                        f"{fasta_prefix}.split.fasta",
                    )
                    
                else:
                    used_fasta = self.fasta

                self.alleles(used_fasta, 
                                                    k=self.alleles_kmer_size,
                                                    w=self.alleles_window_size,
                                                    d=self.alleles_diff_thres,
                                                    m=self.alleles_minimum_similarity,
                                                    tl=self.alleles_trim_length,
                                                    c=c,
                                                    split=self.split,
                                                    threads=self.threads)
                if self.split:
                    self.alleletable = f"{fasta_prefix}.split.split.allele.table"
                else:
                    self.alleletable = f"{fasta_prefix}.allele.table"
            
            if self.alleletable and not self.prunetable:
                self.kprune(self.alleletable, contacts=self.contacts)
        
    
        A = self.HG.clique_expansion_init(self.H, 
                                        #   NW=self.NW,
                                          allelic_factor=self.allelic_factor, 
                                          min_weight=self.min_weight)

        retain_idx1 = HyperGraph.filter_adjacency_matrix(A, self.vertices, self.vertices_idx_sizes, invert=True)
        if retain_idx1 is None:
            retain_idx1 = np.arange(A.shape[0], dtype=np.int64)
        else:
            retain_idx1 = np.asarray(retain_idx1)
            if retain_idx1.dtype == bool:
                retain_idx1 = np.where(retain_idx1)[0].astype(np.int64)
            else:
                retain_idx1 = retain_idx1.astype(np.int64)


        dia = A.diagonal()
        raw_contig_counts = A.shape[0]
        retain_idx = np.where(dia >= self.min_cis_weight)[0]
        retain_idx = np.intersect1d(retain_idx, retain_idx1)
        contig_counts = len(retain_idx)
    
        if len(retain_idx) < raw_contig_counts:
            idx_to_vertices = self.idx_to_vertices
            A = A[retain_idx, :][:, retain_idx]
            self.H, _, _ = extract_incidence_matrix2(self.H, retain_idx)
  
            self.vertices = self.vertices[retain_idx]
            if self.NW is not None:
                self.NW = self.NW[retain_idx, :][:, retain_idx]

            if self.prune_pair_df is not None:
                self.prune_pair_df = self.prune_pair_df[
                    self.prune_pair_df['contig1'].isin(set(list(retain_idx))) &
                    self.prune_pair_df['contig2'].isin(set(list(retain_idx)))
                ]

                self.prune_pair_df['contig1'] = self.prune_pair_df['contig1'].map(idx_to_vertices.get)
                self.prune_pair_df['contig2'] = self.prune_pair_df['contig2'].map(idx_to_vertices.get)
                self.prune_pair_df = self.prune_pair_df.dropna(axis=0).reset_index(drop=True)
                self.prune_pair_df['contig1'] = self.prune_pair_df['contig1'].astype('str')
                self.prune_pair_df['contig2'] = self.prune_pair_df['contig2'].astype('str')
                vertices_idx = self.vertices_idx
                self.prune_pair_df['contig1'] = self.prune_pair_df['contig1'].map(
                    lambda x: vertices_idx.get(x, np.nan))
                self.prune_pair_df['contig2'] = self.prune_pair_df['contig2'].map(
                    lambda x: vertices_idx.get(x, np.nan))
                self.prune_pair_df = self.prune_pair_df.dropna(axis=0).reset_index(drop=True)
                self.P_allelic_idx = [self.prune_pair_df.loc[self.prune_pair_df['type'] == 0, 'contig1'],
                                        self.prune_pair_df.loc[self.prune_pair_df['type'] == 0, 'contig2']]
                self.P_weak_idx = [self.prune_pair_df.loc[self.prune_pair_df['type'] == 1, 'contig1'],
                                    self.prune_pair_df.loc[self.prune_pair_df['type'] == 1, 'contig2']]
            
            logger.info(f"Removed {raw_contig_counts - contig_counts:,} contigs that self edge weight < {self.min_cis_weight} (--min-cis-weight).")

            A = self.HG.clique_expansion_init(self.H, 
                                            # NW=self.NW,
                                            allelic_factor=self.allelic_factor, 
                                            min_weight=self.min_weight)
        vertices_idx = self.vertices_idx
        idx_to_vertices = self.idx_to_vertices

        A = HyperGraph.normalize_adjacency_matrix(A, self.NW)
        # self.H = None 
        raw_A = A
        if self.resolution1 == -1:
            result_K_length = len(self.K)
            tmp_resolution = self.init_resolution1
            max_round = 100 
            round_count = 0
            if not k:
                logger.warning("To automatic search best resolution, the `-n` must be specified.")
            if result_K_length < k:
                while result_K_length < k:
                    if round_count >= max_round:
                        logger.warning("Maximum rounds reached during automatic resolution search.")
                        break
                    # logger.info(f"Automatic search the best resolution ... {tmp_resolution:.1f}")
                   
                    self.cluster_assignments, self.K = HyperGraph.bipartite_leiden_clustering(
                        H=self.H,
                        resolution=tmp_resolution,
                        min_weight=self.min_weight,
                        P_allelic_idx=self.P_allelic_idx,
                        P_weak_idx=self.P_weak_idx,
                        allelic_factor=self.allelic_factor,
                        cross_allelic_factor=self.cross_allelic_factor
                    )
                    self.K = list(map(list, self.K))
                    self.K = self.filter_cluster(verbose=0)
                    result_K_length = len(self.K)
                    logger.info(f"Generated `{result_K_length}` groups at resolution `{tmp_resolution:.2f}`.")
                    tmp_resolution += 0.2
                    round_count += 1

            if result_K_length > k and not self.merge_use_allele and not self.disable_merge_in_first:
                while result_K_length > k:
                    if tmp_resolution <= 0 or round_count >= max_round:
                        logger.warning("Minimum resolution reached during automatic resolution search.")
                        break
                    logger.info(f"Automatic search the best resolution ... {tmp_resolution:.1f}")
                    self.cluster_assignments, self.K = HyperGraph.bipartite_leiden_clustering(
                        H=self.H,
                        resolution=tmp_resolution,
                        min_weight=self.min_weight,
                        P_allelic_idx=self.P_allelic_idx,
                        P_weak_idx=self.P_weak_idx,
                        allelic_factor=self.allelic_factor,
                        cross_allelic_factor=self.cross_allelic_factor
                    )
                    self.K = list(map(list, self.K))
                    self.K = self.filter_cluster(verbose=0)
                    result_K_length = len(self.K)
                    logger.info(f"Generated `{result_K_length}` groups at resolution `{tmp_resolution:.2f}`.")
                    tmp_resolution -= 0.1
                    round_count += 1

        else:
            self.cluster_assignments, self.K = HyperGraph.bipartite_leiden_clustering(
                H=self.H,
                resolution=self.resolution1,
                min_weight=self.min_weight,
                P_allelic_idx=self.P_allelic_idx,
                P_weak_idx=self.P_weak_idx,
                allelic_factor=self.allelic_factor,
                cross_allelic_factor=self.cross_allelic_factor
            )

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
        
            _, _results2 = zip(*_results2)
            _results2 = list_flatten(_results2)
            _results.extend(_results2)

            self.K = _results 
        
        if k and len(self.K) > k:  
            logger.info(f"Merging {len(self.K)} groups into {k} groups ...")
           
            if self.alleletable and self.merge_use_allele:
                at = AlleleTable(self.alleletable, fmt="allele2", sort=False)
                at_df = at.data
                at_df = at_df[at_df['similarity'] >= self.allelic_similarity]
                at_df[1] = at_df[1].map(vertices_idx.get)
                at_df[2] = at_df[2].map(vertices_idx.get)
                at_df.dropna(inplace=True)
                at_df[1] = at_df[1].astype(int)
                at_df[2] = at_df[2].astype(int)
                allelic_idx_set = set(map(tuple, at_df[[1, 2]].values))
                allele_use_method = "merge"
                
            else:
                allelic_idx_set = set()
                allele_use_method = "conflict"
            try:
                self.K = HyperPartition._merge(A, self.K, vertices_idx_sizes, k, 
                                            allelic_idx_set=allelic_idx_set,
                                            min_allelic_overlap=self.min_allelic_overlap, 
                                            allele_use_method=allele_use_method, 
                                            method='mean')
            except AttributeError:
                allelic_idx_set = set()
                self.K = HyperPartition._merge(A, self.K, vertices_idx_sizes, k, 
                                                allelic_idx_set=allelic_idx_set,
                                                min_allelic_overlap=self.min_allelic_overlap, 
                                                allele_use_method=allele_use_method, 
                                                method='sum')

        # self.remove_misassembly()

        if self.is_recluster_contigs:
            self.K = HyperPartition.recluster_contigs(self.K, k, A, raw_A,
                                            vertices_idx_sizes, 
                                            None, 
                                            self.prune_pair_df, self.allelic_similarity)
        
        reclustered_count = 1 if self.is_recluster_by_linkage else 0
        iter_count = 0
        while reclustered_count > 0 and iter_count < 5:
            iter_count += 1
            reclustered_count = self.recluster_by_linkage(A=raw_A, min_ratio=1.5, iter_count=iter_count)
                
            
        if self.refine:
            logger.info("Refining misassembies of allelic contig pairs, which partition into the same group")
            self.K = HyperPartition.refine_allelic_errors(self.K, k, A, 
                                                        vertices_idx_sizes, 
                                                        self.prune_pair_df,
                                                        min_weight=self.min_weight)

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
        total_length = 0
        for i, j in zip(group_info, length_contents):
            table.add_row(str(i), j)
            total_length += int(j.replace(",", ""))
        
        if hasattr(table, "add_section"):
                table.add_section()
        table.add_row("Total", f"{total_length:,}")

        with console_html.capture() as capture:
            console_html.print(table)

        logger.info(capture.get().strip())
        
        
        return self.K  
    
    def single_partition(self, k=None, sort_by_length=True, sort_group=False, is_kprune=False):
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


        if is_kprune:
            if self.fasta is not None and self.alleletable is None and self.prunetable is None:
                c = 100
                fasta_prefix = Path(Path(self.fasta).name)
                while fasta_prefix.suffix in {".fasta", "gz", "fa", '.fa', '.gz', ".split"}:
                    fasta_prefix = fasta_prefix.with_suffix("")
                    
                if self.split:
                    
                    _contigsizes = self.contigsizes[self.contigsizes['length'] > self.min_length]
                    from .utilities import get_fasta_from_split_contig
                    used_fasta = get_fasta_from_split_contig(
                        self.fasta,
                        _contigsizes,
                        f"{fasta_prefix}.split.fasta",
                    )
                    
                else:
                    used_fasta = self.fasta

                self.alleles(used_fasta, 
                                                    k=self.alleles_kmer_size,
                                                    w=self.alleles_window_size,
                                                    d=self.alleles_diff_thres,
                                                    m=self.alleles_minimum_similarity,
                                                    tl=self.alleles_trim_length,
                                                    c=c,
                                                    split=self.split,
                                                    threads=self.threads)
                if self.split:
                    self.alleletable = f"{fasta_prefix}.split.split.allele.table"
                else:
                    self.alleletable = f"{fasta_prefix}.allele.table"
            
            if self.alleletable and not self.prunetable:
                self.kprune(self.alleletable, contacts=self.contacts)
        
    
        A = self.HG.clique_expansion_init(self.H, 
                                        #   NW=self.NW,
                                          allelic_factor=self.allelic_factor, 
                                          min_weight=self.min_weight)

        retain_idx1 = HyperGraph.filter_adjacency_matrix(A, self.vertices, self.vertices_idx_sizes, invert=True)
        if retain_idx1 is None:
            retain_idx1 = np.arange(A.shape[0], dtype=np.int64)
        else:
            retain_idx1 = np.asarray(retain_idx1)
            if retain_idx1.dtype == bool:
                retain_idx1 = np.where(retain_idx1)[0].astype(np.int64)
            else:
                retain_idx1 = retain_idx1.astype(np.int64)


        dia = A.diagonal()
        raw_contig_counts = A.shape[0]
        retain_idx = np.where(dia >= self.min_cis_weight)[0]
        retain_idx = np.intersect1d(retain_idx, retain_idx1)
        contig_counts = len(retain_idx)
    
        if len(retain_idx) < raw_contig_counts:
            idx_to_vertices = self.idx_to_vertices
            A = A[retain_idx, :][:, retain_idx]
            self.H, _, _ = extract_incidence_matrix2(self.H, retain_idx)
  
            self.vertices = self.vertices[retain_idx]
            if self.NW is not None:
                self.NW = self.NW[retain_idx, :][:, retain_idx]

            if self.prune_pair_df is not None:
                self.prune_pair_df = self.prune_pair_df[
                    self.prune_pair_df['contig1'].isin(set(list(retain_idx))) &
                    self.prune_pair_df['contig2'].isin(set(list(retain_idx)))
                ]

                self.prune_pair_df['contig1'] = self.prune_pair_df['contig1'].map(idx_to_vertices.get)
                self.prune_pair_df['contig2'] = self.prune_pair_df['contig2'].map(idx_to_vertices.get)
                self.prune_pair_df = self.prune_pair_df.dropna(axis=0).reset_index(drop=True)
                self.prune_pair_df['contig1'] = self.prune_pair_df['contig1'].astype('str')
                self.prune_pair_df['contig2'] = self.prune_pair_df['contig2'].astype('str')
                vertices_idx = self.vertices_idx
                self.prune_pair_df['contig1'] = self.prune_pair_df['contig1'].map(
                    lambda x: vertices_idx.get(x, np.nan))
                self.prune_pair_df['contig2'] = self.prune_pair_df['contig2'].map(
                    lambda x: vertices_idx.get(x, np.nan))
                self.prune_pair_df = self.prune_pair_df.dropna(axis=0).reset_index(drop=True)
                self.P_allelic_idx = [self.prune_pair_df.loc[self.prune_pair_df['type'] == 0, 'contig1'],
                                        self.prune_pair_df.loc[self.prune_pair_df['type'] == 0, 'contig2']]
                self.P_weak_idx = [self.prune_pair_df.loc[self.prune_pair_df['type'] == 1, 'contig1'],
                                    self.prune_pair_df.loc[self.prune_pair_df['type'] == 1, 'contig2']]
            
            logger.info(f"Removed {raw_contig_counts - contig_counts:,} contigs that self edge weight < {self.min_cis_weight} (--min-cis-weight).")

            A = self.HG.clique_expansion_init(self.H, 
                                            # NW=self.NW,
                                            allelic_factor=self.allelic_factor, 
                                            min_weight=self.min_weight)
        vertices_idx = self.vertices_idx
        idx_to_vertices = self.idx_to_vertices

        A = HyperGraph.normalize_adjacency_matrix(A, self.NW)
        self.H = None 

        if self.resolution1 == -1:
            result_K_length = len(self.K)
            tmp_resolution = self.init_resolution1
            max_round = 100 
            round_count = 0
            if not k:
                logger.warning("To automatic search best resolution, the `-n` must be specified.")
            if result_K_length < k:
                while result_K_length < k:
                    if round_count >= max_round:
                        logger.warning("Maximum rounds reached during automatic resolution search.")
                        break
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
                                                            threads=self.threads,
                                                            method=self.cluster_method)
                    self.K = list(map(list, self.K))
                    self.K = self.filter_cluster(verbose=0)
                    result_K_length = len(self.K)
                    logger.info(f"Generated `{result_K_length}` groups at resolution `{tmp_resolution:.2f}`.")
                    tmp_resolution += 0.2
                    round_count += 1

            if result_K_length > k and not self.merge_use_allele and not self.disable_merge_in_first:
                while result_K_length > k:
                    if tmp_resolution <= 0 or round_count >= max_round:
                        logger.warning("Minimum resolution reached during automatic resolution search.")
                        break
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
                                                            threads=self.threads,
                                                            method=self.cluster_method)
                    self.K = list(map(list, self.K))
                    self.K = self.filter_cluster(verbose=0)
                    result_K_length = len(self.K)
                    logger.info(f"Generated `{result_K_length}` groups at resolution `{tmp_resolution:.2f}`.")
                    tmp_resolution -= 0.1
                    round_count += 1

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
                                                        threads=self.threads,
                                                        method=self.cluster_method)


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
        
            _, _results2 = zip(*_results2)
            _results2 = list_flatten(_results2)
            _results.extend(_results2)

            self.K = _results 
        
        if k and len(self.K) > k:  
            logger.info(f"Merging {len(self.K)} groups into {k} groups ...")
           
            if self.alleletable and self.merge_use_allele:
                at = AlleleTable(self.alleletable, fmt="allele2", sort=False)
                at_df = at.data
                at_df = at_df[at_df['similarity'] >= self.allelic_similarity]
                at_df[1] = at_df[1].map(vertices_idx.get)
                at_df[2] = at_df[2].map(vertices_idx.get)
                at_df.dropna(inplace=True)
                at_df[1] = at_df[1].astype(int)
                at_df[2] = at_df[2].astype(int)
                allelic_idx_set = set(map(tuple, at_df[[1, 2]].values))
                allele_use_method = "merge"
                
            else:
                allelic_idx_set = set()
                allele_use_method = "conflict"
            try:
                self.K = HyperPartition._merge(A, self.K, vertices_idx_sizes, k, 
                                            allelic_idx_set=allelic_idx_set,
                                            min_allelic_overlap=self.min_allelic_overlap, 
                                            allele_use_method=allele_use_method, 
                                            method='mean')
            except AttributeError:
                allelic_idx_set = set()
                self.K = HyperPartition._merge(A, self.K, vertices_idx_sizes, k, 
                                                allelic_idx_set=allelic_idx_set,
                                                min_allelic_overlap=self.min_allelic_overlap, 
                                                allele_use_method=allele_use_method, 
                                                method='sum')

        # self.remove_misassembly()

        if self.is_recluster_contigs:
            self.K = HyperPartition.recluster_contigs(self.K, k, A, raw_A,
                                            vertices_idx_sizes, 
                                            None, 
                                            self.prune_pair_df, self.allelic_similarity)
        
        reclustered_count = 1 if self.is_recluster_by_linkage else 0
        iter_count = 0
        while reclustered_count > 0 and iter_count < 5:
            iter_count += 1
            reclustered_count = self.recluster_by_linkage(A=raw_A, min_ratio=1.5, iter_count=iter_count)
                
            
        if self.refine:
            logger.info("Refining misassembies of allelic contig pairs, which partition into the same group")
            self.K = HyperPartition.refine_allelic_errors(self.K, k, A, 
                                                        vertices_idx_sizes, 
                                                        self.prune_pair_df,
                                                        min_weight=self.min_weight)

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
        total_length = 0
        for i, j in zip(group_info, length_contents):
            table.add_row(str(i), j)
            total_length += int(j.replace(",", ""))
        
        if hasattr(table, "add_section"):
                table.add_section()
        table.add_row("Total", f"{total_length:,}")

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
                                H, vertices_idx_sizes, NW, cis_counts, resolution, 
                                output_dir="./", init_resolution=0.8, min_weight=1, 
                                min_cis_weight = 1.0, allelic_similarity=0.8, 
                                min_allelic_overlap=0.1, allelic_factor=-1, cross_allelic_factor=0.0,
                                is_remove_misassembly=False, is_recluster_contigs=False, refine=False,
                                collapsed_contigs_df=None,
                                min_scaffold_length=10000, threshold=0.01, max_round=1, num=None,
                                kprune_norm_method='cis', count_re=None, threads=1, cluster_method="louvain"):
        """
        single function for incremental_partition.
        """
        for target_logger in ["cphasing.core", "cphasing.utilities"]:
            logging.getLogger(target_logger).setLevel(logging.ERROR)
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
        sub_cis_counts = cis_counts[sub_edge_idx] if cis_counts is not None else None
        ## remove low weigth contigs

        # sub_A = HyperGraph.clique_expansion_init(sub_H, min_weight=min_weight)
        # dia = sub_A.diagonal()
        dia = np.asarray(sub_H.power(2).sum(axis=1)).ravel()
        raw_contig_counts = len(K)
        retain_mask = dia >= min_cis_weight
        K = K[retain_mask]
        contig_counts = len(K)
        if contig_counts == 0:
            return None, None, None
        if (raw_contig_counts - contig_counts) > 0:
            logger.info(f"Removed {raw_contig_counts - contig_counts:,} contigs that self edge weight < {min_cis_weight} (--min-cis-weight).")
            # sub_H, _, sub_edge_idx = extract_incidence_matrix2(H, K)
            # sub_cis_counts = cis_counts[sub_edge_idx] if cis_counts is not None else None
            sub_H = sub_H[retain_mask, :]
            active_edge_mask = np.asarray(sub_H.sum(axis=0)).ravel() > 0
            sub_H = sub_H[:, active_edge_mask]
            if sub_cis_counts is not None:
                sub_cis_counts = sub_cis_counts[active_edge_mask]
    
    
        del H 

        if NW is not None:
            sub_NW = NW[K][:, K]
        else:
            sub_NW = None

        sub_old2new_idx = dict(zip(K, range(len(K))))
        
        sub_vertices_idx_sizes = vertices_idx_sizes.reindex(K)
        sub_vertices_new_idx_sizes = sub_vertices_idx_sizes
        sub_vertices_new_idx_sizes.index = sub_vertices_idx_sizes.index.map(sub_old2new_idx.get)
        
        sub_size_db = sub_vertices_new_idx_sizes.to_dict()
        if 'length' in sub_size_db:
            sub_size_db = sub_size_db['length']

        sub_vertices = np.array(list(idx_to_vertices[i] for i in K))
        sub_vertives_idx = dict(zip(sub_vertices, range(len(sub_vertices))))

        if collapsed_contigs_df is not None:
            sub_collapsed_contigs_df = collapsed_contigs_df[collapsed_contigs_df[0].isin(sub_vertices)]
            sub_collapsed_contigs_df[0] = sub_collapsed_contigs_df[0].map(sub_vertives_idx.get)
       
        if alleletable:
            ## extract sub alleletable
            at = AlleleTable(alleletable, fmt="allele2", sort=False)
            at_df = at.data
            at_df = at_df[
                at_df[1].isin(sub_vertices) &
                at_df[2].isin(sub_vertices)
            ]
            at.data = at_df
            sub_alleletable = str(f"{output_dir}/kprune_workdir/{num}.alleletable")
            at.save(sub_alleletable, index=False)

            sub_output_contacts = f"{output_dir}/kprune_workdir/{num}.hypergraph.expansion.contacts"
            _sub_A = HyperGraph.clique_expansion_init(sub_H, 
                                                    #   cis_count=sub_cis_counts, 
                                                      min_weight=min_weight, use_high_order=False)
            HyperGraph.to_contacts(_sub_A, sub_vertices, NW=sub_NW, min_weight=min_weight,
                                output=sub_output_contacts)

            sub_P_allelic_idx, sub_P_weak_idx, sub_prune_pair_df = HyperPartition.kprune_func(
                                    sub_alleletable, sub_output_contacts,
                                    sub_vertives_idx, 
                                    output_dir=f"{output_dir}/kprune_workdir",
                                    kprune_norm_method=kprune_norm_method,
                                    count_re=count_re,
                                    threads=threads)

        elif  prune_pair_df is not None:      

            # sub_prune_pair_df = prune_pair_df.reindex(list(permutations(K, 2))).dropna().reset_index()
            K_set = set(K)
            sub_prune_pair_df = prune_pair_df[
                prune_pair_df.index.get_level_values(0).isin(K_set) & 
                prune_pair_df.index.get_level_values(1).isin(K_set)
            ].reset_index()

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
        
        # # reset sub_H value to 1 
        # sub_H.data = np.ones_like(sub_H.data)
        sub_A = HyperGraph.clique_expansion_init(sub_H, 
                                                 cis_count=sub_cis_counts,
                                                 min_weight=min_weight)
        if sub_A.dtype != np.float32:
            sub_A = sub_A.astype(np.float32)
        
        sub_A = HyperGraph.normalize_adjacency_matrix(sub_A, sub_NW)


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
                                            method=cluster_method)
        
            new_K = list(filter(
                                    lambda x: sub_vertices_new_idx_sizes.reindex(x).sum().values[0] \
                                        >= min_scaffold_length, new_K)
                    )
            result_K_length = len(new_K)
            auto_round = 2

            if result_K_length < k:
                tmp_resolution += 0.2
                while result_K_length < k:
                    if tmp_resolution > MAX_RESOLUTION or auto_round > 50:
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
                                                        method=cluster_method)
                    new_K = list(filter(
                                    lambda x: sub_vertices_new_idx_sizes.reindex(x).sum().values[0] \
                                        >= min_scaffold_length, new_K)
                    )
                
                    new_K = list(map(list, new_K))
                    tmp_resolution += 0.2  
                    
                    result_K_length = len(new_K)
                    auto_round += 1


            elif result_K_length > k:
                logger.debug(f"Automatic search the best resolution down ...")
                tmp_resolution = init_resolution - 0.1
                auto_round = 1
                while result_K_length > k:
                    if tmp_resolution < 0.01 or auto_round > 50:
                        break
                    # logger.info(f"Automaticly to search best resolution ... {tmp_resolution}")
                    old_K = new_K
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
                                                        method=cluster_method)
                    new_K = list(filter(
                                    lambda x: sub_vertices_new_idx_sizes.reindex(x).sum().values[0] \
                                        >= min_scaffold_length, new_K)
                    )
                    new_K = list(map(list, new_K))
                    if len(new_K) < k:
                        new_K = old_K
                        break

                    tmp_resolution -= 0.1
                    if tmp_resolution < 0.1:
                        tmp_resolution -= 0.02
              
                    result_K_length = len(new_K)
                    auto_round += 1

             
        else:
            logger.debug(f"Using specified resolution: {resolution}")
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
                                                    method=cluster_method)

        
    
        
        ## remove the scaffold that is too short
        # before_filter_K = new_K.copy()
        # new_K = list(map(list, new_K))
        # new_K = list(map(list, filter(
        #                 lambda x: sub_vertices_new_idx_sizes.loc[list(x)].sum().values[0] \
        #                     >= min_scaffold_length, new_K)))
        new_K = list(map(list, filter(
                        lambda x: sum(sub_size_db.get(item, 0) for item in x) >= min_scaffold_length, 
                        new_K)))
        # filtered_K = list(filter(
        #                         lambda x: sub_vertices_new_idx_sizes.reindex(x).sum().values[0] \
        #                             < min_scaffold_length, new_K)
        #         )
        
        
        if k and len(new_K) > k:
            new_K = sorted(new_K, key=lambda x: sub_vertices_new_idx_sizes.loc[x]['length'].sum())
            logger.info(f"FirstGroup{num}: Merging {len(new_K)} groups into {k} groups ...")
            if sub_prune_pair_df is not None:
                allelic_idx_set = set(map(tuple, 
                                        sub_prune_pair_df[(sub_prune_pair_df['type'] == 0) & 
                                                            (sub_prune_pair_df['similarity'] >= allelic_similarity)]
                                        [['contig1', 'contig2']].values)
                                        )

            else:
                allelic_idx_set = set() 

            new_K = HyperPartition._merge(A, new_K, sub_vertices_new_idx_sizes, k, 
                                            allelic_idx_set=allelic_idx_set,
                                            min_allelic_overlap=min_allelic_overlap,
                                            allele_use_method='conflict',
                                            method=merge_method,
                                            )
        

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
            

        sub_new2old_idx = dict(zip(range(len(K)), K))
        new_K = list(map(lambda x: list(map(lambda y: sub_new2old_idx[y], x)), new_K))

        new_K = sorted(new_K, key=lambda x: vertices_idx_sizes.reindex(x).dropna()['length'].sum(), reverse=True)
        logger.info(f"Subgroup {num} result {len(new_K)} groups.")
        
        return A, cluster_assignments, new_K
    

    def kprune(self, alleletable, first_cluster_file=None, contacts=None, is_run=True):
        if is_run:
            if not contacts or not Path(contacts).exists():
                A = HyperGraph.clique_expansion_init(self.H, cis_count=self.cis_counts, use_high_order=False, min_weight=self.min_weight)
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
            if self.count_re:
                cmd.append("--count-re")
                cmd.append(str(self.count_re))
        
            logger.info("Generating the prune table ...")
            flag = run_cmd(cmd, log=f"{self.log_dir}/hyperpartition_kprune.log")
            assert flag == 0, "Failed to execute command, please check log."

            logger.info(f"Prune table stored in `{kprune_output_file}`")
        else:
            logger.warning(f"Load exists prune table `{kprune_output_file}`")
            
        self.prunetable = kprune_output_file
        

        self.P_allelic_idx, self.P_weak_idx, self.prune_pair_df = self.get_prune_pairs()

    @staticmethod
    def alleles(fasta, first_cluster=None, 
                k=19, w=19, m=0.5, d=0.2, c=60, 
                tl=25000, output=None,
                split=False, threads=4):
        from .cli import alleles
        import json 
        import hashlib
        current_params = {
            "fasta": str(Path(fasta).absolute()),
            "k": k, "w": w, "m": m, "d": d, "c": c,
            "tl": tl, "split": split
        }
 
        param_str = json.dumps(current_params, sort_keys=True)
        param_hash = hashlib.md5(param_str.encode()).hexdigest()

        fasta_path = Path(fasta)
        fasta_prefix = Path(fasta_path.name).with_suffix("")
        while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz", ".split", "split"}:
            fasta_prefix = fasta_prefix.with_suffix("")

        if not output:
            output = f"{fasta_prefix}.split.allele.table" if split else f"{fasta_prefix}.allele.table"

        cache_info_file = Path(f".{output}").with_suffix(".params")

        is_cache_valid = False
        if Path(output).exists() and cache_info_file.exists():
            try:
                with open(cache_info_file, 'r') as f:
                    cached_data = json.load(f)
                    if cached_data.get("param_hash") == param_hash:
                        if Path(output).stat().st_mtime > Path(fasta).stat().st_mtime:
                            is_cache_valid = True
            except Exception:
                pass
        if is_cache_valid:
            logger.info(f"Found existing allele table `{output}` with same parameters, skipping.")
            return str(Path(output).absolute())

        with open(cache_info_file, 'w') as f:
            json.dump({"param_hash": param_hash, "params": current_params}, f, indent=4)

        args = [
                    "-f",
                    str(fasta),
                    "-t",
                    str(threads),
                    "-c", str(c),
                    "-k", str(k),
                    "-w", str(w),
                    "-m", str(m),
                    "-d", str(d),
                    "-tl", str(tl),
                ]
        if first_cluster:
            args.append("-fc")
            args.append("first.clusters.txt")

        if split:
            args.append("--split")
        if output:
            args.append("-o")
            args.append(output)
        try:
            alleles.main(
                args=args,
                prog_name='alleles'
            )

        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
        

        if not output:
            if split:
                return f"{fasta_prefix}.split.allele.table"
            else:
                return f"{fasta_prefix}.allele.table"
        else:
            return output

    @staticmethod
    def kprune_func(alleletable, 
                    contacts,
                    vertices_idx,
                    output_dir="./",
                    first_cluster_file=None, 
                    count_re=None,
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
        
        if count_re:
            cmd.append("--count-re")
            cmd.append(str(count_re))

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
        logger.debug(f"Resolution1: {self.resolution1} InitResolution1: {self.init_resolution1} Resolution2: {self.resolution2} InitResolution2: {self.init_resolution2}")
        if not first_cluster:
            logger.info("Starting first round partition ...")
        
        self.P_allelic_idx, self.P_weak_idx = None, None
        A = HyperGraph.clique_expansion_init(self.H, 
                                            min_weight=self.min_weight,
                                            # use_high_order=False
                                            )

        retain_idx1 = HyperGraph.filter_adjacency_matrix(
            A, self.vertices, self.vertices_idx_sizes, invert=True
        )
        if retain_idx1 is None:
            retain_idx1 = np.arange(A.shape[0], dtype=np.int64)
        else:
            retain_idx1 = np.asarray(retain_idx1)
            if retain_idx1.dtype == bool:
                retain_idx1 = np.where(retain_idx1)[0].astype(np.int64)
            else:
                retain_idx1 = retain_idx1.astype(np.int64)

        dia = A.diagonal()

        raw_contig_counts = A.shape[0]
        retain_idx = np.where(dia >= self.min_cis_weight)[0]
        retain_idx = np.intersect1d(retain_idx, retain_idx1)
        contig_counts = len(retain_idx)
        
        if len(retain_idx) < raw_contig_counts:
            A = A[retain_idx, :][:, retain_idx]
            self.H, _, sub_edge_idx = extract_incidence_matrix2(self.H, retain_idx)
            if self.cis_counts is not None:
                try:
                    self.cis_counts = self.cis_counts[sub_edge_idx]
                except IndexError:
                    
                    logger.warning("IndexError when updating cis_counts. Resetting cis_counts to None.")
                    self.cis_counts = None

            if self.NW is not None:
                self.NW = self.NW[retain_idx, :][:, retain_idx]
            self.vertices = self.vertices[retain_idx]
            logger.info(f"Removed {raw_contig_counts - contig_counts:,} contigs that self edge weight < {self.min_cis_weight} (--min-cis-weight).")

        vertices_idx_sizes = self.vertices_idx_sizes
        vertices_idx_sizes = pd.DataFrame(vertices_idx_sizes, index=['length']).T
        vertices_idx = self.vertices_idx
        
        A = HyperGraph.normalize_adjacency_matrix(A, self.NW)
        
        if self.merge_use_allele:
            min_similarity = 0.90
            min_shared = 100
            
            if self.fasta is not None:
                logger.info("Use allelic information for implementing 1st-round merging...")
                if k[1] is not None and k[1] != 0:
                    c = k[1] * 4
                else:
                    c = 60
                fasta_prefix = Path(Path(self.fasta).stem)

                while fasta_prefix.suffix in {".fasta", "gz", "fa", '.fa', '.gz'}:
                    fasta_prefix = fasta_prefix.with_suffix("")
                fasta_stem = fasta_prefix
                _contigsizes = self.contigsizes[self.contigsizes['length'] > self.min_length]
                if self.split:
                    split_fasta_path = get_fasta_from_split_contig(
                        self.fasta,
                        _contigsizes,
                        f"{fasta_stem}.split.fasta",
                    )
                    self.alleles(
                        split_fasta_path,
                        k=self.alleles_kmer_size,
                        w=self.alleles_window_size,
                        d=self.alleles_diff_thres,
                        tl=self.alleles_trim_length,
                        output=f"{fasta_stem}.allele.table",
                        split=True,
                        threads=self.threads,
                        c=c,
                    )
                    fast_alleletable = f"{fasta_stem}.split.allele.table"

                    if Path(split_fasta_path).exists():
                        logger.debug(f"Removing temporary split fasta file: {split_fasta_path}")
                        os.remove(split_fasta_path)
                    if Path(f"{split_fasta_path}.fai").exists():
                        logger.debug(f"Removing temporary split fasta index file: {split_fasta_path}.fai")
                        os.remove(f"{split_fasta_path}.fai")

                else:
                    fast_alleletable = self.alleles(
                        self.fasta,
                        k=self.alleles_kmer_size,
                        w=self.alleles_window_size,
                        d=self.alleles_diff_thres,
                        tl=self.alleles_trim_length,
                        output=f"{fasta_stem}.allele.table",
                        threads=self.threads,
                        c=c,
                    )

                at = AlleleTable(fast_alleletable, fmt="allele2", sort=False)
                at_df = at.data
                if at_df.empty:
                    at_df = None
                    allelic_idx_set = None
                    self.alleletable = None 
                else:
                    at_df = at_df[at_df['similarity'] >= min_similarity]
                    at_df = at_df[(at_df['mzShared'] >= min_shared)]
                    at_df[1] = at_df[1].map(vertices_idx.get)
                    at_df[2] = at_df[2].map(vertices_idx.get)
                    at_df = at_df.dropna()
                    at_df[1] = at_df[1].astype(int)
                    at_df[2] = at_df[2].astype(int)
                    self.alleletable = Path(fast_alleletable).absolute()
                    allelic_idx_set = set(map(tuple, at_df[[1, 2]].values))
            elif self.prune_pair_df is not None:
                logger.info("Use allelic information for implementing 1st-round merging...")
                allelic_idx_set = set(map(tuple, 
                    self.prune_pair_df[(self.prune_pair_df['type'] == 0) & 
                                        (self.prune_pair_df['similarity'] >= min_similarity)]
                    [['contig1', 'contig2']].values)
                    )
                at_df = None
            elif self.alleletable:
                logger.info("Use allelic information for implementing 1st-round merging...")
                at = AlleleTable(self.alleletable, fmt="allele2", sort=False)
                at_df = at.data
                at_df = at_df[at_df['similarity'] >= min_similarity]
                at_df = at_df[(at_df['mzShared'] >= min_shared)]
                at_df[1] = at_df[1].map(vertices_idx.get)
                at_df[2] = at_df[2].map(vertices_idx.get)
                at_df.dropna(inplace=True)
                at_df[1] = at_df[1].astype(int)
                at_df[2] = at_df[2].astype(int)
                allelic_idx_set = set(map(tuple, at_df[[1, 2]].values))
                
            else:
                allelic_idx_set = set()
                at_df = None
        else:
            allelic_idx_set = None
            at_df = None

        raw_A = A.copy()
        if allelic_idx_set and self.merge_use_allele:
            P_allelic_idx = list(zip(*allelic_idx_set))[0], list(zip(*allelic_idx_set))[1]
            if self.allelic_positive_factor != 1.0:
                A = HyperGraph.reweight_adjacency_matrix(A, P_allelic_idx, allelic_factor=self.allelic_positive_factor)
              

        if not first_cluster:
            if k[0] == 1:
                self.K = [list(range(A.shape[0]))]
                _A = A
                logger.info("Only one group specified, all contigs assigned into one group.")
            else:
                if self.resolution1 < 0 :
                    tmp_resolution = self.init_resolution1
                    logger.info(f"Automatic search for best resolution ... {tmp_resolution:.1f}")
                    
                    _, _A, _, self.K = IRMM(self.H, A, self.NW, 
                                        None, None, self.allelic_factor, 
                                        self.cross_allelic_factor, tmp_resolution,
                                        self.min_weight, self.threshold, 
                                        1, threads=self.threads,
                                        method=self.cluster_method)
                    self.K = list(map(list, self.K))
                    self.K = self.filter_cluster(verbose=0)
                    result_K_length = len(self.K)
                    logger.info(f"Generated `{result_K_length}` groups at resolution `{tmp_resolution:.1f}`.")
                    auto_round = 2

                    if result_K_length > k[0] and not self.merge_use_allele and not self.disable_merge_in_first:
                        tmp_resolution = self.init_resolution1 - 0.1
                        auto_round = 1
                        while result_K_length > k[0]:
                            if tmp_resolution < 0.01 or auto_round > 50:
                                break
                            # logger.info(f"Automatic search for best resolution ... {tmp_resolution:.1f}")
                            _, _A, _, new_K = IRMM(self.H, A, self.NW, 
                                    None, None, self.allelic_factor, 
                                        self.cross_allelic_factor, tmp_resolution, 
                                        self.min_weight, self.threshold, 
                                        1, threads=self.threads,
                                        method=self.cluster_method)

                            new_K = list(map(list, new_K))
                            old_K = self.K
                            self.K = new_K
                            self.K = self.filter_cluster(verbose=0)
                            if len(self.K) < k[0]:
                                self.K = old_K
                                break

                            assert len(self.K) > 0, "Couldn't run first cluster."
                            result_K_length = len(self.K)
                            logger.info(f"Generated `{result_K_length}` groups at resolution `{tmp_resolution:.1f}`.")
                            
                            tmp_resolution -= 0.1
                            if tmp_resolution < 0.1:
                                tmp_resolution -= 0.02
                            auto_round += 1
                
                    elif result_K_length < k[0]:
                        while result_K_length < k[0] and k[0] != 0:
                            if tmp_resolution > MAX_RESOLUTION or auto_round > 50:
                                break
                            
                            _, A, _, self.K = IRMM(self.H, A, self.NW, 
                                    None, None, self.allelic_factor, 
                                        self.cross_allelic_factor, tmp_resolution, 
                                        self.min_weight, self.threshold, 
                                        1, threads=self.threads,
                                        method=self.cluster_method)

                            self.K = list(map(list, self.K))
                            self.K = self.filter_cluster(verbose=0)
                            assert len(self.K) > 0, "Couldn't run first cluster."
                            
                            result_K_length = len(self.K)
                            logger.info(f"Generated `{result_K_length}` groups at resolution `{tmp_resolution:.1f}`.")
                            tmp_resolution += 0.2
                            auto_round += 1

                else:
                    _, _A, _, self.K = IRMM(self.H, A, self.NW, 
                                None, None, self.allelic_factor, 
                                    self.cross_allelic_factor, self.resolution1, 
                                    self.min_weight, self.threshold, 
                                    1, threads=self.threads,
                                    method=self.cluster_method)

            A = _A
            self.K = list(map(list, self.K))

            if self.split:
                _, self.K, _ = self.merge_split()
            self.K = self.filter_cluster()

            if k[0] and len(self.K) > k[0]:
                if self.disable_merge_in_first:
                    self.K = sorted(
                        self.K,
                        key=lambda x: vertices_idx_sizes.loc[x]["length"].sum(),
                        reverse=True,
                    )
                    logger.info(
                        f"Discarded `{len(self.K) - k[0]}` smallest groups (--disable-merge-in-first)."
                    )
                    self.K = self.K[: k[0]]
                else:
                    
                    if allelic_idx_set is not None and self.merge_use_allele:
                        self.K = HyperPartition._merge(A, self.K, vertices_idx_sizes, k[0], 
                                                        allelic_idx_set=allelic_idx_set,
                                                        alleletable_df=at_df,
                                                        min_allelic_overlap=0.8,
                                                        allele_use_method='merge',
                                                        method='sum')
                    else:
                        self.K = HyperPartition._merge(A, self.K, vertices_idx_sizes, k[0], 
                                                        
                                                        min_allelic_overlap=self.min_allelic_overlap,

                                                        method='sum')

            self.K = sorted(self.K, key=lambda x: vertices_idx_sizes.loc[x]['length'].sum(), reverse=True)

            reclustered_count = 1 if self.is_recluster_by_linkage else 0
            iter_count = 0
            while reclustered_count > 0 and iter_count < 5:
                iter_count += 1
                reclustered_count = self.recluster_by_linkage(A=raw_A, min_ratio=1.5, iter_count=iter_count)
                
            

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
            if hasattr(table, "add_section"):
                table.add_section()
            table.add_row("Total", f"{total_length:,}")
            
            # first_length_contents = pformat(list(zip(first_group_info, length_contents)))

            logger.info(f"First hyperpartition resulted in {len(self.K)} groups:")
                        # f"{first_length_contents}")
            with console_html.capture() as capture:
                console_html.print(table)
            logger.info(capture.get())   
            mapq_filter = True

            self.to_cluster(f'first.clusters.txt', merge=False)
          
            first_cluster_file = f'first.clusters.txt'
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
        if self.HG.edges.mapq.size and (self.min_quality2 > self.min_quality1) and mapq_filter:
            
            idx_to_vertices = self.idx_to_vertices

            tmp_K = list(map(lambda y: list(
                            map(lambda x: idx_to_vertices[x], y)), 
                            self.K))
            
            self.HG = HyperGraph(self.edges, min_quality=self.min_quality2)

            del self.edges
            gc.collect()
            self.filter_hypergraph()
            self.H, self.vertices = self.get_hypergraph()

            raw_A = HyperGraph.clique_expansion_init(self.H, min_weight=self.min_weight, use_high_order=False)
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

        if self.prunetable:
            self.P_allelic_idx, self.P_weak_idx, self.prune_pair_df = self.get_prune_pairs()
            prune_pair_df = self.prune_pair_df.reset_index().set_index(['contig1', 'contig2'])

        elif self.alleletable:
            # is_run = False if isinstance(first_cluster, ClusterTable) else True
            # is_run = True
            # self.kprune(self.alleletable, first_cluster_file, contacts=self.contacts, is_run=is_run)
            # prune_pair_df = self.prune_pair_df.reset_index().set_index(['contig1', 'contig2'])
            prune_pair_df = None
        else:
            prune_pair_df = None

        total_contig_count = sum(len(sub_k) for sub_k in self.K)
        if total_contig_count == 0:
            total_contig_count = 1

        args = []
        results = []
        self.exclude_groups = []
        idx_to_vertices = self.idx_to_vertices
  
        if self.fasta is not None and self.alleletable is None and self.prunetable is None:
            if k[1] is not None and k[1] != 0:
                c = k[1] * 4
            else:
                c = 100
            
            self.alleletable = str(Path(self.alleles(self.fasta, first_cluster_file, 
                                                    k=self.alleles_kmer_size,
                                                    w=self.alleles_window_size,
                                                    d=self.alleles_diff_thres,
                                                    m=self.alleles_minimum_similarity,
                                                    tl=self.alleles_trim_length,
                                                    c=c,
                                                    split=self.split,
                                                    threads=self.threads)).absolute())


            
        Path("kprune_workdir").mkdir(exist_ok=True)
        current_dir = Path.cwd()

        for num, sub_k in enumerate(self.K, 1):
            sub_threads = max(2, int((len(sub_k) / total_contig_count) * self.threads))
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
                        self.H, vertices_idx_sizes, self.NW, self.cis_counts,
                        self.resolution2, current_dir, self.init_resolution2, 
                        self.min_weight, self.min_cis_weight,
                        self.allelic_similarity,  self.min_allelic_overlap, 
                        self.allelic_factor, self.cross_allelic_factor, self.is_remove_misassembly,
                        self.is_recluster_contigs, self.refine, self.collapsed_contigs_df,
                        self.min_scaffold_length, self.threshold, self.max_round, num,
                        self.kprune_norm_method, self.count_re,
                        sub_threads, self.cluster_method))
            
            # results.append(HyperPartition._incremental_partition(args[-1])

        with silence_logs(["cphasing.algorithms.hypergraph", "cphasing.hyperpartition"], logging.WARNING):
            with parallel_backend('loky', n_jobs=min(self.threads, len(args))):
                try:
                    results = Parallel(n_jobs=min(self.threads, len(args)), return_as="generator")(
                                delayed(HyperPartition._incremental_partition)(*a) for a in args)
                except TypeError:
                    results = Parallel(n_jobs=min(self.threads, len(args)))(
                                delayed(HyperPartition._incremental_partition)(*a) for a in args)
                
                results = list(filter(lambda x: x[2] is not None, results))

            for target_logger in ["cphasing.core", "cphasing.utilities"]:
                logging.getLogger(target_logger).setLevel(logging.INFO)

        if results:
            _, self.cluster_assignments, results = zip(*results)
        else:
            _, self.cluster_assignments, results = [], [], []
        
        if 'args' in locals():
            del args

        gc.collect()

        if self.split:
            new_results = []
          
            for _K in results:
                self.K = _K
                _raw_K, _K, _raw_idx_sizes = self.merge_split()
                new_results.append(_K)

            results = new_results   
          
        self.inc_chr_idx = []
        for i, j in enumerate(results):
            for _ in j:
                self.inc_chr_idx.append(i)
        
 
        second_group_info = [f"{i}g{j}" for i, res in enumerate(results, 1) 
                                            for j, _ in enumerate(res, 1) ]

        if self.exclude_group_to_second:
            exclude_group_init_idx = len(results) + 1


        self.K = list_flatten(results)
        self.K = self.filter_cluster()
        
        if self.exclude_group_to_second:
            self.K.extend(self.exclude_groups)

        reclustered_count = 1 if self.is_recluster_by_linkage else 0
        iter_count = 0
        while reclustered_count > 0 and iter_count < 5:
            iter_count += 1
            _A = HyperGraph.reweight_adjacency_matrix(raw_A, P_allelic_idx=self.P_allelic_idx, 
                                                      P_weak_idx=self.P_weak_idx, allelic_factor=0)
            reclustered_count = self.recluster_by_linkage(A=_A, min_ratio=1.5, iter_count=iter_count)

        second_group_info = [f"{i}g{j}" for i, res in enumerate(results, 1) 
                                        for j, _ in enumerate(res, 1) ]

        if self.split:
            raw_K, _, raw_contig_sizes = self.merge_split()
            length_contents = list(map(
                lambda x: "{:,}".format(int(sum(raw_contig_sizes.get(i, 0) for i in x))), raw_K
            ))
            second_group_info = [f"{self.inc_chr_idx[i] + 1}g{(i + 1)}" for i in range(len(raw_K))]
            chr_hap_counter = defaultdict(int)
            tmp_labels = []
            for chr_i in self.inc_chr_idx[:len(raw_K)]:
                chr_hap_counter[chr_i] += 1
                tmp_labels.append(f"{chr_i + 1}g{chr_hap_counter[chr_i]}")
            second_group_info = tmp_labels
        else:
            vertices_idx_sizes = self.vertices_idx_sizes
            length_contents = list(map(
                lambda x: "{:,}".format(int(sum(vertices_idx_sizes.get(i, 0) for i in x))), self.K))
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

        if hasattr(table, "add_section"):
                table.add_section()

        table.add_row("Total", f"{total_length:,}")
        second_length_contents = pformat(second_length_contents)
        logger.info(f"Second hyperpartition resulted in {len(self.K)} groups:")
                    # f"{second_length_contents}")
        with console_html.capture() as capture:
            console_html.print(table)

        logger.info(capture.get())
    

    @staticmethod
    def _merge(A, K, vertices_idx_sizes, 
                k=None, 
                allelic_idx_set=set(),
                alleletable_df=None,
                min_allelic_overlap=0.1, 
                allele_use_method="conflict", 
                method='mean'):
        if not k:
            return K 
        
        if hasattr(vertices_idx_sizes, 'to_dict'):
            size_db = vertices_idx_sizes.to_dict()
            if 'length' in size_db:
                size_db = size_db['length']
        elif isinstance(vertices_idx_sizes, pd.DataFrame):
            size_db = vertices_idx_sizes['length'].to_dict()
        else:
            size_db = dict(vertices_idx_sizes)

        if alleletable_df is not None:
            mz_df = alleletable_df[[1, 2, 'mzShared']]
            mz_df.set_index([1, 2], inplace=True)
            shared_mz_db = mz_df.to_dict()['mzShared']
            mz_df = alleletable_df[[1, 'mz1']]
            mz_df.drop_duplicates(subset=[1], inplace=True)
            mz_df.set_index([1], inplace=True)
            contig_mz_df = mz_df
            
        else:
            shared_mz_db = None
            contig_mz_df = None

        def group_length(g_tuple):
            if shared_mz_db is not None:
                return float(mz_df.reindex(g_tuple).sum().values[0])
            else:
                # return float(vertices_idx_sizes.reindex(g_tuple).sum().values[0])
                return float(sum(size_db.get(c, 0.0) for c in g_tuple))

        iter_round = 0 
        while k and len(K) > k:
            current_group_number = len(K)
            if iter_round > (current_group_number - k + 100):
                break

            K_tuples = [tuple(sorted(g)) for g in K]
            
            if shared_mz_db is not None:
                group_lengths = np.array([contig_mz_df.reindex(g).sum().values[0] for g in K_tuples])
            else:
                # group_lengths = np.array([vertices_idx_sizes.reindex(g).sum().values[0] for g in K_tuples])
                group_lengths = np.array([sum(size_db.get(c, 0.0) for c in g) for g in K_tuples])
            
            contig_to_group_map = {contig: i for i, group in enumerate(K_tuples) for contig in group}

            value_matrix = np.zeros((current_group_number, current_group_number))
            flag_matrix = np.ones((current_group_number, current_group_number))
            overlap_matrix = np.zeros((current_group_number, current_group_number))
            if allelic_idx_set:
                allelic_len_g1 = defaultdict(float)
                allelic_len_g2 = defaultdict(float)
                allelic_contigs_g1 = defaultdict(set) 
                allelic_contigs_g2 = defaultdict(set)
                contig_partner_groups = defaultdict(set) 
                contig_partner_counts = defaultdict(int) 

                for c1, c2 in allelic_idx_set:
                    g1_idx = contig_to_group_map.get(c1)
                    g2_idx = contig_to_group_map.get(c2)

                    if g1_idx is not None and g2_idx is not None and g1_idx != g2_idx:
                        if g1_idx > g2_idx:
                            g1_idx, g2_idx = g2_idx, g1_idx
                            c1, c2 = c2, c1 

                        if shared_mz_db is not None:
                            allelic_len_g1[(g1_idx, g2_idx)] += shared_mz_db.get((c1, c2), 0.0)
                            allelic_len_g2[(g1_idx, g2_idx)] += shared_mz_db.get((c2, c1), 0.0)
                        else:
                            allelic_contigs_g1[(g1_idx, g2_idx)].add(c1)
                            allelic_contigs_g2[(g1_idx, g2_idx)].add(c2)
                            # allelic_len_g1[(g1_idx, g2_idx)] += vertices_idx_sizes.loc[c1].values[0]
                            # allelic_len_g2[(g1_idx, g2_idx)] += vertices_idx_sizes.loc[c2].values[0]
                            contig_partner_groups[c1].add((g1_idx, g2_idx))
                            contig_partner_groups[c2].add((g1_idx, g2_idx))
                            contig_partner_counts[(c1, g1_idx, g2_idx)] += 1
                            contig_partner_counts[(c2, g1_idx, g2_idx)] += 1

                for (i, j), s1 in allelic_contigs_g1.items():
                    s2 = allelic_contigs_g2.get((i, j), set())
                    # len1 = vertices_idx_sizes.loc[list(s1)]['length'].sum() if len(s1) else 0.0
                    # len2 = vertices_idx_sizes.loc[list(s2)]['length'].sum() if len(s2) else 0.0
                    len1 = sum(size_db.get(c, 0.0) for c in s1) if len(s1) else 0.0
                    len2 = sum(size_db.get(c, 0.0) for c in s2) if len(s2) else 0.0
                    allelic_len_g1[(i, j)] = float(len1)
                    allelic_len_g2[(i, j)] = float(len2)

                for (i, j), g1_allelic_len in allelic_len_g1.items():
                    g2_allelic_len = allelic_len_g2.get((i, j), 0.0)
                    
                    overlap1 = g1_allelic_len / group_lengths[i] if group_lengths[i] > 0 else 0
                    overlap2 = g2_allelic_len / group_lengths[j] if group_lengths[j] > 0 else 0

                    if overlap1 > min_allelic_overlap or overlap2 > min_allelic_overlap:
                        flag_matrix[i, j] = 0
                        overlap_matrix[i, j] = max(overlap1, overlap2)
            
            rows = []
            cols = []
            for gid, group in enumerate(K_tuples):
                rows.extend(group)
                cols.extend([gid] * len(group))
            
            M = csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), 
                           shape=(A.shape[0], current_group_number))
            value_matrix = (M.T.dot(A).dot(M)).toarray()
            np.fill_diagonal(value_matrix, 0)
            if method == "mean":
                group_sizes = np.array([len(g) for g in K_tuples], dtype=np.float32)
                size_outer = np.outer(group_sizes, group_sizes)
                size_outer[size_outer == 0] = 1.0
                value_matrix = value_matrix / size_outer

            # for i in range(current_group_number):
            #     group1 = K_tuples[i]
            #     for j in range(i + 1, current_group_number):
            #         group2 = K_tuples[j]
                    
            #         sub_matrix = A[group1, :][:, group2]
            #         if method == "mean":
            #             value = sub_matrix.mean()
            #         elif method == "median":
            #             if sub_matrix.nnz > 0:
            #                 value = np.median(sub_matrix.toarray())
            #             else:
            #                 value = 0.0
            #         else:
            #             value = sub_matrix.sum()
            #         value_matrix[i, j] = value
   
            # value_matrix += value_matrix.T
            total_value = max(np.triu(value_matrix, 1).sum(), 1e-9)
            res = {}

            for i in range(current_group_number):
                i_value = value_matrix[i].sum()
                for j in range(i+1, current_group_number):
                    j_value = value_matrix[j].sum()
                    value =  value_matrix[i, j]
                    flag = flag_matrix[i, j]
                    raw_Q = value - (j_value * i_value) / total_value
                    # if method == "sum":
                      
                    #     l1_mb = max(group_lengths[i] / 1e6, 0.001)
                    #     l2_mb = max(group_lengths[j] / 1e6, 0.001)
                    #     raw_Q = raw_Q / (l1_mb * l2_mb)

                    if flag:
                        Q = raw_Q
                    else:
                        if allele_use_method == "conflict":
                            Q = - 2 ** 64
                        else:
                            Q = overlap_matrix[i, j] * raw_Q
                            Q = max(Q, 0)
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
    
    @staticmethod
    def _merge_incomplete(
        A,
        K,
        vertices_idx_sizes,
        k=None,
        allelic_idx_set=set(),
        alleletable_df=None,
        min_allelic_overlap=0.1,
        allele_use_method="conflict",
        method="mean",
        max_extra_iter=100,
    ):
        import heapq
        if not k or len(K) <= k:
            return K

        if alleletable_df is not None:
            mz_df = alleletable_df[[1, 'mz1']].drop_duplicates(subset=[1]).set_index(1)
            shared_mz_db = alleletable_df[[1, 2, 'mzShared']].set_index([1, 2])['mzShared'].to_dict()
            use_mz = True
        else:
            shared_mz_db = {}
            mz_df = None
            use_mz = False

        def group_length(g_tuple):
            if use_mz:
                return float(mz_df.reindex(g_tuple).sum().values[0])
            else:
                return float(vertices_idx_sizes.reindex(g_tuple).sum().values[0])

        groups = {i: tuple(sorted(g)) for i, g in enumerate(K)}
        lengths = {gid: group_length(g) for gid, g in groups.items()}

        active = set(groups.keys())
        version = {gid: 0 for gid in groups}

        overlap_info = {}
        conflict_flag = set()

        if allelic_idx_set:
            contig_to_group = {}
            for gid, g in groups.items():
                for c in g:
                    contig_to_group[c] = gid

            allelic_len_g1 = defaultdict(float)
            allelic_len_g2 = defaultdict(float)
            allelic_contigs_g1 = defaultdict(set) 
            allelic_contigs_g2 = defaultdict(set)
            contig_partner_groups = defaultdict(set) 
            contig_partner_counts = defaultdict(int) 

            for c1, c2 in allelic_idx_set:
                g1_idx = contig_to_group.get(c1)
                g2_idx = contig_to_group.get(c2)
                if g1_idx is None or g2_idx is None or g1_idx == g2_idx:
                    continue
                if g1_idx > g2_idx:
                    g1_idx, g2_idx = g2_idx, g1_idx
                    c1, c2 = c2, c1

                if use_mz:
                    allelic_len_g1[(g1_idx, g2_idx)] += shared_mz_db.get((c1, c2), 0.0)
                    allelic_len_g2[(g1_idx, g2_idx)] += shared_mz_db.get((c2, c1), 0.0)
                else:
                    allelic_contigs_g1[(g1_idx, g2_idx)].add(c1)
                    allelic_contigs_g2[(g1_idx, g2_idx)].add(c2)
                    contig_partner_groups[c1].add((g1_idx, g2_idx))
                    contig_partner_groups[c2].add((g1_idx, g2_idx))
                    contig_partner_counts[(c1, g1_idx, g2_idx)] += 1
                    contig_partner_counts[(c2, g1_idx, g2_idx)] += 1

            for (i, j), s1 in allelic_contigs_g1.items():
                s2 = allelic_contigs_g2.get((i, j), set())
                len1 = vertices_idx_sizes.loc[list(s1)]['length'].sum() if len(s1) else 0.0
                len2 = vertices_idx_sizes.loc[list(s2)]['length'].sum() if len(s2) else 0.0
                allelic_len_g1[(i, j)] = float(len1)
                allelic_len_g2[(i, j)] = float(len2)


            for (i, j), len1 in allelic_len_g1.items():
                len2 = allelic_len_g2.get((i, j), 0.0)
                L1 = lengths[i]
                L2 = lengths[j]
                ov1 = len1 / L1 if L1 > 0 else 0.0
                ov2 = len2 / L2 if L2 > 0 else 0.0
                if ov1 > min_allelic_overlap or ov2 > min_allelic_overlap:
                    conflict_flag.add((i, j))
                    overlap_info[(i, j)] = max(ov1, ov2)

        def pair_value(gid1: int, gid2: int) -> float:
            g1 = groups[gid1]
            g2 = groups[gid2]
            sub = A[g1, :][:, g2]
            if method == "mean":
                return float(sub.mean())
            elif method == "median":
                if sub.nnz > 0:
                    return float(np.median(sub.toarray()))
                return 0.0
            else:  # 'sum'
                return float(sub.sum())

        value_cache = {}
        order_ids = sorted(active)
        for i in range(len(order_ids)):
            gi = order_ids[i]
            for j in range(i+1, len(order_ids)):
                gj = order_ids[j]
                v = pair_value(gi, gj)
                value_cache[(gi, gj)] = v

        def total_value():
            return sum(value_cache[k] for k in value_cache)

        total_v = total_value()

        def compute_Q(i: int, j: int) :
            if i > j:
                i, j = j, i
            v = value_cache.get((i, j), 0.0)
            if v <= 0 or total_v <= 0:
                return 0.0

            i_sum = sum(value_cache[(min(i,x), max(i,x))] for x in active if x != i)
            j_sum = sum(value_cache[(min(j,x), max(j,x))] for x in active if x != j)
            raw = v - (j_sum * i_sum) / total_v
            if (i, j) in conflict_flag:
                if allele_use_method == "conflict":
                    return -2**64
                else:
                    ov = overlap_info.get((i, j), 0.0)
                    adj = ov * raw
                    return max(adj, 0.0)
            return raw

        heap = []
        for (i, j), _v in value_cache.items():
            q = compute_Q(i, j)
            heapq.heappush(heap, (-q, version[i], version[j], i, j))

        merges = 0
        max_merges = (len(K) - k) + max_extra_iter

        while len(active) > k and heap and merges < max_merges:
            neg_q, vi, vj, i, j = heapq.heappop(heap)

            if i not in active or j not in active:
                continue
            if vi != version[i] or vj != version[j]:
 
                q_new = compute_Q(i, j)
                heapq.heappush(heap, (-q_new, version[i], version[j], i, j))
                continue

            q = -neg_q
            if q == -2**64:
                continue


            new_group = tuple(sorted(set(groups[i]) | set(groups[j])))
            groups[i] = new_group
            lengths[i] = group_length(new_group)
            version[i] += 1

            active.remove(j)
            version[j] += 1
            del groups[j]

            new_cache = {}
            for (a, b), val in value_cache.items():
                if a == j or b == j:
                    continue
                new_cache[(a, b)] = val
            value_cache = new_cache

            for other in active:
                if other == i:
                    continue
                a, b = (other, i) if other < i else (i, other)
                value_cache[(a, b)] = pair_value(a, b)

            total_v = total_value()

            for other in active:
                if other == i:
                    continue
                a, b = (other, i) if other < i else (i, other)
                q_new = compute_Q(a, b)
                heapq.heappush(heap, (-q_new, version[a], version[b], a, b))

            merges += 1

        final_groups = [list(groups[g]) for g in active]
        final_groups = sorted(final_groups,
                              key=lambda x: vertices_idx_sizes.reindex(x).sum().values[0],
                              reverse=True)
        return final_groups
    
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
        
        try:
            allelic_idx_set = self.allelic_idx_set 
        except AttributeError:
            if self.alleletable is not None:
                at = AlleleTable(self.alleletable, fmt="allele2", sort=False)
                at_df = at.data
                at_df = at_df[at_df['similarity'] >= 0.85]
                at_df[1] = at_df[1].map(vertices_idx.get)
                at_df[2] = at_df[2].map(vertices_idx.get)
                at_df.dropna(inplace=True)
                at_df[1] = at_df[1].astype(int)
                at_df[2] = at_df[2].astype(int)
                allelic_idx_set = set(map(tuple, at_df[[1, 2]].values))
                
            else:
                allelic_idx_set = set()
                at_df = None

        self.K = HyperPartition._merge(A, self.K, vertices_idx_sizes, self.k[0],
                                        allelic_idx_set=allelic_idx_set,
                                        alleletable_df=at_df,
                                        min_allelic_overlap=0.8, 
                                        allele_use_method="merge", method='sum')

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
        A = HyperGraph.clique_expansion_init(sub_H, P_allelic_idx=P_allelic_idx, allelic_factor=0)
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

    def recluster_by_linkage(self, A=None, min_ratio=1.2, iter_count=1):
        if A is None:
            A = HyperGraph.clique_expansion_init(self.H, min_weight=self.min_weight, 
                                                use_high_order=False)
        
        if not isinstance(A, csr_matrix):
            A = A.tocsr()

        num_vertices = A.shape[0]
        num_groups = len(self.K)

        contig_to_group = np.full(num_vertices, -1, dtype=int)
        rows, cols = [], []
        for gid, group in enumerate(self.K):
            g_list = list(group)
            contig_to_group[g_list] = gid
            contig_to_group[g_list] = gid
            rows.extend(g_list)
            cols.extend([gid] * len(g_list))
        
        allelic_lookup = defaultdict(set)
        if hasattr(self, 'allelic_idx_set') and self.allelic_idx_set:
            for c1, c2 in self.allelic_idx_set:
                allelic_lookup[c1].add(c2)
                allelic_lookup[c2].add(c1)
        idx_to_vertices = self.idx_to_vertices
        contig_lens = np.array([self.contig_sizes.get(idx_to_vertices[i], 1) for i in range(num_vertices)], dtype=np.float32)

        M = csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), 
                       shape=(num_vertices, num_groups))

        S = A.dot(M).toarray() 
        length_norm = False
        if length_norm:
            group_total_lengths = np.bincount(cols, weights=contig_lens[rows], minlength=num_groups).astype(np.float32)
            
            dia = A.diagonal()
            valid_mask = contig_to_group != -1
            row_indices = np.where(valid_mask)[0]
            col_indices = contig_to_group[valid_mask]
            S[row_indices, col_indices] -= dia[valid_mask]

            denom_matrix = group_total_lengths + 10000.0 
            norm_matrix = S / denom_matrix
            current_denoms = group_total_lengths[col_indices] - contig_lens[valid_mask] + 10000.0
            current_denoms = np.maximum(current_denoms, 1.0)
            current_norm_scores = S[row_indices, col_indices] / current_denoms

            best_gids = np.argmax(norm_matrix, axis=1)
            best_norm_scores = norm_matrix[np.arange(num_vertices), best_gids]

            total_signals = np.sum(S, axis=1)

            full_current_norm = np.zeros(num_vertices)
            full_current_norm[valid_mask] = current_norm_scores
            specificity = np.divide(S[np.arange(num_vertices), best_gids], total_signals, 
                                    out=np.zeros(num_vertices), where=total_signals>0)
        else:

            best_gids = np.argmax(S, axis=1)
            best_S = S[np.arange(num_vertices), best_gids]
            
            curr_S = np.zeros(num_vertices)
            valid_mask = contig_to_group != -1
            row_indices = np.where(valid_mask)[0]
            col_indices = contig_to_group[valid_mask]
            curr_S[valid_mask] = S[row_indices, col_indices]
            
            total_signals = np.sum(S, axis=1)
            specificity = np.divide(best_S, total_signals, out=np.zeros(num_vertices), where=total_signals>0)
            full_current_norm = curr_S
            best_norm_scores = best_S

        potential_mask = (
            (best_gids != contig_to_group) & 
            (best_norm_scores > (full_current_norm * min_ratio)) & 
            (specificity > 0.9) &
            (total_signals > self.min_weight)
        )
        candidates = np.where(potential_mask)[0]
        if len(candidates) == 0:
            return 0

        candidates = sorted(candidates, key=lambda x: contig_lens[x], reverse=True)

        reclustered_count = 0
        new_K = [set(group) for group in self.K]

        for c_idx in candidates:
            best_gid = best_gids[c_idx]
            current_gid = contig_to_group[c_idx]
            if S[c_idx, best_gid] / (total_signals[c_idx] + 1e-9) < 0.25:
                continue
    
            conflicts = allelic_lookup.get(c_idx, set())
            if conflicts.isdisjoint(new_K[best_gid]):
                if current_gid != -1:
                    new_K[current_gid].discard(c_idx)
                new_K[best_gid].add(c_idx)

                contig_to_group[c_idx] = best_gid
                reclustered_count += 1

            
        if reclustered_count > 0:
            logger.debug(f"Round [{iter_count}]: reclustered {reclustered_count} contigs based on linkage strength.")
            self.K = [sorted(list(g)) for g in new_K]
            
        return reclustered_count
 

    # @staticmethod
    # def refine_allelic_errors(K, k, A, idx_size, prune_pair_df, 
    #                           min_weight, allelic_similarity=0.99, 
    #                           rel_tol=0.2):
    #     if prune_pair_df is None:
    #         return K 
    #     if len(K) == 1:
    #         return K
    
        
    #     tmp_df = prune_pair_df[(prune_pair_df['type'] == 0) &
    #                             (prune_pair_df['similarity'] >= allelic_similarity)]
        
    #     idx_size_db = idx_size.to_dict()['length']

    #     tmp_df['length1'] = tmp_df['contig1'].map(idx_size_db)
    #     tmp_df['length2'] = tmp_df['contig2'].map(idx_size_db)
    
    #     high_similar_pairs = set()
    #     for idx, row in tmp_df.iterrows():
    #         contig1, contig2 = row['contig1'], row['contig2']
    #         mz1, mz2 = row['mz1'], row['mz2']
    #         length1, length2 = row['length1'], row['length2']
    #         if math.isclose(mz1, mz2, rel_tol=rel_tol) or \
    #             math.isclose(length1, length2, rel_tol=rel_tol):
    #             if contig1 > contig2:
    #                 contig1, contig2 = contig2, contig1
    #             high_similar_pairs.add((contig1, contig2))


    #     if not high_similar_pairs:
    #         return K
        
    #     retain_res = OrderedDict()
    #     candicate_error_contig_pairs = OrderedDict()
    #     all_candicate_error_contigs = set()
    #     all_candicate_error_contig_db = defaultdict(set)
    #     for idx, group in enumerate(K):
    #         group_set = set(combinations(sorted(group), 2))
    #         candicate_allelic_error_pairs = group_set.intersection(high_similar_pairs)
    #         candicate_error_contigs = set(list_flatten(candicate_allelic_error_pairs))
    #         all_candicate_error_contigs.update(candicate_error_contigs)
    #         for pair in candicate_allelic_error_pairs:
    #             all_candicate_error_contig_db[pair[0]].add(pair[0])
    #             all_candicate_error_contig_db[pair[0]].add(pair[1]
    #                                                        )

    #         retain_contigs = set(group).difference(candicate_error_contigs)
    #         retain_res[idx] = retain_contigs
    
    #         for pair in candicate_allelic_error_pairs:
    #             if pair not in candicate_error_contig_pairs:
    #                 candicate_error_contig_pairs[pair] = idx
        
    #     for pair in candicate_error_contig_pairs:
    #         idx = candicate_error_contig_pairs[pair]
    #         retain_contigs = retain_res[idx]
    #         other_idxes = sorted(set(list(retain_res)).difference([idx]))

    #         counts_list = []
    #         for contig in pair:
                
    #             interaction_pairs = list(product([int(contig)], retain_contigs))
                
    #             interaction_pairs = [(c1, c2) if c1 < c2 else (c2, c1) for c1, c2 in interaction_pairs] 
    #             contig1_idxes, contig2_idxes = list(zip(*interaction_pairs))
                
    #             sub_A = A[contig1_idxes, :][:, contig2_idxes]

    #             counts_list.append(sub_A.sum())
            
    #         max_idx = pair[np.argmax(counts_list)]
    #         min_idx = pair[np.argmin(counts_list)]

    #         retain_res[idx].add(max_idx)

    #         other_counts_list = []
    #         for other_idx in other_idxes:
    #             other_retain_contigs = retain_res[other_idx]
    #             if len(other_retain_contigs) == 0:
    #                 continue
    #             interaction_pairs = list(product([min_idx], other_retain_contigs))
    #             interaction_pairs = [(c1, c2) if c1 < c2 else (c2, c1) for c1, c2 in interaction_pairs] 
    #             contig1_idxes, contig2_idxes = list(zip(*interaction_pairs))
    #             if len(contig1_idxes) == 0 or len(contig2_idxes) == 0:
    #                 other_counts_list.append(0)
    #             else:
    #                 sub_A = A[contig1_idxes,  contig2_idxes]
    #                 other_counts_list.append(sub_A.sum())
            
    #         other_counts_db = dict(zip(other_idxes, other_counts_list))

    #         other_counts_list = sorted(other_counts_db, 
    #                                    key=lambda x: other_counts_db[x],
    #                                    reverse=True)

    #         for i in range(len(other_counts_list)):
    #             max_chrom = other_counts_list[i]
    #             try:
    #                 second_chrom = other_counts_list[i+1]
    #             except:
    #                 second_chrom = None
                
    #             tmp_error_contigs = all_candicate_error_contig_db[min_idx]
    #             tmp_counts = len(retain_res[max_chrom].intersection(tmp_error_contigs))
    #             if tmp_counts > 0:
    #                 continue 
    #             if second_chrom:
    #                 if (other_counts_db[max_chrom] > min_weight) and (other_counts_db[max_chrom] != other_counts_db[second_chrom]):
    #                     retain_res[max_chrom].add(min_idx)
    #                     break
    #             else:
    #                 if other_counts_db[max_chrom] > min_weight:
    #                     retain_res[max_chrom].add(min_idx)
    #                     break

    #     new_K = []
    #     for idx in sorted(retain_res):
    #         new_K.append(sorted(retain_res[idx]))
        

    #     return new_K 
         
    @staticmethod
    def refine_allelic_errors(K, k, A, idx_size, prune_pair_df, 
                                min_weight, allelic_similarity=0.85, 
                                rel_tol=0.2, min_ratio=1.2):
        """
        Improvements:
        1. Competitive Scoring: Maximize global fitness (Score(u, cur) + Score(v, other)) vs (Score(v, cur) + Score(u, other)).
        2. Strict Constraint Check: Ensure moving a contig doesn't create new conflicts in the target group.
        """
        if prune_pair_df is None:
            return K 
        if len(K) < 2:
            return K
    
        tmp_df = prune_pair_df[(prune_pair_df['type'] == 0) &
                                (prune_pair_df['similarity'] >= allelic_similarity)]
        
        idx_size_db = idx_size.to_dict()['length']
        tmp_df['length1'] = tmp_df['contig1'].map(idx_size_db)
        tmp_df['length2'] = tmp_df['contig2'].map(idx_size_db)
    
        high_similar_pairs = set()
        high_similar_lookup = defaultdict(set)
        
        for idx, row in tmp_df.iterrows():
            contig1, contig2 = int(row['contig1']), int(row['contig2'])
            mz1, mz2 = row['mz1'], row['mz2']
            length1, length2 = row['length1'], row['length2']
            
            if math.isclose(mz1, mz2, rel_tol=rel_tol) or \
                math.isclose(length1, length2, rel_tol=rel_tol):
                if contig1 > contig2:
                    contig1, contig2 = contig2, contig1
                high_similar_pairs.add((contig1, contig2))
                high_similar_lookup[contig1].add(contig2)
                high_similar_lookup[contig2].add(contig1)

        if not high_similar_pairs:
            return K
        
        cluster_sets = {i: set(group) for i, group in enumerate(K)}
        contig_to_group = {c: i for i, g in cluster_sets.items() for c in g}
        group_lengths = {i: sum(idx_size_db.get(c, 0) for c in group) for i, group in enumerate(K)}

        def get_density(contig, group_idx, ignore_contigs=None):
            gset = cluster_sets[group_idx]
            targets = list(gset - (ignore_contigs or set()))
            if not targets:
                return 0.0
            
            strength = A[contig, targets].sum()
            denom = group_lengths[group_idx] + 100000.0
            if ignore_contigs:
                denom -= sum(idx_size_db.get(c, 0) for c in ignore_contigs)
            
            return strength / max(denom, 1.0)

        # corrected_count = 0
        # for c1, c2 in sorted(set(high_similar_pairs), key=lambda x: idx_size_db.get(x[0], 0), reverse=True):
        #     cg1, cg2 = contig_to_group.get(c1), contig_to_group.get(c2)
        #     if cg1 is None or cg2 is None or cg1 != cg2: 
        #         continue
            
        #     origin_gid = cg1
        #     ignore = {c1, c2}
            
        #     best_dest_c1, best_dens_c1 = -1, 0.0
        #     best_dest_c2, best_dens_c2 = -1, 0.0
            
        #     for ogid in cluster_sets:
        #         if ogid == origin_gid: 
        #             continue
        #         if not high_similar_lookup[c1].intersection(cluster_sets[ogid]):
        #             d1 = get_density(c1, ogid)
        #             if d1 > best_dens_c1: 
        #                 best_dest_c1, best_dens_c1 = ogid, d1
        #         if not high_similar_lookup[c2].intersection(cluster_sets[ogid]):
        #             d2 = get_density(c2, ogid)
        #             if d2 > best_dens_c2: 
        #                 best_dest_c2, best_dens_c2 = ogid, d2

        #     curr_dens_c1 = get_density(c1, origin_gid, ignore)
        #     curr_dens_c2 = get_density(c2, origin_gid, ignore)

        #     move_c1, move_c2 = False, False
        #     if best_dest_c1 != -1 and (best_dens_c1 > curr_dens_c1 * min_ratio):
        #         if best_dest_c2 != -1 and (best_dens_c2 > curr_dens_c2 * min_ratio):
        #             if (best_dens_c1 + curr_dens_c2) > (best_dens_c2 + curr_dens_c1): 
        #                 move_c1 = True
        #             else: 
        #                 move_c2 = True
        #         else: 
        #             move_c1 = True
        #     elif best_dest_c2 != -1 and (best_dens_c2 > curr_dens_c2 * min_ratio):
        #         move_c2 = True

        #     target_c = c1 if move_c1 else (c2 if move_c2 else None)
        #     target_g = best_dest_c1 if move_c1 else (best_dest_c2 if move_c2 else None)

        #     if target_c and target_g != -1:
        #         cluster_sets[origin_gid].remove(target_c)
        #         cluster_sets[target_g].add(target_c)
        #         group_lengths[origin_gid] -= idx_size_db[target_c]
        #         group_lengths[target_g] += idx_size_db[target_c]
        #         contig_to_group[target_c] = target_g
        #         corrected_count += 1

        # if corrected_count > 0:
        #     logger.info(f"Refined {corrected_count} allelic errors using normalized linkage density.")
        # return [sorted(list(g)) for g in cluster_sets.values()]

        conflicting_pairs_queue = []
        for c1, c2 in high_similar_pairs:
            g1 = contig_to_group.get(c1)
            g2 = contig_to_group.get(c2)
            if g1 is not None and g2 is not None and g1 == g2:
                combined_len = idx_size_db.get(c1, 0) + idx_size_db.get(c2, 0)
                conflicting_pairs_queue.append((c1, c2, g1, combined_len))
        
        if not conflicting_pairs_queue:
            return K
        
        conflicting_pairs_queue.sort(key=lambda x: x[3], reverse=True)

        def get_connection_strength(contig, target_group_set, ignore_set=None):
            if not target_group_set: 
                return 0.0
            
            if ignore_set:
                targets = list(target_group_set - ignore_set)
            else:
                targets = list(target_group_set)
            
            if not targets: 
                return 0.0

            return A[contig, targets].sum()

        corrected_count = 0
        
        for c1, c2, origin_g, _ in conflicting_pairs_queue:
            curr_g1 = -1
            curr_g2 = -2
            
            for gid, gset in cluster_sets.items():
                if c1 in gset: 
                    curr_g1 = gid
                if c2 in gset: 
                    curr_g2 = gid

            if curr_g1 != curr_g2:
                continue
            
            group_idx = curr_g1
            current_group_set = cluster_sets[group_idx]
            
            other_groups = [g for g in cluster_sets.keys() if g != group_idx]
            
            ignore_pair = {c1, c2}
            score_c1_curr = get_connection_strength(c1, current_group_set, ignore_pair)
            score_c2_curr = get_connection_strength(c2, current_group_set, ignore_pair)

            best_other_g_c1 = None
            best_score_c1_other = 0.0
            
            for og in other_groups:
                target_set = cluster_sets[og]
              
                conflicts = high_similar_lookup[c1]
                if not conflicts.intersection(target_set):
                    s = get_connection_strength(c1, target_set)
                    if s > best_score_c1_other:
                        best_score_c1_other = s
                        best_other_g_c1 = og

            best_other_g_c2 = None
            best_score_c2_other = 0.0
            
            for og in other_groups:
                target_set = cluster_sets[og]
                conflicts = high_similar_lookup[c2]
                if not conflicts.intersection(target_set):
                    s = get_connection_strength(c2, target_set)
                    if s > best_score_c2_other:
                        best_score_c2_other = s
                        best_other_g_c2 = og

            if best_score_c1_other < min_weight: 
                best_other_g_c1 = None
            if best_score_c2_other < min_weight: 
                best_other_g_c2 = None

            score_keep_c1 = score_c1_curr + (best_score_c2_other if best_other_g_c2 is not None else -9999)
            
            score_keep_c2 = score_c2_curr + (best_score_c1_other if best_other_g_c1 is not None else -9999)
            
            if best_other_g_c1 is None and best_other_g_c2 is None:
                len1 = idx_size_db.get(c1, 0)
                len2 = idx_size_db.get(c2, 0)

                if (score_c1_curr < score_c2_curr) or (score_c1_curr == score_c2_curr and len1 <= len2):
                    move_out = c1
                else:
                    move_out = c2

                new_gid = max(cluster_sets.keys()) + 1 if cluster_sets else 0
                cluster_sets[group_idx].remove(move_out)
                cluster_sets[new_gid] = {move_out}
                corrected_count += 1
                continue
            
            move_c1 = False
            move_c2 = False
            
            if best_other_g_c1 is not None and best_other_g_c2 is None:
                move_c1 = True
            elif best_other_g_c1 is None and best_other_g_c2 is not None:
                move_c2 = True
            else:
                if score_keep_c1 >= score_keep_c2:
                    move_c2 = True
                else:
                    move_c1 = True
                    
            if move_c1:
                cluster_sets[group_idx].remove(c1)
                cluster_sets[best_other_g_c1].add(c1)
                corrected_count += 1
            elif move_c2:
                cluster_sets[group_idx].remove(c2)
                cluster_sets[best_other_g_c2].add(c2)
                corrected_count += 1

        if corrected_count > 0:
            logger.info(f"Refined {corrected_count} allelic mis-assignments using competitive scoring.")

        new_K = [sorted(list(cluster_sets[i])) for i in sorted(cluster_sets.keys())]

        return new_K


    def rescue(self, rescue_threshold=0.1, allelic_similarity=0.99):
        """
        rescue unanchor contigs 
        """
        if not self.K:
            return

        logger.info("Start rescuing unanchored contigs...")

        anchored_indices = set(list_flatten(self.K))
        all_indices = set(self.vertices_idx.values())
        unanchored_indices = list(all_indices - anchored_indices)
        
        if not unanchored_indices:
            logger.info("No unanchored contigs found to rescue.")
            return

        allelic_constraints = defaultdict(set)
        if self.prune_pair_df is not None:
                high_sim_df = self.prune_pair_df[
                    (self.prune_pair_df['type'] == 0) & 
                    (self.prune_pair_df['similarity'] >= allelic_similarity)
                ]
                vertices_idx_map = self.vertices_idx

                for _, row in high_sim_df.iterrows():
                    c1_name = row['contig1']
                    c2_name = row['contig2']
                    c1_idx = vertices_idx_map.get(c1_name) if isinstance(c1_name, str) else c1_name
                    c2_idx = vertices_idx_map.get(c2_name) if isinstance(c2_name, str) else c2_name
                    
                    if c1_idx is not None and c2_idx is not None:
                        allelic_constraints[c1_idx].add(c2_idx)
                        allelic_constraints[c2_idx].add(c1_idx)


        logger.debug("Generating adjacency matrix for rescue...")
        A = HyperGraph.clique_expansion_init(self.H, min_weight=self.min_weight)

        rescued_count = 0
        
        vertices_idx_sizes = self.vertices_idx_sizes.to_dict()['length']
        unanchored_indices.sort(key=lambda x: vertices_idx_sizes.get(x, 0), reverse=True)

        for contig_idx in unanchored_indices:
            best_group_idx = -1
            best_score = 0.0
            total_connection = 0.0
            
            group_scores = []
            
            for g_idx, group in enumerate(self.K):
                score = A[contig_idx, group].sum()
                group_scores.append(score)
                total_connection += score

            if total_connection == 0:
                continue

            # Find candidate group
            max_score = max(group_scores)

            ratio = max_score / total_connection
            
            if ratio < rescue_threshold:
                continue

            # Allelic Conflict Check
            target_group = self.K[candidate_g_idx]
            is_conflict = False
            
            if contig_idx in allelic_constraints:
                conflict_partners = allelic_constraints[contig_idx]
                intersection = conflict_partners.intersection(set(target_group))
                if intersection:

                    is_conflict = True
            
            if not is_conflict:
                self.K[candidate_g_idx].append(contig_idx)
                rescued_count += 1
        
        if rescued_count > 0:
            logger.info(f"Rescued {rescued_count} contigs into existing groups.")
            
            # Re-sort groups if needed
            self.K = sorted(self.K, 
                            key=lambda x: self.vertices_idx_sizes.reindex(x)['length'].sum(), 
                            reverse=True)

    @staticmethod
    def _rescue_static(K, H, vertices_idx_sizes, prune_pair_df=None, 
                       min_weight=1.0, rescue_threshold=0.1, 
                       min_contacts=2, ambiguity_threshold=0.85,
                       allelic_similarity=0.99, verbose=False):
        """
        Static method to rescue unanchored contigs.
        Can be used inside _incremental_partition or other static contexts.

        Params:
        -------
        K: list of lists
            Current clustering results [[idx1, idx2...], [idx3...]]
        H: scipy.sparse.csr_matrix or HyperGraph
            Hypergraph or Incidence Matrix.
        vertices_idx_sizes: pd.DataFrame/Series/dict
            Length mapping for vertices.
        prune_pair_df: pd.DataFrame
            Table containing allelic/conflict information.
        min_weight: float
            Minimum weight for clique expansion.
        rescue_threshold: float
            Ratio threshold (score_to_best / total_score).
        min_contacts: int
            Absolute score threshold.
        ambiguity_threshold: float
            Ratio of second_best / best. If > threshold, consider ambiguous and skip.
        
        Returns:
        --------
        new_K: list of lists
            Updated clustering results.
        rescued_count: int
            Number of rescued contigs.
        """
        import numpy as np
        from collections import defaultdict
        
        if not K:
            return K, 0

        total_vertices = H.shape[0]
        anchored_indices = set(list_flatten(K))
        all_indices = set(range(total_vertices))
        unanchored_indices = list(all_indices - anchored_indices)
        
        if not unanchored_indices:
            return K, 0


        if hasattr(H, 'shape'):

             A = HyperGraph.clique_expansion_init(H, min_weight=min_weight)
        else:

             A = H

        allelic_constraints = defaultdict(set)
        if prune_pair_df is not None:

             try:
                 high_sim_rows = prune_pair_df[
                     (prune_pair_df['type'] == 0) & 
                     (prune_pair_df['similarity'] >= allelic_similarity)
                 ][['contig1', 'contig2']].values
                 
                 for c1, c2 in high_sim_rows:
                     allelic_constraints[int(c1)].add(int(c2))
                     allelic_constraints[int(c2)].add(int(c1))
             except KeyError:
                 if verbose: 
                    logger.warning("prune_pair_df format error in rescue.")


        def get_len(idx):
            try:
                if isinstance(vertices_idx_sizes, dict):
                    return vertices_idx_sizes.get(idx, 0)
                res = vertices_idx_sizes.loc[idx]
                if hasattr(res, 'values'): 
                    return res['length'].max() if 'length' in res else res.values[0]
                return res # Series
            except Exception:
                return 0

        unanchored_indices.sort(key=get_len, reverse=True)

        rescued_count = 0
        
        for contig_idx in unanchored_indices:
            group_scores = []
            valid_stats = False
            total_connection = 0.0


            for g_idx, group in enumerate(K):
                if not group: 
                    group_scores.append(0.0)
                    continue

                raw_score = A[contig_idx, group].sum()
                
                group_scores.append(raw_score)
                total_connection += raw_score
            
            if total_connection <= 1e-6:
                continue

            group_scores = np.array(group_scores)
            best_idx = np.argmax(group_scores)
            best_score = group_scores[best_idx]
            

            if best_score < min_contacts:
                continue


            ratio = best_score / total_connection
            if ratio < rescue_threshold:
                continue
            
            sorted_scores = np.sort(group_scores)
            if len(sorted_scores) > 1:
                second_best_score = sorted_scores[-2]
                if second_best_score > 0:
                    ambiguity = second_best_score / best_score
                    if ambiguity > ambiguity_threshold:
                        continue

            target_group = K[best_idx]
            has_conflict = False
            
            if contig_idx in allelic_constraints:
                conflict_partners = allelic_constraints[contig_idx]
                if not conflict_partners.isdisjoint(target_group):
                    has_conflict = True
            
            if not has_conflict:
                K[best_idx].append(contig_idx)

                rescued_count += 1
                if verbose:
                    logger.debug(f"Rescued contig {contig_idx} to group {best_idx} (Score: {best_score:.1f}, Ratio: {ratio:.2f})")

        if verbose and rescued_count > 0:
            logger.info(f"Statically rescued {rescued_count} contigs.")

        return K, rescued_count

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

        def rescue_collapsed(
            self,
            A=None,
            min_signal=50.0,
            ratio_threshold=0.6,
            allelic_similarity=0.99,
            min_allelic_partners=2,
            max_dups_per_contig=2,
            output_collapsed_list="collapsed_list.tsv",
            apply_to_clusters=False,
        ):
            """
            Rescue collapsed contigs:
            - detect contigs that strongly connect to multiple groups (top2/top1 high)
            - and have many allelic partners (suggesting they represent collapsed haplotypes)
            - output a collapsed_list for downstream `pqs dup`
            - optionally (apply_to_clusters=True) duplicate them into multiple groups in self.K (non-strict partition)
    
            Parameters
            ----------
            A : csr_matrix
                adjacency matrix, if None will be computed from self.H.
            min_signal : float
                minimum top1 signal required.
            ratio_threshold : float
                require top2/top1 >= ratio_threshold to call collapsed.
            allelic_similarity : float
                similarity cutoff in prune table for allelic partners.
            min_allelic_partners : int
                minimum number of allelic partners to consider collapsed.
            max_dups_per_contig : int
                maximum number of groups to duplicate into (>=2).
            output_collapsed_list : str
                output file for downstream duplication.
            apply_to_clusters : bool
                if True, will insert contig into multiple groups in self.K (for cluster output only).
            """
            if not self.K or len(self.K) < 2:
                return
    
            if A is None:
                logger.debug("Generating adjacency matrix for collapsed rescue...")
                A = HyperGraph.clique_expansion_init(self.H, min_weight=self.min_weight, use_high_order=False)
                A = HyperGraph.normalize_adjacency_matrix(A, self.NW)
    
            # Build allelic partner counts from prune_pair_df (index-space)
            allelic_partners = defaultdict(set)
            if self.prune_pair_df is not None:
                try:
                    df = self.prune_pair_df[
                        (self.prune_pair_df["type"] == 0) &
                        (self.prune_pair_df["similarity"] >= allelic_similarity)
                    ][["contig1", "contig2"]]
                    for c1, c2 in df.values:
                        c1 = int(c1); c2 = int(c2)
                        allelic_partners[c1].add(c2)
                        allelic_partners[c2].add(c1)
                except Exception as e:
                    logger.warning(f"rescue_collapsed: failed to parse prune_pair_df: {e}")
    
            # Precompute group sets for faster membership changes
            group_sets = [set(g) for g in self.K]
            contig_to_groups = defaultdict(list)
            for gi, g in enumerate(group_sets):
                for c in g:
                    contig_to_groups[int(c)].append(gi)
    
            idx_to_vertices = self.idx_to_vertices
    
            collapsed_records = []  # (contig_name, n_dups, group_ids)
    
            # only consider contigs that are already anchored (collapsed is about missing copies)
            anchored = set(list_flatten(self.K))
    
            for c in list(anchored):
                c = int(c)
    
                # must have enough allelic partners
                if len(allelic_partners.get(c, set())) < min_allelic_partners:
                    continue
    
                # score to each group
                scores = []
                for gi, g in enumerate(group_sets):
                    if not g:
                        scores.append(0.0)
                        continue
                    s = float(A[c, list(g)].sum())
                    scores.append(s)
    
                # get top groups
                order = np.argsort(scores)[::-1]
                top1 = scores[order[0]]
                top2 = scores[order[1]] if len(order) > 1 else 0.0
    
                if top1 < min_signal:
                    continue
                if top1 <= 0:
                    continue
                if (top2 / top1) < ratio_threshold:
                    continue
    
                # decide number of duplicates/groups to place into
                m = min(max_dups_per_contig, len(order))
                # require at least 2 groups
                if m < 2:
                    continue
    
                top_groups = [int(x) for x in order[:m] if scores[int(x)] > 0]
                if len(top_groups) < 2:
                    continue
    
                contig_name = idx_to_vertices.get(c, str(c))
                collapsed_records.append((contig_name, len(top_groups), ",".join(map(str, top_groups))))
    
                if apply_to_clusters:
                    # duplicate contig into these groups (do not remove from original)
                    for gi in top_groups:
                        group_sets[gi].add(c)
    
            if not collapsed_records:
                logger.info("Collapsed rescue: no collapsed contigs detected.")
                return
    
            # Write collapsed list for downstream pqs dup
            with open(output_collapsed_list, "w") as f:
                f.write("contig\tn_dups\tgroups\n")
                for contig_name, ndup, groups in collapsed_records:
                    f.write(f"{contig_name}\t{ndup}\t{groups}\n")
    
            logger.info(f"Collapsed rescue: detected {len(collapsed_records)} collapsed contigs. "
                        f"Written to `{output_collapsed_list}`.")
    
            if apply_to_clusters:
                self.K = [sorted(list(s)) for s in group_sets]


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

        size_dict = self.vertices_idx_sizes
        try: 
            self.inc_chr_idx 
            pop_idx = []
            _K = []
            for i, group in enumerate(self.K):
                _size = sum(size_dict.get(x, 0) for x in group)
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
            _K = []
            for group in self.K:
                _size = sum(size_dict.get(x, 0) for x in group)
                if _size >= self.min_scaffold_length:
                    _K.append(group)
            
        
        logger.info(f"After filtering, {len(_K)} groups remain.")
            
        return _K 
    
    def filter_cluster_ctg(self, verbose=1):
        if verbose == 1: 
            logger.info(f"Removed groups less than {to_humanized2(self.min_scaffold_length)} in length. (--min-scaffold-length)")

        try: 
            self.inc_chr_idx 
            pop_idx = []
            _K = []
            for i, group in enumerate(self.K):
                _size = self.contigsizes.reindex(group).sum().values[0]
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
                filter(lambda x: self.contigsizes.reindex(x).sum().values[0] \
                            >= self.min_scaffold_length, 
                        self.K))
            
        return _K 
    
    def get_raw_contigsizes(self):
        contigsizes = self.contigsizes.copy()
        contigsizes.reset_index(inplace=True)
        # contigsizes[['source_chrom', 'start_end']] = contigsizes['chrom'].str.rsplit("|", n=1, expand=True)
        # contigsizes[['start', 'end']] = contigsizes['start_end'].str.split("_", n=1, expand=True).astype(int)
        
        # return contigsizes.groupby(['source_chrom'])['length'].sum()
        if hasattr(self, '_cached_raw_contigsizes') and self._cached_raw_contigsizes is not None:
            return self._cached_raw_contigsizes
        
        merged_sizes = defaultdict(int)
        for chrom, length in self.contig_sizes.items():
            source_chrom = chrom.rsplit('|', 1)[0] if '|' in chrom else chrom
            merged_sizes[source_chrom] += length
            
        self._cached_raw_contigsizes = dict(merged_sizes)
        return self._cached_raw_contigsizes

    def parse_trimmed_contig(self):
        contigsizes = self.contigsizes.copy()
        contigsizes.reset_index(inplace=True)
        contigsizes[['source_chrom', 'start_end']] = contigsizes['chrom'].str.rsplit("|", n=1, expand=True)
      
        contigsizes[['start', 'end']] = contigsizes['start_end'].str.split("_", n=1, expand=True).dropna().astype(int)


        for i, df in contigsizes.groupby('source_chrom', sort=False):
            df = df.sort_values('start')
            if len(df) > 1:
                first_idx = df.index[0]
                last_idx = df.index[-1] 
                sum_length = df['length'].sum()
                trim_length = df.loc[first_idx, 'start'] 
                if trim_length <= 0:
                    continue 
                # if sum_length <= 3 * trim_length:
                #     continue

                df.loc[first_idx, 'length'] += trim_length

                df.loc[last_idx, 'length'] += trim_length

            elif len(df) == 1:
                sum_length = df['length'].sum()
                first_idx = df.index[0]
                trim_length = df.loc[first_idx, 'start'] 
                if trim_length <= 0:
                    continue 
                # if sum_length <= 3 * trim_length:
                #     continue

                df.loc[first_idx, 'length'] += trim_length * 2

            else:
                pass
            contigsizes.loc[df.index, 'length'] = df['length']
        
        self.contigsizes = contigsizes.set_index('chrom')[['length']]
        # self.contigsizes.to_csv('test.contigsizes', sep='\t', header=None, index=True)
        self.contig_sizes = dict(zip(self.contigsizes.index, self.contigsizes['length']))

    def merge_split(self, secondary=False, confidence_threshold=0.0):
        idx_to_vertices = self.idx_to_vertices
        vertices_to_idx = self.vertices_idx

        clusters = [[idx_to_vertices[x] for x in group] for group in self.K]

        res_db = {}
        for i, group in enumerate(clusters):
            group = list(map(parse_split_contigs, group))
            for contig, start, end in group:
                if contig not in res_db:
                    res_db[contig] = defaultdict(int)
                res_db[contig][i] += (end - start)

        new_res_db = {}
        conflicted_contigs = []
        for contig, dist in res_db.items():
            sorted_dist = sorted(dist.items(), key=lambda x: x[1], reverse=True)
            best_group, best_len = sorted_dist[0]

            total_assigned_len = sum(dist.values())
            
            if total_assigned_len == 0:
                support_ratio = 0.0
            else:
                support_ratio = best_len / total_assigned_len

            if support_ratio < confidence_threshold:
                conflicted_contigs.append((contig, support_ratio))
            
            new_res_db[contig] = (best_group, best_len)
            
        if conflicted_contigs:
            logger.warning(f"Found {len(conflicted_contigs)} contigs with low support ratio (<{confidence_threshold}). "
                           "They might be misassembled or chimeric.")


        merged_contig_sizes =  self.get_raw_contigsizes()
        clusters_db = defaultdict(list)
        for contig, (group_idx, length) in new_res_db.items():
            clusters_db[group_idx].append(contig)
        
        
        if secondary:
            tmp_inc_chr_idx = {i: chr_idx for i, chr_idx in enumerate(self.inc_chr_idx)}
            default_chr_idx = self.inc_chr_idx[-1] if self.inc_chr_idx else 0
            self.inc_chr_idx = []

        tmp_merged_clusters = []
        active_group_mapping = {}
        for group_idx in range(len(clusters)): 
            group = clusters_db.get(group_idx, [])
            if group:
                active_group_mapping[group_idx] = len(tmp_merged_clusters)
                tmp_merged_clusters.append(group)
                if secondary:
                    self.inc_chr_idx.append(tmp_inc_chr_idx.get(group_idx, default_chr_idx))
            
        merged_clusters = tmp_merged_clusters
        
        new_split_clusters = [[] for _ in range(len(merged_clusters))]
        allowed_vertices = set(idx_to_vertices.keys())
        for v, v_idx in vertices_to_idx.items():
            if v_idx not in allowed_vertices:
                continue
            contig, _, _ = parse_split_contigs(v)
            if contig in new_res_db:
                old_group_idx = new_res_db[contig][0]
                if old_group_idx in active_group_mapping:
                    new_idx = active_group_mapping[old_group_idx]
                    new_split_clusters[new_idx].append(v_idx)

        
        keep_mask = [len(g) > 0 for g in merged_clusters]
        if not all(keep_mask):
            old_to_new = {}
            filtered_merged = []
            for old_i, keep in enumerate(keep_mask):
                if keep:
                    old_to_new[old_i] = len(filtered_merged)
                    filtered_merged.append(merged_clusters[old_i])

            merged_clusters = filtered_merged

            filtered_split = [[] for _ in range(len(merged_clusters))]
            for old_i, keep in enumerate(keep_mask):
                if not keep:
                    continue
                new_i = old_to_new[old_i]
                filtered_split[new_i] = new_split_clusters[old_i]
            new_split_clusters = filtered_split

            if secondary:
                self.inc_chr_idx = [x for x, keep in zip(self.inc_chr_idx, keep_mask) if keep]

            logger.debug(f"merge_split: removed {keep_mask.count(False)} empty groups.")



        return merged_clusters, new_split_clusters, merged_contig_sizes
        
    
    def to_cluster(self, output, merge=True):
        try:
            self.inc_chr_idx 
            secondary = True
        except AttributeError:
            secondary = False 
        
        if secondary and self.split and merge:
            clusters, _, _ = self.merge_split(secondary=secondary)
        elif self.split and merge:
            clusters, _, _ = self.merge_split()
        else:
            idx_to_vertices = self.idx_to_vertices

            clusters = list(map(lambda y: list(
                            map(lambda x: idx_to_vertices[x], y)), 
                            self.K))

        if secondary:
            if len(self.inc_chr_idx) != len(clusters):
                logger.warning(
                    f"to_cluster: inc_chr_idx length ({len(self.inc_chr_idx)}) != clusters length ({len(clusters)}); "
                    "will align by truncation + empty filtering."
                )
                n = min(len(self.inc_chr_idx), len(clusters))
                self.inc_chr_idx = self.inc_chr_idx[:n]
                clusters = clusters[:n]

            non_empty = [(c, chr_i) for c, chr_i in zip(clusters, self.inc_chr_idx) if c and len(c) > 0]
            clusters = [c for c, _ in non_empty]
            self.inc_chr_idx = [chr_i for _, chr_i in non_empty]
        else:
            clusters = [c for c in clusters if c and len(c) > 0]

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
