#!/usr/bin/env python

import logging
import os
import dask.array as da
import gc
import numpy as np
import msgspec
import networkx as nx
import pandas as pd 
import igraph as ig

from collections import OrderedDict
# from joblib import Memory
from joblib import Parallel, delayed
from scipy.sparse import (
    identity, 
    dia_matrix, 
    csr_matrix, 
)

logger = logging.getLogger(__name__)

# memory = Memory('./cachedir', verbose=0)

class HyperEdges(msgspec.Struct):
    """
    A type describing the hyper edges

    Params:
    --------
    idx: contig idx 
    row: contigs
    col: hyperedges
    """
    idx: dict 
    row: list
    col: list

class HyperGraph:
    """
    HyperGraph 

    Params:
    --------
    edges: HyperEdges
        HyperEdges
    """
    def __init__(self, edges):
        self.edges = edges 
        self.shape = (len(self.edges.idx), max(self.edges.col) + 1)
        self.nodes = np.array(sorted(self.edges.idx, key=self.edges.idx.get))

    def incidence_matrix(self, remove_zero=True):
        matrix = csr_matrix((np.ones(len(self.edges.row), dtype=np.int8), 
                             (self.edges.row, self.edges.col)
                             ), 
                             shape= self.shape
                             )

        if remove_zero:
            non_zero_contig_idx = np.where(matrix.sum(axis=1).T != 0)[-1]
            matrix = matrix[non_zero_contig_idx]
            self.nodes = self.nodes[non_zero_contig_idx]
            logger.info(f"{self.shape[0] - non_zero_contig_idx.shape[0]} non-singal contigs were removed")
            
            # non_zero_edges_idx = np.where(matrix.sum(axis=0).T != 0)[-1]
            # matrix = matrix[:, non_zero_edges_idx]
            self.shape = matrix.shape 
            
        return matrix 

def calc_new_weight(k, c, e_c, m):
    """
    single function for calculate new weight.

    Params:
    --------
    k: np.array
        overlap of cluster in this edges.
    c: int
        number of clusters.
    e_c: int:
        number of nodes in this hyperedge.
    m: int:
        total number of hyperedges.

    Returns:
    --------
    float: 
        new weight.
    
    Examples:
    --------
    >>> k
    array([3, 4, 1, 4, 5])
    >>> calc_new_weight(k, 5, 17, 1000)
    0.523
    """
    w_new = (1/1+k).sum() * (e_c+c) / m

    return w_new


def reweight(d_array, cluster_assignments, W, n, chunk):
    """
    reweitht hyperedge according by previous cluster results.

    Params:
    --------
    d_array: dask.array
        dask array of hypergraph incidence matrix.
    cluster_assignments: list
        list of cluster assignments from community detection.
    W: dia_matrix
        diagonal matrix of weights.
    n: int
        number of array partitions.
    chunk: int
        size of chunk in dask array.

    Returns:
    --------
    int:
        number of chunk.
    np.array:
        new weights.

    Examples:
    --------
    >>> reweitht(d_array, cluster_assignments, W, 0, 10000)
    """
    array = d_array.partitions[n].compute()
    c = len(cluster_assignments)
    m = W.shape[0]
    W_tmp = []
    idx = 0 + n * chunk

    for e in array:
        e_nodes = set(e.indices)
        e_c = len(e_nodes)
        k = np.array([len(e_nodes & cluster_assignments[i]) for i in range(c)])
        w_prev = W.data[0][idx]
        w_new = calc_new_weight(k, c, e_c, m)
        W_tmp.append(0.5 * (w_new + w_prev))
        idx += 1
    
    return n, np.array(W_tmp)

def extract_incidence_matrix(mat, idx):
    if not isinstance(mat, csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    
    mat = mat[idx].T.tolil()
    ## remove singleton or empty hyperedge
    remove_col_index = np.where(mat.sum(axis=1).T > 1)[1]
    mat = mat[remove_col_index]
    
    return mat.T.tocsr(), remove_col_index

def remove_incidence_matrix(mat, idx):
    if not isinstance(mat, csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    
    idx = np.setdiff1d(np.arange(mat.shape[0]), idx)
    mat = mat[idx].T.tolil()
    ## remove singleton or empty hyperedge
    remove_col_index = np.where(mat.sum(axis=1).T > 1)[1]
    mat = mat[remove_col_index]

    return mat.T.tocsr(), remove_col_index

def IRMM(H, # NW, 
            P_idx=None,
            resolution=1, 
            threshold=0.01, 
            max_round=10, 
            chunk=10000, 
            threads=10, 
            outprefix="all"):
    """
    Iteratively reweight modularity maximization.
    
    Params:
    --------
    H: csr_matrix
        incident matrix of hypergraph
    vertices: list or list-like
        names of vetices.
    P_idx: list, default None
        prune vertices edges
    threshold: float, default 0.01
        threshold of reweight
    max_round: int, default 10
        maximize round of reweight
    chunk: int, default 10000
        number of chunk for parallel
    threads: int, default 10
        number of threads

    Returns:
    --------
    list:
        list of cluster results.
    
    Examples:
    --------
    >>> IRMM(H, vertices)
    """
    m = H.shape[1]
    
    ## diagonal matrix of weights
    W = identity(m, dtype=np.float32)
    
    ## D_e - I 
    D_e_num = (H.sum(axis=0) - 1).astype(np.int8)
    
    # D_e = dia_matrix((D_e_num, np.array([0])), W.shape, dtype=np.int8)
    
    ## inverse diagonal matrix D_e
    D_e_inv = 1/D_e_num
    D_e_inv[D_e_inv == -np.inf] = 0
    D_e_inv = dia_matrix((D_e_inv, np.array([0])), 
                            W.shape, dtype=np.float32)
    
    H_T = H.T
    
    A = H @ W @ D_e_inv @ H_T
    ## normalization
    # A = A.toarray() * NW
    # row, col = np.nonzero(A)
    # values = A[row, col]
    # A = csr_matrix((values, (row, col)), shape=A.shape)

    # A.setdiag(0)
    # A.prune()

    if max_round <= 1:
        del W, D_e_num, D_e_inv, H_T
        gc.collect() 

    if P_idx:
        # P = np.ones((H.shape[0], H.shape[0]), dtype=np.int8)
        # P[np.diag_indices_from(P)] = 0
        # P[P_idx[0], P_idx[1]] = 0
        # A = A * P 
        A = A.tolil()
        A[P_idx[0], P_idx[1]] = 0
        A = A.tocsr()

    try:
        G = ig.Graph.Weighted_Adjacency(A, mode='undirected', loops=False)
    except ValueError:
        return [list(range(A.shape[0]))]

    cluster_assignments = G.community_multilevel(weights='weight', resolution=resolution)
    cluster_assignments = list(cluster_assignments.as_cover())
    cluster_assignments = list(map(set, cluster_assignments))
    # G.write_graphml(f'{outprefix}.out.edgelist')
    
    # G = G.to_networkx()
    # nx.write_weighted_edgelist(G, f'{outprefix}.out.edgelist')
    # cluster_assignments = nx.community.louvain_communities(G, resolution=resolution)
    cluster_stat = list(map(len, cluster_assignments))
    logger.info(f"Cluster Statistics: {cluster_stat}")

    
    del A, G
    gc.collect()
    
    if max_round > 1:
        W_prev = W.copy()

    iter_round = 1
    diff_value = 1

    while diff_value > threshold:

        if iter_round >= max_round:
            break

        a = da.from_array(H_T, chunks=(chunk, H_T.shape[1]))

        args = [(a, cluster_assignments, W, i, chunk) 
                        for i in range(a.npartitions)]

        res = Parallel(n_jobs=min(a.npartitions, threads))(
                        delayed(reweight)(i, j, k, l, m) 
                            for i, j, k, l, m in args)

        res = [i[1] for i in sorted(res, key=lambda x: x[0])]
        W.data[0] = np.concatenate(res)

        del a, args, res
        gc.collect()

        A = H @ W @ D_e_inv @ H_T
        ## normalization
        # A = A.toarray() * NW
        # row, col = np.nonzero(A)
        # values = A[row, col]
        # A = csr_matrix((values, (row, col)), shape=A.shape)
        # A.setdiag(0)
        # A.prune()
        if P_idx:
            A = A.tolil()
            A[P_idx[0], P_idx[1]] = 0
            A = A.tocsr()

        try:
            G = ig.Graph.Weighted_Adjacency(A, mode='undirected', loops=False)
        except ValueError:
            return cluster_assignments
        
        # ## use igraph 
        # cluster_assignments = G.community_multilevel(weights='weight', resolution=resolution)
        # cluster_assignments = list(cluster_assignments.as_cover())
        # cluster_assignments = list(map(set, cluster_assignments))

        ## use networkx 
        G = G.to_networkx()
        cluster_assignments = nx.community.louvain_communities(G, resolution=resolution)
        del A, G
        gc.collect()

        cluster_stat = list(map(len, cluster_assignments))
        logger.info(f"Cluster Statistics: {cluster_stat}")
        

        diff_value = max(abs(W.data[0] - W_prev.data[0]))
        W_prev = W.copy()
        
        iter_round += 1
    
    return cluster_assignments

## old version through hypernetx
def generate_hypergraph(edges):
    """
    generate hypergraph incidence matrix

    Params:
    --------
    edges: dict
        edges for the hypergraph constructions

    Returns:
    --------
    H: csc_matrix
        incidence matrix for hypergraph
    vertices: list
        list of contig names

    Examples:
    --------
    >>> H, vertices = generate_hypergraph(edges)
    """
    pass
    # import hypernetx as hnx
    # HG = hnx.Hypergraph(edges, use_nwhy=True)

    # del edges
    # gc.collect()

    # H = HG.incidence_matrix().astype(np.int8)
    # vertices = list(HG.nodes)

    
    # return H, vertices

# generate_hypergraph = memory.cache(generate_hypergraph)

