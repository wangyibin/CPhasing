#!/usr/bin/env python

import logging
import os
import dask.array as da
import gc
import numpy as np
import pandas as pd
import msgspec
import igraph as ig

from .._config import HYPERGRAPH_ORDER_DTYPE

from joblib import Parallel, delayed
from scipy.sparse import (
    identity, 
    dia_matrix, 
    csr_matrix, 
    coo_matrix,
)

logger = logging.getLogger(__name__)


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
        self.get_data()
       

    def get_data(self):
        """
        initial assign row and edges
        """
        self.nodes = np.array(sorted(self.edges.idx, key=self.edges.idx.get))
        self.row = np.array(self.edges.row)
        self.col = np.array(self.edges.col)
        self.shape = (len(self.nodes), max(self.col) + 1)
        
        self.removed_count = 0

    def remove_rows(self, contigs):
        """
        remove rows by contig list
        """     
        contigs = list(filter(lambda x: x not in set(self.edges.idx.values()), contigs))
        contig_idx = [self.edges.idx[i] for i in contigs]
        
        remove_idx = np.isin(self.row, contig_idx)
        
        self.removed_count += len(contig_idx)
        
        self.row = self.row[~remove_idx]
        self.col = self.col[~remove_idx]


        # self.nodes = np.delete(self.nodes, contig_idx)

        # self.shape = (len(self.nodes), max(self.col) + 1)

    def incidence_matrix(self, min_contacts=3):
        """
        The incidence matrix of hypergraph

        Params:
        --------
        min_contacts: int
            minimum contacts of contig

        """
        matrix = csr_matrix((np.ones(len(self.row), dtype=HYPERGRAPH_ORDER_DTYPE), 
                             (self.row, self.col)
                             ), 
                             shape=self.shape
                             )


        if min_contacts:
            non_zero_contig_idx = np.where(matrix.sum(axis=1).T >= min_contacts)[-1]
            matrix = matrix[non_zero_contig_idx]
            self.nodes = self.nodes[non_zero_contig_idx]
            logger.info(f"Total {self.shape[0] - non_zero_contig_idx.shape[0] - self.removed_count} "
                        f"low-singal (contacts < {min_contacts}) contigs were removed")
            
            non_zero_edges_idx = np.where(matrix.sum(axis=0) >= 2)[-1]
       
            matrix = matrix[:, non_zero_edges_idx]
            self.shape = matrix.shape 

        return matrix 

    @staticmethod
    def clique_expansion_init(H, P_allelic_idx=None, 
                              P_weak_idx=None, 
                              allelic_factor=-1,
                              min_value=1.0):
        """
        clique expansion/reduction
        """
        m = H.shape[1]
        W = identity(m, dtype=np.float32)

        D_e_num = (H.sum(axis=0) - 1).astype(HYPERGRAPH_ORDER_DTYPE)

        ## inverse diagonal matrix D_e
        D_e_inv = 1/D_e_num
        D_e_inv[D_e_inv == -np.inf] = 0
        D_e_inv = dia_matrix((D_e_inv, np.array([0])), 
                                W.shape, dtype=np.float32)
        
        H_T = H.T
        
        A = H @ W @ D_e_inv @ H_T

        mask = A >= min_value
        A.multiply(mask)

        if P_allelic_idx or P_weak_idx:
            # P = np.ones((H.shape[0], H.shape[0]), dtype=np.int8)
            # P[np.diag_indices_from(P)] = 0
            # P[P_idx[0], P_idx[1]] = 0
            # A = A * P 
            A = A.tolil()
            if P_allelic_idx:
                if allelic_factor == 0: 
                    A[P_allelic_idx[0], P_allelic_idx[1]] = 0
                else:
                    A[P_allelic_idx[0], P_allelic_idx[1]] = \
                        allelic_factor * A[P_allelic_idx[0], P_allelic_idx[1]]
            
            if P_weak_idx:
                A[P_weak_idx[0], P_weak_idx[1]] = 0
                
            A = A.tocsr()
        return A



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
        w_prev = W.data[0][idx]
        if e_c <= 2:
            W_tmp.append(w_prev)
            continue
        k = np.array([len(e_nodes & cluster_assignments[i]) for i in range(c)])
        
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

def extract_incidence_matrix2(mat, idx):
    A = coo_matrix(mat, copy=False)
    
    retain = np.isin(A.row, idx)
    row = A.row[retain]
    col = A.col[retain]
    data = A.data[retain]

    A = coo_matrix((data, (row, col)), shape=A.shape).asformat('csr')

    # non_zero_contig_idx = np.where(A.sum(axis=1).T != 0)[-1]
    # zero_contig_idx = np.where(A.sum(axis=1).T == 0)[-1]
    
    A = A[idx]

    if A.shape[0] != 0 and A.shape[1] !=0 : 
        non_zero_edges_idx = np.where(A.sum(axis=0).T > 1)[0]
        remove_edges_idx = np.where(A.sum(axis=0).T <= 1)[0]
        A = A[:, non_zero_edges_idx]
    else:
        remove_edges_idx = np.array([])

    return A, remove_edges_idx

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
            P_allelic_idx=None,
            P_weak_idx=None,
            allelic_factor=-1,
            resolution=1, 
            min_weight=1,
            threshold=0.01, 
            max_round=1, 
            chunk=10000, 
            threads=10, 
            outprefix="all"):
    """
    Iteratively reweight modularity maximization. （Kumar et al. Applied Network Science）
    
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

    # if P_allelic_idx or P_weak_idx:
    #     P_allelic_df = pd.concat(P_allelic_idx, axis=1)
    #     P = csr_matrix((np.ones(len(P_allelic_df), dtype=HYPERGRAPH_ORDER_DTYPE),
    #                     (P_allelic_df['contig1'], P_allelic_df['contig2'])), 
    #                     shape=([H.shape[0], len(P_allelic_df)])) 
        
    #     _P = P.T @ H
    #     remove_idx = (_P.toarray() == 2).any(axis=0)
    #     H = H[:, ~remove_idx]
       
    m = H.shape[1]
    
    ## diagonal matrix of weights
    W = identity(m, dtype=np.float32)
    
    ## D_e - I 
    D_e_num = (H.sum(axis=0) - 1).astype(HYPERGRAPH_ORDER_DTYPE)
    
    ## inverse diagonal matrix D_e
    D_e_inv = 1/D_e_num
    D_e_inv[D_e_inv == -np.inf] = 0
    D_e_inv = dia_matrix((D_e_inv, np.array([0])), 
                            W.shape, dtype=np.float32)
    
    H_T = H.T

  
    A = H @ W @ D_e_inv @ H_T
    
    if min_weight > 0:
        mask = A >= min_weight
        A = A.multiply(mask)
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

    if P_allelic_idx and P_weak_idx:
        P_allelic_df = pd.concat(P_allelic_idx, axis=1)
        # P = csr_matrix((np.ones(len(P_allelic_df), dtype=HYPERGRAPH_ORDER_DTYPE),
        #                 (P_allelic_df['contig1'], P_allelic_df['contig2'])), 
        #                 shape=([A.shape[0], len(P_allelic_df)]))
        
        # print(A.shape, allelic_factor.shape)
        # A = A.multiply(allelic_factor)

        # P = np.ones((H.shape[0], H.shape[0]), dtype=np.int8)
        # P[np.diag_indices_from(P)] = 0
        # P[P_idx[0], P_idx[1]] = 0
        # A = A * P 
        A = A.tolil()
        if P_allelic_idx:
            if allelic_factor == 0: 
                A[P_allelic_idx[0], P_allelic_idx[1]] = 0
            else:
                A[P_allelic_idx[0], P_allelic_idx[1]] = \
                       allelic_factor * A[P_allelic_idx[0], P_allelic_idx[1]]
        
        if P_weak_idx:
            if allelic_factor == 0:
                A[P_weak_idx[0], P_weak_idx[1]] = 0
            else:
                A[P_weak_idx[0], P_weak_idx[1]] = \
                    allelic_factor * A[P_weak_idx[0], P_weak_idx[1]]
            
        A = A.tocsr()
        

    try:
        G = ig.Graph.Weighted_Adjacency(A, mode='undirected', loops=False)
    except ValueError:
        return [list(range(A.shape[0]))]

    cluster_assignments = G.community_multilevel(weights='weight', resolution=resolution)
    cluster_results = list(cluster_assignments.as_cover())
    cluster_results = list(map(set, cluster_results))

  
    cluster_stat = list(map(len, cluster_results))
 
    # del A, G
    # gc.collect()
    
    if max_round > 1:
        W_prev = W.copy()

    iter_round = 1
    diff_value = 1

    while diff_value > threshold:

        if iter_round >= max_round:
            break

        a = da.from_array(H_T, chunks=(chunk, H_T.shape[1]))

        args = [(a, cluster_results, W, i, chunk) 
                        for i in range(a.npartitions)]

        res = Parallel(n_jobs=min(a.npartitions, threads))(
                        delayed(reweight)(i, j, k, l, m) 
                            for i, j, k, l, m in args)

        res = [i[1] for i in sorted(res, key=lambda x: x[0])]
        W.data[0] = np.concatenate(res)

        del a, args, res
        gc.collect()

        A = H @ W @ D_e_inv @ H_T
        if min_weight > 0:
            mask = A >= min_weight
            A = A.multiply(mask)
        ## normalization
        # A = A.toarray() * NW
        # row, col = np.nonzero(A)
        # values = A[row, col]
        # A = csr_matrix((values, (row, col)), shape=A.shape)
        # A.setdiag(0)
        # A.prune()
        if P_allelic_idx or P_weak_idx:
            # P = np.ones((H.shape[0], H.shape[0]), dtype=np.int8)
            # P[np.diag_indices_from(P)] = 0
            # P[P_idx[0], P_idx[1]] = 0
            # A = A * P 
            A = A.tolil()
            if P_allelic_idx:
                if allelic_factor == 0: 
                    A[P_allelic_idx[0], P_allelic_idx[1]] = 0
                else:
                    A[P_allelic_idx[0], P_allelic_idx[1]] = \
                            allelic_factor * A[P_allelic_idx[0], P_allelic_idx[1]]
            if P_weak_idx:
                A[P_weak_idx[0], P_weak_idx[1]] = 0
                
            A = A.tocsr()

        try:
            G = ig.Graph.Weighted_Adjacency(A, mode='undirected', loops=False)
        except ValueError:
            return cluster_assignments
        
        # ## use igraph 
        cluster_assignments = G.community_multilevel(weights='weight', resolution=resolution)
        cluster_results = list(cluster_assignments.as_cover())
        cluster_results = list(map(set, cluster_results))

        ## use networkx 
        # G = G.to_networkx()
        # cluster_assignments = nx.community.louvain_communities(G, resolution=resolution)
        # del A, G
        # gc.collect()

        cluster_stat = list(map(len, cluster_results))
        logger.info(f"Cluster Statistics: {cluster_stat}")
        

        diff_value = max(abs(W.data[0] - W_prev.data[0]))
        W_prev = W.copy()
        
        iter_round += 1
    
    return A, cluster_assignments, cluster_results

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

