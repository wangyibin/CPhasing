#!/usr/bin/env python

import logging
import os
import dask.array as da
import gc
import numpy as np
import pandas as pd 
import igraph as ig
import hypernetx as hnx

from collections import OrderedDict
from joblib import Memory
from joblib import Parallel, delayed
from scipy.sparse import (
    identity, 
    dia_matrix, 
    csr_matrix, 
    lil_matrix
)

logger = logging.getLogger(__name__)

memory = Memory('./cachedir', verbose=0)

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

def IRMM(H, P_idx=None, 
            threshold=0.01, 
            max_round=10, chunk=10000, 
            threads=10):
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
    D_e_num = H.sum(axis=0) - 1
    # D_e = dia_matrix((D_e_num, np.array([0])), W.shape, dtype=np.int8)
    
    ## inverse diagonal matrix D_e
    D_e_inv = 1/D_e_num
    D_e_inv[D_e_inv == -np.inf] = 0
    D_e_inv = dia_matrix((D_e_inv, np.array([0])), 
                            W.shape, dtype=np.float32)
    
    H_T = H.T
    A = H @ W @ D_e_inv @ H_T
    # A.setdiag(0)
    # A.prune()

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

    cluster_assignments = G.community_multilevel(weights='weight')
    cluster_assignments = list(cluster_assignments.as_cover())
    cluster_assignments = list(map(set, cluster_assignments))
    print(list(map(len, cluster_assignments)))
    del A, G
    gc.collect()

    W_prev = W.copy()
    iter_round = 1
    diff_value = 1

    while diff_value > threshold:
        print(diff_value)

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

        cluster_assignments = G.community_multilevel(weights='weight')
        cluster_assignments = list(cluster_assignments.as_cover())
        cluster_assignments = list(map(set, cluster_assignments))
        print(list(map(len, cluster_assignments)))
        del A, G
        gc.collect()

        diff_value = max(abs(W.data[0] - W_prev.data[0]))
        W_prev = W.copy()
        
        iter_round += 1
    
    return cluster_assignments
    # return list(map(lambda y: list(map(lambda x: list(vertices)[x], y)), 
    #             cluster_assignments))

def process_pore_c_table(pore_c_table, npartition, min_order, max_order, min_length):
    
    df = pore_c_table.partitions[npartition]
    df = df.compute()
    
    df = (df.assign(length=lambda x: x.end - x.start)
            .query(f"length >= {min_length}")
            .set_index('read_idx')
    )
    df_grouped = df.groupby('read_idx')['chrom']
    df_grouped_nunique = df_grouped.nunique()
    df = df.loc[(df_grouped_nunique >= min_order) 
                    & (df_grouped_nunique <= max_order)]
    edges = df.groupby('read_idx')['chrom'].unique()

    return edges

process_pore_c_table = memory.cache(process_pore_c_table)

def generate_hypergraph(pore_c_table, 
                            min_order=2, 
                            max_order=15, 
                            min_alignments=500,
                            threads=4):
    """
    generate hypergraph incidence matrix

    Params:
    --------
    pore_c_table: dask.dataframe
        pore_c table, at least have four columns: read_name, chrom, start, end.
    min_order: int, default 2
        minimum contig order of pore-c reads
    max_order: int, default 20
        maximum contig order of pore-c reads
    min_length: int, default 100
        minimum length of alignments
    threads: int, default 10
        number of threads

    Returns:
    --------
    H: csc_matrix
        incidence matrix for hypergraph
    vertices: list
        list of contig names

    Examples:
    --------
    >>> H, vertices = generate_hypergraph("pore_c.pq")
    """
    logger.info("Processing Pore-C table ...")
    
    args = []
    for i in range(pore_c_table.npartitions):
        args.append((pore_c_table, i,
                    min_order, max_order, min_alignments))

    res = Parallel(n_jobs=min(pore_c_table.npartitions, threads))(
                        delayed(process_pore_c_table)(i, j, k, l, m) 
                            for i, j, k, l, m in args)
    
    res_df = pd.concat(res)
    edges = res_df.values

    logger.info(f"Only retained Pore-C concatemer that: \n"
                    f"\talignment length >= {min_alignments}\n"
                    f"\t{min_order} <= contig order <= {max_order}")
    del res, res_df 
    gc.collect()

    edges = dict(enumerate(set(map(tuple, edges))))
    HG = hnx.Hypergraph(edges, use_nwhy=True)

    del edges
    gc.collect()

    H = HG.incidence_matrix().astype(np.int8)
    vertices = list(HG.nodes)

    logger.info(f"Generated hypergraph that containing {H.shape[0]} vertices"    
                    f" and {H.shape[1]} hyperedges.")

    return H, vertices

generate_hypergraph = memory.cache(generate_hypergraph)

