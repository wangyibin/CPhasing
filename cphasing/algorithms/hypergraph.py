#!/usr/bin/env python

import logging
import os
import gc
import numpy as np
import pandas as pd
import msgspec
import igraph as ig

from .._config import HYPERGRAPH_ORDER_DTYPE, HYPERGRAPH_COL_DTYPE

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
    # count: contact of a contig in a read (> 2)
    contigsizes: contig sizes
    mapq: mapq of alignments
    """
    idx: dict 
    row: list
    col: list
    # count: list
    contigsizes: dict 
    mapq: list

    def to_numpy(self):
        self.row = np.array(self.row, dtype=np.int32)
        self.col = np.array(self.col, dtype=HYPERGRAPH_COL_DTYPE)
        # self.count = np.array(self.count, dtype=np.uint32)
        self.mapq = np.array(self.mapq, dtype=np.int8)

    def to_list(self):
        self.row = self.row.tolist()
        self.col = self.col.tolist()
        # self.count = self.count.tolist()
        self.mapq = self.mapq.tolist()
        

def merge_hyperedges(HE_list: list) -> HyperEdges:
    
    if len(HE_list) <= 1:
        return HE_list[0]
    
    init_idx = HE_list[0].idx 
    init_HE = HE_list[0]
    
    for HE in HE_list[1:]:
        assert HE.idx == init_idx, "HyperEdges must build in the same contigsizes"
        init_HE.row += HE.row 
        init_HE.col += HE.col 
        if HE.mapq:
            init_HE.mapq += HE.mapq
        # if HE.count:
        #     init_HE.count += HE.count

    return init_HE


class HyperGraph:
    """
    HyperGraph 

    Params:
    --------
    edges: HyperEdges
        HyperEdges
    """

    def __init__(self, edges, min_quality=1):
        self.edges = edges 
        self.contigsizes = self.edges.contigsizes
        self.min_quality = min_quality
        self.get_data()

    def get_data(self):
        """
        initial assign row and edges
        """
        if len(self.edges.row) == 0:
            raise ValueError("No hyperedges found.")

        if self.min_quality > 1 and len(self.edges.mapq) > 1:
            self.mapq = np.array(self.edges.mapq, dtype=np.int8)

            retain_idx = self.mapq >= self.min_quality
            # self.nodes = np.array(sorted(self.edges.idx, key=self.edges.idx.get))
            self.nodes = np.array([k for k, v in sorted(self.edges.idx.items(), key=lambda item: item[1])])
            self.row = self.edges.row
            total_edge_counts = self.row.shape[0]
            self.row = self.row[retain_idx]
            self.col = self.edges.col[retain_idx]
            # self.count = self.edges.count[retain_idx]
           

            # self.idx = np.sort(np.array(pd.unique(self.row)))
            self.idx = np.unique(self.row)
            remove_idx = np.isin(np.arange(0, len(self.nodes)), self.idx)
            # self.mapq = self.mapq[retain_idx]
            self.remove_contigs = self.nodes[~remove_idx]
            logger.info(f"Removed `{total_edge_counts - np.count_nonzero(retain_idx):,}` "
                            f"hyperedges that mapq < {self.min_quality}.")
        
        else:
            self.mapq = np.array([])
            # self.nodes = np.array(sorted(self.edges.idx, key=self.edges.idx.get))
            self.nodes = np.array([k for k, v in sorted(self.edges.idx.items(), key=lambda item: item[1])])
            self.row = self.edges.row
            self.col = self.edges.col
            # self.count = self.edges.count
            self.remove_contigs = np.array([])

        self.shape = (len(self.nodes), max(self.col) + 1)
        self.removed_count = 0

    def remove_rows(self, contigs):
        """
        remove rows by contig list
        """     
        contigs = [x for x in contigs if x not in self.edges.idx.values()]
        contig_idx = np.array([self.edges.idx[i] for i in contigs], dtype=int)

        remove_idx = np.in1d(self.row, contig_idx)

        self.removed_count += len(contig_idx)

        self.row = self.row[~remove_idx]
        self.col = self.col[~remove_idx]
        # self.count = self.count[~remove_idx]

        # if len(self.mapq) > 1:
        #     self.mapq = self.mapq[~remove_idx]

        # self.nodes = np.delete(self.nodes, contig_idx)

        # self.shape = (len(self.nodes), max(self.col) + 1)

    def extract_rows(self, contigs):
        """
        extract rows by contig list
        """     
        contigs = [x for x in contigs if x not in self.edges.idx.values()]
        contig_idx = [self.edges.idx[i] for i in contigs]

        remove_idx = np.in1d(self.row, contig_idx)

        self.removed_count += len(contig_idx)

        self.row = self.row[remove_idx]
        self.col = self.col[remove_idx]
        # self.count = self.count[remove_idx]
        # if len(self.mapq) > 1:
        #     self.mapq = self.mapq[remove_idx]


    def incidence_matrix(self, min_contacts=1):
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
        # matrix = csr_matrix((self.count, (self.row, self.col)), 
        #                      shape=self.shape)

        if min_contacts:
            non_zero_contig_idx = matrix.sum(axis=1).T.A1 >= min_contacts
            matrix = matrix[non_zero_contig_idx]
            self.remove_contigs = self.nodes[~non_zero_contig_idx]
            self.nodes = self.nodes[non_zero_contig_idx]
           
            logger.info(f"Total {self.shape[0] - matrix.shape[0]:,} "
                        f"low-singal (contacts < {min_contacts}) contigs were removed (--min-contacts). ")
            
            non_zero_edges_idx = matrix.sum(axis=0).A1 >= 2
            matrix = matrix[:, non_zero_edges_idx]
         
            self.shape = matrix.shape 

        else:
            self.remove_contigs = []
       
        del self.row, self.col
        gc.collect()

        return matrix 

    @staticmethod
    def clique_expansion_init(H, P_allelic_idx=None, 
                              P_weak_idx=None, 
                              NW=None,
                              W=None,
                              allelic_factor=-1,
                              min_weight=0.1):
        """
        clique expansion/reduction
        """
        m = H.shape[1]

        # D_e - I
        D_e_num = (H.sum(axis=0) - 1).astype(HYPERGRAPH_ORDER_DTYPE)
        
        # W = identity(m, dtype=np.float32)
        # if W is not None:
        #     W = dia_matrix((W, np.array([0])), shape=(m, m), dtype=np.float32)
        # W = dia_matrix((D_e_num, np.array([0])), (m, m), dtype=np.float32)
        
       
        ## inverse diagonal matrix D_e
        D_e_inv = 1/D_e_num
        D_e_inv[D_e_inv == -np.inf] = 0
        D_e_inv = dia_matrix((D_e_inv, np.array([0])), 
                                (m, m), dtype=np.float32)
    
        if W is not None:
            A = H @ W @ D_e_inv @ H.T
        else:
            A = H.dot(D_e_inv).dot(H.T)

        if NW is not None:
            A = A.toarray() * NW
            row, col = np.nonzero(A)
            values = A[row, col]
            A = csr_matrix((values, (row, col)), shape=A.shape)

        if min_weight > 0:
            mask = A >= min_weight
            # mask = (A >= min_weight).astype(bool)
            # mask += (A <= -min_weight).astype(bool)
            A = A.multiply(mask)
    
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
                    A[P_allelic_idx[0], P_allelic_idx[1]] *= allelic_factor
            
            if P_weak_idx:
                A[P_weak_idx[0], P_weak_idx[1]] = 0 
                
            A = A.tocsr()

        # if P_allelic_idx or P_weak_idx:
        #     A = A.tocoo()
        #     if P_allelic_idx:
        #         if allelic_factor == 0: 
        #             A.data[np.in1d(A.row, P_allelic_idx[0]) & np.in1d(A.col, P_allelic_idx[1])] = 0
        #         else:
        #             A.data[np.in1d(A.row, P_allelic_idx[0]) & np.in1d(A.col, P_allelic_idx[1])] *= allelic_factor
                    
        #     if P_weak_idx:
        #         A.data[np.in1d(A.row, P_weak_idx[0]) & np.in1d(A.col, P_weak_idx[1])] = 0 

        #     A = A.tocsr()


        return A

    
    @staticmethod
    def to_contacts(A, nodes, 
                    NW=None,
                    min_weight=1.0, 
                    output="hypergraph.expansion.contacts"):
        # A = HyperGraph.clique_expansion_init(H, NW=NW, min_weight=min_weight).tocoo()
        A = A.tocoo()
        V = nodes
       
        df = pd.DataFrame({0: V[A.row], 1: V[A.col], 2: A.data})
        df = df[df[0] <= df[1]]
        df = df[(df[2] >= min_weight) | (df[2] <= -min_weight)]
        df = df[abs(df[2]) != np.inf]
        # df = df.query(f"0 <= 1 & 2 >= {min_weight} & abs(2) != inf")
        df.to_csv(output, sep='\t', header=False, index=False)
                
        logger.info(f"Export hypergraph into `{output}`")


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
    # A = coo_matrix(mat, copy=False)
    
    # retain = np.isin(A.row, idx)
    # row, col, data = A.row[retain], A.col[retain], A.data[retain]

    # A = coo_matrix((data, (row, col)), shape=A.shape).asformat('csr')
    # A = A[idx]
    A = mat[idx]
    if A.shape[0] != 0 and A.shape[1] !=0 : 
        A_sum = A.sum(axis=0).T
        non_zero_edges_idx = np.where(A_sum > 1)[0]
        remove_edges_idx = np.where(A_sum <= 1)[0]
        A = A[:, non_zero_edges_idx]
    else:
        remove_edges_idx = np.array([])
        non_zero_edges_idx = np.array([])

    return A, remove_edges_idx, non_zero_edges_idx

def remove_incidence_matrix(mat, idx):   
    if not isinstance(mat, csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    
    idx = np.setdiff1d(np.arange(mat.shape[0]), idx)
    mat = mat[idx].T.tolil()
    ## remove singleton or empty hyperedge
    remove_col_index = np.where(mat.sum(axis=1).T > 1)[1]
    mat = mat[remove_col_index]

    return mat.T.tocsr(), remove_col_index

def IRMM(H, A=None,
            NW=None, 
            P_allelic_idx=None,
            P_weak_idx=None,
            allelic_factor=-1,
            cross_allelic_factor=0.3,
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

    if max_round > 1:
        import dask.array as da 
    # if P_allelic_idx or P_weak_idx:
    #     P_allelic_df = pd.concat(P_allelic_idx, axis=1)
    #     P = csr_matrix((np.ones(len(P_allelic_df), dtype=HYPERGRAPH_ORDER_DTYPE),
    #                     (P_allelic_df['contig1'], P_allelic_df['contig2'])), 
    #                     shape=([H.shape[0], len(P_allelic_df)])) 
        
    #     _P = P.T @ H
    #     remove_idx = (_P.toarray() == 2).any(axis=0)
    #     H = H[:, ~remove_idx]
    
    if A is None or max_round > 1:
        m = H.shape[1]
        
        ## diagonal matrix of weights
        W = identity(m, dtype=np.float32)
        
        ## D_e - I 
        # D_e_num = H.getnnz(axis=0).astype(HYPERGRAPH_ORDER_DTYPE)
    
        D_e_num = (H.sum(axis=0) - 1).astype(HYPERGRAPH_ORDER_DTYPE)
        
        ## inverse diagonal matrix D_e
        D_e_inv = 1/D_e_num
        D_e_inv[D_e_inv == -np.inf] = 0
        D_e_inv = dia_matrix((D_e_inv, np.array([0])), 
                                W.shape, dtype=np.float32)

        if max_round <= 1:
            A = H @ W @ D_e_inv @ H.T

        else:
            H_T = H.T
            A = H @ W @ D_e_inv @ H_T
        
        if min_weight > 0:
            mask = A >= min_weight
            # mask = (A >= min_weight).astype(bool)
            # mask += (A <= -min_weight).astype(bool)
            A = A.multiply(mask)
        # normalization
        if NW is not None:
            A = A.toarray()
            diag_A = np.diagonal(A)
            NW = 1 / np.sqrt(np.outer(diag_A, diag_A))
            np.fill_diagonal(NW, 1)

            A = A * NW
            row, col = np.nonzero(A)
            values = A[row, col]
            A = csr_matrix((values, (row, col)), shape=A.shape, dtype=np.float32)

        A.setdiag(0)
        A.prune()

        if max_round <= 1:
            del W, D_e_num, D_e_inv
            gc.collect() 

        raw_A = A.copy()
        if P_allelic_idx or P_weak_idx:
            A = A.tolil()
            if P_allelic_idx:
        
                if allelic_factor == 0: 
                    A[P_allelic_idx[0], P_allelic_idx[1]] = 0
                else:
                    A[P_allelic_idx[0], P_allelic_idx[1]] *= allelic_factor
            if P_weak_idx:
                if cross_allelic_factor == 0:
                    A[P_weak_idx[0], P_weak_idx[1]] = 0
                else:
                    A[P_weak_idx[0], P_weak_idx[1]] *= cross_allelic_factor
                
            A = A.tocsr()
    else:
        raw_A = A.copy()
        if P_allelic_idx or P_weak_idx:
            A = A.tolil()
            if P_allelic_idx:
        
                if allelic_factor == 0: 
                    A[P_allelic_idx[0], P_allelic_idx[1]] = 0
                else:
                    A[P_allelic_idx[0], P_allelic_idx[1]] *= allelic_factor
            if P_weak_idx:
                if cross_allelic_factor == 0:
                    A[P_weak_idx[0], P_weak_idx[1]] = 0
                else:
                    A[P_weak_idx[0], P_weak_idx[1]] *= cross_allelic_factor
                
            A = A.tocsr()


    try:
        G = ig.Graph.Weighted_Adjacency(A, mode='undirected', loops=False)
    except ValueError:
        return raw_A, A, None, []

    cluster_assignments = G.community_multilevel(weights='weight', resolution=resolution)
    # cluster_assignments = G.community_leiden(weights='weight', resolution=resolution)
    # import leidenalg as la
    # cluster_assignments = la.find_partition(G, la.CPMVertexPartition, resolution_parameter=0.01)
    cluster_results = list(map(set, list(cluster_assignments.as_cover())))
  
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
                    A[P_allelic_idx[0], P_allelic_idx[1]] *= allelic_factor
            if P_weak_idx:
                if cross_allelic_factor == 0:
                    A[P_weak_idx[0], P_weak_idx[1]] = 0
                else:
                    A[P_weak_idx[0], P_weak_idx[1]] *= cross_allelic_factor

            A = A.tocsr()

        try:
            G = ig.Graph.Weighted_Adjacency(A, mode='undirected', loops=False)
        except ValueError:
            return A, None, []
        
        # ## use igraph 
        cluster_assignments = G.community_multilevel(weights='weight', resolution=resolution)
        cluster_results = list(map(set, list(cluster_assignments.as_cover())))

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
    
    return raw_A, A, cluster_assignments, cluster_results
