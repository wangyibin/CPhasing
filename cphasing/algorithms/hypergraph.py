#!/usr/bin/env python

import logging
import os
import gc
import random

import numpy as np
import polars as pl
import pandas as pd
import msgspec
import igraph as ig

from .._config import HYPERGRAPH_ORDER_DTYPE, HYPERGRAPH_COL_DTYPE
from joblib import Parallel, delayed
from scipy.sparse import (
    triu,
    identity, 
    dia_matrix, 
    csr_matrix, 
    coo_matrix,
)

logger = logging.getLogger(__name__)


# class HyperEdges(msgspec.Struct):
#     """
#     A type describing the hyper edges

#     Params:
#     --------
#     idx: contig idx 
#     row: contigs
#     col: hyperedges
#     # count: contact of a contig in a read (> 2)
#     contigsizes: contig sizes
#     mapq: mapq of alignments
#     """
#     idx: dict 
#     row: list
#     col: list
#     count: list
#     contigsizes: dict 
#     mapq: list

#     def to_numpy(self):
#         self.row = np.array(self.row, dtype=np.int32)
#         self.col = np.array(self.col, dtype=HYPERGRAPH_COL_DTYPE)
#         self.count = np.array(self.count, dtype=np.uint32)
#         self.mapq = np.array(self.mapq, dtype=np.int8)

#     def to_list(self):
#         self.row = self.row.tolist()
#         self.col = self.col.tolist()
#         self.count = self.count.tolist()
#         self.mapq = self.mapq.tolist()


def numpy_enc_hook(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist() 
    if isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    raise TypeError(f"Encoding objects of type {type(obj)} is unsupported")

def numpy_dec_hook(type, obj):
    if type is np.ndarray:
        return np.array(obj)
    elif type is list:
        return np.array(obj)
    
    return obj

class HyperEdges(msgspec.Struct):
    idx: dict 
    row: np.ndarray
    col: np.ndarray
    count: np.ndarray
    contigsizes: dict 
    mapq: np.ndarray

    @classmethod
    def from_lists(cls, idx, row, col, count, contigsizes, mapq):
        return cls(
            idx=idx,
            row=np.array(row, dtype=HYPERGRAPH_COL_DTYPE),
            col=np.array(col, dtype=HYPERGRAPH_COL_DTYPE),
            count=np.array(count, dtype=np.uint32),
            contigsizes=contigsizes,
            mapq=np.array(mapq, dtype=np.int8)
        )
    
    def to_numpy(self):
        self.row = np.asarray(self.row, dtype=HYPERGRAPH_COL_DTYPE)
        self.col = np.asarray(self.col, dtype=HYPERGRAPH_COL_DTYPE)
        self.count = np.asarray(self.count, dtype=np.uint32)
        self.mapq = np.asarray(self.mapq, dtype=np.int8)

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
        if HE.count:
            init_HE.count += HE.count

    return init_HE


class HyperGraph:
    """
    HyperGraph 

    Params:
    --------
    edges: HyperEdges
        HyperEdges
    """

    def __init__(self, edges, min_quality=1, mapq_filter=True):
        self.edges = edges 
        self.contigsizes = self.edges.contigsizes
        self.min_quality = min_quality
        self.mapq_filter = mapq_filter
        self.get_data()

    def get_data(self):
        """
        initial assign row and edges
        """
        logger.debug("Starting to process hypergraph data.")
        if len(self.edges.row) == 0:
            raise ValueError("No hyperedges found.")
        
        names = np.array(list(self.edges.idx.keys()))
        indices = np.array(list(self.edges.idx.values()))
        self.nodes = names[np.argsort(indices)]


        if self.min_quality > 0 and len(self.edges.mapq) > 1 and self.mapq_filter:
            self.mapq = np.array(self.edges.mapq, dtype=np.int8)

            retain_idx = self.mapq >= self.min_quality
            # total_edge_counts = self.edges.row.shape[0] 

            self.row = self.edges.row[retain_idx]
            self.col = self.edges.col[retain_idx]
        
            self.count = self.edges.count[retain_idx]

            self.idx = np.unique(self.row)
            remove_idx = np.isin(np.arange(0, len(self.nodes)), self.idx)
            self.remove_contigs = self.nodes[~remove_idx]
            # logger.info(f"Removed `{total_edge_counts - np.count_nonzero(retain_idx):,}` "
            #                 f"hyperedges that mapq < {self.min_quality}.")
            logger.info(f"Total {len(self.nodes) - np.sum(remove_idx):,} "
                        f"low-quality (mapq < {self.min_quality}) contigs were removed. ")
            # logger.info(f"Total {total_edge_counts - np.count_nonzero(retain_idx):,} "
            #             f"hyperedges with low-quality (mapq < {self.min_quality}) contigs were removed. ")
        
        else:
            self.mapq = np.array([])
            self.row = self.edges.row
            self.col = self.edges.col
            self.count = self.edges.count
            self.remove_contigs = np.array([])

        self.shape = (len(self.nodes), self.col.max() + 1)
        self.removed_count = 0
        logger.debug("Finished processing hypergraph data.")

    def remove_rows(self, contigs):
        """
        remove rows by contig list
        """     
        contigs = [x for x in contigs if x not in self.edges.idx.values()]
        contig_idx = np.array([self.edges.idx[i] for i in contigs], dtype=int)

        remove_idx = np.isin(self.row, contig_idx)

        self.removed_count += len(contig_idx)

        self.row = self.row[~remove_idx]
        self.col = self.col[~remove_idx]
        self.count = self.count[~remove_idx]

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

        remove_idx = np.isin(self.row, contig_idx)

        self.removed_count += len(contig_idx)

        self.row = self.row[remove_idx]
        self.col = self.col[remove_idx]
        self.count = self.count[remove_idx]
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

        logger.debug("Starting to construct hypergraph incidence matrix.")
        data = self.count if self.count is not None else np.ones(len(self.row), dtype=np.int32)
       
        # matrix = csr_matrix((data, (self.row, self.col)), shape=self.shape)   

        rows = np.asarray(self.row, dtype=np.int32)
        cols = np.asarray(self.col, dtype=np.int32)
        data = np.asarray(data)

        if not np.all(rows[:-1] <= rows[1:]):
            logger.debug("Sorting hypergraph data for incidence matrix construction.")
            order = np.argsort(rows)
            rows = rows[order]
            cols = cols[order]
            data = data[order]
            del order

        row_counts = np.bincount(rows, minlength=self.shape[0])
        del rows
        indptr = np.concatenate(([0], np.cumsum(row_counts)))
        del row_counts

        matrix = csr_matrix((data, cols, indptr), shape=self.shape)

        logger.debug("Constructed raw hypergraph incidence matrix.")

        if min_contacts:
            non_zero_contig_idx = matrix.sum(axis=1).T.A1 >= min_contacts
            matrix = matrix[non_zero_contig_idx]
            # self.remove_contigs = self.nodes[~non_zero_contig_idx]
            self.remove_contig_idx = np.where(~non_zero_contig_idx)[0]
            self.nodes = self.nodes[non_zero_contig_idx]
           
            logger.info(f"Total {self.shape[0] - matrix.shape[0]:,} "
                        f"low-signal (contacts < {min_contacts}) contigs were removed (--min-contacts). ")
            
            non_zero_edges_idx = matrix.sum(axis=0).A1 >= 2
            matrix = matrix[:, non_zero_edges_idx]
            self.shape = matrix.shape 

        else:
            self.remove_contig_idx = np.array([])
            self.remove_contigs = np.array([])
        
        # node_degrees = matrix.sum(axis=1).A1
        # low_signal_mask = node_degrees >= min_contacts
        # max_contacts_percentile = 99.5
        # if max_contacts_percentile < 100:
        #     upper_bound = np.percentile(node_degrees[node_degrees > 0], max_contacts_percentile)
        #     hub_mask = node_degrees <= upper_bound
        #     logger.info(f"Filtering hub contigs with degree > {upper_bound:.2f} ({max_contacts_percentile}th percentile)")
        # else:
        #     hub_mask = np.ones_like(low_signal_mask, dtype=bool)

        # non_zero_contig_idx = low_signal_mask & hub_mask
        
        # matrix = matrix[non_zero_contig_idx]
        # self.remove_contigs = self.nodes[~non_zero_contig_idx]
        # self.remove_contig_idx = np.where(~non_zero_contig_idx)[0]  
        # self.nodes = self.nodes[non_zero_contig_idx]

        # non_zero_edges_idx = matrix.sum(axis=0).A1 >= 2
        # matrix = matrix[:, non_zero_edges_idx]
    
       
        logger.debug("Finished constructing hypergraph incidence matrix.")
        return matrix 

    def clear(self):
        del self.row, self.col, self.count, self.mapq
        gc.collect()

    def to_adjacency_matrix(self, min_weight=0.1):
        """
        Directly construct adjacency matrix for Hi-C data (k=2).
        This is much faster and memory-efficient than clique_expansion_init.
        """
        logger.info("Directly constructing adjacency matrix for Hi-C data...")
  
        order = np.argsort(self.col)
        row_sorted = self.row[order]
        col_sorted = self.col[order]
        count_sorted = self.count[order] if self.count is not None else np.ones(len(row_sorted))

        col_counts = np.bincount(col_sorted)
        valid_cols_mask = col_counts[col_sorted] == 2
        
      
        r_filtered = row_sorted[valid_cols_mask]
        w_filtered = count_sorted[valid_cols_mask]

        u = r_filtered[0::2]
        v = r_filtered[1::2]
        w = w_filtered[0::2]

        rows = np.concatenate([u, v])
        cols = np.concatenate([v, u])
        weights = np.concatenate([w, w])

        n_nodes = len(self.nodes)
        A = csr_matrix((weights, (rows, cols)), shape=(n_nodes, n_nodes))
        A.sum_duplicates() 

        if min_weight > 0:
            A.data[A.data < min_weight] = 0
            A.eliminate_zeros()

        logger.info(f"Constructed adjacency matrix with {A.nnz} non-zero elements.")
        return A


    @staticmethod
    def clique_expansion_init(H,
                              cis_count=None, 
                              P_allelic_idx=None, 
                              P_weak_idx=None, 
                              NW=None,
                              W=None,
                              allelic_factor=-1,
                              min_weight=0.1,
                              use_high_order=True):
        """
        clique expansion/reduction
        """
        logger.debug("Starting clique expansion/reduction of hypergraph.")
        m = H.shape[1]

        if use_high_order:
            # D_e - I
            # D_e_num = (H.sum(axis=0) - 1).astype(HYPERGRAPH_ORDER_DTYPE)
            if cis_count is None:
                D_e_num = (H.sum(axis=0).A1 - 1).astype(HYPERGRAPH_ORDER_DTYPE)
                # D_e_num = H.getnnz(axis=0).astype(HYPERGRAPH_ORDER_DTYPE)
            else:
                D_e_num = (cis_count - H.sum(axis=0).A1 - 1).astype(HYPERGRAPH_ORDER_DTYPE)
        else:
            D_e_num = np.ones(m, dtype=HYPERGRAPH_ORDER_DTYPE)
        
        # W = identity(m, dtype=np.float32)
        # if W is not None:
        #     W = dia_matrix((W, np.array([0])), shape=(m, m), dtype=np.float32)
        # W = dia_matrix((D_e_num, np.array([0])), (m, m), dtype=np.float32)
        
       
        # inverse diagonal matrix D_e
        D_e_inv = dia_matrix(((1 / np.maximum(D_e_num, 1)).astype(np.float32), [0]), shape=(m, m))

        logger.debug("Calculated inverse diagonal matrix of hyperedges.")

        if W is not None:
            A = H @ W @ D_e_inv @ H.T
        else:
            A = H.dot(D_e_inv).dot(H.T)
        logger.debug("Calculated clique expansion/reduction adjacency matrix.")
     

        if min_weight > 0:
            A.data[A.data < min_weight] = 0
            A.eliminate_zeros()

        if cis_count is not None:
            diag_indices = np.arange(A.shape[0])
            A.setdiag(cis_count[diag_indices])

        if NW is not None:
            A = HyperGraph.normalize_adjacency_matrix(A, NW=NW)

        if P_allelic_idx or P_weak_idx:
            if not isinstance(A, csr_matrix):
                A = A.tocsr()
            
            if P_allelic_idx:
                rows, cols = P_allelic_idx
                if allelic_factor == 0:
                    A[rows, cols] = 0
                else:
                    vals = A[rows, cols]
                    if hasattr(vals, "multiply"):
                        A[rows, cols] = vals.multiply(allelic_factor)
                    else:
                        A[rows, cols] = vals * allelic_factor

            if P_weak_idx:
                rows, cols = P_weak_idx
                A[rows, cols] = 0

        logger.debug("Finished clique expansion/reduction of hypergraph.")

        return A

    @staticmethod
    def normalize_adjacency_matrix(A, NW=None):
        if NW is not None:
            logger.info("Normalizing the weight of edges.")
            A_coo = A.tocoo()
            row, col, data = A_coo.row, A_coo.col, A_coo.data

            scales = NW[row, col]
            data = data * scales

            A = csr_matrix((data, (row, col)), shape=A.shape)

        return A
    
    @staticmethod
    def reweight_adjacency_matrix(A, P_allelic_idx=None, P_weak_idx=None,
                                   allelic_factor=-1, cross_allelic_factor=0):
        if P_allelic_idx or P_weak_idx:
            if not isinstance(A, csr_matrix):
                A = A.tocsr()

            if P_allelic_idx:
                rows, cols = P_allelic_idx
                if allelic_factor == 0:
                    A[rows, cols] = 0
                else:
                    vals = A[rows, cols]
                    
                    mean_val = A.data.mean() if A.nnz > 0 else 1.0
                    
                    if allelic_factor == 'power':
                        if hasattr(vals, "toarray"):
                            d = vals.toarray()
                            d = vals * 2.0 + mean_val * 0.5
                            A[rows, cols] = d
                        else:
                            A[rows, cols] = vals * 2.0 + mean_val * 0.5
                    elif allelic_factor == 'maximum':
                        target_value = mean_val * 2
                        if hasattr(vals, "toarray"):

                            d = vals.toarray()
                            d = np.maximum(d * 10, target_value)
                            A[rows, cols] = d
                        else:
                            A[rows, cols] = np.maximum(vals * 10, target_value)
                    elif allelic_factor == 'dynamic':
                        if hasattr(vals, "toarray"):
                            d = vals.toarray()
                            d = np.power(d, 1.5) * 1.5 
                            A[rows, cols] = d
                        else:
                            A[rows, cols] = np.power(vals, 1.2) * 2.0
                    elif  isinstance(allelic_factor, (int, float)) and allelic_factor > 0:
                        if A.nnz > 0:
                            ref_signal = np.percentile(A.data, 90) 
                        else:
                            ref_signal = 1.0
                        if allelic_factor > 0:
                            baseline = ref_signal * 0.5
                            A[rows, cols] = (vals + baseline) * allelic_factor
                    else:
                        if hasattr(vals, "multiply"):
                            A[rows, cols] = vals.multiply(allelic_factor)
                        else:
                            A[rows, cols] = vals * allelic_factor 

            if P_weak_idx:
                rows, cols = P_weak_idx
                vals = A[rows, cols]
                if hasattr(vals, "multiply"):
                    A[rows, cols] = vals.multiply(cross_allelic_factor)
                else:
                    A[rows, cols] = vals * cross_allelic_factor

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

# def extract_incidence_matrix2(mat, idx):
#     A = mat[idx]
#     if A.shape[0] != 0 and A.shape[1] !=0 : 
#         A_sum = A.sum(axis=0).T
#         non_zero_edges_idx = np.where(A_sum > 1)[0]
#         remove_edges_idx = np.where(A_sum <= 1)[0]
#         A = A[:, non_zero_edges_idx]
#     else:
#         remove_edges_idx = np.array([])
#         non_zero_edges_idx = np.array([])

#     return A, remove_edges_idx, non_zero_edges_idx

def extract_incidence_matrix2(mat, idx):
    A = mat[idx]
    if A.nnz > 0: 
        A_sum = np.bincount(A.indices, weights=A.data, minlength=A.shape[1])
        non_zero_mask = A_sum > 1
        non_zero_edges_idx = np.where(non_zero_mask)[0]
        remove_edges_idx = np.where(~non_zero_mask)[0]

        if len(non_zero_edges_idx) < A.shape[1]:
            A = A.tocsc()[:, non_zero_edges_idx].tocsr()
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

def sparsify_topk(A, topk=100):
    """
    For each vertex keep only top-k incident edges by weight.
    """
    logger.info("Sparsifying matrix by keeping top-k edges per vertex.")

    keep_row = []
    keep_col = []
    for i in range(A.shape[0]):
        row = A.getrow(i).tocoo()
      
        if row.nnz > topk:
            sorted_idx = np.argsort(row.data)[-topk:]
            keep_col.extend(row.col[sorted_idx])
            keep_row.extend([i] * topk)
        else:
            keep_col.extend(row.col)
            keep_row.extend([i] * row.nnz)
    
    A[keep_row, keep_col] = 0
    
    A = A.tocsr()
    A.eliminate_zeros()
    logger.info(f"Sparsified matrix to {A.shape[0]} rows and {A.shape[1]} columns, "
                f"keeping top-{topk} edges per vertex.")

    return A
    


def balance_matrix_iteratively(matrix, max_iter=200, tol=1e-5):
    """
    ICE-like matrix balancing.
    A more robust and efficient implementation.
    """
    if not isinstance(matrix, csr_matrix):
        matrix = csr_matrix(matrix)

  
    good_rows = np.array(matrix.sum(axis=1)).flatten() > 0
    good_cols = np.array(matrix.sum(axis=0)).flatten() > 0
    good_nodes = good_rows & good_cols
    
    if not np.all(good_nodes):
        logger.info(f"Temporarily removing {len(good_nodes) - np.sum(good_nodes)} disconnected nodes for balancing.")
        original_size = matrix.shape[0]
        original_indices = np.arange(original_size)
        kept_indices = original_indices[good_nodes]
        
        sub_matrix = matrix[good_nodes, :][:, good_nodes]
    else:
        sub_matrix = matrix
        kept_indices = None

    m = sub_matrix.copy().astype(float)
    
    for i in range(max_iter):
        row_sums = np.array(m.sum(axis=1)).flatten()
        

        row_sums_non_zero = row_sums[row_sums != 0]
        if len(row_sums_non_zero) > 0 and np.all(np.abs(1 - row_sums_non_zero) < tol):
            logger.info(f"Matrix balancing converged after {i} iterations.")
            break

        s = np.ones(m.shape[0])
        s[row_sums != 0] = 1.0 / row_sums[row_sums != 0]
        
        bias_diag = dia_matrix((s, 0), shape=(m.shape[0], m.shape[0])).tocsr()
        m = bias_diag @ m @ bias_diag
    else:
        logger.warning(f"Matrix balancing did not converge after {max_iter} iterations.")

    if kept_indices is not None:
        final_matrix = csr_matrix((original_size, original_size), dtype=float)
        row_idx, col_idx = m.nonzero()
        final_matrix[kept_indices[row_idx], kept_indices[col_idx]] = m.data
        return final_matrix
    else:
        return m

def balance_matrix_vc(A):
    A.setdiag(0)
    A.prune()

    row_sums = np.array(A.sum(axis=1)).flatten()
    col_sums = np.array(A.sum(axis=0)).flatten()

    row_sums[row_sums == 0] = 1
    col_sums[col_sums == 0] = 1
    A_coo = A.tocoo()
    row, col, data = A_coo.row, A_coo.col, A_coo.data

    normalized_data = data / (row_sums[row] * col_sums[col])

    A = csr_matrix((normalized_data, (row, col)), shape=A.shape)

    return A 


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
            method="louvain",
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
        D_e_num = (H.sum(axis=0) - 1).astype(HYPERGRAPH_ORDER_DTYPE)
        # D_e_num = H.getnnz(axis=0).astype(HYPERGRAPH_ORDER_DTYPE)

        ## inverse diagonal matrix D_e
        D_e_inv = dia_matrix((1 / np.maximum(D_e_num, 1), [0]), shape=(m, m))

        if max_round <= 1:
            A = H @ W @ D_e_inv @ H.T

        else:
            H_T = H.T
            A = H @ W @ D_e_inv @ H_T
        
        if min_weight > 0:
            A.data[A.data < min_weight] = 0
            A.eliminate_zeros()

        # normalization
        if NW is not None:
            A = HyperGraph.normalize_adjacency_matrix(A, NW=NW)

        A.setdiag(0)
        A.eliminate_zeros()

        if max_round <= 1:
            del W, D_e_num, D_e_inv
            gc.collect() 

        raw_A = A.copy()

        if P_allelic_idx or P_weak_idx:
            if not isinstance(A, csr_matrix):
                A = A.tocsr()
            
            if P_allelic_idx:
                rows, cols = P_allelic_idx
                A[rows, cols] *= allelic_factor
            if P_weak_idx:
                rows, cols = P_weak_idx
                A[rows, cols] = 0


    else:
        raw_A = A.copy()

        if P_allelic_idx or P_weak_idx:
            if not isinstance(A, csr_matrix):
                A = A.tocsr()
            if P_allelic_idx:
                rows, cols = P_allelic_idx
                A[rows, cols] *= allelic_factor
            if P_weak_idx:
                rows, cols = P_weak_idx
                A[rows, cols] = 0

    A.setdiag(0)
    A_upper = triu(A, format='csr')
    A_upper.eliminate_zeros() 

    random.seed(os.getenv("CPHASING_RANDOM_SEED", 42))
    ig.set_random_number_generator(random)

    try:
        # G = ig.Graph.Weighted_Adjacency(A_upper, mode='undirected', loops=False)
        A_coo = A_upper.tocoo()
        edges = np.column_stack((A_coo.row, A_coo.col))
        G = ig.Graph(
            n=A.shape[0], 
            edges=edges, 
            edge_attrs={'weight': A_coo.data},
            directed=False
        )
    except ValueError:
        return raw_A, A, None, []
    if allelic_factor < 0 or method == 'leiden':
        cluster_assignments = G.community_leiden(weights='weight', 
                                            resolution_parameter=resolution, 
                                            objective_function='modularity')
    else:
        cluster_assignments = G.community_multilevel(weights='weight', resolution=resolution)

    cluster_results = [set(c) for c in cluster_assignments]
  
    cluster_stat = list(map(len, cluster_results))
    current_membership = cluster_assignments.membership
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

        a, args, res = None, None, None
        gc.collect()

        A = H @ W @ D_e_inv @ H_T
        if min_weight > 0:
            A.data[A.data < min_weight] = 0
            A.eliminate_zeros()
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
            if not isinstance(A, csr_matrix):
                A = A.tocsr()
            if P_allelic_idx:
                rows, cols = P_allelic_idx
                A[rows, cols] *= allelic_factor
            if P_weak_idx:
                rows, cols = P_weak_idx
                A[rows, cols] = 0

        A.setdiag(0)
        A_upper = triu(A, format='csr')
        A_upper.eliminate_zeros()

        try:
            # G = ig.Graph.Weighted_Adjacency(A_upper, mode='undirected', loops=False)
            A_coo = A_upper.tocoo()
            edges = np.column_stack((A_coo.row, A_coo.col))
            G = ig.Graph(
                n=A.shape[0], 
                edges=edges, 
                edge_attrs={'weight': A_coo.data},
                directed=False
            )
        except ValueError:
            return A, None, []
        
        # ## use igraph 
        if allelic_factor < 0 or method == 'leiden':
            cluster_assignments = G.community_leiden(weights='weight', 
                                                resolution_parameter=resolution, 
                                                objective_function='modularity',
                                                initial_membership=current_membership)
        else:
            cluster_assignments = G.community_multilevel(weights='weight', resolution=resolution)
        
        
        cluster_results = [set(c) for c in cluster_assignments]
        current_membership = cluster_assignments.membership
        ## use networkx 
        # G = G.to_networkx()
        # cluster_assignments = nx.community.louvain_communities(G, resolution=resolution)

        cluster_stat = list(map(len, cluster_results))
        logger.info(f"Cluster Statistics: {cluster_stat}")
        

        diff_value = max(abs(W.data[0] - W_prev.data[0]))
        W_prev = W.copy()
        
        iter_round += 1
    
    return raw_A, A, cluster_assignments, cluster_results

