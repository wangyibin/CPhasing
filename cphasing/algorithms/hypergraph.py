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

try:
    from sparse_dot_mkl import dot_product_mkl
    HAS_MKL = True
except ImportError:
    HAS_MKL = False

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

def find_csr_data_offsets(A, rows, cols):
    if len(rows) == 0:
        return np.array([], dtype=np.int64)
    coo = A.tocoo()
    dense_cols = A.shape[1]
    
    csr_keys = coo.row.astype(np.int64) * dense_cols + coo.col.astype(np.int64)
    target_keys = rows.astype(np.int64) * dense_cols + cols.astype(np.int64)
    
    sorter = np.argsort(csr_keys)
    insert_idx = np.searchsorted(csr_keys, target_keys, sorter=sorter)
    insert_idx = np.clip(insert_idx, 0, len(csr_keys) - 1)
    matched_indices = sorter[insert_idx]
    
    valid_mask = (csr_keys[matched_indices] == target_keys)

    return matched_indices[valid_mask]


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
        contig_idx_list = [self.edges.idx.get(c) for c in contigs]
        contig_idx_list = [i for i in contig_idx_list if i is not None]
        contig_idx = np.asarray(contig_idx_list, dtype=int)


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
        data = self.count if self.count is not None else np.ones(len(self.row), dtype=np.float32)
        # data = np.ones(len(self.row), dtype=np.int32)
        # matrix = csr_matrix((data, (self.row, self.col)), shape=self.shape)   

        rows = np.asarray(self.row, dtype=np.int32)
        cols = np.asarray(self.col, dtype=HYPERGRAPH_COL_DTYPE)
        data = np.asarray(data, dtype=np.float32)

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

        matrix = csr_matrix((data, cols, indptr), shape=self.shape, dtype=np.float32)

        logger.debug("Constructed raw hypergraph incidence matrix.")

        if min_contacts:
            non_zero_contig_idx = matrix.sum(axis=1).T.A1 >= min_contacts
            # 
            # D_e = matrix.sum(axis=0).A1  
            # term1 = matrix.dot(D_e)
            # matrix_sq = matrix.copy()
            # matrix_sq.data **= 2
            # term2 = matrix_sq.sum(axis=1).A1
            # vpc_degree = term1 - term2      
            # non_zero_contig_idx = vpc_degree >= min_contacts
            # self.remove_contigs = self.nodes[~non_zero_contig_idx]

            matrix = matrix[non_zero_contig_idx]
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
                ## v0.3.0 
                # D_e_num = (H.sum(axis=0).A1 - 1).astype(HYPERGRAPH_ORDER_DTYPE)
                ## v0.3.1
                # k = H.sum(axis=0).A1
                # D_e_num = (k * (k - 1) / 2.0).astype(HYPERGRAPH_ORDER_DTYPE)
                ## v0.3.2
        
                D_e_num = H.getnnz(axis=0).astype(HYPERGRAPH_ORDER_DTYPE) - 1
                ## v0.3.3
                # D_e_num = H.multiply(H).sum(axis=0).A1 

                ## v0.2.0
                # D_e_num = ((H.sum(axis=0).A1 - 1) * (H.sum(axis=0).A1) / 2).astype(HYPERGRAPH_ORDER_DTYPE)
                # D_e_num = H.getnnz(axis=0).astype(HYPERGRAPH_ORDER_DTYPE)
            else:
                D_e_num = (cis_count - H.sum(axis=0).A1 - 1).astype(HYPERGRAPH_ORDER_DTYPE)
        else:
            D_e_num = np.ones(m, dtype=HYPERGRAPH_ORDER_DTYPE)
        
        # W = identity(m, dtype=np.float32)
        # if W is not None:
        #     W = dia_matrix((W, np.array([0])), shape=(m, m), dtype=np.float32)
        # W = dia_matrix((D_e_num, np.array([0])), (m, m), dtype=np.float32)
      
        # edge_weights = np.ones(H.shape[1], dtype=np.float64)
        # edge_sizes = np.array(H.sum(axis=0)).flatten()
        # norm_factor = np.zeros_like(edge_sizes, dtype=np.float64)
        # valid_edges = edge_sizes >= 2

        # norm_factor[valid_edges] = edge_weights[valid_edges] / (edge_sizes[valid_edges] * (edge_sizes[valid_edges] - 1) / 2.0)

        # Zhou's
        # n = H.shape[0]
        # D_v = np.array(H.sum(axis=1)).flatten()
        # D_v_inv_sqrt = dia_matrix(((1.0 / np.sqrt(np.maximum(D_v, 1))).astype(np.float32), [0]), shape=(n, n))

        # inverse diagonal matrix D_e
        # D_e_inv = dia_matrix(((1 / np.maximum(D_e_num, 1)).astype(np.float32), [0]), shape=(m, m))

        logger.debug("Calculated inverse diagonal matrix of hyperedges.")
        d_e_diag = (1.0 / np.maximum(D_e_num, 1)).astype(np.float32)
        H = H.astype(np.float32, copy=False)

        # from scipy.sparse import csc_matrix
        # if isinstance(H, csr_matrix):
        #     H_T = csc_matrix((H.data, H.indices, H.indptr), shape=(H.shape[1], H.shape[0]))
        # else:
        #     H_T = H.T

        if W is not None:
            scale_vector = W.data[0] * d_e_diag
        else:
            scale_vector = d_e_diag
        
        scale_vector_sqrt = np.sqrt(scale_vector)

        if not isinstance(H, csr_matrix):
            H_scaled = H.tocsr()
        else:
            H_scaled = H.copy()
        H_scaled.data *= scale_vector_sqrt[H_scaled.indices]

        if HAS_MKL:
            try:
                A = dot_product_mkl(H_scaled, H_scaled.T)
            except Exception as e:
                A = H_scaled.dot(H_scaled.T)
        else:
            A = H_scaled.dot(H_scaled.T)
            # A = D_v_inv_sqrt @ A @ D_v_inv_sqrt
            # from scipy.sparse import diags
            # A = H.dot(diags(norm_factor)).dot(D_e_inv).dot(H.T)
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
    def _to_series_contigsizes(contigsizes, nodes):
        """
        Normalize contigsizes input to a pandas Series indexed by contig name.
        Accept:
          - dict {contig: size}
          - pd.Series (index=contig)
          - pd.DataFrame with columns like [chrom,length] or [contig,size]
        """
        if contigsizes is None:
            raise ValueError("contigsizes is required")

        if isinstance(contigsizes, dict):
            return pd.Series(contigsizes, name="size")

        if isinstance(contigsizes, pd.Series):
            if contigsizes.name is None:
                contigsizes = contigsizes.rename("size")
            return contigsizes

        if isinstance(contigsizes, pd.DataFrame):
            df = contigsizes.copy()
            # common: ["chrom","length"]
            cols = {c.lower(): c for c in df.columns}
            if "chrom" in cols and "length" in cols:
                df = df.rename(columns={cols["chrom"]: "contig", cols["length"]: "size"})
            elif "contig" in df.columns and "size" in df.columns:
                pass
            else:
                # fallback: first two columns
                df = df.iloc[:, :2].copy()
                df.columns = ["contig", "size"]

            s = pd.Series(df["size"].values, index=df["contig"].astype(str).values, name="size")
            return s

        raise TypeError(f"Unsupported contigsizes type: {type(contigsizes)}")

    
    @staticmethod
    def filter_adjacency_matrix(A, nodes, contigsizes, invert=False,
                                method="quantile",
                                high_q=0.999,
                                low_q=None,
                                mad_z=3.0):
        """
        Filter contigs by abnormal normalized degree (contacts per length).

        Returns
        -------
        retain_idx : np.ndarray[int]
            indices to retain (0..n-1)
        """
        n = A.shape[0]
        if A is None or n == 0 or A.nnz == 0:
            return np.arange(n, dtype=np.int64)

        A = A.tocsr()
        deg = np.asarray(np.abs(A).sum(axis=1)).ravel().astype(np.float64)

        cs = HyperGraph._to_series_contigsizes(contigsizes, nodes)
        sizes = pd.Index(nodes.astype(str)).map(cs).to_numpy(dtype=np.float64)
        sizes[~np.isfinite(sizes)] = np.nan
        sizes[sizes <= 0] = np.nan

        norm = deg / sizes * 10000.0
        x = norm[np.isfinite(norm)]
        if x.size == 0:
            return np.arange(n, dtype=np.int64)

        keep = np.ones(n, dtype=bool)

        if method == "quantile":
            if high_q is not None:
                hi = float(np.nanquantile(x, float(high_q)))
                keep &= ~(np.isfinite(norm) & (norm > hi))
            if low_q is not None:
                lo = float(np.nanquantile(x, float(low_q)))
                keep &= ~(np.isfinite(norm) & (norm < lo))

        elif method == "mad":
            med = float(np.nanmedian(x))
            mad = float(np.nanmedian(np.abs(x - med)))
            if mad <= 0:
                hi = float(np.nanquantile(x, 0.999))
                keep &= ~(np.isfinite(norm) & (norm > hi))
            else:
                robust_z = 0.6745 * (norm - med) / mad
                keep &= ~(np.isfinite(robust_z) & (robust_z > float(mad_z)))

        else:
            raise ValueError("method must be 'quantile' or 'mad'")

        return np.where(keep)[0].astype(np.int64)

    @staticmethod
    def filter_by_weight(A, min_weight=0.1):
        if min_weight > 0:
            A.data[A.data < min_weight] = 0
            A.eliminate_zeros()

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
    def _normalize_pair_index(P):
        """
        Normalize pair indices into (rows, cols) 1D int arrays.

        Accepts:
          - (rows, cols)
          - list/array of shape (n, 2)
        """
        if P is None:
            return None, None

        if isinstance(P, (tuple, list)) and len(P) == 2 and not (
            isinstance(P[0], (tuple, list, np.ndarray)) and len(np.asarray(P[0]).shape) == 2
        ):
            rows = np.asarray(P[0], dtype=np.int64).ravel()
            cols = np.asarray(P[1], dtype=np.int64).ravel()
            return rows, cols


        arr = np.asarray(P)
        if arr.ndim != 2 or arr.shape[1] != 2:
            raise ValueError(f"Invalid pair index format for P_allelic_idx: shape={arr.shape}")
        rows = arr[:, 0].astype(np.int64, copy=False)
        cols = arr[:, 1].astype(np.int64, copy=False)
        return rows, cols

    
    @staticmethod
    def reweight_adjacency_matrix(A, P_allelic_idx=None, P_weak_idx=None,
                                   allelic_factor=-1, cross_allelic_factor=0):
        if P_allelic_idx or P_weak_idx:
            if not isinstance(A, csr_matrix):
                A = A.tocsr()

            if P_allelic_idx:
                rows, cols = HyperGraph._normalize_pair_index(P_allelic_idx)
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
                    elif allelic_factor == "auto":
                        deg = np.asarray(A.sum(axis=1)).ravel().astype(np.float64)
             
                        deg[deg <= 0] = 1.0

                        vals = A[rows, cols]
                        if hasattr(vals, "toarray"):
                            vals = vals.toarray().ravel()
                        else:
                            vals = np.asarray(vals).ravel()

                        # FIX: rows/cols are guaranteed 1D now
                        ## get row1 * col1 row2 * col2 ...
                        norm = np.sqrt(deg[rows] * deg[cols]) * 0.05

                        A[rows, cols] = norm
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
    def to_contacts(A, nodes=None, 
                    NW=None,
                    min_weight=1.0, 
                    contigsizes=None,
                    output="hypergraph.expansion.contacts"):
        A = A.tocoo()

        mask = (A.row <= A.col) & (np.abs(A.data) >= min_weight) & (np.isfinite(A.data))
        if nodes is not None:
            df = pl.DataFrame({
                "0": nodes[A.row[mask]],
                "1": nodes[A.col[mask]],
                "2": A.data[mask]
            })
            if contigsizes is not None:
                contigsizes_series = pl.from_pandas(contigsizes.reset_index().rename(columns={"index": "contigsizes.0", "size": "contigsizes.1"}))

                contigsizes_dict = {row["chrom"]: row["length"] for row in contigsizes_series.iter_rows(named=True)}

                df = df.with_columns([
                    pl.col("0").map_elements(contigsizes_dict.get).alias("contigsizes.0"),
                    pl.col("1").map_elements(contigsizes_dict.get).alias("contigsizes.1")
                ])
                df = df.with_columns([
                    (pl.col("2") / (pl.col("contigsizes.0") * pl.col("contigsizes.1")).sqrt()).alias("2")
                ])
                df = df.select(["0", "1", "2"])
        else:
            df = pl.DataFrame({
                    "0": A.row[mask],
                    "1": A.col[mask],
                    "2": A.data[mask]
            })

       
        if output:
            df.write_csv(output, separator='\t', include_header=False)    
            logger.info(f"Export hypergraph into `{output}`")
            return 
        else:
            return df

    @staticmethod
    def bipartite_leiden_clustering(H, resolution=1.0, min_weight=0.01,
                                    P_allelic_idx=None, P_weak_idx=None,
                                    allelic_factor=-1.0, cross_allelic_factor=0.0,
                                    max_k=15):
        from scipy.sparse import bmat, diags, triu, csr_matrix
        import igraph as ig

        logger.info("Building Bipartite (Star-Expansion) Graph instead of Clique Expansion...")
        n_contigs, n_hyperedges = H.shape

        k = np.array(H.sum(axis=0)).flatten()
        valid_edges = k <= max_k
        H_sub = H[:, valid_edges]
        k_sub = k[valid_edges]

        weight_decay = 1.0 / np.maximum(k_sub - 1, 1).astype(np.float32)
        H_weighted = H_sub @ diags(weight_decay)

        logger.info("Projecting to Contig-Contig adjacency matrix...")
        A = H_weighted @ H_sub.T
        
        if min_weight > 0:
            A.data[A.data < min_weight] = 0
            A.eliminate_zeros()
            
        A.setdiag(0)
        A.eliminate_zeros()

        if P_allelic_idx or P_weak_idx:
            logger.info("Applying penalizations to adjacency matrix...")
            if P_allelic_idx:
                r, c = HyperGraph._normalize_pair_index(P_allelic_idx)
                if allelic_factor == 0:
                    A[r, c] = 0
                else:
                    A[r, c] *= allelic_factor
                    
            if P_weak_idx:
                r, c = HyperGraph._normalize_pair_index(P_weak_idx)
                A[r, c] *= cross_allelic_factor

        A_upper = triu(A, format='coo')
        A_upper.eliminate_zeros()
        
        edges = np.column_stack((A_upper.row, A_upper.col))
        G = ig.Graph(
            n=n_contigs + n_hyperedges,
            edges=edges.tolist(),
            edge_attrs={'weight': A_upper.data},
            directed=False
        )
        
        G.vs['type'] = [0] * n_contigs + [1] * n_hyperedges
        
        logger.info(f"Running Leiden on Bipartite Graph "
                    f"({n_contigs:,} contigs + {n_hyperedges:,} hyperedges)...")
        
        cluster_assignments = G.community_leiden(
            weights='weight', 
            resolution_parameter=resolution, 
            objective_function='modularity'
        )

        contig_memberships = np.array(cluster_assignments.membership[:n_contigs])
        
        cluster_results = []
        for c_id in np.unique(contig_memberships):
            cluster_results.append(set(np.where(contig_memberships == c_id)[0]))
            
        logger.info(f"Bipartite clustering finished: found {len(cluster_results)} groups.")
        
        return contig_memberships, cluster_results


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
    # w_new = (1/1+k).sum() * (e_c+c) / m

    inv_k = 1.0 / (1.0 + k)
    w_new = inv_k.sum() * (e_c + c) / (e_c * 2.0)

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
        if e_c <= 1:
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
    if mat is None:
        raise ValueError("The incidence matrix (H) is None. Please ensure H isprovided or use the adjacency matrix (A) directly.")
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
        # D_e_num = (H.sum(axis=0) - 1).astype(HYPERGRAPH_ORDER_DTYPE)
        D_e_num = H.getnnz(axis=0).astype(HYPERGRAPH_ORDER_DTYPE) - 1
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
            r, c = HyperGraph._normalize_pair_index(P_allelic_idx)
            offsets = find_csr_data_offsets(A, r, c)
            if len(offsets) > 0:
                A.data[offsets] *= allelic_factor
        if P_weak_idx:
            r, c = HyperGraph._normalize_pair_index(P_weak_idx)
            offsets = find_csr_data_offsets(A, r, c)
            if len(offsets) > 0:
                A.data[offsets] = 0

    A.setdiag(0)
    # if A.nnz > 5000:
    #     threshold_val = np.percentile(A.data, 75)
    #     A.data[A.data < threshold_val] = 0.0
    #     A.eliminate_zeros()
    A_upper = triu(A, format='csr')
    A_upper.eliminate_zeros() 
    
    random_state = 42
    random.seed(os.getenv("CPHASING_RANDOM_SEED", random_state))
    try:
        # raise ImportError("sknetwork not available, falling back to igraph for clustering.")
        # from sknetwork.clustering import Louvain, Leiden
        
        # if allelic_factor < 0 or method == 'leiden':
        #     algo = Leiden(resolution=resolution * 0.9, random_state=random_state)
        # else:
        #     algo = Louvain(resolution=resolution * 0.9, random_state=random_state)
            
        # labels = algo.fit_predict(A)
        
        # unique_labels = np.unique(labels)
        # cluster_results = [set(np.where(labels == lbl)[0]) for lbl in unique_labels]
        # cluster_assignments = type('Dummy', (object,), {
        #     'membership': labels.tolist(), 
        #     '__iter__': lambda self: iter(cluster_results)
        # })()
        
        raise ImportError("hyrex not avaliable.")
        import hyrex as hx 
        # cluster_assignments = hx.graph_louvain(A_upper, resolution=resolution, seed=random_state)
        if P_allelic_idx is not None:
            blocked_pairs = list(zip(P_allelic_idx[0].values.tolist(), P_allelic_idx[1].values.tolist()))
        else:
            blocked_pairs = None
        
        if P_weak_idx is not None:
            if blocked_pairs is None:
                blocked_pairs = list(zip(P_weak_idx[0].values.tolist(), P_weak_idx[1].values.tolist()))
            else:
                blocked_pairs.extend(list(zip(P_weak_idx[0].values.tolist(), P_weak_idx[1].values.tolist())))

        # cluster_assignments = hx.cluster(H, resolution=resolution, seed=random_state, blocked_pairs=blocked_pairs)

        logger.info(f"Running hypergraph denoising with blocked pairs: {len(blocked_pairs) if blocked_pairs else 0}")
        H = hx.denoise_hypergraph(H, 
                                  min_support=2, 
                                  use_weights=False, 
                                  max_motif_order=10,
                                  max_size=16,
                                  mode="filter_original", 
                                  evidence="any_pair",
                                  blocked_pair_strategy="split",
                                  blocked_pairs=blocked_pairs)
        logger.info(f"Hypergraph denoising completed. New hypergraph shape: {H.shape}")
        cluster_assignments = hx.louvain(H, resolution=resolution, 
                                         aggregate_small_edges=False,
                                         init='evidence',
                                         seed=random_state)
        logger.info(f"Hypergraph clustering completed. Found {len(cluster_assignments)} clusters.")

    except ImportError:
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
    previous_memberships = None
    damping = 0.5
    while diff_value > threshold:

        if iter_round >= max_round:
            break

        n_contigs, n_hyperedges = H.shape
        c = len(cluster_results)

        C = csr_matrix(
            (np.ones(n_contigs, dtype=np.float32), (np.arange(n_contigs), current_membership)),
            shape=(n_contigs, c)
        )
        K = (H_T @ C).toarray()
        e_c = np.array(H.sum(axis=0)).flatten() 
        w_prev = W.data[0]

        inv_k = 1.0 / (1.0 + K)
        inv_k_sum = inv_k.sum(axis=1)

        denom = np.maximum(e_c * 2.0, 1.0)
        w_new = inv_k_sum * (e_c + c) / denom
        

        max_cluster_counts = K.max(axis=1)
        purity = max_cluster_counts / np.maximum(e_c, 1.0)
        is_highly_pure_but_has_noise = (purity >= 0.75) & (e_c > max_cluster_counts) & (e_c >= 3)
        if np.any(is_highly_pure_but_has_noise):
            logger.debug(f"[IRMM] Suppressing mismatching monomers in {np.sum(is_highly_pure_but_has_noise):,} structured concatemers")

            pure_read_indices = np.where(is_highly_pure_but_has_noise)[0]
            majority_clusters = K[pure_read_indices].argmax(axis=1)

            read_to_maj_cluster = np.full(n_hyperedges, -1, dtype=np.int32)
            read_to_maj_cluster[pure_read_indices] = majority_clusters
            H_reweighted = H.tocoo()
            
            membership_arr = np.asarray(current_membership, dtype=np.int32)
            monomer_memberships = membership_arr[H_reweighted.row]
            expected_maj_clusters = read_to_maj_cluster[H_reweighted.col]

            noise_mask = (expected_maj_clusters >= 0) & (monomer_memberships != expected_maj_clusters)
            
            if np.any(noise_mask):
                H_reweighted.data[noise_mask] = 0.0
                H_reweighted.eliminate_zeros()

            H = H_reweighted.tocsr()
        low_purity_penalty = np.where((e_c >= 4) & (purity < 0.4), 0.1, 1.0)
        w_new = w_new * low_purity_penalty


        w_updated = 0.5 * (w_new + w_prev)
        W.data[0] = np.where(e_c >= 2, w_updated, w_prev)

        gc.collect()

        A_new = H @ W @ D_e_inv @ H_T
        if min_weight > 0:
            A_new.data[A_new.data < min_weight] = 0
            A_new.eliminate_zeros()

        A = A * (1.0 - damping) + A_new * damping
        if NW is not None:
            A = HyperGraph.normalize_adjacency_matrix(A, NW=NW)

        A.setdiag(0)
        A.eliminate_zeros()
        if P_allelic_idx or P_weak_idx:
            if not isinstance(A, csr_matrix):
                A = A.tocsr()
            if P_allelic_idx:
                # rows, cols = P_allelic_idx
                # A[rows, cols] *= allelic_factor
                r, c = HyperGraph._normalize_pair_index(P_allelic_idx)
                offsets = find_csr_data_offsets(A, r, c)
                if len(offsets) > 0:
                    A.data[offsets] *= allelic_factor
            if P_weak_idx:
                # rows, cols = P_weak_idx
                # A[rows, cols] = 0
                r, c = HyperGraph._normalize_pair_index(P_weak_idx)
                offsets = find_csr_data_offsets(A, r, c)
                if len(offsets) > 0:
                    A.data[offsets] = 0

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
        
        if previous_memberships is not None:
            from sklearn.metrics import adjusted_rand_score
            ari = adjusted_rand_score(previous_memberships, current_membership)
            logger.info(f"[IRMM] Round {iter_round} Similarity to previous round (ARI): {ari:.4f}")
            if ari > 0.99: 
                logger.info("[IRMM] Converged early due to stable partitioning.")
                break
        previous_memberships = current_membership.copy()
        diff_value = max(abs(W.data[0] - W_prev.data[0]))
        W_prev = W.copy()
        
        iter_round += 1
    
    return raw_A, A, cluster_assignments, cluster_results

