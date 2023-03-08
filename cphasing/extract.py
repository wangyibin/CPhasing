#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
extract the hyperedges from pore-c table
"""

import logging
import gc
import msgspec
import os
import os.path as op
import sys

import dask.dataframe as dd

from joblib import Parallel, delayed, Memory
from pathlib import Path

from .core import Pairs
from .utilities import listify, list_flatten 
from ._config import *

logger = logging.getLogger(__name__)
memory = Memory('./cachedir', verbose=0)


class Extractor:
    """
    extract edges from pairs file.

    Params:
    --------
    pairs_pathes: list
        list of pairs file
    threads: int
        number of threads
    
    Examples:
    --------
    >>> extractor = Extractor(pairs_pathes)
 
    """
    def __init__(self, pairs_pathes, threads=4):
        self.pairs_pathes = listify(pairs_pathes)
        self.threads = threads 
        self.edges = self.generate_edges()

        
    def generate_edges(self):
        
        if len(self.pairs_pathes) == 1:
            p = pd.read_csv(self.pairs_pathes[0], sep='\t', comment="#",
                                header=None, index_col=None,
                                usecols=[1, 3])
            
    
        
            res = p.data.values.tolist()
            
            
        else: 
            p_list = self.pairs_pathes
            res_list = Parallel(n_jobs=self.threads)(delayed(
                lambda x: pd.read_csv(x, sep='\t', comment="#",
                                header=None, index_col=None,
                                usecols=[1, 3]))(i) for i in p_list)
            
            res_df = pd.concat(res_list, axis=1)
            res = res_df.values.tolist()

        return res

    def save(self, output):
        with open(output, 'wb') as out:
            out.write(msgspec.msgpack.encode(self.edges))
        logger.info(f"Successful output edges into `{output}`")


def process_pore_c_table(df, npartition, 
                            min_order, max_order, 
                            min_alignments, use_dask=False):
    if use_dask:
        df = df.partitions[npartition]
        df = df.compute()
    df['chrom_idx'] = df['chrom'].cat.codes
    chrom_db = dict(zip(range(len(df.chrom.cat.categories)), 
                        df.chrom.cat.categories))
    df = (df.assign(alignment_length=lambda x: x.end - x.start)
            .query(f"alignment_length >= {min_alignments}")
            .set_index('read_idx')
    )
    df_grouped = df.groupby('read_idx')['chrom_idx']
    df_grouped_nunique = df_grouped.nunique()
    df = df.loc[(df_grouped_nunique >= min_order) 
                    & (df_grouped_nunique <= max_order)]
    edges = df.groupby('read_idx')['chrom_idx'].unique()
    edges = edges.apply(lambda x: list(map(chrom_db.get, x)))

    return edges

process_pore_c_table = memory.cache(process_pore_c_table)

class HyperExtractor:
    """
    Params:
        --------
    pore_c_table_pathes: list
        pore_c table, at least have four columns: read_id x, chrom, start, end.
    min_order: int, default 2
        minimum contig order of pore-c reads
    max_order: int, default 20
        maximum contig order of pore-c reads
    threads: int, default 10
        number of threads
    use_dask: bool, default False
        use dask to parse pore-c table
    """
    def __init__(self, pore_c_table_pathes, 
                            min_order=2, 
                            max_order=15, 
                            min_alignments=500,
                            threads=4,
                            use_dask=False):
        
        self.pore_c_table_pathes = listify(pore_c_table_pathes)
        self.min_order = min_order
        self.max_order = max_order
        self.min_alignments = min_alignments
        self.threads = threads
        self.use_dask = use_dask 

        self.pore_c_tables = self.import_pore_c_table()
        self.edges = self.generate_edges()

        
    def import_pore_c_table(self):
        """
        import pore-c table from pore
        """
        logger.info("Loading Pore-C table ...")
        if len(self.pore_c_table_pathes) == 1:
            if Path(self.pore_c_table_pathes[0]).is_symlink():
                infile = os.readlink(self.pore_c_table_pathes[0])
            else:
                infile = self.pore_c_table_pathes[0]
            
            if self.use_dask:
                df = dd.read_parquet(infile, 
                                    columns=['read_idx', 'chrom',
                                            'start', 'end', 
                                            'pass_filter'],
                                    engine=PQ_ENGINE)
            else:
                df = pd.read_parquet(infile, 
                                        columns=['read_idx', 'chrom',
                                                'start', 'end', 
                                                'pass_filter'],
                                    engine=PQ_ENGINE)
            
            df = df.query("pass_filter == True")
            df_list = [df]
            
        else:
            infiles = []
            for i in self.pore_c_table_pathes:
                if Path(i).is_symlink():
                    infiles.append(os.readlink(i))
                else:
                    infiles.append(i)

            if self.use_dask:
                df_list = list(map(lambda x: dd.read_parquet(
                                    x, columns=['read_idx', 'chrom', 
                                                'start', 'end', 
                                                'pass_filter'], 
                                    engine='pyarrow'), infiles))

                df = dd.concat(df_list, axis=0)
                # df['read_idx'] = df['read_name'].cat.as_known().cat.codes
                df = df.query("pass_filter == True").drop("pass_filter", axis=1)
                df_list = [df]
            else:
                df_list = list(map(lambda x: pd.read_parquet(
                                x, columns=['read_idx', 'chrom', 
                                            'start', 'end', 
                                            'pass_filter'], 
                                    engine='pyarrow',),
                                    #filters=[('pass_filter', '=', True)]),
                                    infiles))
                df_list = Parallel(n_jobs=self.threads)(delayed(
                                    lambda x: x.query("pass_filter == True")
                                                .drop("pass_filter", axis=1)
                                                )(i) for i in df_list)

        return df_list
    

    def generate_edges(self):
        """
        generate hypergraph incidence matrix

       

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
        
        if self.use_dask:
            pore_c_table = self.pore_c_tables[0]
            args = []
            for i in range(pore_c_table.npartitions):
                args.append((pore_c_table, i,
                            self.min_order, self.max_order, 
                            self.min_alignments, self.use_dask))
        else:
            args = []
            for i, pore_c_table in enumerate(self.pore_c_tables):
                args.append((pore_c_table, i, 
                            self.min_order, self.max_order, 
                            self.min_alignments, self.use_dask))
        
        res = Parallel(n_jobs=self.threads)(
                        delayed(process_pore_c_table)(i, j, k, l, m, n) 
                                for i, j, k, l, m, n in args)
        
        res_df = pd.concat(res)
        edges = res_df.values.tolist()
        
        logger.info(f"Only retained Pore-C concatemer that: \n"
                        f"\talignment length >= {self.min_alignments}\n"
                        f"\t{self.min_order} <= contig order <= {self.max_order}")
                        
        del res, res_df, pore_c_table
        gc.collect()

        

        return edges
    
    def save(self, output):
        with open(output, 'wb') as out:
            out.write(msgspec.msgpack.encode(self.edges))

        logger.info(f"Successful output edges into `{output}`")