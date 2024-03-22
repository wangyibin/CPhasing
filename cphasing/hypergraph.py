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


from joblib import Parallel, delayed
from pathlib import Path
from pandarallel import pandarallel 
from pytools import natsorted

from .algorithms.hypergraph import HyperEdges
from .utilities import listify, list_flatten 
from ._config import *

logger = logging.getLogger(__name__)

class Extractor:
    """
    extract edges from pairs file.

    Params:
    --------
    pairs_pathes: list
        list of pairs file
    contig_idx: dict
        dictionary of contig idx
    contigsizes: dict
        dictionary of contig sizes
    threads: int
        number of threads
    
    Examples:
    --------
    >>> extractor = Extractor(pairs_pathes, contig_idx, contigsizes)
 
    """
    def __init__(self, pairs_pathes, contig_idx, contigsizes, threads=4):
        self.pairs_pathes = listify(pairs_pathes)
        self.contig_idx = contig_idx
        self.contigsizes = contigsizes

        self.threads = threads 
        self.edges = self.generate_edges()

    @staticmethod
    def _process_df(df, contig_idx, threads=1):
        pandarallel.initialize(nb_workers=threads, verbose=0)
        df['chrom1'] = df['chrom1'].parallel_map(contig_idx.get)
        df['chrom2'] = df['chrom2'].parallel_map(contig_idx.get)
        df = df.dropna(subset=['chrom1', 'chrom2'], axis=0, how="any")

        return df

    def generate_edges(self):
        """
        """
        logger.info(f"Extract edges from pairs.") 

        if len(self.pairs_pathes) == 1:
            # if self.pairs_pathes[0][-3:] == ".gz":
            #     compression = 'gzip'
            # else:
            #     compression='infer'
            p = pd.read_csv(self.pairs_pathes[0], sep='\t', comment="#",
                                header=None, index_col=None,
                                usecols=[1, 3], names=['chrom1', 'chrom2'], 
                               )
           
            res = Extractor._process_df(p, self.contig_idx, self.threads)    
            res = res.reset_index()
  
            res = pd.concat([res[['chrom1', 'index']].rename(
                                        columns={'chrom1': 'row', 'index': 'col'}),
                              res[['chrom2', 'index']].rename(
                                        columns={'chrom2': 'row', 'index': 'col'})], 
                              axis=1)


        else: 
            p_list = self.pairs_pathes
            # if p_list[0][-3:] == ".gz":
            #     compression = 'gzip'
            # else:
            #     compression='infer'
            
            threads_2 = self.threads // len(p_list) + 1
            threads_1 = int(self.threads / threads_2)
            if threads_1 == 0:
                threads_1 = 1
            
            res = Parallel(n_jobs=min(self.threads, len(p_list)))(delayed(
                lambda x: pd.read_csv(x, sep='\t', comment="#",
                                header=None, index_col=None, 
                                # compression=compression,
                                usecols=[1, 3], names=['chrom1', 'chrom2']))
                                (i) for i in p_list)
            
            args = [ (i, self.contig_idx, threads_2) for i in res ]
        

            res = Parallel(n_jobs=threads_1)(delayed(
                                Extractor._process_df)(i, j, k) for i, j, k in args)
            
            res = pd.concat(res, axis=1).reset_index()
            res = pd.concat([res[['chrom1', 'index']].rename(
                                        columns={'chrom1': 'row', 'index': 'col'}),
                              res[['chrom2', 'index']].rename(
                                        columns={'chrom2': 'row', 'index': 'col'})], 
                              axis=1)
        
        number_of_contigs = len(self.contig_idx)
        length = res['row'].shape[0] * res['row'].shape[1]
        logger.info(f"Result of {length} raw "
                    f"hyperedges of {number_of_contigs} contigs. "
                    "Note: it's not the final statistics for hypergraph.")
        return HyperEdges(idx=self.contig_idx, 
                            row=res['row'].values.flatten().tolist(), 
                            col=res['col'].values.flatten().tolist(),
                            contigsizes=self.contigsizes,
                            mapq=[])

    def save(self, output):
        with open(output, 'wb') as out:
            out.write(msgspec.msgpack.encode(self.edges))
        
        logger.info(f"Successful output graph into `{output}`")
    

def process_pore_c_table(df, contig_idx, threads=1,
                            min_order=2, max_order=50, 
                            min_alignments=50, is_parquet=False):
   
    pandarallel.initialize(nb_workers=threads, verbose=0)
    df['chrom_idx'] = df['chrom'].parallel_map(contig_idx.get).astype('int')
    df = df.dropna(axis=0)
    # chrom_db = dict(zip(range(len(df.chrom.cat.categories)), 
    #                     df.chrom.cat.categories))
    # df = (df.assign(alignment_length=lambda x: x.end - x.start)
    #         .query(f"alignment_length >= {min_alignments}")
    #         .set_index('read_idx'))
    df = df.set_index('read_idx')
    df = df[['chrom_idx', 'mapping_quality']]
    df_grouped = df.groupby('read_idx')['chrom_idx']
    df_grouped_nunique = df_grouped.nunique()
    df = df.loc[(df_grouped_nunique >= min_order) 
                    & (df_grouped_nunique <= max_order)]
    
    df = df[['chrom_idx', 'mapping_quality']].reset_index().drop_duplicates(['read_idx', 'chrom_idx'])
  
    df['read_idx'] = df['read_idx'].astype('category')
    df['read_idx'] = df['read_idx'].cat.codes

    
    return df


class HyperExtractor:
    """
    Params:
        --------
    pore_c_table_pathes: list
        pore_c table, at least have four columns: read_id x, chrom, start, end.
    contig_idx: 
    min_order: int, default 2
        minimum contig order of pore-c reads
    max_order: int, default 50
        maximum contig order of pore-c reads
    min_alignments: int, default 30
        minimum length of alignments
    threads: int, default 10
        number of threads

    """
    HEADER = ["read_idx", "read_length", 
              "read_start", "read_end",  
              "strand", "chrom", "start",
              "end", "mapping_quality", "identity", 
              "filter_reason"]
    def __init__(self, pore_c_table_pathes, 
                            contig_idx,
                            contigsizes,
                            min_order=2, 
                            max_order=50, 
                            min_alignments=30,
                            min_quality=1,
                            threads=4,
                            is_parquet=False):
        
        self.pore_c_table_pathes = listify(pore_c_table_pathes)
        self.contig_idx = contig_idx
        self.contigsizes = contigsizes
        self.min_order = min_order
        self.max_order = max_order
        self.min_alignments = min_alignments
        self.min_quality = min_quality
        self.threads = threads
        self.is_parquet = is_parquet

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
            
            if self.is_parquet:
                df = pd.read_parquet(infile, 
                                        columns=['read_idx', 'chrom',
                                                # 'start', 'end', 
                                                'mapping_quality',
                                                #'filter_reason',
                                                ],
                                    engine=PQ_ENGINE)
            else: 
                try:
                    df = pd.read_csv(infile, 
                                    sep='\t',
                                    header=None,
                                    index_col=None,
                                        usecols=['read_idx', 'chrom',
                                                # 'start', 'end', 
                                                'mapping_quality',
                                                #'filter_reason'
                                                ],
                                        names=self.HEADER)
                except pd.errors.ParserError:
                    logger.error("The input contacts are not the pore-c table format"
                                "do you want to parse Pairs file? please add parameters of `--pairs`")
                    sys.exit(-1)
            # df = df.query("filter_reason == 'pass'")
            if self.min_quality > 1:
                df = df.query(f"mapping_quality >= {self.min_quality}")
            df_list = [df]
        

        else:
            infiles = []
            for i in self.pore_c_table_pathes:
                if Path(i).is_symlink():
                    infiles.append(os.readlink(i))
                else:
                    infiles.append(i)
            
            if self.is_parquet:
                df_list = list(map(lambda x: pd.read_parquet(
                                x, columns=['read_idx', 'chrom', 
                                            # 'start', 'end', 
                                            'mapping_quality',
                                            #'filter_reason'
                                            ], 
                                    engine='pyarrow',),
                                    infiles))
                
                # df_list = Parallel(n_jobs=self.threads)(delayed(
                #                     lambda x: x.query("filter_reason == 'pass'")
                #                                 .drop("filter_reason", axis=1)
                #                                 )(i) for i in df_list)

                
            else:
                df_list = list(map(lambda x: pd.read_csv(
                                x, names=self.HEADER, 
                                    usecols=['read_idx', 'chrom', 
                                            # 'start', 'end', 
                                            'mapping_quality',
                                            #'filter_reason'
                                            ], 
                                    sep='\t',
                                    index_col=None,
                                    header=None,
                                    ),
                                    #filters=[('pass_filter', '=', True)]),
                                    infiles))
                
                if self.min_quality > 1:
                    df_list = Parallel(n_jobs=self.threads)(delayed(
                                        lambda x: x.query(#f"filter_reason == 'pass' &"
                                                          f"'min_quality >= {self.min_quality}")
                                                    # .drop("filter_reason", axis=1)
                                                    )(i) for i in df_list)
                    df = df.query(f"mapping_quality >= {self.min_quality}")
                # else:
                #     df_list = Parallel(n_jobs=self.threads)(delayed(
                #                         lambda x: x.query("filter_reason == 'pass'")
                #                                     .drop("filter_reason", axis=1)
                #                                     )(i) for i in df_list)
            

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
        logger.info(f"Only retained Pore-C concatemer that: \n"
                        # f"\talignment length >= {self.min_alignments} and \n"
                        # f"\t{self.min_order} <= contig order <= {self.max_order} and  \n"
                        f"\tmapping_quality >= {self.min_quality}")

        threads_2 = self.threads // len(self.pore_c_tables) + 1
        threads_1 = int(self.threads / threads_2)
        if threads_1 == 0:
            threads_1 = 1
        args = []
        for i, pore_c_table in enumerate(self.pore_c_tables):
            args.append((pore_c_table, self.contig_idx, threads_2, 
                        self.min_order, self.max_order, 
                        self.min_alignments, self.is_parquet))
        
        res = Parallel(n_jobs=threads_1)(
                        delayed(process_pore_c_table)(i, j, k, l, m, n, o) 
                                for i, j, k, l, m, n, o in args)
        
        idx = 0
        mapping_quality_res = []
        for i, df in enumerate(res):
            mapping_quality_res.append(df['mapping_quality'])
            df['read_idx'] = df['read_idx'] + idx 
            res[i] = df
            idx += len(df)
        
        res_df = pd.concat(res)
        mapping_quality_res = pd.concat(mapping_quality_res)
        
        edges = HyperEdges(idx=self.contig_idx, 
                       row=res_df['chrom_idx'].values.flatten().tolist(),
                       col=res_df['read_idx'].values.flatten().tolist(),
                       mapq=mapping_quality_res.values.flatten().tolist(),
                       contigsizes=self.contigsizes)
        
        number_of_contigs = len(self.contig_idx)
        number_of_hyperedges = res_df['read_idx'].max()
        logger.info(f"Result of {number_of_hyperedges} raw "
                    f"hyperedges of {number_of_contigs} contigs. "
                    "Note: it's not the final statistics for hypergraph.")

        return edges
    
    def save(self, output):
        with open(output, 'wb') as out:
            out.write(msgspec.msgpack.encode(self.edges))

        logger.info(f"Successful output hypergraph into `{output}`")


class HyperExtractorSplit:
    HEADER = ["read_idx", "read_length", 
              "read_start", "read_end",  
              "strand", "chrom", "start",
              "end", "mapping_quality", "identity", 
              "filter_reason"]
    def __init__(self, pore_c_table_pathes, contig_idx, contigsizes, 
                 split=5,
                 min_quality=1,
                 threads=4, 
                is_parquet=False, ):
        self.pore_c_table_pathes = listify(pore_c_table_pathes)
        self.contig_idx = contig_idx
        self.contigsizes = contigsizes 
        self.min_quality = min_quality
        self.threads = threads 
        self.is_parquet = is_parquet

        self.split = split
        self.split_contig_idx = {}
        idx = 0
        for contig in self.contig_idx:
            for i in range(split):
                self.split_contig_idx[f"{contig}_{i}"] = idx
                idx += 1 

        self.pore_c_tables = self.import_pore_c_table()
        self.edges = self.generate_edges()

    @staticmethod
    def process_pore_c_table(df, contig_idx, contig_sizes, 
                             split_contig_idx,
                             split=5, threads=4, 
                             ):
        pandarallel.initialize(nb_workers=threads, verbose=0)
        df = df.set_index('chrom')
        df = df.loc[list(contig_idx.keys())]
        df = df.reset_index().set_index('read_idx')
        df['contigsize'] = df['chrom'].parallel_map(contig_sizes.get)

        df['pos'] = ((df['start'] + df['end']) // 2) // (df['contigsize'] // split)
        df['pos'] = df['pos'].astype(str)
        
        df['chrom'] = df['chrom'].str.cat(df['pos'], sep="_")
        df['chrom_idx'] = df['chrom'].parallel_map(split_contig_idx.get)
        df = df[['chrom_idx', 'mapping_quality', 'chrom']]
        df_grouped = df.groupby('read_idx')['chrom_idx']
        df_grouped_nunique = df_grouped.nunique()
        # df = df.loc[(df_grouped_nunique >= min_order) 
        #             & (df_grouped_nunique <= max_order)]
    
        df = df[['chrom_idx', 'mapping_quality', 'chrom']].reset_index().drop_duplicates(['read_idx', 'chrom_idx'])
        df['read_idx'] = df['read_idx'].astype('category')
        df['read_idx'] = df['read_idx'].cat.codes

        return df         


    def import_pore_c_table(self):
        if len(self.pore_c_table_pathes) == 1:
            if Path(self.pore_c_table_pathes[0]).is_symlink():
                infile = os.readlink(self.pore_c_table_pathes[0])
            else:
                infile = self.pore_c_table_pathes[0]
        
            if self.is_parquet:
                df = pd.read_parquet(infile, 
                                        columns=['read_idx', 'chrom',
                                                'start', 'end', 
                                                'mapping_quality',
                                                #'filter_reason',
                                                ],
                                    engine=PQ_ENGINE) 
            else:
                df = pd.read_csv(infile, 
                                sep='\t',
                                header=None,
                                index_col=None,
                                    usecols=['read_idx', 'chrom',
                                            'start', 'end', 
                                            'mapping_quality',
                                            #'filter_reason'
                                            ],
                                    names=self.HEADER)
                
                if self.min_quality > 1:
                    df = df.query(f"mapping_quality >= {self.min_quality}")
                df_list = [df]

        

        return df_list
    
    def generate_edges(self):
        logger.info("Processing Pore-C table ...")
    
        threads_2 = self.threads // len(self.pore_c_tables) + 1
        threads_1 = int(self.threads / threads_2)
        if threads_1 == 0:
            threads_1 = 1
        args = []
        for i, pore_c_table in enumerate(self.pore_c_tables):
            args.append((pore_c_table, self.contig_idx, self.contigsizes, self.split_contig_idx,
                         self.split, threads_2))
        
        res = Parallel(n_jobs=threads_1)(
            delayed(HyperExtractorSplit.process_pore_c_table)(i, j, k, l, m, n)
                for i, j, k, l, m, n in args
        )

        idx = 0
        mapping_quality_res = []
        for i, df in enumerate(res):
            mapping_quality_res.append(df['mapping_quality'])
            df['read_idx'] = df['read_idx'] + idx 
            res[i] = df
            idx += len(df)
        
        res_df = pd.concat(res)
        mapping_quality_res = pd.concat(mapping_quality_res)


        edges = HyperEdges(idx=self.split_contig_idx, 
                       row=res_df['chrom_idx'].values.flatten().tolist(),
                       col=res_df['read_idx'].values.flatten().tolist(),
                       mapq=mapping_quality_res.values.flatten().tolist(),
                       contigsizes=self.contigsizes)
        

        return edges 

    def save(self, output):
        with open(output, 'wb') as out:
            out.write(msgspec.msgpack.encode(self.edges))

        logger.info(f"Successful output hypergraph into `{output}`")        