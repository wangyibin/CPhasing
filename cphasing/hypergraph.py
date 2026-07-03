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
import shutil
import subprocess

import numpy as np
import pandas as pd
import polars as pl
pl.enable_string_cache()
from joblib import Parallel, delayed
from pathlib import Path
from pandarallel import pandarallel 
from pytools import natsorted
from subprocess import Popen, PIPE

from .pqs import PQS
from .algorithms.hypergraph import HyperEdges, numpy_enc_hook
from .utilities import (listify, 
                        list_flatten, 
                        is_compressed_table_empty, 
                        decompress_cmd)
from ._config import *

logger = logging.getLogger(__name__)


if pd.__version__.split(".")[0] == 1:
    pandas_version = 1
elif pd.__version__.split(".")[0] == 2:
    pandas_version = 2
else:
    pandas_version = 1

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
        
    def __init__(self, pairs_pathes, contig_idx, contigsizes, 
                 min_quality=1, hcr_bed=None, hcr_invert=False,
                 threads=4, edge_length=2e6, split_length=None, 
                 split_contig_boundarys=None, max_q0_ratio=0.0,
                 low_memory=True, log_dir="logs"):
        self.pairs_pathes = listify(pairs_pathes)
        self.contig_idx = contig_idx
        self.contigsizes = contigsizes
        
        max_length = max(self.contigsizes.values())
        if max_length < 2**32 / 2:
            self.pos_dtype = pl.UInt32
        else:
            self.pos_dtype = pl.UInt64
    

        self.min_mapq = min_quality
        logger.debug(f"Minimum mapping quality: {self.min_mapq}")
        self.hcr_bed = hcr_bed
        self.hcr_invert = hcr_invert
        self.threads = threads 

        os.environ["POLARS_MAX_THREADS"] = str(threads)
        
        self.edge_length = edge_length
        self.split_length = split_length
        self.split_contig_boundarys = split_contig_boundarys
        self.max_q0_ratio = max_q0_ratio
        self.low_memory = low_memory

        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self.edges = self.generate_edges()

    @staticmethod
    def _process_df(df, contig_idx, threads=1):
        pandarallel.initialize(nb_workers=min(10, threads), verbose=0)
        df['chrom1'] = df['chrom1'].parallel_map(contig_idx.get)
        df['chrom2'] = df['chrom2'].parallel_map(contig_idx.get)
        df = df.dropna(subset=['chrom1', 'chrom2'], axis=0, how="any")

        return df
    @staticmethod
    def _process_df_pl(df, contig_idx, threads=1):

        mapping_df = pl.DataFrame(
            {
                "chrom": list(contig_idx.keys()),
                "idx": list(contig_idx.values()),
            }
        ).with_columns(
            pl.col("chrom").cast(pl.Categorical)
        )

        df = (
            df
            .join(mapping_df.rename({"chrom": "chrom1"}), on="chrom1", how="inner")
            .drop("chrom1")
            .rename({"idx": "chrom1"})
            .join(mapping_df.rename({"chrom": "chrom2"}), on="chrom2", how="inner")
            .drop("chrom2")
            .rename({"idx": "chrom2"})
        )

        return df
    
    def generate_edges(self):
        """
        """
        logger.info(f"Extract edges from pairs.") 
        if self.low_memory:
            # dtype={'chrom1': 'category', 'chrom2': 'category', 'mapq': 'int8'}
            dtype = {'chrom1': pl.Categorical, 'chrom2': pl.Categorical, 'mapq': pl.Int8}
        else:
            dtype={'mapq': pl.Int8}
        
        if self.edge_length or self.split_length:
            columns = ['chrom1', 'pos1', 'chrom2', 'pos2', 'mapq']
            usecols = [1, 2, 3, 4, 7]
        else:
            columns = ['chrom1', 'chrom2', 'mapq']
            usecols = [1, 3, 7]
        
        if len(self.pairs_pathes) == 1:
            pairs_prefix = Path(Path(self.pairs_pathes[0]).name)
            while pairs_prefix.suffix in {'.pairs', '.gz', '.pqs'}:
                pairs_prefix = pairs_prefix.with_suffix('')
                    

            if Path(self.pairs_pathes[0]).is_dir():
                p = PQS(path=self.pairs_pathes[0], threads=self.threads)
                
                p.init_read()
                if not p.is_pairs():
                    logger.error(f"The input `{self.pairs_pathes[0]}` is not a pairs file, please check it.")
                    sys.exit(-1)
                

                if self.hcr_bed:
                    logger.info(f"Filtering pairs by {self.hcr_bed} ...")
                    cmd = ["cphasing-rs", "pairs-intersect", 
                           self.pairs_pathes[0], self.hcr_bed, 
                           "-q", str(self.min_mapq),
                           "--max-q0-ratio", str(self.max_q0_ratio),
                           "-t", str(self.threads)]

                    if self.hcr_invert:
                        cmd.append("--invert")

                    cmd.extend(["-o", f"{pairs_prefix}.intersect.pqs"])
                    cmd.append(f"2>{self.log_dir}/{pairs_prefix}.pairs.intersect.log")


                    flag = os.system(" ".join(cmd))
                    assert flag == 0, "Failed to execute command, please check log."

                    p = PQS(path=f"{pairs_prefix}.intersect.pqs", threads=self.threads)
                    p.init_read()
                
                chunks = p.read(min_mapq=self.min_mapq, return_as='files')

                res = p.to_hg_df(chunks, self.contig_idx, self.min_mapq, 
                                 edge_length=self.edge_length, 
                                 split_length=self.split_length,
                                    split_contig_boundarys=self.split_contig_boundarys,
                                    max_q0_ratio=self.max_q0_ratio)
                
                
                if Path(f"{pairs_prefix}.intersect.pqs").exists():
                    shutil.rmtree(f"{pairs_prefix}.intersect.pqs")  
                
                res = res.with_row_count("col")
                df1 = res.select([
                    pl.col("chrom1").alias("row"),
                    pl.col("col"),
                    pl.col("mapq")
                ])
            
                df2 = res.select([
                    pl.col("chrom2").alias("row"),
                    pl.col("col"),
                    pl.col("mapq")
                ])

                res = pl.concat([df1, df2]).sort("row")


            else:
                if is_compressed_table_empty(self.pairs_pathes[0]):
                    logger.error(f"The pairs `{self.pairs_pathes[0]}` is empty, can not load anything, please check it.")
                    sys.exit(-1)    

                input_file = ""
                if not self.hcr_bed:
                    input_file = self.pairs_pathes[0]
                else:
                    
                    if str(self.pairs_pathes[0]).endswith(".gz"):
                        cmd0 = decompress_cmd(str(self.pairs_pathes[0]), str(self.threads))
                        cmd = (f"{' '.join(cmd0)} 2>{self.log_dir}/{pairs_prefix}.decompress.hcr.log | "
                                f"cphasing-rs pairs-intersect - {self.hcr_bed} -q {self.min_mapq} -t {self.threads} -o temp.{pairs_prefix}.hcr.pairs.gz")
                    else:
                        cmd = f"cphasing-rs pairs-intersect {self.pairs_pathes[0]} {self.hcr_bed} -q {self.min_mapq} -t {self.threads} -o temp.{pairs_prefix}.hcr.pairs.gz"

                    if self.hcr_invert:
                        cmd += " --invert"
                    cmd += f" 2>{self.log_dir}/{pairs_prefix}.pairs.intersect.log"

                    logger.info(f"Generating hcr pairs by {self.hcr_bed} ...")
                    flag = os.system(cmd)
                    assert flag == 0, "Failed to execute command, please check log."
                    input_file = f"temp.{pairs_prefix}.hcr.pairs.gz"      

                p = pd.read_csv(self.pairs_pathes[0], sep='\t', comment="#", 
                                    header=None, index_col=None, nrows=1)
                if len(p.columns) >= 8  and isinstance(p[7].values[0], np.int64) and p[7].values[0] <= 60:
                        
                    p = pl.read_csv(input_file, separator='\t', has_header=False,
                                    comment_prefix="#", columns=usecols,
                                    new_columns=columns,
                                    dtypes=dtype)
                    if self.hcr_bed:
                        if Path(f"temp.{pairs_prefix}.hcr.pairs.gz").exists():
                            os.remove(f"temp.{pairs_prefix}.hcr.pairs.gz")
                    

                    if self.min_mapq > 0:
                        p = p.filter(pl.col('mapq') >= self.min_mapq)

                    
                    if self.edge_length:
                        edge_length = self.edge_length
                        mapping_df = pl.DataFrame(
                            {
                                "chrom": list(self.contigsizes.keys()),
                                "length": list(self.contigsizes.values()),
                            }
                        ).with_columns(
                            pl.col("chrom").cast(pl.Categorical)
                        )

                        p = (
                            p
                            .join(mapping_df.rename({"chrom": "chrom1"}), on="chrom1", how="inner")
                            .rename({"length": "length1"})
                            .join(mapping_df.rename({"chrom": "chrom2"}), on="chrom2", how="inner")
                            .rename({"length": "length2"})
                            .filter(
                                ((pl.col("pos1") < edge_length)
                                 | (pl.col("pos1") > (pl.col("length1") - edge_length)))
                                & ((pl.col("pos2") < edge_length)
                                   | (pl.col("pos2") > (pl.col("length2") - edge_length)))
                            )
                            .select(columns)
                        )

                    if self.split_length and self.split_contig_boundarys:
                        split_map = pl.DataFrame(
                            {
                                "chrom": list(self.split_contig_boundarys.keys()),
                                "init_idx": list(self.split_contig_boundarys.values()),
                            }
                        ).with_columns(
                            pl.col("chrom").cast(pl.Categorical)
                        )
                       
                        p = (
                            p.with_columns(
                                (pl.col("pos1") // self.split_length).alias("sub_idx1"),
                                (pl.col("pos2") // self.split_length).alias("sub_idx2"),
                            )
                            .join(split_map.rename({"chrom": "chrom1"}), on="chrom1", how="inner")
                            .rename({"init_idx": "init_idx1"})
                            .join(split_map.rename({"chrom": "chrom2"}), on="chrom2", how="inner")
                            .rename({"init_idx": "init_idx2"})
                        )

                        p = (
                            p.with_columns(
                                (pl.col("init_idx1") + pl.col("sub_idx1")).alias("chrom1"),
                                (pl.col("init_idx2") + pl.col("sub_idx2")).alias("chrom2"),
                            )
                            .drop(['init_idx1', 'init_idx2', 'sub_idx1', 'sub_idx2'])
                            .drop_nulls(subset=['chrom1', 'chrom2'])
                            .select(['chrom1', 'chrom2', 'mapq'])
                        )
                        res = p
                    else:
                        res = Extractor._process_df_pl(p, self.contig_idx, self.threads)   
                        
                    # res = res.reset_index(drop=True).reset_index()
            
                    # res = pd.concat([res[['chrom1', 'index', 'mapq']].rename(
                    #                             columns={'chrom1': 'row', 'index': 'col'}),
                    #                 res[['chrom2', 'index', 'mapq']].rename(
                    #                             columns={'chrom2': 'row', 'index': 'col'})], 
                    #                 axis=0)
                    res_list = [res]
                    full_pairs = pl.concat([r for r in res_list if r is not None and len(r) > 0])
                    del res_list
                    full_pairs = full_pairs.with_row_count("col")

                    df1 = full_pairs.select([
                        pl.col("chrom1").alias("row"),
                        pl.col("col"),
                        (pl.col("mapq") if "mapq" in full_pairs.columns else pl.lit(0, dtype=pl.Int8)).alias("mapq")
                    ])
                    df2 = full_pairs.select([
                        pl.col("chrom2").alias("row"),
                        pl.col("col"),
                        (pl.col("mapq") if "mapq" in full_pairs.columns else pl.lit(0, dtype=pl.Int8)).alias("mapq")
                    ])

                    res = pl.concat([df1, df2]).sort("row")

                    del full_pairs, df1, df2
                                
                else:
                    columns.remove('mapq')
                    usecols.remove(7)
                    p = pl.read_csv(input_file, separator='\t', has_header=False,
                                    comment_prefix="#", columns=usecols,
                                    new_columns=columns,
                                    dtypes=dtype)

                    if self.hcr_bed:
                        if Path(f"temp.{pairs_prefix}.hcr.pairs").exists():
                            os.remove(f"temp.{pairs_prefix}.hcr.pairs")
                    
                    if self.edge_length:
                        edge_length = self.edge_length
                        p = (
                            p.with_columns(
                                [
                                    pl.col("chrom1")
                                    .map_elements(
                                        self.contigsizes.get, skip_nulls=False
                                    )
                                    .alias("length1"),
                                    pl.col("chrom2")
                                    .map_elements(
                                        self.contigsizes.get, skip_nulls=False
                                    )
                                    .alias("length2"),
                                ]
                            )
                            .filter(
                                ((pl.col("pos1") < edge_length)
                                    | (pl.col("pos1") > (pl.col("length1") - edge_length)))
                                & ((pl.col("pos2") < edge_length)
                                    | (pl.col("pos2") > (pl.col("length2") - edge_length)))
                            )
                            .select(columns)
                        )

                    if self.split_length and self.split_contig_boundarys:
                        split_map = pl.DataFrame(
                            {
                                "chrom": list(self.split_contig_boundarys.keys()),
                                "init_idx": list(self.split_contig_boundarys.values()),
                            }
                        ).with_columns(
                            pl.col("chrom").cast(pl.Categorical)
                        )

                        p = (
                            p.with_columns(
                                (pl.col("pos1") // self.split_length).alias("sub_idx1"),
                                (pl.col("pos2") // self.split_length).alias("sub_idx2"),
                            )
                            .join(split_map.rename({"chrom": "chrom1"}), on="chrom1", how="inner")
                            .rename({"init_idx": "init_idx1"})
                            .join(split_map.rename({"chrom": "chrom2"}), on="chrom2", how="inner")
                            .rename({"init_idx": "init_idx2"})
                        )

                        p = (
                            p.with_columns(
                                (pl.col("init_idx1") + pl.col("sub_idx1")).alias("chrom1"),
                                (pl.col("init_idx2") + pl.col("sub_idx2")).alias("chrom2"),
                            )
                            .drop(['init_idx1', 'init_idx2', 'sub_idx1', 'sub_idx2'])
                            .drop_nulls(subset=['chrom1', 'chrom2'])
                            .select(['chrom1', 'chrom2', 'mapq'])
                        )
                        res = p
                    else:
                        res = Extractor._process_df_pl(p, self.contig_idx, self.threads)   

                    res_list = [res]
                    full_pairs = pl.concat([r for r in res_list if r is not None and len(r) > 0])
                    del res_list
                    full_pairs = full_pairs.with_row_count("col")
                    df1 = full_pairs.select([
                        pl.col("chrom1").alias("row"),
                        pl.col("col"),
                    ])
                    df2 = full_pairs.select([
                        pl.col("chrom2").alias("row"),
                        pl.col("col"),
                    ])
                    res = pl.concat([df1, df2]).sort("row")
                    del full_pairs, df1, df2


        else: 
            p_list = self.pairs_pathes
            
            threads_2 = self.threads // len(p_list) + 1
            threads_1 = int(self.threads / threads_2)
            if threads_1 == 0:
                threads_1 = 1

            if not self.hcr_bed:
                def get_file(i):
                    if is_compressed_table_empty(i):
                        logger.error(f"The pairs `{i}` is empty, can not load anything, please check it.")
                        return None
                    return i
            else:
                def get_file(i):
                    cmd = f"cphasing-rs pairs-intersect {i} {self.hcr_bed} -q {self.min_mapq} -t {self.threads} -o -"
                    if self.hcr_invert:
                        cmd += " --invert"
                    process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
                    stdout, stderr = process.communicate()

                    return stdout
                


            # dtype={'chrom1': 'category', 'chrom2': 'category', 'mapq': 'int8'}
            p = pd.read_csv(self.pairs_pathes[0], sep='\t', comment="#", 
                                header=None, index_col=None, nrows=1)
            if len(p.columns) >= 8  and isinstance(p[7].values[0], np.int64) and p[7].values[0] <= 60:
                def read_csv(x):
                    df = pl.read_csv(x, separator='\t', has_header=False, 
                                     comment_prefix="#", columns=usecols,
                                    new_columns=columns,
                                    dtypes=dtype)
                    if self.min_mapq > 0:
                        df = df.filter(pl.col('mapq') >= self.min_mapq)

                    return df
                with Parallel(backend="loky", n_jobs=min(self.threads, len(p_list))) as parallel:   
                    res = parallel(delayed(
                                lambda x: read_csv(get_file(x)))(i) for i in p_list)

            else:
                usecols.remove(7)
                columns.remove('mapq')
                def read_csv(x):
                    df = pl.read_csv(x, separator='\t', has_header=False, 
                                     comment_prefix="#", columns=usecols,
                                    new_columns=columns,
                                    dtypes=dtype)
    
                    return df
                
                res = Parallel(n_jobs=min(self.threads, len(p_list)))(delayed(
                                lambda x: read_csv(get_file(x)))(i) for i in p_list)

            if self.edge_length:
                mapping_df = pl.DataFrame(
                    {
                        "chrom": list(self.contigsizes.keys()),
                        "length": list(self.contigsizes.values()),
                    }
                ).with_columns(
                    pl.col("chrom").cast(pl.Categorical)
                )
            if self.split_length and self.split_contig_boundarys:
                split_map = pl.DataFrame(
                        {
                            "chrom": list(self.split_contig_boundarys.keys()),
                            "init_idx": list(self.split_contig_boundarys.values()),
                        }
                    ).with_columns(
                        pl.col("chrom").cast(pl.Categorical)
                    )

            for i, p in enumerate(res):
                if self.edge_length:
                    edge_length = self.edge_length
                    p = (
                        p
                        .join(mapping_df.rename({"chrom": "chrom1"}), on="chrom1", how="inner")
                        .rename({"length": "length1"})
                        .join(mapping_df.rename({"chrom": "chrom2"}), on="chrom2", how="inner")
                        .rename({"length": "length2"})
                        .filter(
                            ((pl.col("pos1") < edge_length)
                             | (pl.col("pos1") > (pl.col("length1") - edge_length)))
                            & ((pl.col("pos2") < edge_length)
                               | (pl.col("pos2") > (pl.col("length2") - edge_length)))
                        )
                        .select(columns)
                    )
                
                if self.split_length and self.split_contig_boundarys:
                    p = (
                        p.with_columns(
                            (pl.col("pos1") // self.split_length).alias("sub_idx1"),
                            (pl.col("pos2") // self.split_length).alias("sub_idx2"),
                        )
                        .join(split_map.rename({"chrom": "chrom1"}), on="chrom1", how="inner")
                        .rename({"init_idx": "init_idx1"})
                        .join(split_map.rename({"chrom": "chrom2"}), on="chrom2", how="inner")
                        .rename({"init_idx": "init_idx2"})
                    )

                    p = (
                        p.with_columns(
                            (pl.col("init_idx1") + pl.col("sub_idx1")).alias("chrom1"),
                            (pl.col("init_idx2") + pl.col("sub_idx2")).alias("chrom2"),
                        )
                        .drop(['init_idx1', 'init_idx2', 'sub_idx1', 'sub_idx2'])
                        .drop_nulls(subset=['chrom1', 'chrom2'])
                        .select(['chrom1', 'chrom2', 'mapq'])
                    )
                    res[i] = p.to_pandas()        
                
            if not self.split_length or not self.split_contig_boundarys:
                args = [ (i, self.contig_idx, threads_2) for i in res ]
                with Parallel(backend="loky", n_jobs=threads_1) as parallel:
                    res = parallel(delayed(
                                    Extractor._process_df_pl)(i, j, k) for i, j, k in args)

                
            if len(p.columns) >= 8  and isinstance(p[7].values[0], np.int64) and p[7].values[0] <= 60:
                res_list = [r for r in res if r is not None and len(r) > 0]
                full_pairs = pl.concat(res_list)
                del res_list
                full_pairs = full_pairs.with_row_count("col")
                df1 = full_pairs.select([
                    pl.col("chrom1").alias("row"),
                    pl.col("col"),
                    (pl.col("mapq") if "mapq" in full_pairs.columns else pl.lit(0, dtype=pl.Int8)).alias("mapq")
                ])
                df2 = full_pairs.select([
                    pl.col("chrom2").alias("row"),
                    pl.col("col"),
                    (pl.col("mapq") if "mapq" in full_pairs.columns else pl.lit(0, dtype=pl.Int8)).alias("mapq")
                ])
                res = pl.concat([df1, df2]).sort("row")
                del full_pairs, df1, df2

            else:
                res_list = [r for r in res if r is not None and len(r) > 0]
                full_pairs = pl.concat(res_list)
                del res_list
                full_pairs = full_pairs.with_row_count("col")
                df1 = full_pairs.select([
                    pl.col("chrom1").alias("row"),
                    pl.col("col"),
                ])
                df2 = full_pairs.select([
                    pl.col("chrom2").alias("row"),
                    pl.col("col"),
                ])
                res = pl.concat([df1, df2]).sort("row")
                del full_pairs, df1, df2
    
                
            
        number_of_contigs = len(self.contig_idx)
        length = res['col'].max() + 1
        logger.info(f"Result of {length:,} raw "
                    f"edges of {number_of_contigs:,} contigs. "
                    "Note: it's not the final statistics for hypergraph.")
        
        logger.debug("Generating hyperedges ...")
        if 'mapq' in res.columns:
            return HyperEdges(idx=self.contig_idx, 
                            row=res['row'].to_numpy(), 
                            col=res['col'].to_numpy(),
                            count=np.ones(len(res['col']), dtype=np.uint32),
                            contigsizes=self.contigsizes,
                            mapq=res['mapq'].to_numpy())

        else:
            return HyperEdges(idx=self.contig_idx, 
                            row=res['row'].to_numpy(), 
                            col=res['col'].to_numpy(),
                            count=np.ones(len(res['col']), dtype=np.uint32),
                            contigsizes=self.contigsizes,
                            mapq=[])

    def save(self, output):
        with open(output, 'wb') as out:
            out.write(msgspec.msgpack.encode(self.edges, enc_hook=numpy_enc_hook))
        
        logger.info(f"Successful output graph into `{output}`")
    
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
        
    def __init__(self, pairs_pathes, contig_idx, contigsizes, 
                 min_quality=1, hcr_bed=None, hcr_invert=False,
                 threads=4, edge_length=2e6, split_length=None, 
                 split_contig_boundarys=None, max_q0_ratio=0.0,
                 low_memory=True, log_dir="logs"):
        self.pairs_pathes = listify(pairs_pathes)
        self.contig_idx = contig_idx
        self.contigsizes = contigsizes
        
        max_length = max(self.contigsizes.values())
        if max_length < 2**32 / 2:
            self.pos_dtype = pl.UInt32
        else:
            self.pos_dtype = pl.UInt64
    

        self.min_mapq = min_quality
        logger.debug(f"Minimum mapping quality: {self.min_mapq}")
        self.hcr_bed = hcr_bed
        self.hcr_invert = hcr_invert
        self.threads = threads 

        self.mapping_df = pl.DataFrame({
            "chrom": list(self.contigsizes.keys()),
            "length": list(self.contigsizes.values()),
        }).with_columns(pl.col("chrom").cast(pl.Categorical))

        self.contig_idx_df = pl.DataFrame({
            "chrom": list(self.contig_idx.keys()),
            "idx": list(self.contig_idx.values()),
        }).with_columns(pl.col("chrom").cast(pl.Categorical))

        self.split_length = split_length
        self.split_contig_boundarys = split_contig_boundarys
        self.split_to_orig_df = None
        self.split_to_orig_idx_df = None
        
        if self.split_contig_boundarys:
            split_to_orig = {}
            split_to_offset = {}
            boundary_keys = set(self.split_contig_boundarys.keys())
            for split_name, idx in self.contig_idx.items():
                matched = False
                base_name = split_name.split('|')[0] if '|' in split_name else split_name

                offset = 0
                if '|' in split_name:
                    try:
                        offset = int(split_name.split('|')[1].split('_')[0])
                    except (IndexError, ValueError):
                        offset = 0

                if base_name in boundary_keys:
                    split_to_orig[split_name] = base_name
                    split_to_offset[split_name] = offset
                    continue
                
                # for sep in ('_', '-', '.', '~'):
                #     if sep in base_name:
                #         base = base_name.rsplit(sep, 1)[0]
                #         if base in boundary_keys:
                #             split_to_orig[split_name] = base
                #             split_to_offset[split_name] = offset
                #             matched = True
                #             break
                if matched:
                    continue
                
                for chrom in boundary_keys:
                    if base_name.startswith(f"{chrom}_") or base_name.startswith(f"{chrom}-") or base_name.startswith(f"{chrom}~") or base_name.startswith(f"{chrom}."):
                        split_to_orig[split_name] = chrom
                        split_to_offset[split_name] = offset
                        matched = True
                        break
                
                if not matched:
                    split_to_orig[split_name] = base_name
                    split_to_offset[split_name] = offset
            
            self.split_to_orig_df = pl.DataFrame({
                "chrom": list(split_to_orig.keys()),
                "orig_chrom": list(split_to_orig.values()),
                "orig_offset": [split_to_offset[k] for k in split_to_orig.keys()]
            }, schema={"chrom": pl.Categorical, "orig_chrom": pl.Categorical, "orig_offset": pl.Int64})

            split_to_orig_idx = {self.contig_idx[k]: (v, split_to_offset[k]) for k, v in split_to_orig.items() if k in self.contig_idx}
            self.split_to_orig_idx_df = pl.DataFrame({
                "idx": list(split_to_orig_idx.keys()),
                "orig_chrom": [x[0] for x in split_to_orig_idx.values()],
                "orig_offset": [x[1] for x in split_to_orig_idx.values()]
            }, schema={"idx": pl.Int64, "orig_chrom": pl.Categorical, "orig_offset": pl.Int64})
        
        
        self.preloaded_hcr_df = None
        if self.hcr_bed:
            logger.info("Pre-loading and sorting HCR BED for pairs filtering...")
            original_val = os.environ.get("POLARS_MAX_THREADS", "1")
            os.environ["POLARS_MAX_THREADS"] = "1"
            try:
                self.preloaded_hcr_df = pl.read_csv(
                    self.hcr_bed,
                    separator="\t",
                    has_header=False,
                    new_columns=["chrom", "pos", "end_bed"],
                    dtypes={"chrom": pl.Categorical, "pos": pl.UInt32, "end_bed": pl.UInt32}
                ).sort("pos")
            finally:
                os.environ["POLARS_MAX_THREADS"] = original_val

        os.environ["POLARS_MAX_THREADS"] = str(threads)
        
        self.edge_length = edge_length
    
        self.max_q0_ratio = max_q0_ratio
        self.low_memory = low_memory

        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self.edges = self.generate_edges()

    @staticmethod
    def _process_df(df, contig_idx, threads=1):
        pandarallel.initialize(nb_workers=min(10, threads), verbose=0)
        df['chrom1'] = df['chrom1'].parallel_map(contig_idx.get)
        df['chrom2'] = df['chrom2'].parallel_map(contig_idx.get)
        df = df.dropna(subset=['chrom1', 'chrom2'], axis=0, how="any")

        return df

    @staticmethod
    def _process_df_pl(df, contig_idx_df, threads=1):
        df = (
            df
            .join(contig_idx_df.rename({"chrom": "chrom1"}), on="chrom1", how="inner")
            .drop("chrom1")
            .rename({"idx": "chrom1"})
            .join(contig_idx_df.rename({"chrom": "chrom2"}), on="chrom2", how="inner")
            .drop("chrom2")
            .rename({"idx": "chrom2"})
        )

        return df

    @staticmethod
    def _filter_pairs_hcr(df, hcr_df, invert=False, contig_idx_df=None, split_to_orig_df=None, split_to_orig_idx_df=None, keep_pos=False):
        if hcr_df is None or df.height == 0:
            return df
        
        orig_cols = df.columns
        if "pos1" not in df.columns or "pos2" not in df.columns:
            return df

        df = df.with_row_count("__orig_idx")
        
        if df["chrom1"].dtype.is_numeric():
            if split_to_orig_idx_df is not None:
                right_idx1 = split_to_orig_idx_df.rename({"idx": "chrom1", "orig_chrom": "join_chrom1", "orig_offset": "orig_offset1"}).with_columns(
                    pl.col("chrom1").cast(df["chrom1"].dtype)
                )
                right_idx2 = split_to_orig_idx_df.rename({"idx": "chrom2", "orig_chrom": "join_chrom2", "orig_offset": "orig_offset2"}).with_columns(
                    pl.col("chrom2").cast(df["chrom2"].dtype)
                )
                df = (
                    df
                    .join(right_idx1, on="chrom1", how="left")
                    .join(right_idx2, on="chrom2", how="left")
                )
            elif contig_idx_df is not None:
                right_idx1 = contig_idx_df.rename({"idx": "chrom1", "chrom": "join_chrom1"}).with_columns(
                    pl.col("chrom1").cast(df["chrom1"].dtype)
                )
                right_idx2 = contig_idx_df.rename({"idx": "chrom2", "chrom": "join_chrom2"}).with_columns(
                    pl.col("chrom2").cast(df["chrom2"].dtype)
                )
                df = (
                    df
                    .join(right_idx1, on="chrom1", how="left")
                    .join(right_idx2, on="chrom2", how="left")
                )
        else:
            if split_to_orig_df is not None:
                right_df1 = split_to_orig_df.rename({"chrom": "chrom1", "orig_chrom": "join_chrom1", "orig_offset": "orig_offset1"}).with_columns(
                    pl.col("chrom1").cast(df["chrom1"].dtype)
                )
                right_df2 = split_to_orig_df.rename({"chrom": "chrom2", "orig_chrom": "join_chrom2", "orig_offset": "orig_offset2"}).with_columns(
                    pl.col("chrom2").cast(df["chrom2"].dtype)
                )
                df = (
                    df
                    .join(right_df1, on="chrom1", how="left")
                    .join(right_df2, on="chrom2", how="left")
                )
            else:
                df = df.with_columns([
                    pl.col("chrom1").alias("join_chrom1"),
                    pl.col("chrom2").alias("join_chrom2")
                ])

        if "orig_offset1" not in df.columns:
            df = df.with_columns(pl.lit(0, dtype=pl.Int64).alias("orig_offset1"))
        if "orig_offset2" not in df.columns:
            df = df.with_columns(pl.lit(0, dtype=pl.Int64).alias("orig_offset2"))

        df = df.with_columns([
            pl.col("orig_offset1").fill_null(0),
            pl.col("orig_offset2").fill_null(0),
            pl.col("join_chrom1").cast(pl.Categorical),
            pl.col("join_chrom2").cast(pl.Categorical)
        ])

        df = df.with_columns([
            (pl.col("pos1").cast(pl.Int64) + pl.col("orig_offset1")).cast(pl.UInt32).alias("orig_pos1"),
            (pl.col("pos2").cast(pl.Int64) + pl.col("orig_offset2")).cast(pl.UInt32).alias("orig_pos2"),
        ])

        sorted_df = df.sort("orig_pos1")
        joined1 = sorted_df.join_asof(
            hcr_df.rename({"chrom": "join_chrom1", "pos": "orig_pos1", "end_bed": "end_bed1"}),
            on="orig_pos1",
            by="join_chrom1",
            strategy="backward"
        )
        del sorted_df, df
        
        sorted_joined1 = joined1.sort("orig_pos2")
        joined2 = sorted_joined1.join_asof(
            hcr_df.rename({"chrom": "join_chrom2", "pos": "orig_pos2", "end_bed": "end_bed2"}),
            on="orig_pos2",
            by="join_chrom2",
            strategy="backward"
        )
        del sorted_joined1, joined1

        is_overlapping1 = (pl.col("orig_pos1") < pl.col("end_bed1")) & pl.col("end_bed1").is_not_null()
        is_overlapping2 = (pl.col("orig_pos2") < pl.col("end_bed2")) & pl.col("end_bed2").is_not_null()
        
        any_in_hcr = is_overlapping1 | is_overlapping2
        
        if invert:
            filtered = joined2.filter(~any_in_hcr)
        else:
            filtered = joined2.filter(any_in_hcr)

        if keep_pos:
            final_cols = [c for c in orig_cols]
        else:
            final_cols = [c for c in orig_cols if c not in ("pos1", "pos2")]
        res = filtered.sort("__orig_idx").select(final_cols)
        del filtered
        
        return res

    def generate_edges(self):
        """
        """
        logger.info(f"Extract edges from pairs.") 
        if self.low_memory:
            dtype = {'chrom1': pl.Categorical, 'chrom2': pl.Categorical, 'mapq': pl.Int8}
        else:
            dtype={'mapq': pl.Int8}
        
        if self.edge_length or self.split_length or (self.preloaded_hcr_df is not None):
            columns = ['chrom1', 'pos1', 'chrom2', 'pos2', 'mapq']
            usecols = [1, 2, 3, 4, 7]
        else:
            columns = ['chrom1', 'chrom2', 'mapq']
            usecols = [1, 3, 7]
        
        if len(self.pairs_pathes) == 1:
            pairs_prefix = Path(Path(self.pairs_pathes[0]).name)
            while pairs_prefix.suffix in {'.pairs', '.gz', '.pqs'}:
                pairs_prefix = pairs_prefix.with_suffix('')
                    

            if Path(self.pairs_pathes[0]).is_dir():
                p = PQS(path=self.pairs_pathes[0], threads=self.threads)
                
                p.init_read()
                if not p.is_pairs():
                    logger.error(f"The input `{self.pairs_pathes[0]}` is not a pairs file, please check it.")
                    sys.exit(-1)
                
                chunks = p.read(min_mapq=self.min_mapq, return_as='files')

                res = p.to_hg_df(chunks, self.contig_idx, self.min_mapq, 
                                 edge_length=self.edge_length, 
                                 split_length=self.split_length,
                                 split_contig_boundarys=self.split_contig_boundarys,
                                 max_q0_ratio=self.max_q0_ratio)
                

                if self.preloaded_hcr_df is not None:
                    keep_pos = bool(self.edge_length or self.split_length)
                    res = Extractor._filter_pairs_hcr(res, self.preloaded_hcr_df, self.hcr_invert, self.contig_idx_df, split_to_orig_df=self.split_to_orig_df, split_to_orig_idx_df=self.split_to_orig_idx_df, keep_pos=keep_pos)

                res = res.with_row_count("col")
                df1 = res.select([
                    pl.col("chrom1").alias("row"),
                    pl.col("col"),
                    pl.col("mapq")
                ])
            
                df2 = res.select([
                    pl.col("chrom2").alias("row"),
                    pl.col("col"),
                    pl.col("mapq")
                ])

                res = pl.concat([df1, df2]).sort("row")


            else:
                if is_compressed_table_empty(self.pairs_pathes[0]):
                    logger.error(f"The pairs `{self.pairs_pathes[0]}` is empty, can not load anything, please check it.")
                    sys.exit(-1)    

                input_file = self.pairs_pathes[0]

                p = pd.read_csv(self.pairs_pathes[0], sep='\t', comment="#", 
                                    header=None, index_col=None, nrows=1)
                if len(p.columns) >= 8  and isinstance(p[7].values[0], np.int64) and p[7].values[0] <= 60:
                        
                    p = pl.read_csv(input_file, separator='\t', has_header=False,
                                    comment_prefix="#", columns=usecols,
                                    new_columns=columns,
                                    dtypes=dtype)

                    if self.min_mapq > 0:
                        p = p.filter(pl.col('mapq') >= self.min_mapq)

                    if self.preloaded_hcr_df is not None:
                        keep_pos = bool(self.edge_length or self.split_length)
                        p = Extractor._filter_pairs_hcr(p, self.preloaded_hcr_df, self.hcr_invert, self.contig_idx_df, split_to_orig_df=self.split_to_orig_df, split_to_orig_idx_df=self.split_to_orig_idx_df, keep_pos=keep_pos)

                    if self.edge_length:
                        edge_length = self.edge_length
                        p = (
                            p
                            .join(self.mapping_df.rename({"chrom": "chrom1"}), on="chrom1", how="inner")
                            .rename({"length": "length1"})
                            .join(self.mapping_df.rename({"chrom": "chrom2"}), on="chrom2", how="inner")
                            .rename({"length": "length2"})
                            .filter(
                                ((pl.col("pos1") < edge_length)
                                 | (pl.col("pos1") > (pl.col("length1") - edge_length)))
                                & ((pl.col("pos2") < edge_length)
                                   | (pl.col("pos2") > (pl.col("length2") - edge_length)))
                            )
                            .select(columns)
                        )

                    if self.split_length and self.split_contig_boundarys:
                        split_map = pl.DataFrame(
                            {
                                "chrom": list(self.split_contig_boundarys.keys()),
                                "init_idx": list(self.split_contig_boundarys.values()),
                            }
                        ).with_columns(
                            pl.col("chrom").cast(pl.Categorical)
                        )
                       
                        p = (
                            p.with_columns(
                                (pl.col("pos1") // self.split_length).alias("sub_idx1"),
                                (pl.col("pos2") // self.split_length).alias("sub_idx2"),
                            )
                            .join(split_map.rename({"chrom": "chrom1"}), on="chrom1", how="inner")
                            .rename({"init_idx": "init_idx1"})
                            .join(split_map.rename({"chrom": "chrom2"}), on="chrom2", how="inner")
                            .rename({"init_idx": "init_idx2"})
                        )

                        p = (
                            p.with_columns(
                                (pl.col("init_idx1") + pl.col("sub_idx1")).alias("chrom1"),
                                (pl.col("init_idx2") + pl.col("sub_idx2")).alias("chrom2"),
                            )
                            .drop(['init_idx1', 'init_idx2', 'sub_idx1', 'sub_idx2'])
                            .drop_nulls(subset=['chrom1', 'chrom2'])
                            .select(['chrom1', 'chrom2', 'mapq'])
                        )
                        res = p
                    else:
                        res = Extractor._process_df_pl(p, self.contig_idx_df, self.threads)   
                        
                    res_list = [res]
                    full_pairs = pl.concat([r for r in res_list if r is not None and len(r) > 0])
                    del res_list
                    full_pairs = full_pairs.with_row_count("col")

                    df1 = full_pairs.select([
                        pl.col("chrom1").alias("row"),
                        pl.col("col"),
                        (pl.col("mapq") if "mapq" in full_pairs.columns else pl.lit(0, dtype=pl.Int8)).alias("mapq")
                    ])
                    df2 = full_pairs.select([
                        pl.col("chrom2").alias("row"),
                        pl.col("col"),
                        (pl.col("mapq") if "mapq" in full_pairs.columns else pl.lit(0, dtype=pl.Int8)).alias("mapq")
                    ])

                    res = pl.concat([df1, df2]).sort("row")

                    del full_pairs, df1, df2
                                
                else:
                    columns.remove('mapq')
                    usecols.remove(7)
                    p = pl.read_csv(input_file, separator='\t', has_header=False,
                                    comment_prefix="#", columns=usecols,
                                    new_columns=columns,
                                    dtypes=dtype)

                    if self.preloaded_hcr_df is not None:
                        keep_pos = bool(self.edge_length or self.split_length)
                        p = Extractor._filter_pairs_hcr(p, self.preloaded_hcr_df, self.hcr_invert, self.contig_idx_df, split_to_orig_df=self.split_to_orig_df, split_to_orig_idx_df=self.split_to_orig_idx_df, keep_pos=keep_pos)
                        
                          
                    if self.edge_length:
                        edge_length = self.edge_length
                        p = (
                            p
                            .join(self.mapping_df.rename({"chrom": "chrom1"}), on="chrom1", how="inner")
                            .rename({"length": "length1"})
                            .join(self.mapping_df.rename({"chrom": "chrom2"}), on="chrom2", how="inner")
                            .rename({"length": "length2"})
                            .filter(
                                ((pl.col("pos1") < edge_length)
                                    | (pl.col("pos1") > (pl.col("length1") - edge_length)))
                                & ((pl.col("pos2") < edge_length)
                                    | (pl.col("pos2") > (pl.col("length2") - edge_length)))
                            )
                            .select(columns)
                        )

                    if self.split_length and self.split_contig_boundarys:
                        split_map = pl.DataFrame(
                            {
                                "chrom": list(self.split_contig_boundarys.keys()),
                                "init_idx": list(self.split_contig_boundarys.values()),
                            }
                        ).with_columns(
                            pl.col("chrom").cast(pl.Categorical)
                        )

                        p = (
                            p.with_columns(
                                (pl.col("pos1") // self.split_length).alias("sub_idx1"),
                                (pl.col("pos2") // self.split_length).alias("sub_idx2"),
                            )
                            .join(split_map.rename({"chrom": "chrom1"}), on="chrom1", how="inner")
                            .rename({"init_idx": "init_idx1"})
                            .join(split_map.rename({"chrom": "chrom2"}), on="chrom2", how="inner")
                            .rename({"init_idx": "init_idx2"})
                        )

                        p = (
                            p.with_columns(
                                (pl.col("init_idx1") + pl.col("sub_idx1")).alias("chrom1"),
                                (pl.col("init_idx2") + pl.col("sub_idx2")).alias("chrom2"),
                            )
                            .drop(['init_idx1', 'init_idx2', 'sub_idx1', 'sub_idx2'])
                            .drop_nulls(subset=['chrom1', 'chrom2'])
                            .select(['chrom1', 'chrom2', 'mapq'])
                        )
                        res = p
                    else:
                        res = Extractor._process_df_pl(p, self.contig_idx_df, self.threads)   

                    res_list = [res]
                    full_pairs = pl.concat([r for r in res_list if r is not None and len(r) > 0])
                    del res_list
                    full_pairs = full_pairs.with_row_count("col")
                    df1 = full_pairs.select([
                        pl.col("chrom1").alias("row"),
                        pl.col("col"),
                    ])
                    df2 = full_pairs.select([
                        pl.col("chrom2").alias("row"),
                        pl.col("col"),
                    ])
                    res = pl.concat([df1, df2]).sort("row")
                    del full_pairs, df1, df2


        else: 
            p_list = self.pairs_pathes
            
            threads_2 = self.threads // len(p_list) + 1
            threads_1 = int(self.threads / threads_2)
            if threads_1 == 0:
                threads_1 = 1

            def get_file(i):
                if is_compressed_table_empty(i):
                    logger.error(f"The pairs `{i}` is empty, can not load anything, please check it.")
                    return None
                return i

            p = pd.read_csv(self.pairs_pathes[0], sep='\t', comment="#", 
                                header=None, index_col=None, nrows=1)
            
            if len(p.columns) >= 8  and isinstance(p[7].values[0], np.int64) and p[7].values[0] <= 60:
                def read_csv(x):
                    df = pl.read_csv(x, separator='\t', has_header=False, 
                                     comment_prefix="#", columns=usecols,
                                    new_columns=columns,
                                    dtypes=dtype)
                    if self.min_mapq > 0:
                        df = df.filter(pl.col('mapq') >= self.min_mapq)

                    if self.preloaded_hcr_df is not None:
                        keep_pos = bool(self.edge_length or self.split_length)
                        df = Extractor._filter_pairs_hcr(df, self.preloaded_hcr_df, self.hcr_invert, self.contig_idx_df, split_to_orig_df=self.split_to_orig_df, split_to_orig_idx_df=self.split_to_orig_idx_df, keep_pos=keep_pos)
                        

                    return df
                with Parallel(backend="threading", n_jobs=min(self.threads, len(p_list))) as parallel:   
                    res = parallel(delayed(
                                lambda x: read_csv(get_file(x)))(i) for i in p_list)

            else:
                usecols.remove(7)
                columns.remove('mapq')
                def read_csv(x):
                    df = pl.read_csv(x, separator='\t', has_header=False, 
                                     comment_prefix="#", columns=usecols,
                                    new_columns=columns,
                                    dtypes=dtype)
                    
                    if self.preloaded_hcr_df is not None:
                        keep_pos = bool(self.edge_length or self.split_length)
                        df = Extractor._filter_pairs_hcr(df, self.preloaded_hcr_df, self.hcr_invert, self.contig_idx_df, split_to_orig_df=self.split_to_orig_df, split_to_orig_idx_df=self.split_to_orig_idx_df, keep_pos=keep_pos)
    
                    return df
                
                with Parallel(backend="threading", n_jobs=min(self.threads, len(p_list))) as parallel:
                    res = parallel(delayed(
                                lambda x: read_csv(get_file(x)))(i) for i in p_list)

            if self.split_length and self.split_contig_boundarys:
                split_map = pl.DataFrame(
                        {
                            "chrom": list(self.split_contig_boundarys.keys()),
                            "init_idx": list(self.split_contig_boundarys.values()),
                        }
                    ).with_columns(
                        pl.col("chrom").cast(pl.Categorical)
                    )

            for i, p in enumerate(res):
                if p is None:
                    continue
                if self.edge_length:
                    edge_length = self.edge_length
                    p = (
                        p
                        .join(self.mapping_df.rename({"chrom": "chrom1"}), on="chrom1", how="inner")
                        .rename({"length": "length1"})
                        .join(self.mapping_df.rename({"chrom": "chrom2"}), on="chrom2", how="inner")
                        .rename({"length": "length2"})
                        .filter(
                            ((pl.col("pos1") < edge_length)
                             | (pl.col("pos1") > (pl.col("length1") - edge_length)))
                            & ((pl.col("pos2") < edge_length)
                               | (pl.col("pos2") > (pl.col("length2") - edge_length)))
                        )
                        .select(columns)
                    )
                
                if self.split_length and self.split_contig_boundarys:
                    p = (
                        p.with_columns(
                            (pl.col("pos1") // self.split_length).alias("sub_idx1"),
                            (pl.col("pos2") // self.split_length).alias("sub_idx2"),
                        )
                        .join(split_map.rename({"chrom": "chrom1"}), on="chrom1", how="inner")
                        .rename({"init_idx": "init_idx1"})
                        .join(split_map.rename({"chrom": "chrom2"}), on="chrom2", how="inner")
                        .rename({"init_idx": "init_idx2"})
                    )

                    p = (
                        p.with_columns(
                            (pl.col("init_idx1") + pl.col("sub_idx1")).alias("chrom1"),
                            (pl.col("init_idx2") + pl.col("sub_idx2")).alias("chrom2"),
                        )
                        .drop(['init_idx1', 'init_idx2', 'sub_idx1', 'sub_idx2'])
                        .drop_nulls(subset=['chrom1', 'chrom2'])
                        .select(['chrom1', 'chrom2', 'mapq'])
                    )
                    res[i] = p      
                
            if not self.split_length or not self.split_contig_boundarys:
                args = [ (i, self.contig_idx_df, threads_2) for i in res if i is not None ]
                with Parallel(backend="threading", n_jobs=threads_1) as parallel:
                    res = parallel(delayed(
                                    Extractor._process_df_pl)(i, j, k) for i, j, k in args)

                
            if len(p.columns) >= 8  and isinstance(p[7].values[0], np.int64) and p[7].values[0] <= 60:
                res_list = [r for r in res if r is not None and len(r) > 0]
                full_pairs = pl.concat(res_list)
                del res_list
                full_pairs = full_pairs.with_row_count("col")
                df1 = full_pairs.select([
                    pl.col("chrom1").alias("row"),
                    pl.col("col"),
                    (pl.col("mapq") if "mapq" in full_pairs.columns else pl.lit(0, dtype=pl.Int8)).alias("mapq")
                ])
                df2 = full_pairs.select([
                    pl.col("chrom2").alias("row"),
                    pl.col("col"),
                    (pl.col("mapq") if "mapq" in full_pairs.columns else pl.lit(0, dtype=pl.Int8)).alias("mapq")
                ])
                res = pl.concat([df1, df2]).sort("row")
                del full_pairs, df1, df2

            else:
                res_list = [r for r in res if r is not None and len(r) > 0]
                full_pairs = pl.concat(res_list)
                del res_list
                full_pairs = full_pairs.with_row_count("col")
                df1 = full_pairs.select([
                    pl.col("chrom1").alias("row"),
                    pl.col("col"),
                ])
                df2 = full_pairs.select([
                    pl.col("chrom2").alias("row"),
                    pl.col("col"),
                ])
                res = pl.concat([df1, df2]).sort("row")
                del full_pairs, df1, df2
        
        number_of_contigs = len(self.contig_idx)
        max_col = res['col'].max()
        length = max_col + 1 if max_col is not None else 0
        
        logger.info(f"Result of {length:,} raw "
                    f"edges of {number_of_contigs:,} contigs. "
                    "Note: it's not the final statistics for hypergraph.")
        
        logger.debug("Generating hyperedges ...")
        if 'mapq' in res.columns:
            return HyperEdges(idx=self.contig_idx, 
                            row=res['row'].to_numpy(), 
                            col=res['col'].to_numpy(),
                            count=np.ones(len(res['col']), dtype=np.uint32),
                            contigsizes=self.contigsizes,
                            mapq=res['mapq'].to_numpy())

        else:
            return HyperEdges(idx=self.contig_idx, 
                            row=res['row'].to_numpy(), 
                            col=res['col'].to_numpy(),
                            count=np.ones(len(res['col']), dtype=np.uint32),
                            contigsizes=self.contigsizes,
                            mapq=[])

    def save(self, output):
        with open(output, 'wb') as out:
            out.write(msgspec.msgpack.encode(self.edges, enc_hook=numpy_enc_hook))
        
        logger.info(f"Successful output graph into `{output}`")


def process_pore_c_table_slow(df, contig_idx, contigsizes, threads=1,
                            min_order=2, max_order=10, 
                            min_alignments=50, is_parquet=False,
                            edge_length=2e6, split_length=None,
                            split_contig_boundarys=None,
                          ):
    
    if isinstance(df, pl.DataFrame):
        df = df.lazy()

    if split_length and split_contig_boundarys is not None:
        mapping_df = pl.DataFrame({
            "chrom": list(split_contig_boundarys.keys()),
            "init_idx": list(split_contig_boundarys.values()),
        }).with_columns(
            pl.col("chrom").cast(pl.Categorical)
        ).lazy()

        df = (
            df
            .join(mapping_df, on="chrom", how="inner")
            .with_columns([
                ((pl.col('start') + (pl.col('end') - pl.col('start')) // 2) // split_length).alias('sub_idx')
            ])
            .with_columns([
                (pl.col('init_idx') + pl.col('sub_idx')).cast(pl.UInt32).alias('chrom_idx')
            ])
        )
       
    else:
        mapping_df = pl.DataFrame({
            "chrom": list(contig_idx.keys()),
            "chrom_idx": list(contig_idx.values()),
        }).with_columns(
            pl.col("chrom").cast(pl.Categorical)
        ).lazy()
        
        df = df.join(mapping_df, on="chrom", how="inner")
        df = df.with_columns(pl.col('chrom_idx').cast(pl.UInt32))
        
        if edge_length:
            mapping_df_len = pl.DataFrame({
                "chrom": list(contigsizes.keys()),
                "length": list(contigsizes.values()),
            }).with_columns(
                pl.col("chrom").cast(pl.Categorical)
            ).lazy()
            
            df = (
                df
                .join(mapping_df_len, on="chrom", how="inner")
                .with_columns([
                    pl.col('length').cast(pl.UInt32),
                    (pl.col('start') + (pl.col('end') - pl.col('start')) // 2).alias('pos')
                ])
                .filter(
                    (pl.col('pos') < edge_length) | ((pl.col('length') - pl.col('pos')) < edge_length)
                )
            )
        
    df = df.select(['read_idx', 'chrom_idx', 'mapping_quality'])

    df = df.with_columns([
        pl.col('chrom_idx').n_unique().over('read_idx').alias('chrom_idx_nunique')
    ])
    
    df = df.filter(
        (pl.col('chrom_idx_nunique') >= min_order) & 
        (pl.col('chrom_idx_nunique') <= max_order)
    )


    df = (
        df
        .group_by(['read_idx', 'chrom_idx'])
        .agg([
            pl.col('chrom_idx').count().alias('chrom_idx_count'),
            pl.col('mapping_quality').max() 
        ])
    )


    df = df.select(['read_idx', 'chrom_idx', 'mapping_quality', 'chrom_idx_count'])

    return df.collect()

def process_pore_c_table(df, contig_idx, contigsizes, threads=1, 
                            min_order=2, max_order=10, 
                            min_alignments=50, is_parquet=False,
                            edge_length=2e6, split_length=None,
                            split_contig_boundarys=None,
                            mapping_df=None, mapping_df_len=None, split_map=None
                          ):
    
    if isinstance(df, pl.DataFrame):
        df = df.lazy()

    if split_length and split_contig_boundarys is not None:
        if split_map is None:
            split_map = pl.DataFrame({
                "chrom": list(split_contig_boundarys.keys()),
                "init_idx": list(split_contig_boundarys.values()),
            }).with_columns(pl.col("chrom").cast(pl.Categorical))
        
        s_map = split_map.lazy() if isinstance(split_map, pl.DataFrame) else split_map

        df = (
            df
            .join(s_map, on="chrom", how="inner")
            .with_columns([
                ((pl.col('start') + (pl.col('end') - pl.col('start')) // 2) // split_length).alias('sub_idx')
            ])
            .with_columns([
                (pl.col('init_idx') + pl.col('sub_idx')).cast(pl.UInt32).alias('chrom_idx')
            ])
        )
       
    else:
        if mapping_df is None:
            mapping_df = pl.DataFrame({
                "chrom": list(contig_idx.keys()),
                "chrom_idx": list(contig_idx.values()),
            }).with_columns(pl.col("chrom").cast(pl.Categorical))
            
        m_df = mapping_df.lazy() if isinstance(mapping_df, pl.DataFrame) else mapping_df
        
        df = df.join(m_df, on="chrom", how="inner")
        df = df.with_columns(pl.col('chrom_idx').cast(pl.UInt32))
        
        if edge_length:
            if mapping_df_len is None:
                mapping_df_len = pl.DataFrame({
                    "chrom": list(contigsizes.keys()),
                    "length": list(contigsizes.values()),
                }).with_columns(pl.col("chrom").cast(pl.Categorical))
                
            m_len = mapping_df_len.lazy() if isinstance(mapping_df_len, pl.DataFrame) else mapping_df_len
            
            df = (
                df
                .join(m_len, on="chrom", how="inner")
                .with_columns([
                    pl.col('length').cast(pl.UInt32),
                    (pl.col('start') + (pl.col('end') - pl.col('start')) // 2).alias('pos')
                ])
                .filter(
                    (pl.col('pos') < edge_length) | ((pl.col('length') - pl.col('pos')) < edge_length)
                )
            )
        
    df = df.select(['read_idx', 'chrom_idx', 'mapping_quality'])

    df = df.with_columns([
        pl.col('chrom_idx').n_unique().over('read_idx').alias('chrom_idx_nunique')
    ])
    
    df = df.filter(
        (pl.col('chrom_idx_nunique') >= min_order) & 
        (pl.col('chrom_idx_nunique') <= max_order)
    )

    df = (
        df
        .group_by(['read_idx', 'chrom_idx'])
        .agg([
            pl.col('chrom_idx').count().alias('chrom_idx_count'),
            pl.col('mapping_quality').max() 
        ])
    )

    df = df.select(['read_idx', 'chrom_idx', 'mapping_quality', 'chrom_idx_count'])

    return df.collect()

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
                            edge_length=2e6,
                            split_length=None,
                            split_contig_boundarys=None,
                            hcr_bed=None,
                            hcr_invert=False,
                            threads=4,
                            is_parquet=False,
                            log_dir="logs"):
    
        self.pore_c_table_pathes = listify(pore_c_table_pathes)
        self.contig_idx = contig_idx
        self.contigsizes = contigsizes
        self.min_order = min_order
        self.max_order = max_order
        self.min_alignments = min_alignments
        self.min_quality = min_quality
        self.edge_length = edge_length
        self.split_length = split_length
        self.split_contig_boundarys = split_contig_boundarys
        self.hcr_bed = hcr_bed
        self.hcr_invert = hcr_invert
        self.threads = threads
        self.mapping_df = pl.DataFrame({
            "chrom": list(self.contig_idx.keys()),
            "chrom_idx": list(self.contig_idx.values()),
        }).with_columns(pl.col("chrom").cast(pl.Categorical))

        self.mapping_df_len = None
        if self.edge_length:
            self.mapping_df_len = pl.DataFrame({
                "chrom": list(self.contigsizes.keys()),
                "length": list(self.contigsizes.values()),
            }).with_columns(pl.col("chrom").cast(pl.Categorical))

        self.split_map = None
        if self.split_length and self.split_contig_boundarys:
            self.split_map = pl.DataFrame({
                "chrom": list(self.split_contig_boundarys.keys()),
                "init_idx": list(self.split_contig_boundarys.values()),
            }).with_columns(pl.col("chrom").cast(pl.Categorical))


        os.environ["POLARS_MAX_THREADS"] = str(threads)

        self.is_parquet = is_parquet

        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        
        self.pore_c_tables = self.import_pore_c_table()
        self.edges = self.generate_edges()

        
    def import_pore_c_table_slow(self):
        """
        import pore-c table from pore
        """
        logger.info("Loading Pore-C table ...")
        pl.enable_string_cache()
        if self.edge_length or self.split_length or self.hcr_bed:
            columns = ['read_idx', 'chrom', 'start', 'end', 'mapping_quality']
            schema = {'read_idx': pl.UInt32, 
                        'start': pl.UInt32,
                        'end': pl.UInt32,
                      'mapping_quality': pl.Int8, 
                      'chrom': pl.Categorical}
            usecols = [0, 5, 6, 7, 8]
        else:
            columns = ['read_idx', 'chrom', 'mapping_quality']
            schema = {'read_idx': pl.UInt32, 'mapping_quality': pl.Int8, 'chrom': pl.Categorical}
            usecols = [0, 5, 8]

        if len(self.pore_c_table_pathes) == 1:
            if self.hcr_bed is None:
                infile = self.pore_c_table_pathes[0]
            else:
                porec_prefix = str(Path(self.pore_c_table_pathes[0]).name).replace(".gz", "").rsplit(".", 1)[0]
                if str(self.pore_c_table_pathes[0]).endswith(".gz"):
                    cmd0 = decompress_cmd(str(self.pore_c_table_pathes[0]), str(self.threads))
                    cmd = f"{' '.join(cmd0)} 2>{self.log_dir}/{porec_prefix}.decompress.hcr.log | cphasing-rs porec-intersect - {self.hcr_bed} -o temp.{porec_prefix}.hcr.porec.gz"
                else:
                    cmd = f"cphasing-rs porec-intersect {self.pore_c_table_pathes[0]} {self.hcr_bed} -o temp.{porec_prefix}.hcr.porec.gz"
                if self.hcr_invert:
                    cmd += " --invert"

                cmd += f" 2>{self.log_dir}/{porec_prefix}.porec.intersect.log"
                logger.info(f"Generating hcr porec table by {self.hcr_bed} ...")
                
                flag = subprocess.run(cmd, shell=True)
                assert flag.returncode == 0, "Failed to execute command, please check log."

                infile = f"temp.{porec_prefix}.hcr.porec.gz"

            if self.is_parquet:
                df = pd.read_parquet(infile, 
                                        columns=['read_idx', 'chrom',
                                                # 'start', 'end', 
                                                'mapping_quality',
                                                #'filter_reason',
                                                ],
                                    engine=PQ_ENGINE)
                
            else: 
                logger.debug("Start to load one porec table.")
                try:
                    df = pl.read_csv(infile, separator='\t', has_header=False,
                                    columns=usecols,
                                    dtypes=schema,
                                    new_columns=columns)
                except pl.exceptions.NoDataError:
                    logger.error(f"The pore-c table `{infile}` is empty, can not load anything, please check it.")
                    sys.exit(-1)
                except IsADirectoryError:
                    logger.error(f"The pore-c table `{infile}` is a directory, may be you want to load pairs.pqs file please specified `--pairs` parameter.")
                    sys.exit(-1)
                except pl.exceptions.OutOfBoundsError:
                    logger.error(f"The pore-c table `{infile}` is incorrect, it's not a pore-c table, "
                                 f"may be it is a pairs, if that you should change `-pct` to `-prs`, or the pore-c table is incomplete.")
                    sys.exit(-1)
                
            if self.hcr_bed:
                if Path(f"temp.{porec_prefix}.hcr.porec.gz").exists():
                    os.remove(f"temp.{porec_prefix}.hcr.porec.gz")

            if self.min_quality > 0:
                df = df.filter(pl.col('mapping_quality') >= self.min_quality)
            
           
            df_list = [df]
            logger.debug("Successful load porec table.")

        else:
            infiles = []
            for i in self.pore_c_table_pathes:
                infiles.append(i)

            if not self.hcr_bed:
                def get_file(x):
                    return x 
            else:
                def get_file(x):
                    porec_prefix = str(Path(x).name).replace(".gz", "").rsplit(".", 1)[0]
                    cmd = f"cphasing-rs porec-intersect {x} {self.hcr_bed} -o temp.{porec_prefix}.hcr.porec.gz"
                    if self.hcr_invert:
                        cmd += " --invert"
                    process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
                    stdout, stderr = process.communicate()

                    return f"temp.{porec_prefix}.hcr.porec.gz"

            if self.is_parquet:
                df_list = list(map(lambda x: pd.read_parquet(
                                x, columns=['read_idx', 'chrom', 
                                            # 'start', 'end', 
                                            'mapping_quality',
                                            #'filter_reason'
                                            ], 
                                    engine=PQ_ENGINE,),
                                    infiles))
                
                # df_list = Parallel(n_jobs=self.threads)(delayed(
                #                     lambda x: x.query("filter_reason == 'pass'")
                #                                 .drop("filter_reason", axis=1)
                #                                 )(i) for i in df_list)

                
            else:
                df_list = list(map(lambda x: pl.read_csv(
                    get_file(x), separator='\t', has_header=False,
                    columns=usecols, dtypes={'mapping_quality': pl.Int8, 'chrom': pl.Categorical},
                    new_columns=columns,
                ), infiles))
                
                if self.min_quality > 0:
                    df_list = Parallel(n_jobs=self.threads)(delayed(
                        lambda x: x.filter(pl.col('mapping_quality') >= self.min_quality))(i) for i in df_list)
                
                # if self.edge_length:
                #     df_list = Parallel(n_jobs=self.threads)(delayed(
                #         lambda x: x.with_columns([
                #             pl.col('chrom').map_elements(self.contigsizes.get).alias('length'),
                #         ]).filter(
                #             (pl.col('start') < self.edge_length) | (pl.col('end') > pl.col('length') - self.edge_length)
                #         ).select(['read_idx', 'chrom', 'mapping_quality'])
                #     )(i) for i in df_list)
               
        return df_list
    
    def import_pore_c_table(self):
        """
        import pore-c table from pore
        """
        logger.info("Loading Pore-C table ...")
        pl.enable_string_cache()
        if self.edge_length or self.split_length or self.hcr_bed:
            columns = ['read_idx', 'chrom', 'start', 'end', 'mapping_quality']
            schema = {'read_idx': pl.UInt32, 
                        'start': pl.UInt32,
                        'end': pl.UInt32,
                      'mapping_quality': pl.Int8, 
                      'chrom': pl.Categorical}
            usecols = [0, 5, 6, 7, 8]
        else:
            columns = ['read_idx', 'chrom', 'mapping_quality']
            schema = {'read_idx': pl.UInt32, 'mapping_quality': pl.Int8, 'chrom': pl.Categorical}
            usecols = [0, 5, 8]

        hcr_df = None
        if self.hcr_bed:
            logger.info("Pre-loading and sorting HCR BED...")
            hcr_df = pl.read_csv(
                self.hcr_bed,
                separator="\t",
                has_header=False,
                new_columns=["chrom", "start_bed", "end_bed"],
                dtypes={"chrom": pl.Categorical, "start_bed": pl.UInt32, "end_bed": pl.UInt32}
            ).sort("start_bed").rename({"start_bed": "pos"}) 

        def filter_hcr_in_memory(df, preloaded_hcr_df, invert=False):
            if preloaded_hcr_df is None or df.height == 0:
                return df
                
            orig_order = df.with_row_count("__orig_idx")
            df_with_pos = orig_order.with_columns(
                (pl.col("start") + (pl.col("end") - pl.col("start")) // 2).alias("pos")
            ).sort("pos")
            
            joined = df_with_pos.join_asof(
                preloaded_hcr_df,
                on="pos",
                by="chrom",
                strategy="backward"
            )
            
            is_overlapping = (pl.col("pos") < pl.col("end_bed")) & pl.col("end_bed").is_not_null()
            if invert:
                filtered = joined.filter(~is_overlapping)
            else:
                filtered = joined.filter(is_overlapping)
                
            return filtered.sort("__orig_idx").drop(["__orig_idx", "pos", "end_bed"])

        def _load_and_filter_single(infile):
            if self.is_parquet:
                df = pl.read_parquet(infile, columns=columns)
            else:
                logger.debug(f"Start to load and parse porec table: {infile}")
                try:
                    df = pl.read_csv(infile, separator='\t', has_header=False,
                                    columns=usecols,
                                    dtypes=schema,
                                    new_columns=columns)
                except pl.exceptions.NoDataError:
                    logger.error(f"The pore-c table `{infile}` is empty, cannot load anything, please check it.")
                    sys.exit(-1)
                except IsADirectoryError:
                    logger.error(f"The pore-c table `{infile}` is a directory, may need --pairs specified.")
                    sys.exit(-1)
                except pl.exceptions.OutOfBoundsError:
                    logger.error(f"The pore-c table `{infile}` format error, please check core columns.")
                    sys.exit(-1)

            if self.min_quality > 0:
                df = df.filter(pl.col('mapping_quality') >= self.min_quality)
            
            if hcr_df is not None:
                df = filter_hcr_in_memory(df, hcr_df, self.hcr_invert)
                
            if not (self.edge_length or self.split_length) and 'start' in df.columns:
                df = df.drop(['start', 'end'])

            return df

        df_list = Parallel(n_jobs=min(self.threads, len(self.pore_c_table_pathes)), backend="threading")(
            delayed(_load_and_filter_single)(infile) for infile in self.pore_c_table_pathes
        )
               
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

        logger.debug("Processing Pore-C table ...")
        if len(self.pore_c_tables) > 1:
            for i, pore_c_table in enumerate(self.pore_c_tables):
                args.append((pore_c_table, self.contig_idx, 
                             self.contigsizes, threads_2, 
                            self.min_order, self.max_order, 
                            self.min_alignments, self.is_parquet, 
                            self.edge_length, self.split_length,
                            self.split_contig_boundarys,
                            self.mapping_df, self.mapping_df_len, self.split_map))
           
            res = Parallel(n_jobs=threads_1)(
                            delayed(process_pore_c_table)(*a) for a in args)
        else:
            res = [process_pore_c_table(self.pore_c_tables[0], self.contig_idx, 
                                        self.contigsizes, threads_2, 
                                        self.min_order, self.max_order, 
                                        self.min_alignments, self.is_parquet,
                                        self.edge_length, self.split_length,
                                        self.split_contig_boundarys,
                                        self.mapping_df, self.mapping_df_len, self.split_map)]

        if len(res) < 1:
            raise ValueError("No pore-c reads are retained, please check the pore-c table.")

        # rows = []
        # cols = []
        # counts = []
        # mapqs = []
        
        # current_read_idx_offset = 0
        # for df in res:
        #     rows.append(df['chrom_idx'].to_numpy())
        #     cols.append(df['read_idx'].to_numpy() + current_read_idx_offset)
        #     counts.append(df['chrom_idx_count'].to_numpy())
        #     mapqs.append(df['mapping_quality'].to_numpy())

        #     if len(df) > 0:
        #         current_read_idx_offset = cols[-1].max() + 1

        #     del df 

        # del res 
        
        # edges = HyperEdges(idx=self.contig_idx, 
        #                row=np.concatenate(rows).astype(np.uint32),
        #                col=np.concatenate(cols).astype(np.uint32),
        #                count=np.concatenate(counts).astype(np.uint32),
        #                mapq=np.concatenate(mapqs).astype(np.int8),
        #                contigsizes=self.contigsizes)
        # del rows, cols, counts, mapqs

        current_read_idx_offset = 0
        for i in range(len(res)):
            if len(res[i]) > 0:
                res[i] = res[i].with_columns(pl.col("read_idx") + current_read_idx_offset)
                current_read_idx_offset = res[i]["read_idx"].max() + 1

        full_df = pl.concat(res).sort("chrom_idx")
        
        del res
        gc.collect()

        edges = HyperEdges(
            idx=self.contig_idx, 
            row=full_df['chrom_idx'].to_numpy().astype(np.int32),
            col=full_df['read_idx'].to_numpy().astype(np.int32),
            count=full_df['chrom_idx_count'].to_numpy().astype(np.uint32),
            mapq=full_df['mapping_quality'].to_numpy().astype(np.int8),
            contigsizes=self.contigsizes
        )
        del full_df
        gc.collect()

        
        number_of_contigs = len(self.contig_idx)
        if len(edges.col) == 0:
            logger.error("No valid hyperedges found after all filters. Check your min_order and mapping_quality.")
            number_of_hyperedges = 0
        else:
            number_of_hyperedges = edges.col.max() + 1
            
        logger.info(f"Result of {number_of_hyperedges:,} raw "
                    f"hyperedges of {number_of_contigs:,} contigs. "
                    "Note: it's not the final statistics for hypergraph.")

        return edges
    
    def save(self, output):
        with open(output, 'wb') as out:
            out.write(msgspec.msgpack.encode(self.edges, enc_hook=numpy_enc_hook))

        logger.info(f"Successful output hypergraph into `{output}`")
