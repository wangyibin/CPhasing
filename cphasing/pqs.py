#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Library of .pqs file format, which is dictory contain several parquet files and a _contigsizes file.
"""

import argparse
import logging
import datetime
import os
import os.path as op
import sys
import io 
import shutil
import pprint 

import pandas as pd
import numpy as np

import bisect
import polars as pl 
import glob 

from collections import defaultdict
from multiprocessing import Pool    
from pytools import natsorted
from joblib import Parallel, delayed
from pathlib import Path
from tempfile import TemporaryDirectory
from polars import String, Categorical, UInt32, UInt64, Utf8, UInt8

from .__init__ import __url__
from ._config import *
from .core import Pairs2
from .utilities import read_chrom_sizes, list_flatten, run_cmd
from .utilities import xopen, binnify, decompress_cmd


logger = logging.getLogger(__name__)


__PQS_VERSION__ = "0.1.0"

PAIRS_HEADER = [
    "read_idx",
    "chrom1",
    "pos1",
    "chrom2",
    "pos2",
    "strand1",
    "strand2",
    "mapq",
]
PAIRS_HEADER_WITHOUT_MAPQ = [
    "read_idx",
    "chrom1",
    "pos1",
    "chrom2",
    "pos2",
    "strand1",
    "strand2",
]
PAIRS_SCHEMA_32 = {
    "read_idx": pl.Utf8,
    "chrom1": pl.Categorical,
    "pos1": pl.UInt32,
    "chrom2": pl.Categorical,
    "pos2": pl.UInt32,
    "strand1": pl.Categorical,
    "strand2": pl.Categorical,
    "mapq": pl.UInt8
}

PAIRS_SCHEMA_64 = {
    "read_idx": pl.Utf8,
   "chrom1": pl.Categorical,
    "pos1": pl.UInt64,
    "chrom2": pl.Categorical,
    "pos2": pl.UInt64,
    "strand1": pl.Categorical,
    "strand2": pl.Categorical,
    "mapq": pl.UInt8
}

PAIRS_DTYPES_32 = {
    "read_idx": "U",
    "chrom1": "U",
    "pos1": "uint32",
    "chrom2": "U",
    "pos2": "uint32",
    "strand1": "U",
    "strand2": "U",
    "mapq": "uint8"
}

PAIRS_DTYPES_64 = {
    "read_idx": "U",
    "chrom1": "U",
    "pos1": "uint64",
    "chrom2": "U",
    "pos2": "uint64",
    "strand1": "U",
    "strand2": "U",
    "mapq": "uint8"
}


_README = """
# .pqs file format
The _contigsizes file contains the size of each contig in the .pqs file.
The _metadata file contains the metadata of the .pqs file.
The q0, q1 directory contains the parquet files of the data.

q0 mean the mapping quality of data >= 0.
q1 mean the mapping quality of data >= 1.

.pqs/  
|-- _contigsizes 
|-- _metadata
|-- _readme 
|-- q0/
|   |-- 0.parquet
|   |-- 1.parquet
|   |-- ...
|-- q1/
|   |-- 0.parquet
|   |-- 1.parquet
|   |-- ...
"""


class PQS:
    """
    A class to represent a .pqs file.
    """
    def __init__(self, path=None, threads=10):
        self._readme = _README
        self.path = path 
        self.threads = threads

        os.environ["POLARS_MAX_THREADS"] = str(self.threads)


    @property
    def contigsizes(self):
        """
        Get the contig sizes of the .pqs file.
        """

        return self._contigsizes

    @property
    def contigsizes_db(self):
        """
        Get the contig sizes of the .pqs file as dict.
        """

        return self.contigsizes.to_dict()['length']

    @property
    def chromsizes(self):
        """
        Get the contig sizes of the .pqs file.
        """

        return self._contigsizes

    @property
    def metadata(self):
        """
        Get the metadata of the .pqs file.
        """

        return self._metadata
    
    @property
    def readme(self):
        """
        Get the readme of the .pqs file.
        """
        return self._readme
    
    @property
    def q0(self):
        """
        Get the q0 directory of the .pqs file.
        """
        return self._q0
    
    @property
    def q1(self):
        """
        Get the q1 directory of the .pqs file.
        """
        return self._q1
    
    @property
    def schema(self):
        """
        Get the schema of the .pqs file.
        """
        return self._schema

    @staticmethod
    def is_changed(path):
        """
        Check if the .pqs file is changed.
        """
        assert PQS.is_pqs(path), "The given path is not a .pqs file."
        with open(Path(path) / "_metadata", "r") as f:
            metadata = eval(f.read())
        
        input_file_path = Path(path).absolute()
        file_name = f"_metadata"
        prefix = f"_metadata"

        flag = True
      
        if Path(f"{input_file_path}/._metadata.md5.txt").exists():
            text = os.popen(f"md5sum -c {input_file_path}/.{prefix}.md5.txt 2>/dev/null").read()
            if not text.strip():
                flag =  True
            if text.strip().split()[-1] == "OK":
                flag = False
            else:
                os.system(f"md5sum {input_file_path}/{file_name} > {input_file_path}/.{prefix}.md5.txt")
                flag = True
        else:
            os.system(f"md5sum {input_file_path}/{file_name} > {input_file_path}/.{prefix}.md5.txt")
            flag =  True
                
        return flag

    @staticmethod
    def is_pqs(path):
        """
        Check if the given path is a .pqs file.
        """
        path = Path(path)
        if not path.exists():
            return False
        
        if path.is_dir():
            if (path / "_metadata").exists() and (path / "_contigsizes").exists():
                if not os.listdir(path / "q0") and not os.listdir(path / "q1"):
                    return False
                _metadata = eval(open(path / "_metadata", "r").read())
                try:
                    if _metadata["is_pqs"]:
                        return True 
                    
                except KeyError:
                    return False
            else:
                return False
    
    def is_pairs(self, path=None):
        """
        Check if the given path is a .pairs file.
        """
        if path is None:
            path = self.path

        path = Path(path)
        if not path.exists():
            return False
        
        if not self.is_pqs(path):
            return False
        
        
        if self._metadata["format"] == "pairs":
            return True
        else:
            return False
        

    def init_read(self, path=None):
        """
        load chromsizes and metadata
        """
        assert path is not None or self.path is not None, "Please provide a path of .pqs."
        if path is None:
            path = self.path
        else:
            self.path = path
        
        assert self.is_pqs(
            path
        ), "The given path is not a .pqs file or could not found this path."

        self._contigsizes = read_chrom_sizes(Path(path) / "_contigsizes").sort_index()
    
        with open(Path(path) / "_metadata", "r") as f:
            self._metadata = eval(f.read())

        self._schema = self._metadata["schema"]
    
    def read(self, path=None, column_names=None,
                columns=None, min_mapq=0, 
                return_as="generator", low_memory=True):
        """
        Read the .pqs file from the given path.
        """
        assert path is not None or self.path is not None, "Please provide a path of .pqs."
        if path is None:
            path = self.path

        # assert self.is_pqs(
        #     path
        # ), "The given path is not a .pqs file or could not found this path."


        if columns is None and column_names is None:
            columns = list(range(len(self._metadata["columns"])))
        elif columns is None and column_names is not None:
            if 'mapq' in column_names:
                if not self._metadata["is_with_mapq"]:
                    column_names.remove('mapq')

            columns = [self._metadata["columns"].index(col) for col in column_names]
        elif columns is not None and column_names is not None:
            logger.warning(
                "Both columns and column_names are provided, column_names will be used."
            )
            columns = [self._metadata["columns"].index(col) for col in column_names]
        
        self._q0 = Path(path) / "q0"
        self._q1 = Path(path) / "q1"
        
        if self._metadata["is_with_mapq"]:
            if min_mapq == 0:
                target_dir = self._q0
                is_filter = False
            else:
                target_dir = self._q1
                if min_mapq > 1:
                    is_filter = True
                else:
                    is_filter = False
        else:
            target_dir = self._q0
            is_filter = False

        if min_mapq > 0:
            logger.info(f"Filtered the data with mapq >= {min_mapq}.")
        
        if return_as == "generator":
            if low_memory:
                for file in target_dir.iterdir():
                    df = pl.read_parquet(file, columns=columns)
                    if is_filter:
                        df = df.filter(pl.col("mapq") >= min_mapq)
                    else:
                        if self._metadata['is_with_mapq']:
                            if column_names:
                                if 'mapq' not in column_names:
                                    df = df.drop('mapq')
                
                    yield df
            else:
                df = pl.read_parquet(target_dir/"*.parquet", columns=columns)
                if is_filter:
                    df = df.filter(pl.col("mapq") >= min_mapq).drop('mapq')
                else:
                    if self._metadata['is_with_mapq']:
                            if column_names:
                                if 'mapq' not in column_names:
                                    df = df.drop('mapq')

                yield df

        elif return_as == "files":
            for file in target_dir.iterdir():
                yield file 
        
        else:
            raise ValueError("The return_as must be 'generator' or 'files'.")
        

    def to_depth(self, chunks, output, min_mapq=0, binsize=10000):
        """
        Convert the .pqs file to depth file.
        """
        is_with_mapq = self._metadata["is_with_mapq"]
        chromsizes_db = self.contigsizes_db
        chrom_list = []
        start_list = []
        end_list = []
        
        for chrom, clen in chromsizes_db.items():
            n_bins = int(np.ceil(clen / binsize))
            binedges = np.arange(0, (n_bins + 1)) * binsize
            binedges[-1] = clen
            chrom_list.extend([chrom] * n_bins)
            start_list.extend(binedges[:-1])
            end_list.extend(binedges[1:])

        bins = pd.DataFrame({
            "chrom": chrom_list,
            "start": start_list,
            "end": end_list
        })

        bins['count'] = 0
        
        def process_chunk(chunk):
            chunk = chunk.with_columns(
                ((pl.col('pos1')) // binsize).alias('pos1'),
                ((pl.col('pos2')) // binsize).alias('pos2')
            )
            chunk1 = (
                chunk.select(["chrom1", "pos1", "mapq"])
                .group_by(["chrom1", "pos1"])
                .agg([pl.col("mapq").count().alias("count")])
            )
            chunk2 = (
                chunk.select(["chrom2", "pos2", "mapq"])
                .group_by(["chrom2", "pos2"])
                .agg([pl.col("mapq").count().alias("count")])
            )

            # chunk = chunk1.join(
            #     chunk2,
            #     left_on=["chrom1", "pos1"],
            #     right_on=["chrom2", "pos2"],
            #     how="outer",

            # )

            # chunk = chunk.drop_nulls(subset=["chrom1", "pos1", "chrom2", "pos2" ])
        
            # chunk = chunk.with_columns(
            #     (pl.col('count') + pl.col('count_right')).alias('count')
            # )
  
            # chunk = chunk.drop(['count_right'])
            # chunk = chunk.with_columns(
            #     ((pl.col('pos1') + 1) * binsize - binsize).alias('start'),
            #     ((pl.col('pos1') + 1) * binsize).alias('end')
            # ).drop(['pos1'])

            # chunk = chunk.select(['chrom1', 'start', 'end', 'count'])
            # chunk.columns = ['chrom', 'start', 'end', 'count']
            # chunk = chunk.to_pandas()
     
            # chunk['chrom'] = chunk['chrom'].astype(str)
            # chunk['length'] = chunk['chrom'].map(chromsizes_db.get)
            # chunk["end"] = np.where(
            #     chunk["end"] > chunk["length"], chunk["length"], chunk["end"]
            # )
            # chunk.drop(['length'], axis=1, inplace=True)
            # chunk.fillna(0, inplace=True)
            # chunk['end']  = chunk['end'].astype(int)
            # chunk.set_index(['chrom', 'start', 'end'], inplace=True)

            # return chunk

    
            chunk1 = chunk1.with_columns(
                ((pl.col('pos1') + 1) * binsize - binsize).alias('start'),
                ((pl.col('pos1') + 1) * binsize).alias('end')
            ).drop(['pos1'])

            chunk2 = chunk2.with_columns(
                ((pl.col('pos2') + 1) * binsize - binsize).alias('start'),
                ((pl.col('pos2') + 1) * binsize).alias('end')
            ).drop(['pos2'])

            chunk1 = chunk1.select(['chrom1', 'start', 'end', 'count'])
            chunk2 = chunk2.select(['chrom2', 'start', 'end', 'count'])

            chunk1.columns = ['chrom', 'start', 'end', 'count']
            chunk2.columns = ['chrom', 'start', 'end', 'count']
            

            chunk1 = chunk1.with_columns([
                pl.col('chrom').map_elements(chromsizes_db.get).alias('length')
            ])

            chunk1 = chunk1.with_columns([
                pl.when(pl.col('end') > pl.col('length'))
                .then(pl.col('length'))
                .otherwise(pl.col('end'))
                .alias('end')
            ])
            chunk1 = chunk1.drop(['length'])

            chunk2 = chunk2.with_columns([
                pl.col('chrom').map_elements(chromsizes_db.get).alias('length')
            ])

            chunk2 = chunk2.with_columns([
                pl.when(pl.col('end') > pl.col('length'))
                .then(pl.col('length'))
                .otherwise(pl.col('end'))
                .alias('end')
            ])
            chunk2 = chunk2.drop(['length'])

            return chunk1, chunk2
        
        def process_chunk_to_depth(chunks):
            for chunk in chunks:
                yield process_chunk(chunk)
        
       
        # iterator = process_chunk_to_depth(chunks)
        
        iterator = Parallel(n_jobs=self.threads, return_as='generator')(
            delayed(process_chunk_to_depth_global)(*arg) for arg in ((chunk, binsize, chromsizes_db, min_mapq, is_with_mapq) for chunk in chunks)
        )
        
        bins = pl.DataFrame(bins)
        bins = bins.cast(
            {
                "chrom": pl.Utf8,
                "start": self._schema["pos1"],
                "end": self._schema["pos1"],
                "count": pl.UInt32,
            }
        )

        for i, results in enumerate(iterator):
            for result in results:
                result = result.cast(
                    {
                        "chrom": pl.Utf8,
                        "start": self._schema["pos1"],
                        "end": self._schema["pos1"],
                        "count": pl.UInt32,
                    }
                )
              
                bins = bins.join(
                    result, how="left", on=["chrom", "start", "end"], suffix="_right"
                )

                bins = bins.with_columns([pl.col("count").fill_null(0) + pl.col("count_right").fill_null(0)])
                
                bins = bins.drop(["count_right"])
                bins = bins.cast(
                    {
                        "chrom": pl.Utf8,
                        "start": self._schema["pos1"],
                        "end": self._schema["pos1"],
                        "count": pl.UInt32,
                    }
                )
        
        os.environ["POLARS_MAX_THREADS"] = str(self.threads)

        bins.write_csv(output, separator='\t', include_header=False)
        logger.info(f"Successful output depth file into `{output}`.")

    def to_cool(self, chunks, output, 
                binsize=10000, min_mapq=1, low_memory=False):
        """
        Convert the .pqs file to .cool file.
        """
        import cooler 
        from cooler.create import create, create_cooler

        is_with_mapq = self._metadata["is_with_mapq"]

        bins = binnify(self.contigsizes['length'], binsize=binsize)
        if len(bins) > 2**32:
            self.bin_id_dtype = pl.UInt64
        elif len(bins) > 2**16:
            self.bin_id_dtype = pl.UInt32
        elif len(bins) > 2**8:
            self.bin_id_dtype = pl.UInt16
        else:
            self.bin_id_dtype = pl.UInt8

        bin_offset = bins.groupby('chrom').size().cumsum()
        bin_offset = bin_offset.shift(1).fillna(0).astype(int)
        bin_offset_db = bin_offset.to_dict()

        def process_chunk(chunk):
            bin1_id = (chunk["pos1"] // binsize) + chunk["chrom1"].map_elements(
                bin_offset_db.get
            )
            bin2_id = (chunk["pos2"] // binsize) + chunk["chrom2"].map_elements(
                bin_offset_db.get
            )
            chunk = (
                chunk.with_columns([bin1_id.alias("bin1_id"), bin2_id.alias("bin2_id")])
                .drop(["chrom1", "chrom2", "pos1", "pos2"])
                .filter(pl.col("bin1_id") <= pl.col("bin2_id"))
            )

            chunk = chunk.group_by(["bin1_id", "bin2_id"], maintain_order=True).agg(
                pl.count("bin1_id").alias("count")
            )

            chunk = chunk.sort(['bin1_id', 'bin2_id'])
     
            return chunk
        
        def process_chunk_to_cool(chunks):
            for chunk in chunks:
                yield process_chunk(chunk).to_pandas()
        
        if low_memory:
            iterator = process_chunk_to_cool(chunks)

            create_cooler(output, bins, iterator, 
                      triucheck=False, dupcheck=False, boundscheck=False,)
        
        else:
            chunks = list(map(lambda x: Path(x).absolute(), chunks))
            results = Parallel(n_jobs=self.threads)(
                delayed(process_chunk_to_cool_global)(*arg)
                for arg in ((chunk, binsize, bin_offset_db, self._schema,
                             min_mapq, is_with_mapq) for chunk in chunks)
            )

            os.environ["POLARS_MAX_THREADS"] = str(self.threads)

            pixels = pl.concat(results)

            pixels = pixels.group_by(["bin1_id", "bin2_id"], maintain_order=True).agg(
                pl.sum("count").alias("count")
            )
            pixels = pixels.sort(["bin1_id", "bin2_id"])  
            create(output, bins, pixels.to_pandas(), dtypes={'count': np.int32},
                    triucheck=False, dupcheck=False, boundscheck=False)

        logger.info(f"Successful output cooler file into `{output}`.")
    
    def write(self):
        """
        Write the .pqs file to the given path.
        """
        pass
    
    def to_hg_df(self, chunks, contig_idx, 
                 min_mapq=1, edge_length=0,
                 split_length=None, split_contig_boundarys=None,
                 bed=None, hcr_binsize=10000):
        from intervaltree import IntervalTree
        pl.enable_string_cache()
        contigsizes = self.contigsizes_db

        if bed is None:
            bed_dict = None
        else:
            bed_df = pl.read_csv(bed, separator="\t", has_header=False,
                                columns=[0, 1, 2], new_columns=['chrom', 'start', 'end'])

            if edge_length > 0:
                bed_df = bed_df.with_columns([
                    pl.col('chrom').map_elements(contigsizes.get).alias('length')
                ]).filter(
                    (pl.col('start') < edge_length) | (pl.col('end') > pl.col('length') - edge_length)
                )
 
            bed_df = bed_df.to_pandas()
            bed_df.set_index(['chrom'], inplace=True)
            bed_df['region'] = bed_df.apply(lambda x: (x['start'], x['end']), axis=1)
            bed_dict = defaultdict(list)

            for chrom, region in bed_df.groupby(level=0):
                bed_dict[chrom] = region['region'].values.tolist()
                # bed_dict[chrom] = (region['start'].values, region['end'].values)
                # bed_dict[chrom] = IntervalTree.from_tuples(region['region'])

        logger.info("Parsing pqs ...")

        args = []
        for chunk in chunks:
            args.append((Path(chunk).absolute(), bed_dict, contigsizes, contig_idx, min_mapq, edge_length,
                         self.schema, split_length, split_contig_boundarys))
      
        results = Parallel(n_jobs=self.threads)(
                    delayed(process_chunk_hg)(*arg) for arg in args
                )
        
        results = list(filter(lambda x: x is not None, results))
        results = list(filter(lambda x: len(x) > 0, results))
        if len(results) == 0:
            logger.warning("No data found in the given region.")
            return
        
        hg_df = pl.concat(results)

        return hg_df

    def to_cis_depth_df(self, chunks, window_size=500, min_mapq=0):
        pl.enable_string_cache()
        is_with_mapq = self._metadata["is_with_mapq"]   
        args = []
        for chunk in chunks:
            args.append((chunk, window_size, min_mapq, is_with_mapq))
        
        results = Parallel(n_jobs=self.threads)(
                    delayed(process_chunk_cis_depth)(*arg) for arg in args
                )
        
        cis_depth_df = pl.concat(results).collect() 

        return cis_depth_df.to_pandas()
    
    def to_clm(self, chunks, output,
               min_mapq=0,
               min_count=1):
        """
        Convert the .pqs file to .clm file.
        """
        pl.enable_string_cache()

        is_with_mapq = self._metadata["is_with_mapq"]    

        args = []
        

        with TemporaryDirectory(suffix="_clm", delete=False, dir="./") as tmpdir:
            for chunk in chunks:
                args.append((chunk, self.contigsizes_db, self._schema, min_mapq, is_with_mapq, tmpdir))
            Parallel(n_jobs=self.threads, backend="multiprocessing")(
                delayed(process_chunk_clm)(*arg)
                (*arg) for arg in args
            )
        
        cmd = ["cphasing-rs", "clm-merge", tmpdir, "-o", output]
        
        flag = run_cmd(cmd, log=os.devnull)
        assert flag == 0, "Failed to merge the clm files."

    
    def intersect(self, chunks, bed, output, min_mapq=1):
        from intervaltree import IntervalTree
        is_with_mapq = self._metadata["is_with_mapq"]
        bed_df = pl.read_csv(bed, separator="\t", has_header=False,
                             columns=[0, 1, 2], new_columns=['chrom', 'start', 'end'])
        bed_df = bed_df.to_pandas()
        bed_df.set_index(['chrom'], inplace=True)
        bed_df['region'] = bed_df.apply(lambda x: (x['start'], x['end']), axis=1)
        bed_dict = defaultdict(list)
        for chrom, region in bed_df.groupby(level=0):
            # bed_dict[chrom] = IntervalTree.from_tuples(region['region'])
            bed_dict[chrom] = region['region'].values.tolist()

        Path(output).mkdir(exist_ok=True)
        Path(f"{output}/q0").mkdir(exist_ok=True)
        Path(f"{output}/q1").mkdir(exist_ok=True)
      
        Parallel(n_jobs=self.threads)(
            delayed(process_chunk_intersect_global)(chunk, bed_dict, is_with_mapq, min_mapq, output) for chunk in chunks
        )

        shutil.copy(f"{self.path}/_contigsizes", output)
        shutil.copy(f"{self.path}/_metadata", output)
        shutil.copy(f"{self.path}/_readme", output)
    
        logger.info(f"Successful output intersect file into `{output}`.")


    def to_pairs(self, chunks, output):
        """
        Convert the .pqs file to .pairs file.
        """
        pass
    
    
    def from_pairs(self, pairs,
                   path=None,
                   chunksize=1e5):
        """
        Create a .pqs file from the given pairs.
        """
        assert path is not None or self.path is not None, "Please provide a path to save the .pqs file."
        if path is None:
            path = self.path
        
        chunksize = int(chunksize)

        if Path(path).exists():
            logger.warning(f"The path `{path}` already exists, the content will be overwritten.")
            try:
                exist_chunksize = int(eval(open(Path(path) / "_metadata", "r").read())["chunksize"])
                if exist_chunksize > chunksize:
                    logger.warning(f"The chunksize of the existing .pqs file is {exist_chunksize}, "
                                    f"which is larger than the given chunksize {chunksize}.")
                    logger.warning(f"Remove the existing .pqs file.")

                    shutil.rmtree(f"{path}/q0")
                    shutil.rmtree(f"{path}/q1")
            except:
                pass 

        p = Pairs2(pairs)
        is_pairs_with_mapq = p.is_pairs_with_mapq()
        logger.debug(f"pairs with mapq: {is_pairs_with_mapq}")
        contigsizes = p.chromsizes
        max_size = max(contigsizes["length"])
        if is_pairs_with_mapq:
            columns = list(range(8))
            if max_size < 2**32:
                schema = PAIRS_SCHEMA_32
                dtypes = PAIRS_DTYPES_32
                HEADER = PAIRS_HEADER
            else:
                schema = PAIRS_SCHEMA_64
                dtypes = PAIRS_DTYPES_64
                HEADER = PAIRS_HEADER
        else:
            columns = list(range(7))
            if max_size < 2**32:
                schema = PAIRS_SCHEMA_32
                schema.pop("mapq")
                dtypes = PAIRS_DTYPES_32
                HEADER = PAIRS_HEADER_WITHOUT_MAPQ
            else:
                schema = PAIRS_SCHEMA_64
                schema.pop("mapq")
                dtypes = PAIRS_DTYPES_64
                HEADER = PAIRS_HEADER_WITHOUT_MAPQ
        

        Path(path).mkdir(exist_ok=True)
        pairs_path = Path(pairs).absolute()
        q0 = Path(path).absolute() / "q0"
        q1 = Path(path).absolute() / "q1"
        q0.mkdir(exist_ok=True)
        q1.mkdir(exist_ok=True)
        with TemporaryDirectory(suffix="_split_pairs", delete=True, dir="./") as tmpdir:
            os.chdir(tmpdir)

            cmd0 = decompress_cmd(str(pairs_path))

            cmd = [
                "cphasing-rs",
                "pairs-split",
                "-",
                "--chunksize",
                str(int(chunksize)),
                " 2>/dev/null",
            ] 
            os.system(" ".join(cmd0) + " 2>/dev/null" + " | " + " ".join(cmd))

            args = []
            for i, chunk in enumerate(natsorted(os.listdir())):
                args.append((i, chunk, q0, q1, columns, schema, is_pairs_with_mapq))


            Parallel(n_jobs=self.threads)(
                delayed(process_chunk_split_pairs)(*arg) for arg in args
            )

            metadata = {
                "format-version": __PQS_VERSION__,
                "format-url": __url__,
                "creation_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "is_pqs": True,
                "format": "pairs",
                "chunksize": chunksize,
                "is_with_mapq": is_pairs_with_mapq,
                "columns": HEADER,
                "dtypes": dtypes,
                "schema": schema
            }

            os.chdir('..')
        
        contigsizes.to_csv(Path(path) / "_contigsizes", sep="\t", header=False, index=False)
        with open(Path(path) / "_metadata", "w") as f:
            pprint.pprint(metadata, f, sort_dicts=False)
        with open(Path(path) / "_readme", "w") as f:
            f.write(self._readme)

        logger.info(f"Successful output .pqs file into `{path}`.")

    def from_porec_table(self, porec):
        """
        Create a .pqs file from the given porec table.
        """
        pass


def process_chunk_to_depth_global(chunk, binsize, chromsizes_db,
                                 min_mapq, is_with_mapq):

    os.environ["POLARS_MAX_THREADS"] = "1"

    columns = (
        ["chrom1", "pos1", "chrom2", "pos2", "mapq"]
        if is_with_mapq
        else ["chrom1", "pos1", "chrom2", "pos2"]
    )

    chunk = pl.read_parquet(chunk,
                            columns=columns)
    
    if min_mapq > 1:
        if is_with_mapq:
            chunk = chunk.filter(pl.col("mapq") >= min_mapq)

    chunk = chunk.with_columns(
        ((pl.col('pos1')) // binsize).alias('pos1'),
        ((pl.col('pos2')) // binsize).alias('pos2')
    )
    chunk1 = (
        chunk.select(["chrom1", "pos1", "mapq"])
        .group_by(["chrom1", "pos1"])
        .agg([pl.col("mapq").count().alias("count")])
    )
    chunk2 = (
        chunk.select(["chrom2", "pos2", "mapq"])
        .group_by(["chrom2", "pos2"])
        .agg([pl.col("mapq").count().alias("count")])
    )

    chunk1 = chunk1.with_columns(
        ((pl.col('pos1') + 1) * binsize - binsize).alias('start'),
        ((pl.col('pos1') + 1) * binsize).alias('end')
    ).drop(['pos1'])

    chunk2 = chunk2.with_columns(
        ((pl.col('pos2') + 1) * binsize - binsize).alias('start'),
        ((pl.col('pos2') + 1) * binsize).alias('end')
    ).drop(['pos2'])

    chunk1 = chunk1.select(['chrom1', 'start', 'end', 'count'])
    chunk2 = chunk2.select(['chrom2', 'start', 'end', 'count'])

    chunk1.columns = ['chrom', 'start', 'end', 'count']
    chunk2.columns = ['chrom', 'start', 'end', 'count']
    

    chunk1 = chunk1.with_columns([
        pl.col('chrom').map_elements(chromsizes_db.get).alias('length')
    ])

    chunk1 = chunk1.with_columns([
        pl.when(pl.col('end') > pl.col('length'))
        .then(pl.col('length'))
        .otherwise(pl.col('end'))
        .alias('end')
    ])
    chunk1 = chunk1.drop(['length'])

    chunk2 = chunk2.with_columns([
        pl.col('chrom').map_elements(chromsizes_db.get).alias('length')
    ])

    chunk2 = chunk2.with_columns([
        pl.when(pl.col('end') > pl.col('length'))
        .then(pl.col('length'))
        .otherwise(pl.col('end'))
        .alias('end')
    ])
    chunk2 = chunk2.drop(['length'])

    return chunk1, chunk2

def process_chunk_intersect_global(chunk_name, bed_dict, is_with_mapq, min_mapq, output):

    chunk = pl.read_parquet(chunk_name)
    chunk_name = Path(chunk_name).stem
    if min_mapq > 1:
        chunk = chunk.filter(pl.col("mapq") >= min_mapq)

    chunk = chunk.filter(
        pl.struct(["chrom1", "pos1"]).apply(lambda x: is_in_regions(x["chrom1"], x["pos1"], bed_dict)), 
        pl.struct(["chrom2", "pos2"]).apply(lambda x: is_in_regions(x["chrom2"], x["pos2"], bed_dict))
    )

    chunk.write_parquet(f"{output}/q0/{chunk_name}.parquet")
    if is_with_mapq:
        chunk = chunk.filter(pl.col("mapq") >= 1)
    chunk.write_parquet(f"{output}/q1/{chunk_name}.parquet")



def process_chunk_to_cool_global(chunk, binsize, 
                                 bin_offset_db, 
                                 schema,
                                 min_mapq,
                                 is_with_mapq,
                                ):

    os.environ["POLARS_MAX_THREADS"] = "1"

    columns = ["chrom1", "pos1", "chrom2", "pos2", "mapq"] if is_with_mapq else ["chrom1", "pos1", "chrom2", "pos2"]
    chunk = pl.read_parquet(chunk, columns=columns)
    if min_mapq > 1:
        if is_with_mapq:
            chunk = chunk.filter(pl.col("mapq") >= min_mapq).drop("mapq")
    
    bin1_id = (chunk["pos1"] // binsize) + chunk["chrom1"].map_elements(
            bin_offset_db.get, skip_nulls=False
    ).cast(schema["pos1"])
    bin2_id = (chunk["pos2"] // binsize) + chunk["chrom2"].map_elements(
        bin_offset_db.get, skip_nulls=False
    ).cast(schema["pos2"])
    chunk = (
        chunk.with_columns([bin1_id.alias("bin1_id"), bin2_id.alias("bin2_id")])
        .drop(["chrom1", "chrom2", "pos1", "pos2"])
    )
    chunk = chunk.with_columns([
                pl.when(chunk['bin1_id'] > chunk['bin2_id'])
                .then(chunk['bin2_id'])
                .otherwise(chunk['bin1_id'])
                .alias('bin1_id'),
                pl.when(chunk['bin1_id'] > chunk['bin2_id'])
                .then(chunk['bin1_id'])
                .otherwise(chunk['bin2_id'])
                .alias('bin2_id')
            ])
        # .filter(pl.col("bin1_id") <= pl.col("bin2_id"))


    chunk = chunk.group_by(["bin1_id", "bin2_id"], maintain_order=True).agg(
        pl.count("bin1_id").alias("count")
    )

    chunk = chunk.sort(['bin1_id', 'bin2_id'])
    
    return chunk


def is_in_regions2(chrom, pos, bed_dict):
    if chrom not in bed_dict:
        return False
    
    return bed_dict[chrom].overlaps(pos)

def is_in_regions(chrom, pos, bed_dict):
    if chrom not in bed_dict:
        return False

    intervals = bed_dict[chrom]
    if len(intervals) == 0:
        return False 
  
    idx = bisect.bisect_right(intervals, (pos, float('inf')))  

    if idx == 0:
        return False
    
    start, end = intervals[idx - 1]

    return start <= pos < end


def process_chunk_hg(chunk_name, bed_dict, contigsizes,
                     contig_idx, min_mapq, edge_length, schema,
                     split_length=None, split_contig_boundarys=None):
    pl.enable_string_cache()
    os.environ["POLARS_MAX_THREADS"] = "1"
    if not Path(chunk_name).exists():
        return None

    columns = ["chrom1", "pos1", "chrom2", "pos2", "mapq"]
    chunk = pl.scan_parquet(chunk_name).select(columns)
    chunk_name = Path(chunk_name).stem

    
    if min_mapq > 1:
        chunk = chunk.filter(pl.col("mapq") >= min_mapq)

    # chunk = chunk.filter(pl.col("chrom1") != pl.col("chrom2"))
    chunk = chunk.collect()
    if split_length and split_contig_boundarys:
        edge_length = 0

    if edge_length > 0:
        mapping_df = pl.DataFrame(
            {
                "chrom": list(contigsizes.keys()),
                "length": list(contigsizes.values()),
            }
        ).with_columns(
            pl.col("chrom").cast(pl.Categorical)
        )
        
        mapping1 = mapping_df.rename({"chrom": "chrom1"})
        chunk = (
            chunk
            .join(mapping1, on="chrom1", how="inner")
            .rename({"length": "length1"})
        )

        mapping2 = mapping_df.rename({"chrom": "chrom2"})
        chunk = (
            chunk
            .join(mapping2, on="chrom2", how="inner")
            .rename({"length": "length2"})
        )
        chunk = (
            chunk
            .filter(
                (
                    ((pl.col("pos1") < edge_length)
                     | (pl.col("pos1") > (pl.col("length1") - edge_length)))
                    & ((pl.col("pos2") < edge_length)
                       | (pl.col("pos2") > (pl.col("length2") - edge_length)))
                )
            ).select(["chrom1", "pos1", "chrom2", "pos2", "mapq"])
        )

    if bed_dict:
        # chunk = chunk.filter(
        #     pl.struct(["chrom1", "pos1"]).apply(
        #         lambda x: is_in_regions(x["chrom1"], x["pos1"], bed_dict),
                
                
        #     ),
        #     pl.struct(["chrom2", "pos2"]).apply(
        #         lambda x: is_in_regions(x["chrom2"], x["pos2"], bed_dict),
                
        #     ),
        # )
        logger.warning("incomplete, skip the filter of bed regions.")
        pass 
    

    if split_length and split_contig_boundarys:

        chunk = (
            chunk.with_columns((pl.col("pos1") // split_length).alias('sub_idx1'),
                               (pl.col("pos2") // split_length).alias('sub_idx2'),
            )   
        )
        split_map = pl.DataFrame(
            {
                "chrom": list(split_contig_boundarys.keys()),
                "init_idx": list(split_contig_boundarys.values()),
            }
        ).with_columns(
            pl.col("chrom").cast(pl.Categorical)
        )

        chunk = (
            chunk
            .join(split_map, left_on="chrom1", right_on="chrom", how="inner")
            .rename({"init_idx": "init_idx1"})
        )
        chunk = (
            chunk
            .join(split_map, left_on="chrom2", right_on="chrom", how="inner")
            .rename({"init_idx": "init_idx2"})
        )
   

        chunk = chunk.with_columns(
            (pl.col("init_idx1") + pl.col("sub_idx1")).alias("chrom1"),
            (pl.col("init_idx2") + pl.col("sub_idx2")).alias("chrom2"),
        ).drop(["sub_idx1", "sub_idx2", "init_idx1", "init_idx2"])

        chunk = chunk.drop_nulls(subset=["chrom1", "chrom2"])

    else: 
        mapping_df = pl.DataFrame(
            {
                "chrom": list(contig_idx.keys()),
                "chrom_idx": list(contig_idx.values()),
            }
        ).with_columns(
            pl.col("chrom").cast(pl.Categorical)
        )
     
        mapping1 = mapping_df.rename({"chrom": "chrom1"})
        chunk = chunk.join(mapping1, on="chrom1", how="inner")
        chunk = chunk.drop("chrom1").rename({"chrom_idx": "chrom1"})
        
        mapping2 = mapping_df.rename({"chrom": "chrom2"})
        chunk = chunk.join(mapping2, on="chrom2", how="inner")
        chunk = chunk.drop("chrom2").rename({"chrom_idx": "chrom2"})

        chunk = chunk.drop_nulls(subset=["chrom1", "chrom2"])

        # chunk = chunk.with_columns(
        #     pl.col("chrom1").map_elements(contig_idx.get, skip_nulls=False).alias("chrom1"),
        #     pl.col("chrom2").map_elements(contig_idx.get, skip_nulls=False).alias("chrom2"),
        # ).drop_nulls(subset=["chrom1", "chrom2"])

    return chunk



def process_chunk_clm(chunk, contigsizes_db, schema,
                      min_mapq=0, is_with_mapq=True,
                      tmp_dir="."):
    pl.enable_string_cache()
    os.environ["POLARS_MAX_THREADS"] = "1"

    columns = ["chrom1", "pos1", "chrom2", "pos2", "mapq"] if is_with_mapq else ["chrom1", "pos1", "chrom2", "pos2"]
    
    chunk_name = Path(chunk).stem

    chunk = pl.read_parquet(chunk, columns=columns)
    if min_mapq > 1:
        if is_with_mapq:
            chunk = chunk.filter(pl.col("mapq") >= min_mapq).drop("mapq")   

    chunk = chunk.filter(pl.col("chrom1") != pl.col("chrom2"))
    chunk = chunk.with_columns(
        pl.col("chrom1")
        .map_elements(contigsizes_db.get)
        .cast(schema["pos1"])
        .alias("length1"),
        pl.col("chrom2")
        .map_elements(contigsizes_db.get)
        .cast(schema["pos2"])
        .alias("length2"),
    )

    chunk = chunk.with_columns(
        ((pl.col("length1") - pl.col("pos1")) + pl.col("pos2") )
        .cast(pl.String)
        .alias("++"),
        (
            (pl.col("length1") - pl.col("pos1")) 
            + (pl.col("length2") - pl.col("pos2")) 
        )
        .cast(pl.String)
        .alias("+-"),
        (pl.col("pos1") + pl.col("pos2"))
        .cast(pl.String)
        .alias("-+"),
        (pl.col("pos1") + (pl.col("length2") - pl.col("pos2")))
        .cast(pl.String)
        .alias("--"),
    )

    chunk = chunk.drop(['length1', 'length2', 'pos1', 'pos2'])
    try:
        chunk = chunk.unpivot(["++", "+-", "-+", "--"], index=["chrom1", "chrom2"])
    except AttributeError:
        logger.error("The version of polars is too low, please update to the version of >= 1.*")
        sys.exit(1)

    chunk = (
        chunk.group_by(["chrom1", "chrom2", "variable"])
        .agg(pl.col("value"))
        .sort(["chrom1", "chrom2", "variable"])
    )

    chunk = chunk.with_columns(
        (
            pl.col("chrom1")
            + pl.col("variable").str.slice(0, 1)
            + pl.lit(" ")
            + pl.col("chrom2")
            + pl.col("variable").str.slice(1, 1)
        ).cast(Categorical).alias("chrom"),
    ).drop(["chrom1", "chrom2", "variable"])

    chunk = chunk.with_columns(
        # (pl.col("chrom1") + pl.col("variable").str.slice(0, 1) + pl.lit(" ") + pl.col("chrom2") + pl.col("variable").str.slice(1, 1)).alias("chrom"),
        pl.col("value").list.eval(pl.element().len()).alias("count").map_elements(lambda x: x[0]),
        pl.col("value").list.eval(pl.element().str.join(" ")).alias("value").map_elements(lambda x: x[0]),

    ).select(["chrom", "count", "value"])


    chunk.write_csv(f"{tmp_dir}/{chunk_name}.clm", include_header=False, separator="\t")
    
    # return chunk

def process_chunk_cis_depth(chunk, window_size, min_mapq, is_with_mapq):
    os.environ["POLARS_MAX_THREADS"] = "1"
    chunk = pl.scan_parquet(chunk).select(["chrom1", "pos1", "chrom2", "pos2", "mapq"])
    if is_with_mapq:
        if min_mapq > 1:
            chunk = chunk.filter(pl.col("mapq") >= min_mapq)
        else:
            chunk = chunk.drop("mapq")

    chunk = chunk.filter(pl.col("chrom1") == pl.col("chrom2")).drop("chrom2")

    chunk = chunk.with_columns(
        (pl.when(pl.col('pos1') > pl.col('pos2'))
            .then(pl.col('pos2'))
            .otherwise(pl.col('pos1'))).alias('pos1'),
        (pl.when(pl.col('pos1') > pl.col('pos2'))
            .then(pl.col('pos1'))
            .otherwise(pl.col('pos2'))).alias('pos2')
    )

    chunk = chunk.with_columns(
        (pl.col('pos1') - 1) // window_size,
        (pl.col('pos2') - 1) // window_size,
    )

    return chunk


def process_chunk_split_pairs(i, chunk, q0, q1, columns, schema, is_pairs_with_mapq):
    os.environ["POLARS_MAX_THREADS"] = "1"
    chunk = pl.read_csv(
        chunk,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        new_columns=PAIRS_HEADER,
        schema=schema,
        columns=columns,
        truncate_ragged_lines=True,
    )
    chunk = chunk.cast(schema)
    chunk.write_parquet(q0 / f"{i}.parquet")
    if is_pairs_with_mapq:
        chunk = chunk.filter(pl.col("mapq") >= 1)
    chunk.write_parquet(q1 / f"{i}.parquet")