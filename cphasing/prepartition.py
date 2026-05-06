#!/usr/bin/env python
# -*- coding:utf-8 -*-


import argparse
import logging
import os
import os.path as op
import sys
import subprocess


import pandas as pd
import polars as pl 
import numpy as np

from collections import defaultdict
from pathlib import Path
from subprocess import PIPE, Popen

from .utilities import cmd_exists, run_cmd

logger = logging.getLogger(__name__)


class PrePartition:
    PAF_HADER = ["contig1", "length1", "start1", "end1", "strand",
                 "contig2", "length2", "start2", "end2", "mismatch",
                 "matches", "mapq", "identity"]
    META = {
        "contig1": "object",
        "length1": int,
        "start1": int,
        "end1": int,
        "strand": "category",
        "contig2": "object",
        "length2": int,
        "start2": int,
        "end2": int,
        "mismatch": int,
        "matches": int,
        "mapq": int,
    }
    def __init__(self, ref, query, aligner='wfmash', 
                 min_mapq=1, min_identity=0.90,
                 force=False, 
                 output="prepartition.clusters.txt",
                 threads=10,
                 log_dir="logs"):
        self.ref = Path(ref)
        self.query = Path(query)

        self.prefix_ref = Path(ref).stem
        self.prefix_query = Path(query).stem 

        self.min_mapq = min_mapq 
        self.min_identity = min_identity

        self.aligner = aligner
        self.threads = threads
        self.force = force
        self.output = output
        self.paf = f"{self.prefix_ref}.{self.prefix_query}.paf"

        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        if not cmd_exists(self.aligner):
            logger.error(f'No such command of `{self.aligner}`')
            sys.exit()

        


    def mapping(self):

        if self.force is False:
            logger.debug("Force is False.")
            if Path(self.paf).exists() and Path(self.paf).stat().st_size > 0:
                logger.warning(f"Using existing mapping results: `{self.paf}`")
                return self.paf 
        
        logger.info("Mapping ...")

        if self.aligner == 'wfmash':
            cmd = ["wfmash", str(self.ref), str(self.query), 
                        "-N", 
                        "-n", '1',
                        "-m",
                        "-t", str(self.threads)]

        else:
            logger.error(f"Unsupported aligner: `{self.aligner}`")
            sys.exit()

        pipelines = []
        try:
            pipelines.append(
                Popen(cmd, stdout=open(self.paf, "w"),
                      stderr=open(f"{self.log_dir}/ref.align.log", "w"),
                      bufsize=-1)
            )
            pipelines[-1].wait()
        except:
            raise Exception('Failed to execute command:' 
                                f'\t{cmd}')
        finally:
            for p in pipelines:
                if p.poll() is None:
                    p.terminate()
                else:
                    assert pipelines != [], \
                        "Failed to execute command, please check log."
                if p.returncode != 0:
                    raise Exception('Failed to execute command:' 
                                        f'\t{" ".join(cmd)}')

        return self.paf
    

    def read_paf(self):
        logger.info(f"Load alignments results `{self.paf}`")
        
        try:
            df = pd.read_csv(self.paf, sep='\s+', header=None, usecols=list(range(13)),
                             index_col=None, names=self.PAF_HADER[:13])
            
            df['identity'] = df['identity'].map(lambda x: x.replace("id:f:", "") if isinstance(x, str) and x.startswith("id:f:") else x)
            df['identity'] = pd.to_numeric(df['identity'], errors='coerce')
            
            df = df.dropna().astype(self.META)
            
        except pd.errors.ParserError:
            df = pd.read_csv(self.paf, sep='\s+', header=None, usecols=list(range(13)),
                             index_col=None, names=self.PAF_HADER[:13])
            
            df['identity'] = pd.to_numeric(df['identity'], errors='coerce')
            df = df.dropna().astype(self.META)

        df = df[df['identity'] > self.min_identity]
        df = df[df['mapq'] >= self.min_mapq]

        self.paf_df = df

        return df
    
    def cluster(self):
        logger.info("Resolving multimapping and clustering contigs...")
        
        if self.paf_df.empty:
            logger.warning("No alignments passed the filtering criteria.")
            return None

        best_idx = self.paf_df.groupby('contig1')['matches'].idxmax()
        best_df = self.paf_df.loc[best_idx]
        
        cluster_dict = {}
        idx = 1
        for group, df_group in best_df.groupby("contig2"):
            cluster_dict[idx] = df_group['contig1'].tolist()
            idx += 1
        from .core import ClusterTable
        clustertable = ClusterTable.from_dict(cluster_dict, self.output)

        return clustertable
    

    def run(self):
        self.mapping()
        self.read_paf()
        clustertable = self.cluster()
        clustertable.save(self.output)