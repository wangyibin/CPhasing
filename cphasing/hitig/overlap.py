#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
get overlap regions 
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

from pathlib import Path 
from pyranges import PyRanges
from cphasing.utilities import ( cmd_exists,
                                 get_contig_length,
                                 get_contig_size_from_fasta,
                                 xopen
                                 )

from line_profiler import profile

logger = logging.getLogger(__name__)

class OverlapFinder:
    """
    Get overlap regions from gfa
    """
    def __init__(self, gfa):
        self.gfa = gfa

    def parse_gfa(self):
        length_db = {}
        rd_db = {}
        overlap_db = []
        with xopen(self.gfa) as f:
            for line in f:
                if line.startswith("L"):
                    parts = line.strip().split("\t")
                    contig1, starnd1, contig2, strand2, overlap = parts[1], parts[2], parts[3], parts[4], parts[5]
                    overlap = overlap.strip("M")
                    overlap = int(overlap)

                    overlap_db.append((contig1, starnd1, contig2, strand2, overlap))
                elif line.startswith("S"):
                    parts = line.strip().split("\t")
                    contig, _, length, rd = parts[1], parts[2], parts[3], parts[4]
                    length = int(length.split(":")[-1])

                    rd = int(rd.split(":")[-1])
                    length_db[contig] = length
                    rd_db[contig] = rd
    
   

class OverlapFinder2:
    """
    Get overlap regions from selfalign
    """
    PAF_HADER = ["contig1", "length1", "start1", "end1", "strand",
                 "contig2", "length2", "start2", "end2", 
                 "matches", "aligns", "mapq", "dv"]
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
        "matches": int,
        "aligns": int,
        "mapq": int,
        "dv": "object"
    }
    def __init__(self, fasta, 
                 aligner="minigraph",
                 log_dir="log",
                 threads=10):
        self.fasta = Path(fasta )

        self.prefix = Path(fasta).stem
        self.aligner = aligner
        self.threads = threads
        self.contigsizes = get_contig_length(fasta)
        self.contigsizes = get_contig_size_from_fasta(self.fasta)

        self.paf = f"{self.prefix}.self.mapping.paf"
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        if not cmd_exists(self.aligner):
            logger.error(f'No such command of `{self.aligner}`')
            sys.exit()
    
    def mapping(self):
        
        if self.aligner == "minigraph":
            cmd = ["minigraph", "-xasm", "-N", str(100), 
                        "--secondary=yes",
                        str(self.fasta), str(self.fasta),
                        "-t", str(self.threads)]
            
    def read_paf(self):
        logger.info(f"Load alignments results `{self.paf}`")

        try:
            df = pd.read_csv(self.paf, sep='\s+', header=None, usecols=list(range(12)) + [16],
                             index_col=None, names=self.PAF_HADER)
            
            df = df.dropna().astype(self.META)

        except pd.errors.ParserError:
            df = pd.read_csv(self.paf, sep='\s+', header=None, usecols=list(range(12))  + [16],
                             index_col=None, names=self.PAF_HADER)
            
            df = df.dropna().astype(self.META)

        # df.drop(['length1', 'length2'], inplace=True, axis=1)
        df['dv'] = df['dv'].str.split(":").map(lambda x: x[-1]).astype(np.float64)
        df['identity'] = df['matches'] / df['aligns']
        
        df.query('dv < 0.01', inplace=True)
        df.query("identity > 0.85", inplace=True)

        return df
    
    @profile
    def get_overlap_regions(self):  
        min_length = 5000

        def func(row):
            return (((row['start1'] - min_length) < 0 or (row['length1'] - row['end1'] < min_length))
                     and 
                     ((row['start2'] - min_length) < 0 or (row['length2'] - row['end2'] < min_length))
                    )
        
        self.paf_df = self.paf_df[self.paf_df.apply(func, axis=1)]
        

        self.paf_df.loc[(self.paf_df['start1'] - min_length) < 0, 'start1'] = 0
        self.paf_df.loc[(self.paf_df['start2'] - min_length) < 0, 'start2'] = 0

        idx = self.paf_df.loc[(self.paf_df['length1'] - self.paf_df['end1']) < min_length].index
        self.paf_df.loc[idx, 'end1'] = self.paf_df.loc[idx, 'length1']

        idx = self.paf_df.loc[(self.paf_df['length2'] - self.paf_df['end2']) < min_length].index
        self.paf_df.loc[idx, 'end2'] = self.paf_df.loc[idx, 'length2']

        df1 = self.paf_df[['contig1', 'start1', 'end1']]
        df1.columns = ["Chromosome", "Start", "End"]
        df2 = self.paf_df[['contig2', 'start2', 'end2']]
        df2.columns = ["Chromosome", "Start", "End"]
        df = pd.concat([df1, df2], axis=0)
        df.to_csv(f"{self.prefix}.tmp.overlap.bed", sep="\t", index=False, header=False)

        cmd = (f"sort -k1,1V -k 2,2n {self.prefix}.tmp.overlap.bed | "
                f"bedtools merge -i - 2>/dev/null > {self.prefix}.tmp2.overlap.bed")

        os.system(cmd)
        
        cmd = (f"awk '{{print $1,0,$2}}' OFS='\\t' {self.contigsizes} | "
               f"bedtools subtract -a - -b {self.prefix}.tmp2.overlap.bed 2>/dev/null > {self.prefix}.overlap.bed")

        os.system(cmd)

        if Path(f"{self.prefix}.tmp.overlap.bed").exists():
            os.remove(f"{self.prefix}.tmp.overlap.bed")
        if Path(f"{self.prefix}.tmp2.overlap.bed").exists():
            os.remove(f"{self.prefix}.tmp2.overlap.bed")
       
        return f"{self.prefix}.overlap.bed"


    def run(self):
        self.paf_df = self.read_paf()
        self.overlap_df = self.get_overlap_regions()
        
