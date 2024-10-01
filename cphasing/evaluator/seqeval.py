#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
evaluate the sequences compontents
"""


import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

from collections import OrderedDict, defaultdict
from pandarallel import pandarallel
from pathlib import Path
from pyranges import PyRanges
from subprocess import PIPE, Popen

from ..utilities import (cmd_exists,
                         get_contig_size_from_fasta,
                         read_chrom_sizes)


logger = logging.getLogger(__name__)


class HomoAnalysis:
    PAF_HADER = ["contig1", "length1", "start1", "end1", "strand",
                 "contig2", "length2", "start2", "end2", "mismatch",
                 "matches", "mapq", "identity"]
    def __init__(self, fasta, 
                min_similarity=0.998,
                aligner='wfmash',
                threads=10,
                log_dir="logs",
                force=False):

        self.fasta = fasta 
        self.fasta_prefix = Path(fasta).stem
        self.paf = f"{self.fasta_prefix}.selfalign.paf"
        self.contigsizes = read_chrom_sizes(str(get_contig_size_from_fasta(self.fasta)))
        self.total_length = self.contigsizes['length'].sum()
        
        self.min_similarity = min_similarity

        self.threads = threads
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.aligner = aligner
        if not cmd_exists(self.aligner):
            logger.error(f'No such command of `{self.aligner}`.')
            sys.exit()

        self.force = force 

    def mapping(self):

        if self.force is False:
            logger.debug("Force is False")
            if Path(self.paf).exists():
                logger.warning(f"Using existing mapping results: `{self.paf}`")
                return self.paf

        logger.info("Mapping ...")
        cmd = ["wfmash", "-t", str(self.threads), str(self.fasta), "-m"]

        pipelines = []
        try:
            pipelines.append(
                Popen(cmd, stdout=open(self.paf, "w"),
                      stderr=open(f"{self.log_dir}/seqeval.self.align.log", "w"),
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
    

    def mapping2(self):

        if self.force is False:
            logger.debug("Force is False")
            if Path(self.paf).exists():
                logger.warning(f"Using existing mapping results: `{self.paf}`")
                return self.paf

        logger.info("Mapping ...")
        cmd = ["minigraph", "-t", str(self.threads), "-DP", "-cxasm", str(self.fasta), str(self.fasta)]

        pipelines = []
        try:
            pipelines.append(
                Popen(cmd, stdout=open(self.paf, "w"),
                      stderr=open(f"{self.log_dir}/seqeval.self.align.log", "w"),
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
        df = pd.read_csv(self.paf, sep='\t', header=None, usecols=range(13),
                         names=self.PAF_HADER, index_col=None)
        df['identity'] = df['identity'].map(lambda x: x.replace("id:f:", "")).astype('float64')
        df = df[df['contig1'] != df['contig2']]
    
        # df = df.sort_values(['contig2', 'start2'])

        self.paf_df = df[df['identity'] >= self.min_similarity]

        return df 
    
    def read_paf2(self):
        logger.info(f"Load alignments results `{self.paf}`")
        df = pd.read_csv(self.paf, sep='\t', header=None, usecols=range(13),
                         names=self.PAF_HADER, index_col=None)
        df['identity'] = df['mismatch'] / df['matches']


        # df['identity'].map(lambda x: x.replace("id:f:", "")).astype('float64')
        df = df[df['contig1'] != df['contig2']]
    
        # df = df.sort_values(['contig2', 'start2'])

        self.paf_df = df[df['identity'] >= self.min_similarity]
   
        return self.paf_df
    
    def calcualte_homo(self):
        
        df1 = self.paf_df[['contig1', 'start1', 'end1']]
        df1.columns = ['Chromosome', 'Start', 'End']
        df2 = self.paf_df[['contig2', 'start2', 'end2']]
        df2.columns = ['Chromosome', 'Start', 'End']

        df = pd.concat([df1, df2], axis=0)
        df.sort_values(by=['Chromosome', 'Start'], inplace=True)

        pandarallel.initialize(nb_workers=self.threads, progress_bar=False, verbose=0)

        def func(df):
            if len(df) <= 1:
                return df 
            
            return PyRanges(df).merge().df

        df = df.groupby(['Chromosome'], sort=False).parallel_apply(func).reset_index(drop=True)

        contigsizes = self.contigsizes.to_dict()['length']

        df['region_length'] = df['End'] - df['Start'] + 1

        df = df.groupby('Chromosome', as_index=False).agg({'region_length': 'sum'})

        df['contig_length'] = df['Chromosome'].map(contigsizes.get)
        df['coverage'] = df['region_length'] / df['contig_length']

        
        self.total_highly_homologous = df['region_length'].sum()
         
        print(f"all: {self.total_highly_homologous / self.total_length:.2%}")

        for i in range(1, 11, 1):
            print(f"coverage ({i/10.0}):{df[df['coverage'] >= i/10.0]['region_length'].sum() / self.total_length:.2%}")
      
    def run(self):

        if self.aligner == 'minigraph':
            self.mapping2()
            self.read_paf2()
        else:
            self.mapping()
            self.read_paf()

        self.calcualte_homo()
