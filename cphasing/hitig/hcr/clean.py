#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
clean low coverage contigs of a fasta
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np



from collections import OrderedDict
from pathlib import Path
from pyranges import PyRanges
from subprocess import PIPE, Popen

from ...utilities import (
    cmd_exists,
    read_fasta, 
    get_contig_size_from_fasta, 
    read_chrom_sizes
    )

logger = logging.getLogger(__name__)


def clean(fasta, low_coverage_df, output, break_pos=None):
    fa_db = read_fasta(fasta)
    low_coverage_contigs = set(low_coverage_df.index.tolist())

    if break_pos:
        break_pos_db = OrderedDict()
        with open(break_pos, 'r') as fp:
            for line in fp:
                if line.strip():
                    line_list = line.strip().split()
                    contig, positions = line_list 
                    positions = positions.split(",")
                    positions = list(map(int, positions))
                    if contig not in break_pos_db:
                        break_pos_db[contig] = []
                    for pos in positions:
                        break_pos_db[contig].append(pos)
    with open(output, 'w') as out:
        for contig in fa_db:
            if contig in low_coverage_contigs:
                if break_pos:
                    if contig.rsplit(":", 1)[0] in break_pos_db:

                        out.write(f">{contig}\n{fa_db[contig]}\n")
                continue

            out.write(f">{contig}\n{fa_db[contig]}\n")

    logger.info(f"Successful output low coverage contigs removed fasta in `{output}`")


class Clean:
    """
    remove false assemblies or false duplication, which false dupliaction will only reatin one copy.
    """
    PAF_HADER = ["contig1", "length1", "start1", "end1", "strand",
                 "contig2", "length2", "start2", "end2", "mismatch",
                 "matches", "mapq", "identity"]
    
    DEPTH_HEADER = ["chrom", "start", "end", "count"]
    def __init__(self, fasta, depth, 
                    max_cn=2,
                    min_coverage=0.5,
                    min_length=50000, ## minimum length of collapsed contigs
                    output=None, 
                    skip_remove=False,
                    threads=10, 
                    log_dir="logs",
                    force=False):
        self.fasta = fasta
        self.fasta_prefix = Path(fasta).stem
        self.paf = f"{self.fasta_prefix}.selfalign.paf"
        self.depth = depth 
        self.output = output
        self.force = force
        self.max_cn = max_cn 
        self.min_coverage = min_coverage 
        self.min_length = min_length

        self.contigsizes = read_chrom_sizes(str(get_contig_size_from_fasta(self.fasta)))

        self.skip_remove = skip_remove

        if not self.skip_remove:
            self.mapping()
            self.paf_df = self.read_paf()

        self.depth_df = self.read_depth()
        self.peak = self.get_main_peak()
        self.depth_df['CN'] = self.depth_df['count'] / self.peak 
        self.category_df = self.category()

        self.junk_contigs = set()
        self.remove_dup_contigs = set()

        self.junk_contig_length = 0
        self.remove_dup_contig_length = 0
        self.collapsed_contig_length = 0
        
        self.threads = threads
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        if not cmd_exists("wfmash"):
            logger.error(f'No such command of `wfmash`.')
            sys.exit()

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
                      stderr=open(f"{self.log_dir}/clean.self.align.log", "w"),
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
    
        df = df.sort_values(['contig2', 'start2'])

        self.paf_df = df 

        return df 

    def get_coverage(self):

        tmp_df = self.paf_df[['contig1', 'start1', 'end1']]
        tmp_df.columns = ['Chromosome', 'Start', 'End']

        tmp_df = PyRanges(tmp_df).merge().df

        tmp_df['Length'] = tmp_df['End'] - tmp_df['Start']
        tmp_df = tmp_df.groupby('Chromosome')['Length'].sum()

        tmp_db = tmp_df.to_dict()
        self.paf_df['coverage1'] = self.paf_df['contig1'].map(tmp_db.get) / self.paf_df['length1']


    def read_depth(self):
        logger.info(f"Load depth file: `{self.depth}`")
        df = pd.read_csv(self.depth, sep='\t', header=None, 
                            names=self.DEPTH_HEADER, index_col=None)
        df.columns = ['chrom', 'start', 'end', 'count']

        return df 

    def get_main_peak(self):
        
        self.depth_df['count'] = self.depth_df['count'].map(np.round)
        depthHist = self.depth_df.groupby(['count'])['chrom'].count()
        depthHist = depthHist[depthHist.index > 0]
        max_values = round(depthHist.argmax()  * 4)
        depthHist = depthHist[depthHist.index < max_values]
        depthHist = depthHist.to_dict()
        ## find wave vally
        xxxyyy = list(depthHist.items())
        xxxyyy.sort(key = lambda x:x[0])

        xxx = [x[0] for x in xxxyyy]
        yyy = [y[1] for y in xxxyyy]
        xxx = xxx[:max_values]
        yyy = yyy[:max_values]

        peak_value = xxx[np.argmax(yyy)]

        return peak_value

    def get_collapsed(self):
        min_coverage = self.min_coverage
        tmp_df = self.depth_df[(self.depth_df['CN'] > 1.5) & (self.depth_df['CN'] <= self.max_cn + 0.5)][['chrom', 'start', 'end']]
        tmp_df.columns = ['Chromosome', 'Start', 'End']
        tmp_df = PyRanges(tmp_df).merge().df
        contig_length_db = self.contigsizes.to_dict()['length']
        
        tmp_df['Length'] = tmp_df['Chromosome'].map(contig_length_db.get)
        tmp_df['coverage'] = (tmp_df['End'] - tmp_df['Start']) / tmp_df['Length']

        self.collapsed_contigs = set(tmp_df[tmp_df['coverage'] > min_coverage]['Chromosome'].values.tolist())

        self.collapsed_contig_length = self.contigsizes.reindex(self.collapsed_contigs).sum()

        logger.info(f"Total `{self.collapsed_contig_length}` collapsed contigs be found.")

        return self.collapsed_contigs 
        
    
    def remove_junk(self):

        junk_df = self.category_df[self.category_df['type'] == 0]

        self.junk_contigs = set(junk_df['chrom'].values.tolist())

        self.junk_contig_length += self.contigsizes.reindex(self.junk_contigs).sum()

        return self.junk_contigs
        
    
    def remove_dup(self):
        min_coverage = self.min_coverage 
        false_df = self.category_df[self.category_df['type'] == 1] 


        alignment_df = self.paf_df 
        alignment_df['coverage2'] = (alignment_df['end2'] - alignment_df['start2']) / alignment_df['length2']
        alignment_df.set_index('contig1', inplace=True)
        alignment_cache = {row.chrom: alignment_df.loc[row.chrom] 
                           for idx, row in false_df.iterrows() if row.chrom in alignment_df.index}
        for idx, row in false_df.iterrows():
            try:
                tmp_df = alignment_cache[row.chrom]
            except KeyError:
                continue 
            
            contig1 = row.chrom

            if isinstance(tmp_df, pd.DataFrame):
                if (tmp_df['coverage1'] <= min_coverage).any():
                    continue

                if (tmp_df['length1'] < tmp_df['length2']).any():
                    remove_contig = contig1 
                    remove_length = tmp_df['length1'].values[0]
                    if remove_contig not in self.remove_dup_contigs:
                        self.remove_dup_contig_length += remove_length
                        self.remove_dup_contigs.add(remove_contig)
                    
                else:
                    tmp_df2 = tmp_df[(tmp_df["length2"] < tmp_df["length1"]) 
                                        & (tmp_df['coverage2'] > min_coverage)]
                    for idx, row in tmp_df2.iterrows():
                        remove_contig = row['contig2']
                        remove_length = row['length2']
                        if remove_contig not in self.remove_dup_contigs:
                            self.remove_dup_contig_length += remove_length
                            self.remove_dup_contigs.add(remove_contig)
                           
                
            else:
                if tmp_df['coverage1'] <= min_coverage and tmp_df['coverage2'] <= min_coverage:
                    continue
                contig2 = tmp_df['contig2']
                if tmp_df['length1'] < tmp_df['length2']:
                    remove_contig = contig1
                    remove_length = tmp_df['length1']
                else:
                    remove_contig = contig2
                    remove_length = tmp_df['length2']

                if remove_contig not in self.remove_dup_contigs:
                    self.remove_dup_contig_length += remove_length
                    self.remove_dup_contigs.add(remove_contig)
               

        logger.info(f"Total `{self.remove_dup_contig_length}` false duplicated.")

        return self.remove_dup_contigs

    def category(self):
        
        def func(cn):
            
            if cn <= 0.25:
                return 0
            elif 0.25 < cn <= 0.75:
                return 1
            else:
                return 2
        
        self.cn_df = self.depth_df.groupby('chrom')['count'].mean()
        cn_db = self.cn_df.to_dict()
        self.cn_df = self.cn_df.to_frame()
        self.cn_df['CN'] = self.cn_df['count'] / self.peak
        
        self.depth_df['type'] = self.depth_df['CN'].map(func)

        tmp_df = self.depth_df.copy()
        tmp_df['length'] = (tmp_df['end'] - tmp_df['start'])
        length_db = tmp_df.groupby('chrom')['length'].sum().to_dict()
        
        tmp_df = tmp_df.groupby(['chrom', 'type'])['length'].sum().reset_index()
        tmp_df = tmp_df.loc[tmp_df.groupby('chrom')['length'].idxmax()].reset_index(drop=True)


        # res_df = self.depth_df.groupby(['chrom', 'type'])['CN'].count().reset_index()
        # self.category_df = res_df.loc[res_df.groupby('chrom')['CN'].idxmax()].reset_index(drop=True)
        self.category_df = tmp_df.drop(['length'], axis=1)
        self.category_df['CN'] = self.category_df['chrom'].map(cn_db.get)
        
        return self.category_df 
    
    def output_fasta(self):
        
        fasta_db = read_fasta(self.fasta)
        output = f"{self.fasta_prefix}.cleaned.fasta" if self.output is None else self.output
        with open(output, 'w') as out:
            for contig in fasta_db:
                if contig in self.remove_dup_contigs or contig in self.junk_contigs:
                    continue 
                if contig in self.collapsed_contigs:
                    out.write(f">{contig}_d2\n{fasta_db[contig]}")

                out.write(f">{contig}\n{fasta_db[contig]}")
            
        logger.debug("Output new fasta `{output}`")

    def run(self):
        
        self.get_coverage()
        if not self.skip_remove:
            self.remove_junk()
            self.remove_dup()
        self.get_collapsed()
        
        output_clean = f"{self.fasta_prefix}.remove.txt"
        
        if not self.skip_remove: 
            contigs = self.junk_contigs | self.remove_dup_contigs
            self.cn_df.loc[list(contigs)].to_csv(output_clean, sep='\t', header=None, index=True)
            logger.info(f"Output removed contigs information to `{output_clean}`")
        output_collapsed = f"{self.fasta_prefix}.collapsed.txt"
        self.cn_df.loc[list(self.collapsed_contigs)].to_csv(output_collapsed, sep='\t', header=None, index=True)
        logger.info(f"Output collapsed contigs information to `{output_collapsed}`")

        self.output_fasta()