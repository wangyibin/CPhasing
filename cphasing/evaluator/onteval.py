#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
evaluate the asssembly by ultra-long reads
"""

import argparse
import logging
import gc
import os
import os.path as op
import sys


import bisect
import pandas as pd
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, MaxNLocator
import seaborn as sns

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from joblib import Parallel, delayed
from pandarallel import pandarallel
from pathlib import Path
from pyranges import PyRanges

# from line_profiler import LineProfiler

from ..agp import import_agp
from ..utilities import listify, xopen, list_flatten

logger = logging.getLogger(__name__)


def get_file_type(input_file):
    with xopen(input_file) as fp:
        for line in fp:
            if line.startswith(">"):
                return "fasta"
            elif line.startswith("@"):
                return "fastq"
            else:
                return

class GapEvaluator:
    """
    Evaluate assembly by ont data at gap region (developing)
    """

    def __init__(self, agp):
        self.agp = agp
        self.agp_df, self.gap_df = import_agp(self.agp)
        self.gap_num = len(self.gap_df)

    def get_gap_bed(self):
        self.gap_df
        

    def run(self):
        pass 


@dataclass
class LIS:
    # chrom: str = ""
    start: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.int32))
    # end: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.int32))
    split_idx: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.int32))
    mapq: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.int8))
    nm: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.int32))
    fragment_length: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.int32))

    def __len__(self):
        return len(self.start)
    
    def pop(self, index):
        start_val = self.start[index]
        end_val = self.end[index]
        split_idx_val = self.split_idx[index]
        mapq_val = self.mapq[index]

        self.start = np.delete(self.start, index)
        # self.end = np.delete(self.end, index)
        self.split_idx = np.delete(self.split_idx, index)
        self.mapq = np.delete(self.mapq, index)
        self.nm = np.delete(self.nm, index)
        self.fragment_length = np.delete(self.fragment_length, index)

        return (start_val, end_val, split_idx_val, mapq_val)

    def append(self, start, #end, 
               split_idx, mapq, nm, fragment_length):
        self.start = np.append(self.start, start)
        # self.end = np.append(self.end, end)
        self.split_idx = np.append(self.split_idx, split_idx)
        self.mapq = np.append(self.mapq, mapq)
        self.nm = np.append(self.nm, nm)
        self.fragment_length = np.append(self.fragment_length, fragment_length)

    def count_by_mapq(self, min_mapq=2):
        return np.sum(self.mapq >= min_mapq)

class LISContainer(list):
    raw_read_name = ""
    def __init__(self):
        super().__init__()
    
    def idxmax(self):

        # max_count_lis = self[0]
        max_idx = 0
        for idx, lis in enumerate(self[1:], 1):
            if self[max_idx].count_by_mapq() < lis.count_by_mapq():
                # max_count_lis = lis
                max_idx = idx
        
        return max_idx
    
    def join(self, min_mapq=2):
        max_idx = self.idxmax()
        max_lis = self[max_idx]
        
        split_idx = np.array(max_lis.split_idx)
        init_array = np.zeros(max(split_idx) + 1)
        
        non_dup = pd.Index(split_idx).duplicated(keep='first')
        split_idx = split_idx[~non_dup]
      
        init_array[split_idx] = np.array(max_lis.nm)[~non_dup]

        
        df = np.array(max_lis.fragment_length)[~non_dup]
        
        total_windows = len(split_idx)
        non_max_idx_list = []

        for idx, lis in enumerate(self):
            if idx == max_idx:
                continue
            non_max_idx_list.append(idx)

            split_idx2 = np.array(lis.split_idx)
            non_dup2 = pd.Index(split_idx2).duplicated(keep='first')
            split_idx2 = split_idx2[~non_dup2]
            _filter_split_idx2 = np.isin(split_idx2, split_idx)
            if not any(_filter_split_idx2):
                continue

            filter_split_idx2 = split_idx2[_filter_split_idx2]
            mapq2 = np.array(lis.mapq)[~non_dup2]

            filter_mapq2 = mapq2 >= min_mapq
            if not any(filter_mapq2):
                continue

            nm2 = np.zeros(max(split_idx2) + 1)
            nm2[filter_split_idx2] = np.array(lis.nm)[~non_dup2][_filter_split_idx2]
            
            filter_idx = split_idx2[_filter_split_idx2 & filter_mapq2]
            init_array[filter_idx] = - init_array[filter_idx] + nm2[filter_idx]

        error_idx = init_array < 0
        error_array = init_array[error_idx]
    
        if len(error_array) == 0:
            return total_windows, [np.nan]
        
        true_idx = np.where(error_idx)[0]
        true_idx = np.where(np.isin(true_idx, split_idx))[0]
      
        window_length = df[true_idx]
    
        error_rate = np.abs(error_array) / window_length
        
        return total_windows, error_rate

class SwitchError:
    PAF_HEADER = [
        'read_name', 'read_length', 
        'read_start', 'read_end', 
        'strand', 'chrom', 
        'chrom_length', 'start', 
        'end', 'matches', 
        'align_length', 'mapping_quality',
        "NM"]

    PAF_HEADER2 = [
        'read_name', 
        'read_start', 'read_end', 
        'strand', 'chrom', 
        'start', 'end',
        'mapping_quality',
        "NM"]
    def __init__(self, fasta, reads, outprefix="output",
                 window=5000, min_windows=2, 
                 maximum_gap=1000000,
                 chunksize=1e6,
                 threads=10, is_hifi=False):
        
        self.fasta = fasta 
        self.reads = reads 

        self.outprefix = outprefix
        paf = f"{self.outprefix}.paf"
        self.paf_path = Path(paf)
        self.window = window
        self.min_windows = min_windows
        self.maximum_gap = maximum_gap
        self.chunksize = chunksize
        self.threads = threads
        self.is_hifi = is_hifi
    
    def mapping(self):
        checkMinimap2CMD = "minimap2 --version"
        minimapVer = "".join(os.popen(checkMinimap2CMD).readlines())
        minimapVer = float(minimapVer.split('-')[0].split('.')[1])
        if minimapVer < float(24):
            print("Warnning: The minimap2 version should be 2.24 or higher.")
            sys.exit()
        preset = "map-ont" if not self.is_hifi else "map-hifi"

        if ',' in self.reads:
            readsLst = self.reads.split(',')
        else:
            readsLst = [self.reads]

        file_type = get_file_type(readsLst[0])
        if not file_type:
            logger.error("Input error: please input fasta or fastq.")
            sys.exit()

        if readsLst[0][-3:] == ".gz":
            cat_cmd = "pigz -p 4 -dc"
        else:
            cat_cmd = "cat"

        min_length = self.window * self.min_windows

        if len(readsLst) > 1:
            minimap2CMD = "{} {} | cphasing-rs slidefq - -w {} -l {} -f {}\
                                | minimap2 -I 16g -t {} --qstrand -cx {} \
                                -p.3 {} - > {}.paf".format(
                                    cat_cmd, ' '.join(readsLst), self.window, min_length, 
                                            file_type, self.threads, preset, self.fasta,  self.outprefix)
        else:
            minimap2CMD = "cphasing-rs slidefq {} -w {} -l {} -f {} \
                            | minimap2 -I 16g -t {} --qstrand -cx {} \
                            -p.3 {} - > {}.paf".format(
                                readsLst[0], self.window, min_length, file_type, self.threads, preset, self.fasta,  self.outprefix)
    
        logger.info("Running Command:")
        logger.info(f"\t\t{minimap2CMD}")
        flag = os.system(minimap2CMD)
        assert flag == 0, "Failed to run the `minimap2`"

        return self.outprefix + ".paf"

    def read_paf(self):
        self.paf_df = pd.read_csv(self.paf_path, sep='\t', header=None, index_col=None,
                            names=self.PAF_HEADER2, usecols=[0, 2, 3, 4, 5, 7, 8, 11, 12])
        
        self.paf_df['NM'] = self.paf_df['NM'].str.replace("NM:i:", "").astype(int)
        
        raw_read_names = self.paf_df['read_name'].str.rsplit("_", n=1, expand=True)
        self.paf_df.drop(['read_name'], axis=1, inplace=True)
        self.paf_df['raw_read_name'] = raw_read_names[0].astype('category').cat.codes
        
        self.paf_df['split_idx'] = raw_read_names[1].astype(int)
        self.paf_df['fragment_length'] = self.paf_df['end'] - self.paf_df['start']
    
        self.paf_df = self.paf_df[self.paf_df['read_end'] - self.paf_df['read_start'] > self.window * 0.9]

        return self.paf_df

    def read_paf_by_chunk(self):
        chunks = pd.read_csv(self.paf_path, sep='\t', header=None, index_col=None,
                            names=self.PAF_HEADER2, usecols=[0, 2, 3, 4, 5, 7, 8, 11, 12],
                            chunksize=self.chunksize)
        
        df1 = chunks.__next__()

        df1['NM'] = df1['NM'].str.replace("NM:i:", "").astype(int)
        
        raw_read_names = df1['read_name'].str.rsplit("_", n=1, expand=True)
        df1.drop(['read_name'], axis=1, inplace=True)
        df1['raw_read_name'] = raw_read_names[0] # .astype('category').cat.codes
        
        df1['split_idx'] = raw_read_names[1].astype(int)
        df1['fragment_length'] = df1['end'] - df1['start']
    
        df1 = df1[df1['read_end'] - df1['read_start'] > self.window * 0.9]

        for df2 in chunks:
            df2['NM'] = df2['NM'].str.replace("NM:i:", "").astype(int)
        
            raw_read_names = df2['read_name'].str.rsplit("_", n=1, expand=True)
            df2.drop(['read_name'], axis=1, inplace=True)
            df2['raw_read_name'] = raw_read_names[0] #.astype('category').cat.codes
            
            df2['split_idx'] = raw_read_names[1].astype(int)
            df2['fragment_length'] = df2['end'] - df2['start']
        
            df2 = df2[df2['read_end'] - df2['read_start'] > self.window * 0.9]
            another_df = df2.loc[df2['raw_read_name'] == df1.tail(1)['raw_read_name'].values[0]]

            df1 = pd.concat([df1, another_df], axis=0)
            df1['raw_read_name'] = df1['raw_read_name'].astype('category').cat.codes
            yield df1 
            df1 = df2
        else:
            df1['raw_read_name'] = df1['raw_read_name'].astype('category').cat.codes
            yield df1 

    
    def find_lis(self, input_df):
        def sub_func(tmp_df):
            lis = LIS()
            # chrom = tmp_df['chrom'].values[0]
            popular_strand = tmp_df.groupby('strand')['raw_read_name'].count().idxmax()

            ascending = True if popular_strand == '+' else False
            tmp_df.sort_values(by=['split_idx'], inplace=True, ascending=ascending)

            # lis.chrom = chrom
            starts = tmp_df['start'].values
            # ends = tmp_df['end'].values
            split_idxs = tmp_df['split_idx'].values
            mapqs = tmp_df['mapping_quality'].values
            nms = tmp_df['NM'].values
            fragment_lengths = tmp_df['fragment_length']
    
        
            for start, split_idx, mapq, nm, fragment_length in zip(
                starts, split_idxs, mapqs, nms, fragment_lengths):

                pos = bisect.bisect_left(lis.start, start)
                if pos == len(lis.start):
                    if lis.start.size > 0 and (start - lis.start[(pos - 1)]) > self.maximum_gap:
                        continue

                    lis.append(start, split_idx, mapq, nm, fragment_length)
            
            return lis
    
        def func(df):
            res = LISContainer()
            res.raw_read_name = df['raw_read_name'].values[0]
            res.extend(df.groupby('chrom', sort=False).apply(sub_func).values.tolist())
            
            return res

        
        pandarallel.initialize(nb_workers=self.threads, verbose=0)
        res = input_df.groupby('raw_read_name', sort=False).parallel_apply(func)


        def process(lis_list):
            window_num, error_rate = lis_list.join()
            return window_num, error_rate
        
        args = []
        window_nums = []
        error_rates = []
        for read, lis_list in res.items():
            if isinstance(lis_list, pd.Series):
                continue
            if len(lis_list) == 1:
                window_nums.append(len(lis_list[0]))
                continue 
            # args.append(lis_list)
            window_num, error_rate = process(lis_list)
            window_nums.append(window_num)
            error_rates.append(error_rate)

        # res = Parallel(n_jobs=self.threads)(
        #     delayed(process)(i) for i in args
        # )
        # window_nums, error_rates = list(zip(*res))
        
        error_rates = list(filter(lambda x: ~np.isnan(x), list_flatten(error_rates)))

        return window_nums, error_rates 
    
    def calculate(self, window_nums, error_rates):
        total_window_num = sum(list_flatten(window_nums))
        # error_rates = list(filter(lambda x: ~np.isnan(x), list_flatten(error_rates)))
        error_rates = list_flatten(error_rates)
        error_rate = np.median(error_rates)
        error_rate = error_rate if not np.isnan(error_rate) else 0

        print(f"Total window counts: {total_window_num}", file=sys.stdout)
        print(f"Error window counts: {len(error_rates)}", file=sys.stdout)
        print(f"Hamming error rate: {len(error_rates)/total_window_num:.4%}", file=sys.stdout)
        print(f"Switch error rate: {error_rate:.4%}", file=sys.stdout)

    def run(self):

        if not Path(self.paf_path).exists():
            self.mapping()
        else:
            logger.warning(f"Using existed mapping results: `{self.paf_path}`")
        # self.read_paf()
        window_nums = []
        error_rates = []
        # for df in self.read_paf_by_chunk():

        #     _window_nums, _error_rates = self.find_lis(df)
        #     window_nums.extend(_window_nums)
        #     error_rates.extend(_error_rates)

        res = Parallel(n_jobs=self.threads)(delayed(self.find_lis)(df) for df in self.read_paf_by_chunk())

        
        window_nums, error_rates = list(zip(*res))
        
        self.calculate(window_nums, error_rates)
        

def get_mode(arr):
    values, counts = np.unique(arr, return_counts=True)
    mode_index = np.argmax(counts)
    return values[mode_index]

class ComponentAnalysis:
    """
    Analysis the components of contigs, including junk, false duplication, haploid, collapsed.
    """
    DEPTH_HEADER = ["chrom", "start", "end", "count"]
    def __init__(self, depth, outprefix):
        self.depth_path = depth
        self.outprefix = outprefix

    def read_depth(self):
        logger.info(f"Load depth file: `{self.depth_path}`")
        df = pd.read_csv(self.depth_path, sep='\t', header=None, 
                            names=self.DEPTH_HEADER, index_col=None)
        df.columns = ['chrom', 'start', 'end', 'count']

        popular_window_size = get_mode(df['end'] - df['start'])
        
        df = df[(df['end'] - df['start']) == popular_window_size]

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

        self.xxx = xxx 
        self.yyy = yyy 

        peak_value = xxx[np.argmax(yyy)]

        return peak_value

    def category(self):
        def func(cn):
            
            if cn <= 0.25:
                return 0
            elif 0.25 < cn <= 0.5:
                return 1
            elif 0.5 < cn <= 1.5:
                return 2
            else:
                return 3
        
        self.depth_df['type'] = self.depth_df['CN'].map(func)
    
    def analysis(self):

        # for type, tmp_df in self.depth_df.groupby('type'):
            # tmp_df2 = tmp_df[['chrom', 'start', 'end']]
            # tmp_df2.columns = ['Chromosome', 'Start', 'End']

            # tmp_df2 = PyRanges(tmp_df2).merge().df
        
        res_df = self.depth_df.groupby('type')['CN'].count()
        total_count = res_df.sum()
        
        res_df = res_df / total_count 

        res_df = res_df.rename(index={0: 'Erroneous', 1: 'Duplicated',
                             2: 'Haploid', 3: 'Collapsed'})
        
        return res_df

    def plot(self, output):

        fig, ax = plt.subplots(figsize=(5.7, 5))
        plt.rcParams['font.family'] = 'Arial'
        plt.plot(self.xxx, self.yyy, '-', label='Depth', color="#209093")
        plt.ylim(0)
        plt.xlim(0)
        plt.axvline(x=self.peak, linewidth=1, color='#b02418', linestyle='--')
        plt.text(self.peak, ax.get_ylim()[1] * 0.98 , f"{str(self.peak)}",
                fontsize=14, color='k', )
        plt.axvline(x=self.peak / 2, linewidth=1, color='#253761', linestyle='--')
        plt.text(self.peak / 2, ax.get_ylim()[1] * 0.98 , f"{self.peak/ 2:.1f}",
                fontsize=14, color='k', )
        plt.axvline(x=self.peak * 2, linewidth=1, color='#253761', linestyle='--')
        plt.text(self.peak * 2, ax.get_ylim()[1] * 0.98 , f"{str(self.peak * 2)}",
                fontsize=14, color='k', )  

        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        plt.tick_params(which='both', width=1.5, length=5)
        sns.despine()

        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.ylabel("Frequency", fontsize=20)
        plt.xlabel("Coverage", fontsize=20)

       
        plt.savefig(f'{output}.hist.plot.png', dpi=600, bbox_inches='tight')
        plt.savefig(f'{output}.hist.plot.pdf', dpi=600, bbox_inches='tight')

    def run(self):
        self.depth_df = self.read_depth()
        self.peak = self.get_main_peak()
        self.depth_df['CN'] = self.depth_df['count'] / self.peak 

        self.category()
        component_df = self.analysis()
        component_df.to_csv(f"{self.outprefix}.component.tsv", sep='\t', header=None, index=True)

        self.plot(self.outprefix)
