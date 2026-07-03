#!/usr/bin/env python
# -*- coding:utf-8 -*-


import logging
import gc
import os
import os.path as op
import shutil
import sys
import warnings
import io

import cooler 
import numpy as np
import pandas as pd
import polars as pl

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MaxNLocator
from joblib import Parallel, delayed
from rich.console import Console
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from sklearn.mixture import GaussianMixture
from collections import Counter, defaultdict, OrderedDict
from itertools import product
from scipy.sparse import triu
from pandarallel import pandarallel 
from pathlib import Path

try:
    from .algorithms.hypergraph import extract_incidence_matrix2, HyperGraph
    from .core import Tour
    from .utilities import (run_cmd, list_flatten, read_fasta,
                            get_contig_size_from_fasta, 
                            read_chrom_sizes, xopen)
    from .agp import agp2cluster, import_agp, agp2tour
    from .cli import build
    from .cli import utils
    agp2fasta = utils.get_command(utils, cmd_name='agp2fasta')
    agp_dup = utils.get_command(utils, cmd_name='agp-dup')
except ImportError:
    from cphasing.utilities import run_cmd, list_flatten


logger = logging.getLogger(__name__)

def find_haploid_peak(raw_counts,
                      left_trim_quantile=0.005,
                      right_trim_quantile=0.995,
                      max_multiplier=1.8,
                      smooth_sigma=1.2,
                      min_prom_frac=0.05,
                      min_width=2,
                      debug=False):

    data = pd.Series(raw_counts).dropna()

    data = data[data > 0]

    hi = data.quantile(right_trim_quantile)
    data = data[data <= hi]

    lo = data.quantile(left_trim_quantile)
    data = data[data >= lo]

    if len(data) == 0:
        raise ValueError("No data left after trimming")


    q75, q25 = np.percentile(data, [75, 25])
    iqr = q75 - q25 if (q75 - q25) > 0 else data.std()
    bin_width = 2 * iqr / (len(data) ** (1/3))
    if not np.isfinite(bin_width) or bin_width <= 0:
        bin_width = max(1, data.mean() / 50)
    bins = int((data.max() - data.min()) / bin_width) + 30
    bins = max(60, min(bins, 600))

    hist_counts, bin_edges = np.histogram(data, bins=bins)
    centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    smooth = gaussian_filter1d(hist_counts.astype(float), smooth_sigma)


    min_prominence = smooth.max() * min_prom_frac
    peaks, props = find_peaks(smooth,
                              prominence=min_prominence,
                              width=min_width)
    if len(peaks) == 0:

        mode_idx = np.argmax(smooth)
        peak_cov = centers[mode_idx]
    else:

        prom = props['prominences']
        widths = props['widths']
        score = prom * widths

        order = np.argsort(score)[::-1]
        candidate_indices = peaks[order]

        sorted_by_pos = np.sort(peaks)
        peak_cov = centers[candidate_indices[0]]
        if len(sorted_by_pos) >= 2:
            left_cov = centers[sorted_by_pos[0]]
            second_cov = centers[sorted_by_pos[1]]
            if left_cov < second_cov * 0.7:
                left_prom = prom[np.where(peaks == sorted_by_pos[0])[0][0]]
                if left_prom >= 0.5 * prom.max():
                    peak_cov = left_cov

        mean_cov = data.mean()
        if peak_cov > mean_cov * 1.2 and len(peaks) > 1:
            smaller = centers[peaks][centers[peaks] < peak_cov]
            if len(smaller):
                best_small = smaller[np.argmax([prom[np.where(peaks == p)[0][0]] for p in peaks if centers[p] in smaller])]
                peak_cov = best_small

    peak_cov = float(peak_cov)

    lower = peak_cov * 0.10
    upper = peak_cov * max_multiplier

    if debug:
        print(f"[DEBUG] peak={peak_cov:.2f} lower={lower:.2f} upper={upper:.2f} bins={bins}")

    return peak_cov, lower, upper

class CollapseFromDepth:
    """
    get collapsed contigs from read depth
    """
    def __init__(self, depth, cn_offset=0, peak_value=None):
        self.depth = depth
        self.depth_df = self.read_depth()
        self.cn_offset = cn_offset
        self.contigsizes = self.depth_df.groupby('chrom').agg({'end': 'max'})
        self.contigsizes.columns = ['length']
        self.contigsizes_db = self.contigsizes.to_dict()['length']
        self.peak_value = peak_value

    def read_depth(self):
        df = pd.read_csv(self.depth, sep='\t', header=None)
        df.columns = ['chrom', 'start', 'end', 'count']
        df['length'] = df['end'] - df['start']
        return df

    def plot_distribution(self, output="plot", lower_value=0.1, upper_value=1.75 ):
        data = self.depth_df['count']
        fig, ax = plt.subplots(figsize=(5.7, 5))
        plt.rcParams['font.family'] = 'Arial'

        data = data[data <= np.percentile(data, 98) * 1.5]
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))


        kdelines = sns.kdeplot(data, ax=ax, color='#209093', linewidth=2)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        formatter = plt.gca().get_yaxis().get_major_formatter()
        plt.gca().yaxis.set_major_formatter(formatter)
        plt.gca().yaxis.get_offset_text().set_fontsize(14)
        # plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        formatter = plt.gca().get_xaxis().get_major_formatter()
        plt.gca().xaxis.set_major_formatter(formatter)
        plt.gca().xaxis.get_offset_text().set_fontsize(14)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel("Coverage", fontsize=20)
        plt.ylabel("Density", fontsize=20)

        x = kdelines.lines[0].get_xdata()
        y = kdelines.lines[0].get_ydata()

        if self.peak_value is not None:
            peak_val = float(self.peak_value)
        else:
            peak_ind = find_peaks(y, distance=10)[0]
            median_value = np.quantile(data, .3)
            peak_ind = list(filter(lambda j: x[j] > median_value, peak_ind))
            if len(peak_ind) == 0:
                max_idx = np.argsort(x)[len(x)//2]
            else:
                max_idx = peak_ind[np.argmax(y[peak_ind])]
            peak_val = float(x[max_idx])


        # ax.fill_between((x[max_idx] * lower_value, x[max_idx] * upper_value), 
        #             0, ax.get_ylim()[1], alpha=0.5 , color='#bcbcbc')
       
        ax.axvline(peak_val, linestyle='--', color='#b02418')
        ax.text(int(peak_val), ax.get_ylim()[1] * 0.98, str(int(peak_val)), fontsize=14, color='#cb6e7f')    
        ax.axvline(x=peak_val * 2, linewidth=1, color='#253761', linestyle='--')
        ax.text(peak_val * 2, ax.get_ylim()[1] * 0.98, f"{peak_val * 2:.1f}", fontsize=14, color='#253761')

        sns.despine()
        plt.legend([], frameon=False)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        plt.tick_params(which='both', width=1.5, length=5)

        plt.savefig(f'{output}.kde.plot.png', dpi=600, bbox_inches='tight')
        plt.savefig(f'{output}.kde.plot.pdf', dpi=600, bbox_inches='tight')

        return peak_val, peak_val * lower_value, peak_val * upper_value
        
    def run(self):
        peak_depth, lower_depth, upper_depth = self.plot_distribution()

        rd_df = self.depth_df
        rd_df['CN'] = rd_df['count'] / peak_depth

        high_depth_threshold = peak_depth * 1.5 
        rd_df['is_high'] = rd_df['count'] > high_depth_threshold
        contig_stats = rd_df.groupby('chrom').agg({
            'count': 'median',
            'is_high': 'mean',
            'end': 'max'
        }).rename(columns={'is_high': 'high_coverage_ratio', 'end': 'length'})
        contig_stats['CN'] = contig_stats['count'] / peak_depth
        _contig_stats = contig_stats.copy()
        _contig_stats['CN'] = _contig_stats['CN'] - self.cn_offset
        collapsed_mask = (_contig_stats['CN'] > 1.5) & (_contig_stats['high_coverage_ratio'] > 0.5)
        collapsed_df = _contig_stats[collapsed_mask]
        # collapsed_df['CN'] = round(collapsed_df['CN'])


        contig_depth = rd_df.groupby('chrom')['count'].median().to_frame()
        contig_depth['CN'] = contig_depth['count'] / peak_depth


        # collapsed_df = contig_depth[contig_depth['CN'] > 1.5]
        
        collapsed_df['length'] = collapsed_df.index.map(self.contigsizes_db.get)

        
        collapsed_df = collapsed_df.reset_index()
        logger.info(f"Identified {collapsed_df['length'].sum():,} bp contigs with read depth >= {upper_depth:.2f}")
        collapsed_df.drop(columns=["length", "high_coverage_ratio"], inplace=True)

        collapsed_df.to_csv("contigs.collapsed.contig.list", sep='\t', header=False, index=False)
        collapsed_df['CN'] = round(collapsed_df['CN'])
        with open("contigs.collapsed.dup.contig.list", 'w') as out:
            for idx, row in collapsed_df.iterrows():
                contig, cn = row['chrom'], row['CN']
                for i in range(2, int(cn) + 1):
                    print(contig, f"{contig}_d{i}", sep='\t', file=out)

        low_df = contig_depth[contig_depth['CN'] < 0.1]
        low_df['length'] = low_df.index.map(self.contigsizes_db.get)
        low_df = low_df.reset_index()
        logger.info(f"Identified {low_df['length'].sum():,} bp contigs with read depth < {lower_depth:.2f}")
        low_df.drop(columns=["length"], inplace=True)
        low_df.to_csv("contigs.low.contig.list", sep='\t', header=False, index=False)
    

    def run_gmm(self):
        peak_depth, lower_depth, upper_depth = self.plot_distribution()
        
        rd_df = self.depth_df
        
        contig_depth = rd_df.groupby('chrom')['count'].mean().to_frame()
        
  
        contig_depth['length'] = contig_depth.index.map(self.contigsizes_db.get)
        valid_contigs = contig_depth[contig_depth['length'] >= 10000].copy() # 仅使用 >10kb 的 contig 进行训练

        X = valid_contigs['count'].values.reshape(-1, 1)
        
  
        gmm = GaussianMixture(n_components=2, covariance_type='full', random_state=42)
        gmm.fit(X)
        
        means = gmm.means_.flatten()
        labels = gmm.predict(X)
        
        collapsed_label = np.argmax(means)
        haploid_label = np.argmin(means)

        probs = gmm.predict_proba(X)
        
        valid_contigs['cluster'] = labels
        valid_contigs['prob_collapsed'] = probs[:, collapsed_label]
        

        valid_contigs['CN'] = valid_contigs['count'] / peak_depth
        

        is_collapsed = (valid_contigs['cluster'] == collapsed_label) & \
                        (valid_contigs['prob_collapsed'] > 0.8) & \
                        (valid_contigs['CN'] > 1.35) 

        collapsed_df = valid_contigs[is_collapsed].copy()
        
        collapsed_df = collapsed_df.reset_index()
        logger.info(f"Identified {collapsed_df['length'].sum():,} bp contigs as collapsed (GMM method)")

        output_df = collapsed_df[['chrom', 'count']].copy()
        output_df.to_csv("contigs.collapsed.contig.list", sep='\t', header=False, index=False)

        low_df = contig_depth[contig_depth['count'] / peak_depth < 0.1].copy()
        low_df = low_df.reset_index()
        logger.info(f"Identified {low_df['length'].sum():,} bp contigs with read depth < {lower_depth:.2f}")
        low_df[['chrom', 'count']].to_csv("contigs.low.contig.list", sep='\t', header=False, index=False)
  

class CollapseFromGfa:
    """
    get collapsed contigs from gfa 
    """
    def __init__(self, gfa, peak_value=None):
        self.gfa = gfa
        self.peak_value = peak_value

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
    
        self.length_db = length_db 
        self.rd_db = rd_db 
        self.overlap_db = overlap_db

    def split_depth_by_bin(self, binsize=10000):
        length_df = pd.DataFrame(self.length_db.items(), columns=["chrom", "length"])
        length_df.set_index('chrom', inplace=True)
        
        self.bins = cooler.util.binnify(length_df.iloc[:, 0], binsize)
    
    def get_depth_from_rd(self):

        self.depth_df = self.bins.copy()
        self.depth_df['count'] = self.depth_df['chrom'].map(self.rd_db.get)
        self.gfa.replace(".gfa", ".depth")
        self.depth_df.to_csv(self.gfa.replace(".gfa", ".depth"), sep='\t', header=False, index=False)
       
    def plot_distribution(self, output="plot", lower_value=0.1, upper_value=1.5 ):
        data = self.depth_df['count']

        if self.peak_value is not None:
            peak_cov = float(self.peak_value)
            lower_thr = peak_cov * lower_value
            upper_thr = peak_cov * upper_value
        else:
            peak_cov, lower_thr, upper_thr = find_haploid_peak(
                data,
                max_multiplier=upper_value
            )

        fig, ax = plt.subplots(figsize=(5.7, 5))
        plt.rcParams['font.family'] = 'Arial'

        data = data[data <= np.percentile(data, 98) * 1.5]
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))


        kdelines = sns.kdeplot(data, ax=ax, color='#209093', linewidth=2)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        formatter = plt.gca().get_yaxis().get_major_formatter()
        plt.gca().yaxis.set_major_formatter(formatter)
        plt.gca().yaxis.get_offset_text().set_fontsize(14)
        # plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        formatter = plt.gca().get_xaxis().get_major_formatter()
        plt.gca().xaxis.set_major_formatter(formatter)
        plt.gca().xaxis.get_offset_text().set_fontsize(14)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel("Coverage", fontsize=20)
        plt.ylabel("Density", fontsize=20)

        x = kdelines.lines[0].get_xdata()
        y = kdelines.lines[0].get_ydata()

        ax.axvline(peak_cov, ls='--', color='#b02418')
        ax.text(peak_cov, ax.get_ylim()[1]*0.95, f"{peak_cov:.1f}",
                color='#b02418', ha='center')
        # ax.axvspan(lower_thr, upper_thr, color='#ccc', alpha=0.25)
        ax.axvline(upper_thr, ls='--', color='#253761')
        ax.text(upper_thr, ax.get_ylim()[1]*0.95, f"{upper_thr:.1f}",
                color='#253761', ha='center')
        
        sns.despine()
        plt.legend([], frameon=False)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        plt.tick_params(which='both', width=1.5, length=5)

        plt.savefig(f'{output}.kde.plot.png', dpi=600, bbox_inches='tight')
        plt.savefig(f'{output}.kde.plot.pdf', dpi=600, bbox_inches='tight')

        
        # return int(x[max_idx]), x[max_idx] * lower_value, x[max_idx] * upper_value
       
        return peak_cov, lower_thr, upper_thr
        

    def run(self):
        self.parse_gfa()
        self.split_depth_by_bin()
        self.get_depth_from_rd()
        min_length = 50000
        peak_depth, lower_depth, upper_depth = self.plot_distribution()
        rd_df = pd.DataFrame(self.rd_db.items(), columns=['contig', 'count'])

        rd_df['length'] = rd_df['contig'].map(self.length_db.get)
        # rd_df = rd_df[rd_df['length'] >= min_length]

        rd_df['CN'] = rd_df['count'] / peak_depth

        collapsed_df = rd_df[rd_df['count'] >= upper_depth]
        logger.info(f"Identified {collapsed_df['length'].sum():,} bp contigs with read depth >= {upper_depth:.2f}")
        collapsed_df.drop(columns=["length"], inplace=True)
        collapsed_df.to_csv("contigs.collapsed.contig.list", sep='\t', header=False, index=False)
        ## generate duplication contig list for collapsed contigs
        collapsed_df['CN'] = round(collapsed_df['CN'])
        with open("contigs.collapsed.dup.contig.list", 'w') as out:
            for idx, row in collapsed_df.iterrows():
                contig, cn = row['contig'], row['CN']
                for i in range(2, int(cn) + 1):
                    print(contig, f"{contig}_d{i}", sep='\t', file=out)
        low_df = rd_df[rd_df['count'] < lower_depth]
        logger.info(f"Identified {low_df['length'].sum():,} bp contigs with read depth < {lower_depth:.2f}")
        low_df.drop(columns=["length"], inplace=True)
        low_df.to_csv("contigs.low.contig.list", sep='\t', header=False, index=False)


class CollapseContigs:
    def __init__(self, cool_path):
        self.cool_path = cool_path 
        self.cool = cooler.Cooler(self.cool_path)
        
        pass 

class CollapsedRescue:
    """
    Rescue the collapsed contigs into a group
    """
    def __init__(self, HG, clustertable, alleletable, contigsizes,
                    collapsed_contigs, allelic_similarity: float=.85,
                    min_contacts=5, min_cis_weight=20, min_weight=2,
                    only_unplaced_rescue=False, threads=10,
                    disable_conflict_check=True):
        
        self.HG = HG 
        self.clustertable = clustertable

        self.alleletable = alleletable 
        if self.alleletable is not None:
            self.alleletable.data = self.alleletable.data[
                self.alleletable.data['similarity'] >= allelic_similarity]

        self.contigsizes = contigsizes
        self.collapsed_contigs = collapsed_contigs
        self.allelic_similarity = allelic_similarity
        self.min_contacts = min_contacts
        self.min_cis_weight = min_cis_weight
        self.min_weight = min_weight
        self.only_unplaced_rescue = only_unplaced_rescue
        self.threads = threads
        self.disable_conflict_check = disable_conflict_check

        self.log_dir = "logs"
        Path(self.log_dir).mkdir(exist_ok=True)

        self.hap_groups = OrderedDict()
        for group_name, contigs in self.clustertable.data.items():
            if "g" in group_name and group_name.rsplit("g", 1)[-1].isdigit():
                hap = group_name.rsplit("g", 1)[0]
            else:
                hap = group_name
            if hap not in self.hap_groups:
                self.hap_groups[hap] = []
            self.hap_groups[hap].append(contigs)

        self.H = HG.incidence_matrix(min_contacts=self.min_contacts)
        
        self.vertices = self.HG.nodes 

    @property
    def vertices_idx(self):
        return dict(zip(self.vertices, 
                        range(len(self.vertices))))
    
    @property
    def idx_to_vertices(self):
        idx_to_vertices = dict(zip(range(len(self.vertices)), self.vertices))
        
        return idx_to_vertices

    @staticmethod
    def _rescue(A, ):
        pass

    def filter_hypergraph(self):
        
        A = HyperGraph.clique_expansion_init(self.H, P_allelic_idx=None, min_weight=self.min_weight)
        lengths = self.contigsizes.reindex(self.vertices).values.flatten()
        lengths[np.isnan(lengths)] = 0
        lengths_series = pd.Series(lengths, index=self.vertices)
        
        retain_idx1 = None # HyperGraph.filter_adjacency_matrix(A, self.vertices, lengths_series, invert=True)
        if retain_idx1 is None:
            retain_idx1 = np.arange(A.shape[0], dtype=np.int64)
        else:
            retain_idx1 = np.asarray(retain_idx1)
            if retain_idx1.dtype == bool:
                retain_idx1 = np.where(retain_idx1)[0].astype(np.int64)
            else:
                retain_idx1 = retain_idx1.astype(np.int64)

        dia = A.diagonal()
        raw_contig_counts = A.shape[0]
        retain_idx = np.where(dia > self.min_cis_weight)[0]
        dia = A.diagonal()
        raw_contig_counts = A.shape[0]
        retain_idx = np.where(dia > self.min_cis_weight)[0]
        retain_idx = np.intersect1d(retain_idx, retain_idx1)
        contig_counts = len(retain_idx)

        if len(retain_idx) < raw_contig_counts:
            A = A[retain_idx, :][:, retain_idx]
            self.H, _, _ = extract_incidence_matrix2(self.H, retain_idx)
            self.vertices = self.vertices[retain_idx]
            logger.info(f"Filtered {raw_contig_counts - contig_counts} contigs with less than {self.min_cis_weight} cis weight")

    def rescue1(self):
        vertices_idx = self.vertices_idx
        idx_to_vertices = self.idx_to_vertices
        new_cluster_data = {}
        for k, hap_group in enumerate(self.hap_groups):
            groups = self.hap_groups[hap_group]
            if len(groups) < 2:
                continue 
            contigs = list_flatten(groups)
            contigs_idx = list(filter(lambda x: x is not None, map(vertices_idx.get, contigs)))
            
            sub_old2new_idx = dict(zip(contigs_idx, range(len(contigs_idx))))
            if self.alleletable is None:
                sub_alleletable = None 
                P_allelic_idx = None
                P_allelic_sim = None
            else:
                sub_alleletable = self.alleletable.data[
                    self.alleletable.data[1].isin(contigs) & self.alleletable.data[2].isin(contigs)]
                sub_alleletable[1] = sub_alleletable[1].map(vertices_idx.get).map(sub_old2new_idx.get)
                sub_alleletable[2] = sub_alleletable[2].map(vertices_idx.get).map(sub_old2new_idx.get)
                sub_alleletable.dropna(axis=0, inplace=True)
                sub_alleletable[1] = sub_alleletable[1].astype(int)
                sub_alleletable[2] = sub_alleletable[2].astype(int)
                P_allelic_idx = [sub_alleletable[1], sub_alleletable[2]]
                sub_alleletable.set_index([1], inplace=True)

            sub_H, _, _ = extract_incidence_matrix2(self.H, contigs_idx)
            sub_A = self.HG.clique_expansion_init(sub_H, P_allelic_idx=P_allelic_idx, allelic_factor=0)
           
            # sub_allelic = set(map(tuple, sub_alleletable[[1, 2]].values.tolist()))

            groups_idx = list(map(lambda x: list(filter(lambda x: x is not None, map(vertices_idx.get, x))), groups))

            groups_new_idx = list(map(lambda x: list(map(sub_old2new_idx.get, x)), groups_idx))
            groups_new_idx_db = {}
            for i in range(len(groups_new_idx)):
                groups_new_idx_db.update(dict(zip(groups_new_idx[i], [i] * len(groups_new_idx[i]))))
            
            sub_collapsed_contigs = self.collapsed_contigs.reindex(contigs).dropna(axis=0)
            sub_collapsed_contigs.index = sub_collapsed_contigs.index.map(vertices_idx.get)
            sub_collapsed_contigs_new = sub_collapsed_contigs.copy()
            sub_collapsed_contigs_new.index = sub_collapsed_contigs_new.index.map(sub_old2new_idx.get)
            sub_collapsed_contigs_idx = sub_collapsed_contigs.index.tolist()
            sub_collapsed_contigs_idx_new = sub_collapsed_contigs_new.index.tolist()
            sub_collapsed_contigs_idx_new = [int(idx) for idx in sub_collapsed_contigs_idx_new if pd.notnull(idx)]
            sub_collapsed_contigs_cn_dict = sub_collapsed_contigs_new['CN'].to_dict()
            
            res = []
           
            # for j in sub_collapsed_contigs_idx_new:
            #     if j is None:
            #         continue
            #     if sub_alleletable is None:
            #         tmp_allelic_table = []
            #     else:
            #         try:
            #             tmp_allelic_table = sub_alleletable.loc[j]
            #             try:
            #                 tmp_allelic_table.set_index([2], inplace=True)
            #             except:
            #                 tmp_allelic_table.to_frame().set_index([2], inplace=True)

            #         except KeyError or AttributeError:
            #             tmp_allelic_table = []
            #     shared_similarity = np.zeros(len(groups))
            #     tmp_res = []
            #     for i in range(len(groups)):
            #         if i == groups_new_idx_db[j]:
            #             tmp_res.append(0)
            #         else:
            #             if len(tmp_allelic_table) > 1:
            #                 tmp = tmp_allelic_table.reindex(groups_new_idx[i]).dropna()
            #                 if len(tmp) > 1:
            #                     shared_similarity[i] = tmp['mzShared'].sum()

            #             tmp_res.append(sub_A[j, groups_new_idx[i]].mean())
                
            #     cn = int(sub_collapsed_contigs_cn_dict[j])
                
            #     # tmp_res_db = dict(zip(range(len(tmp_res)), tmp_res))
                
            #     # tmp_res_idx = sorted(tmp_res_db, key=lambda x: tmp_res_db[x], reverse=True)
            #     # for max_idx in tmp_res_idx[:(cn-1)]:
            #     #     if tmp_res_db[max_idx] > self.min_weight:
            #     #         groups_new_idx[max_idx].append(j)
            #     #         logger.debug(f"Rescued contig {idx_to_vertices[max_idx]} into group {max_idx}")

            #     # res.append(tmp_res)
            #     tmp_res_db = dict(zip(range(len(tmp_res)), tmp_res))
                
            #     tmp_res_idx = sorted(tmp_res_db, key=lambda x: tmp_res_db[x], reverse=True)
                
            #     added_count = 0
            #     original_group_idx = groups_new_idx_db.get(j) # Get original group index

            #     for max_idx in tmp_res_idx:
            #         if added_count >= (cn - 1):
            #             break
                    
            #         if tmp_res_db[max_idx] <= self.min_weight:
            #             continue
                        
            #         if max_idx == original_group_idx:
            #             continue
    
            #         groups_new_idx[max_idx].append(j)
            #         logger.debug(f"Rescued contig {idx_to_vertices[max_idx]} into group {max_idx}")
            #         added_count += 1

            #     res.append(tmp_res)
            for j in sub_collapsed_contigs_idx_new:
                if j is None: 
                    continue
                
                if sub_alleletable is None:
                    tmp_allelic_table = []
                else:
                    try:
                        tmp_allelic_table = sub_alleletable.loc[j]
                        try:
                            tmp_allelic_table.set_index([2], inplace=True)
                        except:
                            tmp_allelic_table.to_frame().set_index([2], inplace=True)

                    except KeyError or AttributeError:
                        tmp_allelic_table = []
                shared_similarity = np.zeros(len(groups))
                tmp_res = []
                for i in range(len(groups)):
                    if i == groups_new_idx_db[j]:
                        tmp_res.append(0)
                    else:
                        if len(tmp_allelic_table) > 1:
                            tmp = tmp_allelic_table.reindex(groups_new_idx[i]).dropna()
                            if len(tmp) > 1:
                                shared_similarity[i] = tmp['mzShared'].sum()

                        tmp_res.append(sub_A[j, groups_new_idx[i]].mean())
                
                cn = int(sub_collapsed_contigs_cn_dict[j])

                group_scores = np.array(tmp_res)
                
                original_g_idx = groups_new_idx_db.get(j)

                valid_indices = [idx for idx, s in enumerate(group_scores) 
                                 if idx != original_g_idx and s > self.min_weight]
                
                if not valid_indices:
                    continue

                valid_indices.sort(key=lambda idx: group_scores[idx], reverse=True)
                
                targets = []
                primary_score = group_scores[valid_indices[0]]
                
                for candidate_idx in valid_indices:
                    if len(targets) >= (cn - 1):
                        break
                    
                    score = group_scores[candidate_idx]
                    
                    if score < primary_score * 0.25: 
                        continue
                    
                    targets.append(candidate_idx)
                
                for t_idx in targets:
                    groups_new_idx[t_idx].append(j)
                    logger.debug(f"Rescued contig {idx_to_vertices[t_idx]} into group {t_idx} (Score: {group_scores[t_idx]} vs Best: {primary_score})")

                res.append(tmp_res)            

            new_groups = []
            for group_idx in groups_new_idx:
                tmp = list(map(lambda x: contigs_idx[x], group_idx))
                tmp = list(map(idx_to_vertices.get, tmp))
                new_groups.append(tmp)

            for k, group in enumerate(new_groups):
                new_cluster_data[f'{hap_group}g{k+1}'] = group 
        
        self.clustertable.data = new_cluster_data
       
        new_contigs = self.clustertable.contigs
        count = Counter(new_contigs)
        rescued_contigs = [k for k, v in count.items() if v > 1]
        
        with open(f"collapsed.rescue.contigs.list", 'w') as out:
            for contig in rescued_contigs:
                for i in range(2, count[contig] + 1):
                    print(f"{contig}\t{contig}_d{i}", file=out)
        rescued_length = self.contigsizes.reindex(rescued_contigs).sum().values[0]

        logger.info(f"Rescued {len(rescued_contigs)} contigs, total length: {rescued_length:,} bp")
        self.clustertable.save("collapsed.rescue.clusters.txt")

    @staticmethod
    def _run_single_kprune_static(hap_group, groups, vertices_idx, idx_to_vertices, min_weight, 
                                  alleletable_data, threads, sub_H, all_collapsed, kprune_workdir_str, log_dir):
        group_contigs = list_flatten(groups)
        involved_contigs = [c for c in group_contigs if c in vertices_idx]
        if len(involved_contigs) < 2:
            return None

        involved_idx = [vertices_idx[c] for c in involved_contigs]
        sub_old2new = dict(zip(involved_idx, range(len(involved_idx))))
        sub_new2old = {v: k for k, v in sub_old2new.items()}

       
        sub_raw_A = HyperGraph.clique_expansion_init(sub_H, min_weight=min_weight)

        sub_alleletable = alleletable_data[
            alleletable_data[1].isin(involved_contigs)
            & alleletable_data[2].isin(involved_contigs)
        ].copy()

        if sub_alleletable.empty:
            return None

        k_workdir = Path(kprune_workdir_str)
        tmp_allele_path = k_workdir / f"allele.{hap_group}.table"
        sub_alleletable.reset_index().reset_index().to_csv(tmp_allele_path, sep='\t', header=False, index=False)

        sub_involved_vertices = np.array([idx_to_vertices[sub_new2old[i]] for i in range(len(involved_idx))])
        df = HyperGraph.to_contacts(sub_raw_A, sub_involved_vertices, min_weight=min_weight, output=None)
        df = df.to_pandas()
        
        tmp_contacts_path = k_workdir / f"contacts.{hap_group}.table"
        df.to_csv(tmp_contacts_path, sep='\t', header=None, index=None)

        tmp_collapsed_list = k_workdir / f"collapsed.{hap_group}.list"
        with open(tmp_collapsed_list, 'w') as out:
            for contig in all_collapsed.intersection(involved_contigs):
                print(contig, file=out)

        prune_table_path = k_workdir / f"prune.{hap_group}.table"
        cmd = ["cphasing-rs", "kprune", str(tmp_allele_path), str(tmp_contacts_path), 
               str(prune_table_path), "-t", str(threads), "--whitelist", str(tmp_collapsed_list), "-p"]
        
        run_cmd(cmd, log=f"{log_dir}/collapsed.kprune.{hap_group}.log")

        return prune_table_path
    
    def rescue(self):
        """
        Rescue the collapsed contigs into a group with Re-assignment logic.
        """
        vertices_idx = self.vertices_idx
        idx_to_vertices = self.idx_to_vertices

        median_len = np.nanmedian(self.contigsizes.values)
        if np.isnan(median_len) or median_len == 0:
            median_len = 1000000

        placed_contigs = set()
        for groups in self.hap_groups.values():
            placed_contigs.update(list_flatten(groups))
        
        all_collapsed = set(self.collapsed_contigs.index)
        unplaced_collapsed = [
            c for c in self.collapsed_contigs.index 
            if c not in placed_contigs and c in vertices_idx
        ]

        P_allelic_idx_list = []
        P_weak_idx_list = []
        P_allelic_sim_list = []
        total_hap_groups = len(self.hap_groups)
        logger.info(f"Segmenting rescue into {total_hap_groups} homological groups to parallelize kprune...")
        kprune_workdir = Path("collapsed_kprune_workdir")
        kprune_workdir.mkdir(exist_ok=True)
    

        max_workers = min(self.threads, total_hap_groups, 8) 
        kprune_threads = max(1, self.threads // max_workers)
        if self.alleletable is not None:
            full_alleletable_data = self.alleletable.data.copy()
            results = Parallel(n_jobs=max_workers, backend="loky", require=None)(
                delayed(self._run_single_kprune_static)(
                    hap_group, 
                    groups, 
                    vertices_idx, 
                    idx_to_vertices,
                    self.min_weight, 
                    full_alleletable_data[
                        full_alleletable_data[1].isin(set(list_flatten(groups))) & 
                        full_alleletable_data[2].isin(set(list_flatten(groups)))
                    ].copy(),
                    kprune_threads, 
                    extract_incidence_matrix2(self.H, [vertices_idx[c] for c in list_flatten(groups) if c in vertices_idx])[0], 
                    all_collapsed, 
                    str(kprune_workdir), 
                    self.log_dir
                )
                for hap_group, groups in self.hap_groups.items()
            )
            for prune_table_path in results:
                if prune_table_path and prune_table_path.exists() and prune_table_path.stat().st_size > 0:
                    prune_table = pd.read_csv(prune_table_path, sep='\t', header=None)
                    prune_table[0] = prune_table[0].map(vertices_idx.get)
                    prune_table[1] = prune_table[1].map(vertices_idx.get)
                    
                    prune_table.dropna(subset=[0, 1], inplace=True)
                    prune_table[0] = prune_table[0].astype(int)
                    prune_table[1] = prune_table[1].astype(int)

                    allelic_table = prune_table[prune_table[6] == 0]
                    cross_allelic_table = prune_table[prune_table[6] == 1]

                    P_allelic_idx_list.append([allelic_table[0].values, allelic_table[1].values])
                    P_weak_idx_list.append([cross_allelic_table[0].values, cross_allelic_table[1].values])
                    P_allelic_sim_list.append(allelic_table[5].values)

        if kprune_workdir.exists():
            shutil.rmtree(kprune_workdir)

        if P_allelic_idx_list:
            global_allelic_rows = np.concatenate([x[0] for x in P_allelic_idx_list])
            global_allelic_cols = np.concatenate([x[1] for x in P_allelic_idx_list])
            P_allelic_idx = [global_allelic_rows, global_allelic_cols]
            
            global_weak_rows = np.concatenate([x[0] for x in P_weak_idx_list])
            global_weak_cols = np.concatenate([x[1] for x in P_weak_idx_list])
            P_weak_idx = [global_weak_rows, global_weak_cols]

            P_allelic_sim = np.concatenate(P_allelic_sim_list)
        else:
            P_allelic_idx = None
            P_weak_idx = None
            P_allelic_sim = None


        valid_placed = [c for c in placed_contigs if c in vertices_idx]
        all_involved = valid_placed + unplaced_collapsed
        involved_idx = [vertices_idx[c] for c in all_involved]
        
        sub_old2new = dict(zip(involved_idx, range(len(involved_idx))))
        sub_new2old = {v: k for k, v in sub_old2new.items()} 
        
        sub_H, _, _ = extract_incidence_matrix2(self.H, involved_idx)
        sub_raw_A = HyperGraph.clique_expansion_init(
            sub_H,
            min_weight=self.min_weight,
        )

        # if self.alleletable is None:
        #     P_allelic_idx = None
        #     P_weak_idx = None
        #     P_allelic_sim = None

        # else:
        #     sub_alleletable = self.alleletable.data[
        #         self.alleletable.data[1].isin(all_involved)
        #         & self.alleletable.data[2].isin(all_involved)
        #     ].copy()
        #     sub_alleletable[1] = (
        #         sub_alleletable[1].map(vertices_idx.get).map(sub_old2new.get)
        #     )
        #     sub_alleletable[2] = (
        #         sub_alleletable[2].map(vertices_idx.get).map(sub_old2new.get)
        #     )
        #     sub_alleletable.dropna(axis=0, inplace=True)
        #     P_allelic_idx = [
        #         sub_alleletable[1].astype(int),
        #         sub_alleletable[2].astype(int),
        #     ]
        #     P_allelic_sim = sub_alleletable['similarity'].astype(float).values
        #     P_weak_idx = None


        #     sub_involved_vertices = np.array([idx_to_vertices[sub_new2old[i]] for i in range(len(involved_idx))])

        #     df = HyperGraph.to_contacts(sub_raw_A, sub_involved_vertices, min_weight=self.min_weight, output=None)

        
        #     df = df.to_pandas()

        #     df.to_csv("collapsed.hypergraph.expansion.contacts", sep='\t', header=None, index=None)
        #     tmp_collapsed_contigs_list = "tmp.collapsed.contigs.list"
        #     with open(tmp_collapsed_contigs_list, 'w') as out:
        #         for contig in all_collapsed:
        #             print(contig, file=out)

        #     cmd = ["cphasing-rs", "kprune", self.alleletable.filename, "collapsed.hypergraph.expansion.contacts", 
        #            "collapsed.prune.table", "-t", str(self.threads), "--whitelist", tmp_collapsed_contigs_list, "-p"]
        #     flag = run_cmd(cmd, log=f"{self.log_dir}/collapsed.hyperpartition_kprune.log")
        #     assert flag == 0, "Failed to execute command, please check log."

        #     if Path(tmp_collapsed_contigs_list).exists():
        #         Path(tmp_collapsed_contigs_list).unlink()

        #     if not Path("collapsed.prune.table").exists() or Path("collapsed.prune.table").stat().st_size == 0:
        #         logger.warning("No pruned results from kprune, skipping rescue.")
        #         P_allelic_idx = None
        #         P_weak_idx = None
        #         P_allelic_sim = None
        #         P_mz_shared = None
        #     else:
        #         prune_table = pd.read_csv("collapsed.prune.table", sep='\t', header=None)
        #         prune_table[0] = prune_table[0].map(vertices_idx.get).map(sub_old2new.get)
        #         prune_table[1] = prune_table[1].map(vertices_idx.get).map(sub_old2new.get)

        #         allelic_table = prune_table[prune_table[6] == 0]
        #         cross_allelic_table = prune_table[prune_table[6] == 1]
        #         P_allelic_idx = [allelic_table[0].astype(int), allelic_table[1].astype(int)]
        #         P_weak_idx = [cross_allelic_table[0].astype(int), cross_allelic_table[1].astype(int)]
        #         P_allelic_sim = allelic_table[5].astype(float).values
        #         P_mz_shared = allelic_table[4].astype(float).values
        
        # sub_A = HyperGraph.reweight_adjacency_matrix(sub_raw_A, P_allelic_idx=P_allelic_idx, P_weak_idx=P_weak_idx, allelic_factor=0)
        if P_allelic_idx is not None:
            keep_mask = np.isin(P_allelic_idx[0], involved_idx) & np.isin(P_allelic_idx[1], involved_idx)
            rows = np.vectorize(sub_old2new.get)(P_allelic_idx[0][keep_mask]).astype(int)
            cols = np.vectorize(sub_old2new.get)(P_allelic_idx[1][keep_mask]).astype(int)
            
            unique_edges = set()
            for r, c in zip(rows, cols):
                if r != c:
                    unique_edges.add((min(r, c), max(r, c)))
            
            if unique_edges:
                sub_P_allelic_idx = [
                    np.array([e[0] for e in unique_edges], dtype=int),
                    np.array([e[1] for e in unique_edges], dtype=int)
                ]
                sub_P_allelic_sim = np.ones(len(unique_edges), dtype=float)
            else:
                sub_P_allelic_idx = None
                sub_P_allelic_sim = None
            
            if P_weak_idx is not None:
                keep_mask_weak = np.isin(P_weak_idx[0], involved_idx) & np.isin(P_weak_idx[1], involved_idx)
                rows_weak = np.vectorize(sub_old2new.get)(P_weak_idx[0][keep_mask_weak]).astype(int)
                cols_weak = np.vectorize(sub_old2new.get)(P_weak_idx[1][keep_mask_weak]).astype(int)
                
                unique_edges_weak = set()
                for r, c in zip(rows_weak, cols_weak):
                    if r != c:
                        unique_edges_weak.add((min(r, c), max(r, c)))
                
                if unique_edges_weak:
                    sub_P_weak_idx = [
                        np.array([e[0] for e in unique_edges_weak], dtype=int),
                        np.array([e[1] for e in unique_edges_weak], dtype=int)
                    ]
                else:
                    sub_P_weak_idx = None
            else:
                sub_P_weak_idx = None
        else:
            sub_P_allelic_idx = None
            sub_P_weak_idx = None
            sub_P_allelic_sim = None
            
        sub_A = HyperGraph.reweight_adjacency_matrix(
            sub_raw_A, 
            P_allelic_idx=sub_P_allelic_idx, 
            P_weak_idx=sub_P_weak_idx, 
            allelic_factor=0
        )
        group_to_idx = {}
        for hap_group, groups in self.hap_groups.items():
            for i, g_contigs in enumerate(groups):
                g_name = f"{hap_group}g{i+1}"
                group_to_idx[g_name] = [sub_old2new[vertices_idx[c]] for c in g_contigs if c in vertices_idx]

        allele_db = defaultdict(dict)
        if P_allelic_idx is not None:
            for u, v, sim in zip(P_allelic_idx[0], P_allelic_idx[1], P_allelic_sim):
                allele_db[u][v] = sim
                allele_db[v][u] = sim

        if unplaced_collapsed:
            u_lengths = self.contigsizes.reindex(unplaced_collapsed).values.flatten()
            u_lengths[np.isnan(u_lengths)] = 0
            unplaced_collapsed = [
                c for _, c in sorted(zip(u_lengths, unplaced_collapsed), reverse=True)
            ]
            
            logger.info(f"Attempting to pre-rescue {len(unplaced_collapsed)} unplaced collapsed contigs...")
          
            for u_contig in unplaced_collapsed:
                u_idx_new = sub_old2new[vertices_idx[u_contig]]
                
                hap_group_scores = {}
                for hap_group, groups in self.hap_groups.items():
                    hg_idx_new = []
                    for i in range(len(groups)):
                        g_name = f"{hap_group}g{i+1}"
                        if g_name in group_to_idx:
                            hg_idx_new.extend(group_to_idx[g_name])
                            
                    if not hg_idx_new:
                        continue
                        
                    interactions = sub_raw_A[u_idx_new, hg_idx_new]
                    if hasattr(interactions, "toarray"):
                        interactions = interactions.toarray()
                    elif hasattr(interactions, "A"):
                        interactions = interactions.A
                    interactions = np.array(interactions, dtype=float).flatten()
                    
                    target_contigs = [idx_to_vertices[sub_new2old[idx]] for idx in hg_idx_new]
                    lengths = self.contigsizes.reindex(target_contigs).values.flatten()
                    lengths[np.isnan(lengths) | (lengths == 0)] = median_len
                    rel_lengths = lengths / median_len
                    
                    norm_interactions = interactions / rel_lengths
                    
                    score = norm_interactions.mean()
                    hap_group_scores[hap_group] = score
                    
                if not hap_group_scores:
                    continue
                    
                best_hap_group = max(hap_group_scores, key=hap_group_scores.get)
                best_hg_score = hap_group_scores[best_hap_group]
      
                if best_hg_score <= 0.01:
                    continue
        
                best_sub_score = 0
                best_sub_idx = None
                
                groups = self.hap_groups[best_hap_group]
                for i in range(len(groups)):
                    g_name = f"{best_hap_group}g{i+1}"
                    g_idx_new = group_to_idx.get(g_name)
                    if not g_idx_new:
                        continue

                    interactions = sub_A[u_idx_new, g_idx_new]
                    if hasattr(interactions, "toarray"):
                        interactions = interactions.toarray()
                    elif hasattr(interactions, "A"):
                        interactions = interactions.A
                    interactions = np.array(interactions, dtype=float).flatten()
                    
                    target_contigs = [idx_to_vertices[sub_new2old[idx]] for idx in g_idx_new]
                    lengths = self.contigsizes.reindex(target_contigs).values.flatten()
                    lengths[np.isnan(lengths) | (lengths == 0)] = median_len
                    rel_lengths = lengths / median_len
                    
                    norm_interactions = interactions / rel_lengths

                    score = norm_interactions.mean()

                    if score > best_sub_score:
                        best_sub_score = score
                        best_sub_idx = i
                
                if best_sub_idx is not None and best_sub_score > 0.01:
                    self.hap_groups[best_hap_group][best_sub_idx].append(u_contig)
                    best_g_name = f"{best_hap_group}g{best_sub_idx+1}"
                    group_to_idx[best_g_name].append(u_idx_new)
                    logger.debug(f"Pre-rescued unplaced {u_contig} to {best_g_name} (HG Score: {best_hg_score:.2f}, Sub Score: {best_sub_score:.2f})")
                # elif best_hg_score > 0.01:
                #     self.hap_groups[best_hap_group].append([u_contig])
                #     new_idx = len(self.hap_groups[best_hap_group])
                #     best_g_name = f"{best_hap_group}g{new_idx}"
                #     group_to_idx[best_g_name] = [u_idx_new]
                #     logger.debug(f"Pre-rescued unplaced {u_contig} to NEW group {best_g_name} (HG Score: {best_hg_score:.2f})")

            current_placed = set()
            for groups in self.hap_groups.values():
                current_placed.update(list_flatten(groups))
                
            pre_rescued_contigs = [c for c in unplaced_collapsed if c in current_placed]
            pre_rescued_count = len(pre_rescued_contigs)
            
            if pre_rescued_count > 0:
                pre_rescued_length = int(np.nansum(self.contigsizes.reindex(pre_rescued_contigs).values))
            else:
                pre_rescued_length = 0

            logger.info(
                f"Pre-rescued {pre_rescued_count} unplaced collapsed contigs, total length: {pre_rescued_length:,} bp"
            )

        
        new_cluster_data = {}
        is_unphased = len(self.hap_groups) > 0 and all(len(g) == 1 for g in self.hap_groups.values())
        if is_unphased:
            groups_list = []
            original_names = []
            for hg, gs in self.hap_groups.items():
                groups_list.append(gs[0])
                original_names.append(hg)
            working_hap_groups = {"_GLOBAL_": groups_list}
            global_name_map = {"_GLOBAL_": original_names}
        else:
            working_hap_groups = self.hap_groups
            global_name_map = None

        for k, hap_group in enumerate(working_hap_groups):
            groups = working_hap_groups[hap_group]
            if len(groups) < 2:
                for i, g_contigs in enumerate(groups):
                    if global_name_map:
                        g_name = global_name_map[hap_group][i]
                    else:
                        g_name = f"{hap_group}g{i+1}"
                    new_cluster_data[g_name] = g_contigs
                continue

            contigs = list_flatten(groups)

            contigs_idx = list(
                filter(lambda x: x is not None, map(vertices_idx.get, contigs))
            )
            sub_old2new_idx = dict(zip(contigs_idx, range(len(contigs_idx))))

            if self.alleletable is None:
                sub_alleletable = None
                P_allelic_idx = None
                P_allelic_sim = None
            else:
                sub_alleletable = self.alleletable.data[
                    self.alleletable.data[1].isin(contigs)
                    & self.alleletable.data[2].isin(contigs)
                ].copy()
                sub_alleletable[1] = (
                    sub_alleletable[1].map(vertices_idx.get).map(sub_old2new_idx.get)
                )
                sub_alleletable[2] = (
                    sub_alleletable[2].map(vertices_idx.get).map(sub_old2new_idx.get)
                )
                sub_alleletable.dropna(axis=0, inplace=True)
                P_allelic_idx = [
                    sub_alleletable[1].astype(int),
                    sub_alleletable[2].astype(int),
                ]
                P_allelic_sim = sub_alleletable['similarity'].astype(float).values
                sub_alleletable.set_index([1], inplace=True)

            sub_H, _, _ = extract_incidence_matrix2(self.H, contigs_idx)
            sub_A = self.HG.clique_expansion_init(
                sub_H, P_allelic_idx=P_allelic_idx, allelic_factor=0
            )

            current_groups_idx = list(
                map(
                    lambda x: list(
                        filter(lambda x: x is not None, map(vertices_idx.get, x))
                    ),
                    groups,
                )
            )
            groups_new_idx = list(
                map(lambda x: list(map(sub_old2new_idx.get, x)), current_groups_idx)
            )

            groups_new_idx_db = {}
            for i in range(len(groups_new_idx)):
                for idx in groups_new_idx[i]:
                    groups_new_idx_db[idx] = i

            sub_collapsed_contigs = self.collapsed_contigs.reindex(contigs).dropna(
                axis=0
            )
            sub_collapsed_contigs_idx_new = [
                sub_old2new_idx.get(vertices_idx.get(idx))
                for idx in sub_collapsed_contigs.index
            ]
            sub_collapsed_contigs_idx_new = [
                int(idx) for idx in sub_collapsed_contigs_idx_new if idx is not None
            ]
            sub_collapsed_contigs_cn_dict = dict(
                zip(sub_collapsed_contigs_idx_new, sub_collapsed_contigs["CN"])
            )
            allele_db = defaultdict(dict)
            if P_allelic_idx is not None and P_allelic_sim is not None:
                for u, v, sim in zip(P_allelic_idx[0], P_allelic_idx[1], P_allelic_sim):
                    allele_db[u][v] = sim
                    allele_db[v][u] = sim
            original_groups_new_idx = [list(g) for g in groups_new_idx]
            if isinstance(self.contigsizes, pd.Series):
                contig_len_db = self.contigsizes.to_dict()
            elif isinstance(self.contigsizes, pd.DataFrame):
                contig_len_db = self.contigsizes.iloc[:, 0].to_dict()
            else:
                contig_len_db = {}
            if not self.only_unplaced_rescue:

                for j in sub_collapsed_contigs_idx_new:
                    group_scores = []
                    for i in range(len(groups_new_idx)):
                        if len(groups_new_idx[i]) > 0:
                            interactions = sub_A[j, groups_new_idx[i]]
                            if hasattr(interactions, "toarray"):
                                interactions = interactions.toarray()
                            elif hasattr(interactions, "A"): 
                                interactions = interactions.A
                            interactions = np.asarray(interactions).flatten()
            
                            target_contigs = [idx_to_vertices[contigs_idx[idx]] for idx in groups_new_idx[i]]
                            lengths = self.contigsizes.reindex(target_contigs).values.flatten()
                            lengths[np.isnan(lengths) | (lengths == 0)] = median_len
                            rel_lengths = lengths / median_len
                            
                            norm_interactions = interactions / rel_lengths
                            
                            score = norm_interactions.mean()
                        else:
                            score = 0
                     
                        group_scores.append(score)

                    group_scores = np.array(group_scores)

                    best_group_idx = np.argmax(group_scores)
                    primary_score = group_scores[best_group_idx]

                
                    original_g_idx = groups_new_idx_db.get(j)

                    if primary_score <= 0.01:
                        continue

                    if original_g_idx is not None and best_group_idx != original_g_idx:
                        logger.debug(
                            f"Contig {idx_to_vertices[contigs_idx[j]]} re-assigned: Group {original_g_idx} -> {best_group_idx}"
                        )

                    cn = int(round(sub_collapsed_contigs_cn_dict[j]))
                    if cn < 2:
                        continue

                    candidates = []
                    for idx, s in enumerate(group_scores):
                        
                        if idx == best_group_idx:
                            continue
                        if s > 0.01:
                            candidates.append((idx, s))

                    candidates.sort(key=lambda x: x[1], reverse=True)

                    rescued_count = 0
                    j_name = idx_to_vertices[contigs_idx[j]]
                    j_len = contig_len_db.get(j_name, median_len)
                    if np.isnan(j_len) or j_len == 0:
                        j_len = median_len
                    j_alleles = allele_db.get(j, {})
                    conflicted_candidates = []

                    for candidate_idx, score in candidates:
                        logger.debug(
                            f"Evaluating candidate group {candidate_idx} for contig {idx_to_vertices[
                                contigs_idx[j]
                            ]} (Score: {score:.2f}, Primary: {primary_score:.2f})"
                        )
                        if rescued_count >= (cn - 1):
                            break

                        if score < (primary_score * 0.025):
                            continue
                        
                        has_conflict = False
                        if not self.disable_conflict_check:
                            candidate_members = groups_new_idx[candidate_idx]
                            conflict_members = set(candidate_members).intersection(j_alleles.keys())
                            
                            if conflict_members:
                                total_overlap_for_j = 0.0
                                for member in conflict_members:
                                    sim = float(j_alleles[member])
                                    m_name = idx_to_vertices[contigs_idx[member]]
                                    m_len = contig_len_db.get(m_name, median_len)
                                    if np.isnan(m_len) or m_len == 0:
                                        m_len = median_len

                                    overlap_bp = min(j_len, m_len) * sim
                                    # conflict_m = m_len > 0 and (overlap_bp / m_len) > 0.50
                                    # conflict_j = j_len > 0 and (overlap_bp / j_len) > 0.50
                                    
                                    # if conflict_m or conflict_j:
                                    #     has_conflict = True
                                    #     break
                                    if m_len > 0 and (overlap_bp / m_len) > 0.50:
                                        has_conflict = True
                                        break
                            
                                if not has_conflict and j_len > 0 and (total_overlap_for_j / j_len) > 0.50:
                                    has_conflict = True
                                
                        if has_conflict:
                            conflicted_candidates.append((candidate_idx, score))
                            continue 

                        if j not in groups_new_idx[candidate_idx]:
                            groups_new_idx[candidate_idx].append(j)
                            logger.debug(
                                f"Rescued {idx_to_vertices[contigs_idx[j]]} to group {candidate_idx} (Score: {score:.2f})"
                            )
                            rescued_count += 1
                    
                    if rescued_count < (cn - 1) and conflicted_candidates:
                        def get_overlap(c_idx):
                            overlap = 0.0
                            conflict_members = set(groups_new_idx[c_idx]).intersection(j_alleles.keys())
                            for m in conflict_members:
                                m_name = idx_to_vertices[contigs_idx[m]]
                                m_len = contig_len_db.get(m_name, median_len)
                                if np.isnan(m_len) or m_len == 0:
                                    m_len = median_len
                                overlap += min(j_len, m_len) * float(j_alleles[m])
                            return overlap
                        conflicted_candidates.sort(key=lambda x: (get_overlap(x[0]), -x[1]))

                        for candidate_idx, score in conflicted_candidates:
                            if rescued_count >= (cn - 1):
                                break
                            
                            if j not in groups_new_idx[candidate_idx]:
                                groups_new_idx[candidate_idx].append(j)
                                logger.debug(
                                    f"[Fallback] Rescued {j_name} to group {candidate_idx} bypassing conflicts (Score: {score:.2f})"
                                )
                                rescued_count += 1

                    # while rescued_count < (cn - 1):
                    #         groups_new_idx.append([j])
                    #         new_idx = len(groups_new_idx) - 1
                    #         logger.debug(
                    #             f"[New Group] Rescued {j_name} by creating a new group {new_idx} due to lack of candidates"
                    #         )
                    #         rescued_count += 1
            for k, idx_list in enumerate(groups_new_idx):
                converted_names = [idx_to_vertices[contigs_idx[i]] for i in idx_list]
                
              
                if global_name_map and k < len(global_name_map[hap_group]):
                    g_name = global_name_map[hap_group][k]
                else:
                    g_name = f"{hap_group}g{k+1}"
                jrrrt6iol=6
                new_cluster_data[g_name] = list(
                    OrderedDict.fromkeys(converted_names)
                ) 


        self.clustertable.data = new_cluster_data

        new_contigs_all = self.clustertable.contigs
        count_stats = Counter(new_contigs_all)
        rescued_contigs = [c for c, count in count_stats.items() if count > 1]

        with open("collapsed.rescue.contigs.list", "w") as out:
            for contig in rescued_contigs:
                for i in range(2, count_stats[contig] + 1):
                    print(f"{contig}\t{contig}_d{i}", file=out)

        if rescued_contigs:
            lengths = self.contigsizes.reindex(rescued_contigs).values.flatten()
            lengths[np.isnan(lengths)] = 0
            extra_copies = np.array([count_stats[c] - 1 for c in rescued_contigs])
            
            rescued_length = int(np.sum(lengths * extra_copies))
            rescued_copies_count = int(np.sum(extra_copies))
        else:
            rescued_length = 0
            rescued_copies_count = 0

        logger.info(
            f"Rescued {rescued_copies_count} copies from {len(rescued_contigs)} unique contigs, total rescued length: {rescued_length:,} bp"
        )
        stat_table = self.clustertable.get_stat_table(self.contigsizes)
        _console = Console(file=io.StringIO(), force_terminal=False)
        _console.print(stat_table)
        
        logger.info(
            f"Cluster stats after rescue:\n{_console.file.getvalue()}"
        )
        self.clustertable.save("collapsed.rescue.clusters.txt")


class CollapsedRescue2:
    """
    Rescue the collapsed contigs into a ordered and oriented group
    """

    def __init__(
        self,
        HG,
        agp,
        fasta,
        alleletable,
        split_contacts,
        collapsed_contigs,
        allelic_similarity: float = 0.85,
    ):
        self.HG = HG
        self.agp_file = agp

        self.cluster_df = agp2cluster(self.agp_file, store=False)
        self.cluster_df = self.cluster_df.reset_index()
        self.cluster_df = self.cluster_df[
            self.cluster_df["id"].apply(lambda x: len(x) > 0)
        ]
        self.cluster_df.set_index("chrom", inplace=True)
        self.cluster_data = self.cluster_df.to_dict()["id"]

        self.fasta = fasta
        self.fasta_path = Path(self.fasta).absolute()

        self.contigsizes = read_chrom_sizes(
            str(get_contig_size_from_fasta(self.fasta_path))
        )
        self.alleletable = alleletable
        self.alleletable.data = self.alleletable.data[
            self.alleletable.data["similarity"] >= allelic_similarity
        ]
        # self.alleletable.data = self.alleletable.data[
        #     self.alleletable.data["mzShared"] > 200
        # ]

        self.split_contacts = split_contacts

        if self.split_contacts:
            self.split_link_df = pl.read_csv(
                self.split_contacts,
                separator="\t",
                has_header=False,
                dtypes={
                    "column_1": pl.Categorical,
                    "column_2": pl.Categorical,
                    "column_3": pl.UInt32,
                },
            ).to_pandas()
            self.split_link_df.columns = [0, 1, 2]

        else:
            self.split_link_df = None

        self.collapsed_contigs = collapsed_contigs
        self.allelic_similarity = allelic_similarity

        self.H = HG.incidence_matrix()
        self.vertices = self.HG.nodes

    @property
    def hap_groups(self):
        db = OrderedDict()
        for group in self.cluster_data:
            try:
                hap, idx = group.rsplit("g", 1)
            except:
                if group not in db:
                    db[group] = []
                db[group].append(self.cluster_data[group])
                continue
            finally:
                try:
                    if hap not in db:
                        db[hap] = []

                    db[hap].append(self.cluster_data[group])
                except UnboundLocalError:
                    logger.warning(
                        "Unexpect group name `{group}`, must be Chr[\d+]g[\d+]"
                    )

        return db

    @property
    def vertices_idx(self):
        return dict(zip(self.vertices, range(len(self.vertices))))

    @property
    def idx_to_vertices(self):
        idx_to_vertices = dict(zip(range(len(self.vertices)), self.vertices))

        return idx_to_vertices

    def get_group_with_orientation(self):
        db = agp2tour(self.agp_file, store=False)
        for group in db:
            db[group] = OrderedDict(db[group])

        return db

    @staticmethod
    def _rescue(
        A,
    ):
        pass

    def rescue(self):
        vertices_idx = self.vertices_idx
        idx_to_vertices = self.idx_to_vertices
        new_cluster_data = {}
        rescued_data = OrderedDict()
        for k, hap_group in enumerate(self.hap_groups):
            groups = self.hap_groups[hap_group]
            if len(groups) < 2:
                continue
            contigs = list_flatten(groups)
            contigs_idx = list(
                filter(lambda x: x is not None, map(vertices_idx.get, contigs))
            )

            sub_old2new_idx = dict(zip(contigs_idx, range(len(contigs_idx))))
            sub_alleletable = self.alleletable.data[
                self.alleletable.data[1].isin(contigs)
                & self.alleletable.data[2].isin(contigs)
            ]
            sub_alleletable[1] = (
                sub_alleletable[1].map(vertices_idx.get).map(sub_old2new_idx.get)
            )
            sub_alleletable[2] = (
                sub_alleletable[2].map(vertices_idx.get).map(sub_old2new_idx.get)
            )
            P_allelic_idx = [sub_alleletable[1], sub_alleletable[2]]
            sub_alleletable.set_index([1], inplace=True)

            sub_H, _ = extract_incidence_matrix2(self.H, contigs_idx)
            sub_A = self.HG.clique_expansion_init(
                sub_H, P_allelic_idx=P_allelic_idx, allelic_factor=0
            )

            # sub_allelic = set(map(tuple, sub_alleletable[[1, 2]].values.tolist()))

            groups_idx = list(
                map(
                    lambda x: list(
                        filter(lambda x: x is not None, map(vertices_idx.get, x))
                    ),
                    groups,
                )
            )

            groups_new_idx = list(
                map(lambda x: list(map(sub_old2new_idx.get, x)), groups_idx)
            )
            groups_new_idx_db = {}
            rescued_new_idx_db = OrderedDict()
            for i in range(len(groups_new_idx)):
                groups_new_idx_db.update(
                    dict(zip(groups_new_idx[i], [i] * len(groups_new_idx[i])))
                )

            sub_collapsed_contigs = self.collapsed_contigs.reindex(contigs).dropna(
                axis=0
            )
            sub_collapsed_contigs_idx = list(
                map(vertices_idx.get, sub_collapsed_contigs.index.tolist())
            )
            sub_collapsed_contigs_idx_new = list(
                map(sub_old2new_idx.get, sub_collapsed_contigs_idx)
            )

            res = []

            for j in sub_collapsed_contigs_idx_new:
                try:
                    tmp_allelic_table = sub_alleletable.loc[j]
                    try:
                        tmp_allelic_table.set_index([2], inplace=True)
                    except:
                        tmp_allelic_table.to_frame().set_index([2], inplace=True)

                except KeyError:
                    tmp_allelic_table = []
                shared_similarity = np.zeros(len(groups))
                tmp_res = []
                for i in range(len(groups)):
                    if i == groups_new_idx_db[j]:
                        tmp_res.append(0)
                    else:
                        if len(tmp_allelic_table) > 1:
                            tmp = tmp_allelic_table.reindex(groups_new_idx[i]).dropna()
                            if len(tmp) > 1:
                                shared_similarity[i] = tmp["mzShared"].sum()

                        tmp_res.append(sub_A[j, groups_new_idx[i]].mean())

                max_idx = np.argmax(tmp_res)
                groups_new_idx[max_idx].append(j)
                if max_idx not in rescued_new_idx_db:
                    rescued_new_idx_db[max_idx] = []
                rescued_new_idx_db[max_idx].append(j)
                res.append(tmp_res)

            new_rescued_groups = []
            new_groups = []
            for n, group_idx in enumerate(groups_new_idx):
                tmp = list(map(lambda x: contigs_idx[x], group_idx))
                tmp = list(map(idx_to_vertices.get, tmp))
                new_groups.append(tmp)

                tmp = list(map(lambda x: contigs_idx[x], rescued_new_idx_db[n]))
                tmp = list(map(idx_to_vertices.get, tmp))
                new_rescued_groups.append(tmp)

            for k, group in enumerate(new_rescued_groups):
                rescued_data[f"{hap_group}g{k+1}"] = group

            for k, group in enumerate(new_groups):
                new_cluster_data[f"{hap_group}g{k+1}"] = group

        self.rescued_data = rescued_data

    def assign(self):
        """
        assign rescue collapsed contigs into correct position and orientation
        """
        split_contig_size_df = self.contigsizes // 2
        split_contig_size_db = split_contig_size_df.to_dict()['length']

        
        self.split_link_df['contig1'] = self.split_link_df[0].str.rsplit('_', n=1).str[0]
        self.split_link_df['contig2'] = self.split_link_df[1].str.rsplit('_', n=1).str[0]

        self.split_link_df['length1'] = self.split_link_df['contig1'].map(split_contig_size_db.get)
        self.split_link_df['length2'] = self.split_link_df['contig2'].map(split_contig_size_db.get)

        self.split_link_df[2] = self.split_link_df[2] / (self.split_link_df['length1'] * self.split_link_df['length2'])
        
        self.split_link_df.drop(['contig1', 'contig2', 'length1', 'length2'], inplace=True, axis=1)
        split_link_df2 = self.split_link_df.copy()
        split_link_df2.columns = [1, 0, 2]
        self.split_link_df = pd.concat([self.split_link_df, split_link_df2], axis=0)
        self.split_link_df = self.split_link_df.drop_duplicates(subset=[0, 1], keep='first')
        self.split_link_df.set_index([0, 1], inplace=True)


        assign_result = OrderedDict()
        for group in self.rescued_data:
            contigs = self.cluster_data[group]
            split_contigs = list_flatten([list(map(lambda x: f"{x}_0", contigs)), 
                                          list(map(lambda x: f"{x}_1", contigs))])
            res = []
            for contig in self.rescued_data[group]:
                contig1 = f"{contig}_0"
                contig2 = f"{contig}_1"
                tmp_data = self.split_link_df.reindex(list(set(list(product([contig1, contig2], split_contigs))))).dropna()
                tmp_data = tmp_data.reset_index()
                tmp_data = tmp_data.loc[tmp_data[0].str.rsplit("_", n=1).str[0] != tmp_data[1].str.rsplit("_", n=1).str[0]]
                tmp_data = tmp_data.sort_values(2, ascending=False)
                
                if tmp_data.empty:
                    continue
                split_contigs.append(contig1)
                split_contigs.append(contig2)

               
                res.append(tmp_data.loc[tmp_data[2].idxmax()][[0, 1]].values.tolist())
                
            assign_result[group] = res
        
        
        db = self.get_group_with_orientation()
        for group in assign_result:
            tmp_orientations = db[group]
    
            for pair in assign_result[group]:
                contig1, contig2 = pair
                contig1, suffix1 = contig1.rsplit("_", 1)
                contig2, suffix2 = contig2.rsplit("_", 1)

                
                strand2 = tmp_orientations[contig2]
                tmp_orientations = list(tmp_orientations.items())
                idx = tmp_orientations.index((contig2, strand2))
                if strand2 == "+":
                    if suffix2 == "1":
                        idx = idx + 1
                        if suffix1 == "0":
                            strand1 = "+"
                        else:
                            strand1 = "-"
                    else:
                        idx = idx
                        if suffix1 == "0":
                            strand1 = "-"
                        else:
                            strand1 = "+"
                else:
                    if suffix2 == "1":
                        idx = idx
                        if suffix1 == "0":
                            strand1 = "-"
                        else:
                            strand1 = "+"
                    else:
                        idx = idx + 1 
                        if suffix1 == "0":
                            strand1 = "+"
                        else:
                            strand1 = "-"

                
                
                tmp_orientations.insert(idx, (contig1, strand1))
                tmp_orientations = OrderedDict(tmp_orientations)
            db[group] = tmp_orientations

        outdir = "tour"
        Path(outdir).mkdir(exist_ok=True)
        tour = Tour
        for group in db:
            db[group] = list(db[group].items())

            tour.from_tuples(db[group])
            with open(f'{str(outdir)}/{group}.tour', 'w') as out:
                print(' '.join(map(str, tour.data)), file=out)
                logger.debug(f'Output tour: `{out.name}`')
        

        os.chdir(outdir)
        try:
            build.main(args=[str(self.fasta_path),
                            "-oa", "groups.dup.raw.agp",
                            "--only-agp",], 
                       prog_name="build")

        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
            
        
        try:
            agp2fasta.main(args=["groups.dup.raw.agp", 
                                 str(self.fasta_path), 
                                 "--contigs", "-o", 
                                 "contigs.fasta"],
                            prog_name="agp2fasta")
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e

        try:
            agp_dup.main(args=["groups.dup.raw.agp", 
                                "-o",
                                "groups.dup.renamed.agp"
                                 ],
                            prog_name="agp2fasta")
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
            
        shutil.copy(f"groups.dup.raw.agp", "../")
        shutil.copy(f"groups.dup.renamed.agp", "../")
        shutil.copy(f"contigs.collapsed.contig.list", "../")
        shutil.copy(f"contigs.fasta", "../contigs.dup.fasta")

        os.chdir("..")
        
        shutil.rmtree(outdir)

    def run(self):
        self.rescue()
        if self.split_link_df is not None:
            self.assign()


def convert_matrix_with_dup_contigs(cool, dup_contig_path, output, threads=4):
    """
    Params:
    --------
    cool: str
        Path to cool file
    dup_contigs: str
        Path of two columns of `raw_contig dup_contig`
    output: str
        Path of output cool
    """
    pandarallel.initialize(nb_workers=threads, verbose=0)

    cool = cooler.Cooler(cool)
    bins = cool.bins()[:]

    bins_chrom_idx_db = bins.reset_index().groupby('chrom')['index'].apply(list).to_dict()
    bins_chrom = bins.reset_index().set_index('chrom')
    max_bin_id = len(bins)

    matrix = cool.matrix(balance=False, sparse=True)[:].tocsr().astype('float64')
 
    dup_contigs_db = defaultdict(list)
    with open(dup_contig_path) as fp:
        for line in fp:
            line_list = line.strip().split()[:2]
            dup_contigs_db[line_list[0]].append(line_list[1])

   
    bin_res = [bins.reset_index()]
    res = []

    for contig in dup_contigs_db:
        dup_contigs = dup_contigs_db[contig]
        cn = len(dup_contigs) + 1
        
        idxes = bins_chrom_idx_db[contig]
        tmp_bins = bins_chrom.loc[contig]
        tmp_raw_matrix = matrix[idxes]
        
       

        for i, dup_contig in enumerate(dup_contigs):
            
            length_of_tmp_bins = len(tmp_bins)
            tmp_bins_2 = tmp_bins.copy()
            tmp_bins_2['index'] = range(max_bin_id, max_bin_id + length_of_tmp_bins)
            
            tmp_bins_2['chrom'] = dup_contig 
            row_bin_idx_db = dict(zip(range(0, length_of_tmp_bins), 
                                      range(max_bin_id, max_bin_id + length_of_tmp_bins)))
            col_bin_idx_db = dict(zip(tmp_bins['index'].values.tolist(), 
                                    tmp_bins_2['index'].values.tolist()))
            bin_res.append(tmp_bins_2)
            max_bin_id += length_of_tmp_bins
            

            tmp_matrix = (tmp_raw_matrix / cn ).tocoo()
            
            tmp_pixels = pd.Series.sparse.from_coo(tmp_matrix)
            tmp_pixels = (
                            tmp_pixels.reset_index()
                                    .rename(columns={'level_0': 'bin1_id', 
                                                    'level_1': 'bin2_id',
                                                    0: 'count'})
                        )

            def get_col_idx(x):
                try:
                    return col_bin_idx_db[x]
                except KeyError:
                    return x
            
            tmp_pixels['bin1_id'] = tmp_pixels['bin1_id'].parallel_apply(lambda x: row_bin_idx_db[x])
            tmp_pixels['bin2_id'] = tmp_pixels['bin2_id'].parallel_apply(get_col_idx)
            res.append(tmp_pixels) 
        
        rows, cols = triu(matrix).tocsr()[idxes].nonzero()
        matrix[rows, cols] /= float(cn)
        rows, cols = triu(matrix).tocsr()[:, idxes].nonzero()
        matrix[rows, cols] /= float(cn)

    pixels = pd.Series.sparse.from_coo(triu(matrix).tocoo())
    
    pixels = pixels.reset_index()
    pixels.columns = ['bin1_id', 'bin2_id', 'count']
    
    del matrix
    gc.collect()
    new_pixels = pd.concat(res, axis=0)
    
    new_pixels_up = new_pixels[new_pixels['bin1_id'] <= new_pixels['bin2_id']]
    new_pixels_lower = new_pixels[new_pixels['bin1_id'] > new_pixels['bin2_id']]
    new_pixels_lower.columns = ['bin2_id', 'bin1_id', 'count']
    
    new_pixels = pd.concat([pixels, new_pixels_up, new_pixels_lower], axis=0)

    new_pixels = new_pixels.drop_duplicates(subset=['bin1_id', 'bin2_id'])

    bins = pd.concat(bin_res, axis=0)
    bins = bins.drop(['index'], axis=1)

    cooler.create_cooler(output, bins, new_pixels, dtypes={'count': 'float64'})

    
def get_dup_prune_pairs():
    pass



def fasta_dup(fasta, collapsed_list, output):
    fasta = read_fasta(fasta)
    collapsed_db = {}
    with open(collapsed_list) as fp:
        for line in fp:
            line = line.strip().split()
            if line[0] not in collapsed_db:
                collapsed_db[line[0]] = []
            collapsed_db[line[0]].append(line[1])
    
    
    for contig in fasta:
        if contig in collapsed_db:
            for dup_contig in collapsed_db[contig]:
                print(f">{dup_contig}\n{fasta[contig]}", file=output)            
        print(f">{contig}\n{fasta[contig]}", file=output)

    
if __name__ == "__main__":  
    cool = sys.argv[1]
    dup_contigs_path = sys.argv[2]
    output = sys.argv[3]
    convert_matrix_with_dup_contigs(cool, dup_contigs_path, output)
        
        
        