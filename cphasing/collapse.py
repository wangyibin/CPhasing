#!/usr/bin/env python
# -*- coding:utf-8 -*-


import logging
import gc
import os
import os.path as op
import shutil
import sys
import warnings

import cooler 
import numpy as np
import pandas as pd
import polars as pl

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MaxNLocator
from scipy.signal import find_peaks

from collections import Counter, defaultdict, OrderedDict
from itertools import product
from scipy.sparse import triu
from pandarallel import pandarallel 
from pathlib import Path

try:
    from .algorithms.hypergraph import extract_incidence_matrix2, HyperGraph
    from .core import Tour
    from .utilities import (run_cmd, list_flatten, 
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

class CollapseFromDepth:
    """
    get collapsed contigs from read depth
    """
    def __init__(self, depth):
        self.depth = depth
        self.depth_df = self.read_depth()

        self.contigsizes = self.depth_df.groupby('chrom').agg({'end': 'max'})
        self.contigsizes.columns = ['length']
        self.contigsizes_db = self.contigsizes.to_dict()['length']

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

        peak_ind = find_peaks(y, distance=10)[0]
        median_value = np.quantile(data, .3)
        peak_ind = list(filter(lambda j: x[j] > median_value, peak_ind))
        if len(peak_ind) == 0:
            max_idx = np.argsort(x)[len(x)//2]
        else:
            max_idx = peak_ind[np.argmax(y[peak_ind])]

        # ax.fill_between((x[max_idx] * lower_value, x[max_idx] * upper_value), 
        #             0, ax.get_ylim()[1], alpha=0.5 , color='#bcbcbc')
       
        ax.axvline(x[max_idx], linestyle='--', color='#b02418')
        ax.text(int(x[max_idx]), ax.get_ylim()[1] * 0.98, str(int(x[max_idx])), fontsize=14, color='#cb6e7f')    
        ax.axvline(x=x[max_idx] * 2, linewidth=1, color='#253761', linestyle='--')
        ax.text(x[max_idx] * 2, ax.get_ylim()[1] * 0.98, f"{x[max_idx] * 2:.1f}", fontsize=14, color='#253761')

        sns.despine()
        plt.legend([], frameon=False)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        plt.tick_params(which='both', width=1.5, length=5)

        plt.savefig(f'{output}.kde.plot.png', dpi=600, bbox_inches='tight')
        plt.savefig(f'{output}.kde.plot.pdf', dpi=600, bbox_inches='tight')

        return int(x[max_idx]), x[max_idx] * lower_value, x[max_idx] * upper_value
        
    def run(self):
        peak_depth, lower_depth, upper_depth = self.plot_distribution()

        rd_df = self.depth_df
        rd_df['CN'] = rd_df['count'] / peak_depth
        
        contig_depth = rd_df.groupby('chrom')['count'].mean().to_frame()
        contig_depth['CN'] = contig_depth['count'] / peak_depth

        collapsed_df = contig_depth[contig_depth['CN'] > 1.5]
        
        collapsed_df['length'] = collapsed_df.index.map(self.contigsizes_db.get)

        collapsed_df = collapsed_df.reset_index()
        logger.info(f"Identified {collapsed_df['length'].sum():,} bp contigs with read depth >= {upper_depth:.2f}")
        collapsed_df.drop(columns=["length"], inplace=True)
        collapsed_df.to_csv("contigs.collapsed.contig.list", sep='\t', header=False, index=False)

        low_df = contig_depth[contig_depth['CN'] < 0.1]
        low_df['length'] = low_df.index.map(self.contigsizes_db.get)
        low_df = low_df.reset_index()
        logger.info(f"Identified {low_df['length'].sum():,} bp contigs with read depth < {lower_depth:.2f}")
        low_df.drop(columns=["length"], inplace=True)
        low_df.to_csv("contigs.low.contig.list", sep='\t', header=False, index=False)
        

class CollapseFromGfa:
    """
    get collapsed contigs from gfa 
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

        peak_ind = find_peaks(y, distance=10)[0]
        median_value = np.quantile(data, .3)
        peak_ind = list(filter(lambda j: x[j] > median_value, peak_ind))
        if len(peak_ind) == 0:
            max_idx = np.argsort(x)[len(x)//2]
        else:
            max_idx = peak_ind[np.argmax(y[peak_ind])]

        ax.axvline(x[max_idx], linestyle='--', color='#b02418')
        ax.text(int(x[max_idx]), ax.get_ylim()[1] * 0.98, str(int(x[max_idx])), fontsize=14, color='#cb6e7f')    
        ax.axvline(x=x[max_idx] * 2, linewidth=1, color='#253761', linestyle='--')
        ax.text(x[max_idx] * 2, ax.get_ylim()[1] * 0.98, f"{x[max_idx] * 2:.1f}", fontsize=14, color='#253761')

        sns.despine()
        plt.legend([], frameon=False)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        plt.tick_params(which='both', width=1.5, length=5)

        plt.savefig(f'{output}.kde.plot.png', dpi=600, bbox_inches='tight')
        plt.savefig(f'{output}.kde.plot.pdf', dpi=600, bbox_inches='tight')

        return int(x[max_idx]), x[max_idx] * lower_value, x[max_idx] * upper_value
        

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
                    collapsed_contigs, allelic_similarity: float=.85):
        
        self.HG = HG 
        self.clustertable = clustertable
        self.alleletable = alleletable 
        self.alleletable.data = self.alleletable.data[
            self.alleletable.data['similarity'] >= allelic_similarity]

        self.contigsizes = contigsizes
        self.collapsed_contigs = collapsed_contigs
        self.allelic_similarity = allelic_similarity
        
        self.hap_groups = self.clustertable.hap_groups

        self.H = HG.incidence_matrix(min_contacts=5)
        
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
        
        A = HyperGraph.clique_expansion_init(self.H, P_allelic_idx=None, min_weight=2)
        dia = A.diagonal()
        raw_contig_counts = A.shape[0]
        retain_idx = np.where(dia > 20)[0]
        contig_counts = len(retain_idx)

        if len(retain_idx) < raw_contig_counts:
            A = A[retain_idx, :][:, retain_idx]
            self.H, _, _ = extract_incidence_matrix2(self.H, retain_idx)
            self.vertices = self.vertices[retain_idx]
            logger.info(f"Filtered {raw_contig_counts - contig_counts} contigs with less than 20 cis weight")

    def rescue(self):
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
            sub_collapsed_contigs_idx = list(map(vertices_idx.get, sub_collapsed_contigs.index.tolist()))
            sub_collapsed_contigs_idx_new = list(map(sub_old2new_idx.get, sub_collapsed_contigs_idx))

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
                                shared_similarity[i] = tmp['mzShared'].sum()

                        tmp_res.append(sub_A[j, groups_new_idx[i]].mean())
                     

                groups_new_idx[np.argmax(tmp_res)].append(j)
                logger.debug(f"Rescued contig {idx_to_vertices[j]} into group {np.argmax(tmp_res)}")
                res.append(tmp_res)
            
            new_groups = []
            for group_idx in groups_new_idx:
                tmp = list(map(lambda x: contigs_idx[x], group_idx))
                tmp = list(map(idx_to_vertices.get, tmp))
                new_groups.append(tmp)

            for k, group in enumerate(new_groups):
                new_cluster_data[f'{hap_group}g{k+1}'] = group 
        
    
        self.clustertable.data = new_cluster_data
        ## duplicated
       
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


class CollapsedRescue2:
    """
    Rescue the collapsed contigs into a ordered and oriented group
    """
    def __init__(self, HG, agp, fasta,
                 alleletable, split_contacts, 
                    collapsed_contigs, allelic_similarity: float=.85):
        
        self.HG = HG 
        self.agp_file = agp

        self.cluster_df = agp2cluster(self.agp_file, store=False)
        self.cluster_df = self.cluster_df.reset_index()
        self.cluster_df = self.cluster_df[self.cluster_df['id'].apply(lambda x: len(x) > 0)]
        self.cluster_df.set_index('chrom', inplace=True)
        self.cluster_data = self.cluster_df.to_dict()['id']

        self.fasta = fasta 
        self.fasta_path = Path(self.fasta).absolute()
        
        self.contigsizes = read_chrom_sizes(str(get_contig_size_from_fasta(self.fasta_path)))
        self.alleletable = alleletable 
        self.alleletable.data = self.alleletable.data[
            self.alleletable.data['similarity'] >= allelic_similarity]

        self.split_contacts = split_contacts

        if self.split_contacts:
            self.split_link_df  = pl.read_csv(self.split_contacts, separator='\t', has_header=False,
                            dtypes={"column_1": pl.Categorical, 
                                    "column_2": pl.Categorical,
                                    "column_3": pl.UInt32}).to_pandas()
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
                    logger.warning("Unexpect group name `{group}`, must be Chr[\d+]g[\d+]")
            
        return db

    @property
    def vertices_idx(self):
        return dict(zip(self.vertices, 
                        range(len(self.vertices))))
    
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
    def _rescue(A, ):
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
            contigs_idx = list(filter(lambda x: x is not None, map(vertices_idx.get, contigs)))
            
            sub_old2new_idx = dict(zip(contigs_idx, range(len(contigs_idx))))
            sub_alleletable = self.alleletable.data[
                self.alleletable.data[1].isin(contigs) & self.alleletable.data[2].isin(contigs)]
            sub_alleletable[1] = sub_alleletable[1].map(vertices_idx.get).map(sub_old2new_idx.get)
            sub_alleletable[2] = sub_alleletable[2].map(vertices_idx.get).map(sub_old2new_idx.get)
            P_allelic_idx = [sub_alleletable[1], sub_alleletable[2]]
            sub_alleletable.set_index([1], inplace=True)
            
            sub_H, _ = extract_incidence_matrix2(self.H, contigs_idx)
            sub_A = self.HG.clique_expansion_init(sub_H, P_allelic_idx=P_allelic_idx, allelic_factor=0)
           
            # sub_allelic = set(map(tuple, sub_alleletable[[1, 2]].values.tolist()))

            groups_idx = list(map(lambda x: list(filter(lambda x: x is not None, map(vertices_idx.get, x))), groups))

            groups_new_idx = list(map(lambda x: list(map(sub_old2new_idx.get, x)), groups_idx))
            groups_new_idx_db = {}
            rescued_new_idx_db = OrderedDict()
            for i in range(len(groups_new_idx)):
                groups_new_idx_db.update(dict(zip(groups_new_idx[i], [i] * len(groups_new_idx[i]))))
                
            sub_collapsed_contigs = self.collapsed_contigs.reindex(contigs).dropna(axis=0)
            sub_collapsed_contigs_idx = list(map(vertices_idx.get, sub_collapsed_contigs.index.tolist()))
            sub_collapsed_contigs_idx_new = list(map(sub_old2new_idx.get, sub_collapsed_contigs_idx))

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
                                shared_similarity[i] = tmp['mzShared'].sum()

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
                rescued_data[f'{hap_group}g{k+1}'] = group 


            for k, group in enumerate(new_groups):
                new_cluster_data[f'{hap_group}g{k+1}'] = group 
            
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
        
        ## replace 
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
    

    
if __name__ == "__main__":  
    cool = sys.argv[1]
    dup_contigs_path = sys.argv[2]
    output = sys.argv[3]
    convert_matrix_with_dup_contigs(cool, dup_contigs_path, output)
        
        
        