#!/usr/bin/env python
# -*- coding:utf-8 -*-


import logging
import os
import os.path as op
import sys
import shutil

import dask.dataframe as dd
import glob
import numpy as np
import tempfile
import pandas as pd

from collections import OrderedDict 
from itertools import permutations
from joblib import Parallel, delayed
from pathlib import Path
from pyfaidx import Fasta

from .algorithms.hypergraph import (
    IRMM,
    extract_incidence_matrix, 
    remove_incidence_matrix
    )
from ._config import *
from .core import CountRE
from .utilities import run_cmd, listify, list_flatten

logger = logging.getLogger(__name__)


class Partitioner:
    """
    Single partiiton function.

    Params:
    --------
    count_re: str
        CountRE file.
    pairtable: str
        Pairs table file.
    k: int
        Groups of partition.
    maxLinkDensity: int, default 2
        Density threshold before marking contig as repetive.
    minREs: int, default 10
        Minimum number of RE sites in a contig to be clustered.
    nonInfomativeRation: int default 3
        Cutoff for recovering skipped contigs back into the clusters.
    quiet: bool, default False.
        Turn off all logging.
    
    Returns:
    --------
    object:
        object of Partitioner
    
    Examples:
    --------
    >>> p = Partitioner('Chr01.counts_AAGCTT.txt', 'Chr01.pairs.prune.txt', 4, quiet=True)
    >>> p.run()
    """
    def __init__(self, count_re, pairtable, k, 
                maxLinkDensity=2, minREs=10, 
                nonInformativeRatio=3, 
                quiet=False):
        
        self.count_re = count_re
        self.pairtable = pairtable
        self.k = k
        self.maxLinkDensity = maxLinkDensity
        self.minREs = minREs
        self.nonInformativeRatio = nonInformativeRatio

        self.quiet = quiet

    @staticmethod
    def get_partition_res(counts_file, k, targetDir):
        """
        get partition result from target directory.

        Params:
        --------
        counts_file: str
            countRE file.
        k: int
            Groups of partition.
        targetDir: str
            Target directory for collections.
        
        Returns:
        --------
        list:
            list of partition results
        
        Examples:
        --------
        >>> Partitioner.get_partition_res("Chr01.counts_AAGCTT.txt", 4, "./")
        ['./Chr01.counts_AAGCTT.4g1.txt',
         './Chr01.counts_AAGCTT.4g2.txt',
         './Chr01.counts_AAGCTT.4g3.txt',
         './Chr01.counts_AAGCTT.4g4.txt']
        """
        prefix = counts_file.replace('.txt', '')
        targets = glob.glob('{}/{}.{}g*txt'.format(targetDir, prefix, k))
        targets = sorted(targets)
        return targets

    @staticmethod 
    def partition(count_re, pairtable, k, 
                maxLinkDensity=2, minREs=10, 
                nonInformativeRatio=3, log=sys.stderr):

        cmd = ['allhic', 'partition', count_re, 
                pairtable, str(k), 
                '--maxLinkDensity', str(maxLinkDensity),
                '--minREs', str(minREs),
                '--nonInformativeRatio', str(nonInformativeRatio)]
        
        run_cmd(cmd, log, out2err=True)

        targets = Partitioner.get_partition_res(count_re, k, "./")

        return targets 

    @property
    def is_complete_partition(self):
        """
        complete partition is that partition k groups, instead of k+.

        Returns:
        --------
        bool,
            True or False of complete partititon
        
        Examples:
        --------
        >>> p.is_complete_partition
        True
        """
        
        if len(self.partition_res) == self.k:
            return True
        else:
            return False

    @property
    def lengths(self):
        """
        The length of each group.

        Returns:
        --------
        int:
            The length of each group.
        Examples:
        --------
        >>> p.lengths
        [119558114, 102471817, 102184770, 84028906]
        """
        return list(map(lambda x: CountRE(x).length, self.partition_res))
    
    @property
    def range_value(self):
        """
        Range value of partition results.

        Returns:
        --------
        int:
            Range value of partition results.

        Examples:
        --------
        >>> p.range_value
        35529208   
        """
        return max(self.lengths) - min(self.lengths)

    def run(self):
        log = os.devnull if self.quiet else sys.stderr
        self.partition_res = Partitioner.partition(self.count_re, 
                                    self.pairtable, 
                                    self.k,
                                    self.maxLinkDensity,
                                    self.minREs,
                                    self.nonInformativeRatio,
                                    log)
       

class AdaptivePartitioner(Partitioner):

    def __init__(self, count_re, pairtable, k, 
                maxLinkDensity_range=(2, 10),
                maxLinkDensity_step=2,
                minREs_range=(50, 300), 
                minREs_step=5,
                nonInformativeRatio=3,
                tmp_dir='partition_tmp',
                threads=10
                ):

        self.count_re = count_re
        self.prefix = self.count_re.replace(".txt", "")
        self.pairtable = pairtable
        self.k = k
        self.maxLinkDensity_tuple = maxLinkDensity_range
        self.maxLinkDensity_step = maxLinkDensity_step
        self.minREs_tuple = minREs_range
        self.minREs_step = minREs_step
        self.nonInformativeRatio = nonInformativeRatio
        
        self.tmp_dir = tmp_dir
        self.threads = threads

        logger.info(f"AdaptivePartition Params:\n"
                       f"\tmaxLinkDensity_tuple={str(self.maxLinkDensity_tuple)}\n "
                       f"       maxLinkDensity_step={self.maxLinkDensity_step}\n"
                       f"\tminREs_tuple={str(self.minREs_tuple)}\n"
                       f"\tminREs_step={self.minREs_step}\n"
                       f"\tnonInfomativeRatio={self.nonInformativeRatio}" )

    @staticmethod
    def find_best_partition(counts_file, pairs_file, k, minREs, maxLinkDensity):
        workdir = "{}_{}".format(minREs, maxLinkDensity)
        os.makedirs(workdir)
        os.chdir(workdir)
        os.link("../../" + counts_file, "./" + counts_file)
        os.link("../../" + pairs_file, "./" + pairs_file)
        p = Partitioner(counts_file, pairs_file, k, 
                            maxLinkDensity, minREs, quiet=True)        
        p.run()

        if not p.is_complete_partition:
            os.chdir("../")
            shutil.rmtree(workdir)
            return None

        lengths = p.lengths
        range_value = p.range_value

        result = (minREs, 
                    maxLinkDensity, 
                    ",".join(map(str, lengths)), 
                    range_value, 
                    range_value/sum(lengths))
 
        os.chdir("../")

        return result
    
    def run(self):

        with tempfile.TemporaryDirectory(prefix=self.tmp_dir, dir='./') as tmpDir:
            logger.info('Working on temporary directory: {}'.format(tmpDir))
            os.chdir(tmpDir)
            
            params = []
            for minREs in range(self.minREs_tuple[0], 
                                self.minREs_tuple[1], self.minREs_step):
                for maxLinkDensity in range(self.maxLinkDensity_tuple[0], 
                                            self.maxLinkDensity_tuple[1], 
                                            self.maxLinkDensity_step):
                    params.append((self.count_re, self.pairtable, 
                                    self.k, minREs, maxLinkDensity))
            
            logger.info("Finding the best partition result ...")
            results = Parallel(n_jobs=self.threads)(delayed(
                                            self.find_best_partition)(c, p, g, r, l)
                                                for c, p, g, r, l in params)

            results = filter(lambda x: x is not None, results)
            self.results = sorted(results, key=lambda x: x[4])
            try:
                best = self.results[0]
            except IndexError:
                logger.error("Couldn't found best results")
                return
            logger.info('Best result is [{}]'.format(best))

            os.system("cp {}_{}/*{}g*.txt ../".format(best[0], best[1], self.k))
            
            os.chdir("../")
            logger.info("Removed temporary directory.")


            self.save_partition_table()
            logger.info("Written all partition parameter results.")

    def save_partition_table(self):
        output_table = "{}.partition.table".format(self.prefix)
        with open(output_table, 'w') as out:
            print("\t".join(['#minREs', 'maxLinkDensity', 'groups',
                                'range', 'value']), file=out)
            for res in self.results:
                print("\t".join(map(str, res)), file=out)

    
class HyperPartition:
    """
    Method of contigs partition based on hypergraph partition.

    Params:
    --------
    pore_c_tables: list
        Pathes of pore_c_table.
    k: str
        Number of groups. Set to k1:k2 to perform multipartition.

    """
    def __init__(self, pore_c_tables, 
                    k,
                    fasta,
                    prune=None,
                    min_order=2, max_order=15,
                    min_alignments=500,
                    min_length=10000, threshold=0.01,
                    max_round=10, threads=4,
                    use_dask=False):
        
        self.pore_c_tables = listify(pore_c_tables)
        self.k = k

        self.fasta = fasta
        self.prune = prune 
        self.min_order = min_order
        self.max_order = max_order
        self.min_alignments = min_alignments
        self.min_length = min_length
        self.threshold = threshold
        self.max_round = max_round
        self.threads = threads
        self.use_dask = use_dask

        self.data = self.import_pore_c_table()
        self.contigs = self.get_contigs()
        self.H, self.vertices = self.get_hypergraph()
        
        self.filter_hypergraph()

        if prune:
            self.P_idx, self.prune_pair_df = self.get_prune_pairs()
        else:
            self.P_idx, self.prune_pair_df = None, None

    @property
    def vertices_idx(self):
        return dict(zip(self.vertices, 
                        range(len(self.vertices))))
    @property
    def contig_sizes(self):
        fasta = Fasta(self.fasta)
        contig_size_db = OrderedDict(list(map(lambda x: (x.name, len(x)), list(fasta))))
        
        return contig_size_db 

    def get_contigs(self):
        fasta = Fasta(self.fasta)
        contigs = list(map(lambda x: x.name, fasta))
        lengths = list(map(len, list(fasta)))
        length_db = dict(zip(contigs, lengths))
        contigs = sorted(contigs, key=lambda x: length_db[x], reverse=True)

        return contigs

    def get_prune_pairs(self):
        
        vertices_idx = self.vertices_idx

        pair_df = pd.read_csv(self.prune, sep='\t', 
                                header=None, index_col=None)
        pair_df[0] = pair_df[0].map(lambda x: vertices_idx.get(x, np.nan))
        pair_df[1] = pair_df[1].map(lambda x: vertices_idx.get(x, np.nan))
        pair_df = pair_df.dropna(axis=0)
        pair_df2 = pair_df.reindex(columns=[1, 0])
        pair_df2.columns = [0, 1]
        pair_df = pd.concat([pair_df, pair_df2], axis=0)
        pair_df = pair_df.drop_duplicates(subset=[0, 1])
  
        P_idx = [pair_df[0], pair_df[1]]
       
        return P_idx, pair_df

    def import_pore_c_table(self):
        logger.info("Loading Pore-C table ...")
        if len(self.pore_c_tables) == 1:
            if Path(self.pore_c_tables[0]).is_symlink():
                infile = os.readlink(self.pore_c_tables[0])
            else:
                infile = self.pore_c_tables[0]
            
            if self.use_dask:
                df = dd.read_parquet(infile, 
                                    columns=['read_idx', 'chrom',
                                            'start', 'end', 
                                            'pass_filter'],
                                    engine=PQ_ENGINE)
            else:
                df = pd.read_parquet(infile, 
                                     columns=['read_idx', 'chrom',
                                              'start', 'end', 
                                              'pass_filter'],
                                    engine=PQ_ENGINE)
            
            df = df.query("pass_filter == True")
            df_list = [df]
        else:
            infiles = []
            for i in self.pore_c_tables:
                if Path(i).is_symlink():
                    infiles.append(os.readlink(i))
                else:
                    infiles.append(i)

            if self.use_dask:
                df_list = list(map(lambda x: dd.read_parquet(
                                    x, columns=['read_idx', 'chrom', 
                                                'start', 'end', 
                                                'pass_filter'], 
                                    engine='pyarrow'), infiles))

                df = dd.concat(df_list, axis=0)
                # df['read_idx'] = df['read_name'].cat.as_known().cat.codes
                df = df.query("pass_filter == True").drop("pass_filter", axis=1)
                df_list = [df]
            else:
                df_list = list(map(lambda x: pd.read_parquet(
                                x, columns=['read_idx', 'chrom', 
                                            'start', 'end', 
                                            'pass_filter'], 
                                    engine='pyarrow',),
                                    #filters=[('pass_filter', '=', True)]),
                                    infiles))
                df_list = Parallel(n_jobs=self.threads)(delayed(
                                    lambda x: x.query("pass_filter == True")
                                                .drop("pass_filter", axis=1)
                                                )(i) for i in df_list)

        return df_list
    
    def get_hypergraph(self):
        from .algorithms.hypergraph import generate_hypergraph
        H, vertices = generate_hypergraph(self.data,
                                            self.min_order,
                                            self.max_order,
                                            self.min_alignments,
                                            self.threads,
                                            self.use_dask)

        return H, vertices

    def filter_hypergraph(self):
        ## remove too short contigs
        contig_sizes = self.contig_sizes
        vertices_idx = self.vertices_idx
        short_contigs = [i for i in contig_sizes
                                if contig_sizes[i] < self.min_length]

        short_contig_idx = []
        for i in short_contigs:
            try:
                short_contig_idx.append(vertices_idx[i])
            except KeyError:
                continue
        
        if len(short_contig_idx) == 0:
            return

        logger.info(f"Total {len(short_contig_idx)} contigs were removed, "
                        "because it's length too short.")
        self.H, _ = remove_incidence_matrix(self.H, short_contig_idx)
        self.vertices = np.delete(self.vertices, short_contig_idx)

    def merge_group(self):
        """
        merge group by signal.
        """
        self.k 
        self.K
        pass
    
    @classmethod
    def _multi_partition(self, k, prune_pair_df, H, threshold, max_round):
        """
        single function for multi_partition.
        """
        k = np.array(list(k))
        sub_H, _ = extract_incidence_matrix(H, k)

        sub_old2new_idx = dict(zip(k, range(len(k))))
        sub_new2old_idx = dict(zip(range(len(k)), k))

        # sub_vertices = list(np.array(self.vertices)[k])
        # sub_vertives_idx = dict(zip(sub_vertices, range(len(sub_vertices))))

        sub_prune_pair_df = prune_pair_df.reindex(list(permutations(k, 2)))
        sub_prune_pair_df = sub_prune_pair_df.dropna().reset_index()
    
        sub_prune_pair_df[0] = sub_prune_pair_df[0].map(lambda x: sub_old2new_idx[x])
        sub_prune_pair_df[1] = sub_prune_pair_df[1].map(lambda x: sub_old2new_idx[x])

        sub_P_idx = [sub_prune_pair_df[0], sub_prune_pair_df[1]]
        
        new_K = IRMM(sub_H, sub_P_idx, threshold, 
                        max_round, threads=1)
        
        ## remove single group
        new_K = filter(lambda x: len(x) > 1, new_K)
        new_K = list(map(lambda x: list(map(lambda y: sub_new2old_idx[y], x)), new_K))
        
        return new_K

    def multi_partition(self):
        """
        multiple partition for autopolyploid.
        """
        prune_pair_df = self.prune_pair_df.reset_index().set_index([0, 1])

        self.K = IRMM(self.H, None, self.threshold, 
                        self.max_round, threads=self.threads)
        self.K = filter(lambda x: len(x) > 1, self.K)

        # results = []
        args = []
        for k in self.K:
            args.append((k, prune_pair_df, self.H, self.threshold, self.max_round))
           
        results = Parallel(n_jobs=self.threads)(
                        delayed(self._multi_partition)
                                (i, j, k, l, m) for i, j, k, l, m in args)

        self.K = list_flatten(results)

    def single_partition(self):
        
        logger.info("Starting to cluster ...")
        self.K = IRMM(self.H, self.P_idx,
                        self.threshold, self.max_round,
                        threads=self.threads)
        logger.info("Cluster done.")


        self.K = filter(lambda x: len(x) > 1, self.K)
        self.K = sorted(self.K, key=lambda x: len(x), reverse=True)
        
        return self.K

    def to_cluster(self, output):
        idx_to_vertices = dict(zip(range(len(self.vertices)), self.vertices))
        clusters = list(map(lambda y: list(
                        map(lambda x: idx_to_vertices[x], y)), 
                        self.K))
        with open(output, 'w') as out:
            for i, group in enumerate(clusters, 1):
                print(f'group{i}\t{len(group)}\t{" ".join(group)}', 
                        file=out)

        logger.info(f"Successful output cluster results in `{output}`.")