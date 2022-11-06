#!/usr/bin/env python
# -*- coding:utf-8 -*-


import logging
import os
import os.path as op
import sys
import shutil

import glob
import tempfile
import pandas as pd

from collections import OrderedDict 
from joblib import Parallel, delayed

from .core import CountRE
from .utilities import run_cmd, listify

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
        p = Partitioner(counts_file, pairs_file, k, maxLinkDensity, minREs, quiet=True)        
        p.run()

        if not p.is_complete_partition:
            os.chdir("../")
            shutil.rmtree(workdir)
            return None

        lengths = p.lengths
        range_value = p.range_value

        result = (minREs, maxLinkDensity, ",".join(map(str, lengths)), range_value, 
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
    def __init__(self, pore_c_tables):
        self.pore_c_tables = listify(pore_c_tables)
        self.data = self.import_pore_c_table()
        self.edges = self.get_hyperedges()
        self.H = self.get_hypergraph()
    
    def get_hypergraph(self):
        import hypernetx as hnx
        H = hnx.Hypergraph(dict(enumerate(self.edges)), use_nwhy=True)
        return H

    def import_pore_c_table(self):
        if len(self.pore_c_tables) == 1:
            df = pd.read_parquet(self.pore_c_tables[0])
        else:
            df_list = list(map(pd.read_parquet, self.pore_c_tables))
            df = pd.concat(df_list, axis=0)
        
        return df 
    
    def partition(self):
        import hypernetx.algorithms.hypergraph_modularity as hmod 
        HG = hmod.precompute_attributes(self.H)
        self.K = hmod.kumar(HG)

        return self.K

    def get_hyperedges(self):
        df2 = self.data.set_index('read_name')
        df2 = df2.loc[(self.data.groupby('read_name')['chrom'].nunique() >= 2)]
        edges = df2.groupby('read_name')['chrom'].unique().values.tolist()
        edges = list(set(map(tuple, map(sorted, edges))))

        return edges

    def to_cluster(self, output):
        with open(output, 'w') as out:
            for i, group in enumerate(self.K):
                print(f'{i}\t{len(group)}\t{" ".join(group)}', 
                        file=out)

