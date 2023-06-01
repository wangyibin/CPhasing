#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
generate prune contig pair list
"""


import logging
import gc
import os
import os.path as op
import sys

import cooler
import igraph as ig
import numpy as np
import pandas as pd 

from collections import defaultdict
from joblib import Parallel, delayed
from scipy.sparse import triu

from .core import AlleleTable 

logger = logging.getLogger(__name__)


class KPruner:
    """
    generate a prune list by allele table
        rewrite of allhic prune to adapt the cool file.

    Params:
    --------
    alleletable: str
        alleletable file
    coolfile: str
        path to whole contig contacts cool file
    threads: int
        number of threads

    Examples:
    --------
    >>> kp = KPruner('allele.table', 'sample.whole.contig.cool')
    run allelic and weak with allelic contacts removing
    >>> kp.run()
    save the prune list
    >>> kp.save_prune_list('prune.contig.list')
    
    """
    def __init__(self, alleletable, coolfile, threads=4):
        self.alleletable = AlleleTable(alleletable, sort=False, fmt='allele2')
        self.coolfile = coolfile 
        self.threads = threads

        self.allele_group_db = self.get_allele_group()

        self.cool = cooler.Cooler(self.coolfile)
        self.pixels1 = self.cool.matrix(balance=False, join=True, 
                                       as_pixels=True)[:][['chrom1', 'chrom2', 'count']]
        
        self.contig_pairs = self.pixels1[['chrom1', 'chrom2']].values.tolist()
        self.contig_pairs = set(map(tuple, self.contig_pairs))

        self.pixels2 = self.pixels1.copy()
        self.pixels2.columns = ['chrom2', 'chrom1', 'count']
     
        self.pixels = pd.concat([self.pixels1, self.pixels2], axis=0)
        self.pixels.set_index(['chrom1', 'chrom2'], inplace=True)
        self.pixels.replace([np.inf, -np.inf], 0, inplace=True)
        del self.pixels1, self.pixels2
        gc.collect()

        self.score_db = self.get_score_db()

        self.allelic_prune_list = []
        self.weak_prune_list = []

    def get_allele_group(self):
        tmp_df = self.alleletable.data.drop_duplicates(1)
        #return tmp_df.groupby(1)[2].apply(list).to_dict()
      
        return tmp_df.set_index(1).to_dict()[2]
    
    def get_score_db(self):
        db = self.pixels.to_dict()['count']
        db = defaultdict(lambda :0, db)

        return db
        
    def remove_allelic(self):
        _allelic = self.alleletable.data[
                        self.alleletable.data[1] < self.alleletable.data[2]
                        ][[1, 2]].values.tolist()
        _allelic = list(map(tuple, _allelic))
        _allelic = list(set(_allelic).intersection(self.contig_pairs))
        logger.info(f"Removed {len(_allelic)} allelic contacts.")
        
        self.allelic_prune_list.extend(_allelic)

        self.pixels.loc[self.pixels.reindex(_allelic).dropna().index] = 0
        self.pixels = self.pixels.loc[self.pixels['count'] > 0]
        self.contig_pairs = self.contig_pairs - set(_allelic)
    
    @staticmethod
    def is_strong_contact(allele_pair1, allele_pair2, score_db):
        
        edges = [(a, b) for a in allele_pair1 for b in allele_pair2]
        scores = [score_db[i] for i in edges]

        edges = [(0, 2), (0, 3), (1, 2), (1, 3)]

        ## igraph quickly (0.000331s)
        g = ig.Graph.Bipartite([0, 0, 1, 1], edges)
        g.es['weight'] = scores

        matching = g.maximum_bipartite_matching(weights='weight')
        
        if matching.match_of(0) == 2:
            return True
        else:
            return False

    @staticmethod
    def _remove_weak_with_allelic(ctg1, ctg2, allele_group_db, score_db):
        try:
            allele1 = allele_group_db[ctg1]
            allele2 = allele_group_db[ctg2]
        except KeyError:
            return None
        
        flag = KPruner.is_strong_contact((ctg1, allele1), (ctg2, allele2), score_db)
        if not flag:
            return (ctg1, ctg2)
        else:
            return None

    def remove_weak_with_allelic(self):
        args = []
        res = []
        score_db = self.score_db 
        allele_group_db = self.allele_group_db
        for ctg1, ctg2 in self.contig_pairs:
            # args.append((ctg1, ctg2, allele_group_db, score_db))
            res.append(KPruner._remove_weak_with_allelic(ctg1, ctg2, allele_group_db, score_db))
        # res = Parallel(n_jobs=self.threads)(
        #         delayed(KPruner._remove_weak_with_allelic)(i, j, k, l)
        #             for i, j, k, l in args )
        
        res = list(filter(lambda x: x is not None, res))
        weak_contacts = len(res)
        logger.info(f"Removed {weak_contacts} weak contacts with allelic.")
        self.weak_prune_list.extend(res)
    
    def save_prune_list(self, output, symmetric=False):
        allele_data = self.alleletable.data.set_index([1, 2])
        
        allelic_res = allele_data.reindex(self.allelic_prune_list) 
        allelic_res = allelic_res.reset_index()
        allelic_res.columns = self.alleletable.AlleleHeader2[2:]
        allelic_res['type'] = 0 

        weak_res = pd.DataFrame(self.weak_prune_list)
        weak_res.columns = self.alleletable.AlleleHeader2[2:4]
        weak_res['type'] = 1 
        weak_res['mz1'] = 0
        weak_res['mz2'] = 0
        weak_res['mzShared'] = 0
        weak_res['similarity'] = 0
        
        if symmetric:
            allelic_res2 = allelic_res.rename(columns={'contig1': 'contig2',
                                                        'conrig2': 'conrig1'})
            weak_res2 = weak_res.rename(columns={'contig1': 'contig2',
                                                        'conrig2': 'conrig1'})
            res = pd.concat([allelic_res, allelic_res2, 
                                weak_res, weak_res2], axis=0)    
        else:
            res = pd.concat([allelic_res, weak_res], axis=0)
        res.to_csv(output, sep='\t', header=None, index=None)
        
      

    def run(self):
        self.remove_allelic()
        self.remove_weak_with_allelic()
    
    

    