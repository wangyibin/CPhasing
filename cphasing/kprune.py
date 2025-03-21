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
import tempfile

import cooler
import igraph as ig
import numpy as np
import msgspec
import pandas as pd 

from collections import defaultdict
from copy import deepcopy
# from deprecated.sphinx import deprecated
from itertools import product
from joblib import Parallel, delayed
from multiprocessing import Manager
from pathlib import Path
from scipy.sparse import triu



from .core import AlleleTable, CountRE
from .utilities import list_flatten, run_cmd

logger = logging.getLogger(__name__)


class ContigPairs(msgspec.Struct):
    ctg1: int 
    ctg2: int
    l1: int
    l2: int
    scores: list

    
class CandidateContacts(msgspec.Struct):
    data: list 


class KPrunerRust:
    def __init__(self, alleletable, contacts,
                    #count_re,
                    output="prune.contig.table", 
                    method="greedy", 
                    first_cluster=None,
                    whitelist=None,
                    norm_method="none",
                    threads=4, 
                    log_dir="logs"):
        self.alleletable = alleletable
        self.contacts = contacts
        # self.count_re = count_re
        self.output = output
        self.method = method 
        self.first_cluster = "none" if not first_cluster else first_cluster
        self.whitelist = "none" if not whitelist else whitelist
        self.norm_method = norm_method
        self.threads = threads
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

    def run(self):
        cmd = ["cphasing-rs", "kprune", self.alleletable,
               self.contacts, self.output, 
               "-n", self.norm_method, "-f", self.first_cluster,
                "-w", self.whitelist,
                "-t", str(self.threads), "-m", self.method]
        flag = run_cmd(cmd, log=f"{self.log_dir}/kprune.log")
        assert flag == 0, "Failed to execute command, please check log."
        logger.info(f"Successful output the allelic and cross-allelic information into `{self.output}`")


# @deprecated(version="v0.0.56")
class KPruner:
    """
    generate a prune list by kmer similarity allele table
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
    def __init__(self, alleletable, coolfile, 
                sort_by_similarity=True, count_re=None, 
                whitelist=None,
                chunksize=100000,
                threads=4):
        self.alleletable = AlleleTable(alleletable, sort=False, fmt='allele2')
        if sort_by_similarity:
            _alleletable = self.alleletable.data.sort_values(['similarity'], ascending=False)
            self.alleletable.data = _alleletable.query('similarity < 1.0')

        if whitelist:
            whitelist = set(whitelist)
            
            _alleletable = self.alleletable.data[
                self.alleletable.data[1].isin(whitelist) & 
                self.alleletable.data[2].isin(whitelist)]

            self.alleletable.data = _alleletable
             
        self.coolfile = coolfile 
        
        self.chunksize = chunksize
        self.threads = threads

        self.cool = cooler.Cooler(self.coolfile)
        self.contig2idx = dict(zip(self.cool.chromnames, range(len(self.cool.chromnames))))
        self.idx2contig = dict(zip(range(len(self.cool.chromnames)), self.cool.chromnames))
        self.pixels1 = self.cool.matrix(balance=False, join=True, 
                                       as_pixels=True)[:][['chrom1', 'chrom2', 'count']]
        
        
        self.contig_pairs = self.pixels1[['chrom1', 'chrom2']].applymap(self.contig2idx.get).values.tolist()
        self.contig_pairs = set(map(tuple, self.contig_pairs))

        self.count_re = CountRE(count_re, minRE=1) if count_re else None
        if count_re:
            self.re_count_db = self.count_re.data.to_dict()['RECounts']
            self.re_count_db = defaultdict(lambda :0, self.re_count_db)
            self.normalize_score(self.pixels1)

        self.pixels2 = self.pixels1.copy()
        self.pixels2.columns = ['chrom2', 'chrom1', 'count']
     
        self.pixels = pd.concat([self.pixels1, self.pixels2], axis=0)
        self.pixels['chrom1'] = self.pixels['chrom1'].map(self.contig2idx.get)
        self.pixels['chrom2'] = self.pixels['chrom2'].map(self.contig2idx.get)
        self.pixels.set_index(['chrom1', 'chrom2'], inplace=True)
        self.pixels.replace([np.inf, -np.inf], 0, inplace=True)
        del self.pixels1, self.pixels2
        gc.collect()

        self.allele_group_db = self.get_allele_group()
        self.score_db = self.get_score_db()

        self.allelic_prune_list = []
        self.weak_prune_list = []

    def get_allele_group(self):
        tmp_df = self.alleletable.data #.drop_duplicates(1)
        tmp_df = tmp_df.applymap(self.contig2idx.get)
        tmp_df = tmp_df.groupby(1, sort=False)[2].apply(lambda x: list(x))
        
        return tmp_df.to_dict()
        # return tmp_df.set_index(1).to_dict()[2]
    
    def normalize_score(self, pixel):

        pixel['RE1'] = pixel['chrom1'].map(self.re_count_db.get)
        pixel['RE2'] = pixel['chrom2'].map(self.re_count_db.get)
        pixel['count'] = pixel['count'] / ( pixel['RE1'] *  pixel['RE2'])
        pixel.drop(['RE1', 'RE2'], axis=1, inplace=True)

    def get_score_db(self):
        db = self.pixels.to_dict()['count']
        db = defaultdict(lambda :0, db)

        return db

    def remove_allelic(self):
        _allelic = self.alleletable.data[
                        self.alleletable.data[1] < self.alleletable.data[2]
                        ][[1, 2]]
        _allelic = _allelic.applymap(self.contig2idx.get).values.tolist()

        _allelic = list(map(tuple, _allelic))
        _allelic = list(set(_allelic).intersection(self.contig_pairs))
        logger.info(f"Removed {len(_allelic)} allelic contacts.")
        
        self.allelic_prune_list.extend(_allelic)

        # self.pixels.loc[self.pixels.reindex(_allelic).dropna().index] = 0
        # self.pixels = self.pixels.loc[self.pixels['count'] > 0]

        self.pixels = self.pixels.drop(_allelic, axis=0)
        self.contig_pairs = self.contig_pairs - set(_allelic)
    
    @staticmethod
    def is_strong_contact(allele_pair1, allele_pair2, score_db):
        
        contig_edges = [(a, b) for a in allele_pair1 for b in allele_pair2]
        scores = [score_db.get(i, 0) for i in contig_edges]
        
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
    def _remove_weak_with_allelic(ctg1, ctg2, l1, l2, scores):#ctg1, ctg2, allele_group_db, score_db):
        # try:
        #     allele1 = allele_group_db[ctg1]
        #     allele2 = allele_group_db[ctg2]
        # except KeyError:
        #     return None
        edges = list(product(range(l1), range(l1, l1 + l2)))
        # flag = KPruner.is_strong_contact2((ctg1, allele1), (ctg2, allele2), score_db)
        flag = KPruner.is_strong_contact2(l1, l2, edges, scores)
        if not flag:
            return (ctg1, ctg2)
        else:
            return None

    @staticmethod
    def _remove_weak_with_allelic_chunk(candidate_contact_file):
        with open(candidate_contact_file, 'rb') as fp:
            cc = msgspec.msgpack.decode(fp.read(), type=CandidateContacts)
        
        res = []
        for cp in cc.data:
            data = [cp["ctg1"], cp["ctg2"], cp["l1"], cp["l2"], cp["scores"]]
            res.append(KPruner._remove_weak_with_allelic(*data))

        return res 
        

    def remove_weak_with_allelic(self):
        args = []
        res = []
        score_db = self.score_db 
        allele_group_db = self.allele_group_db
        
        line_counts = 0
        tmp_file_idx = 0
        tmp_file_list = []

        if len(self.contig_pairs) > self.chunksize:
            with tempfile.TemporaryDirectory(prefix="kprune_tmp", dir="./") as tmpDir:
                logger.info(f"Working on a temporary directory `{tmpDir}`.")
                for ctg1, ctg2 in self.contig_pairs:
                    line_counts += 1
                    try:
                        allele1 = allele_group_db[ctg1]
                        allele2 = allele_group_db[ctg2]
                    except KeyError:
                        continue
                    
                    tmp_allele_pair1 = [ctg1, *allele1]
                    tmp_allele_pair2 = [ctg2, *allele2]
                    l1 = len(tmp_allele_pair1)
                    l2 = len(tmp_allele_pair2)
            
                    contig_edges = list(product(tmp_allele_pair1, tmp_allele_pair2))
                    scores = [score_db.get(i, .0) for i in contig_edges]

                    cp = ContigPairs(ctg1=ctg1, ctg2=ctg2, l1=l1, l2=l2, scores=scores)
                    args.append(cp)

                    if line_counts % self.chunksize == 0:
                        tmp_file = f"{tmpDir}/{tmp_file_idx}.bin"
                        with open(tmp_file, 'wb') as out:
                            out.write(msgspec.msgpack.encode(CandidateContacts(data=args)))
                        
                        tmp_file_list.append(tmp_file)
                        tmp_file_idx += 1
                        args=[]
                else:
                    tmp_file = f"{tmpDir}/{tmp_file_idx}.bin"
                    with open(tmp_file, 'wb') as out:
                        out.write(msgspec.msgpack.encode(CandidateContacts(data=args)))
                    
                    tmp_file_list.append(tmp_file)
                    args=[]


                res = Parallel(n_jobs=self.threads)(
                        delayed(KPruner._remove_weak_with_allelic_chunk)(i,)
                            for i in tmp_file_list)
        
                res = list_flatten(res)

            logger.info(f"Removed temporary directory: `{tmpDir}`.")
        else:
            for ctg1, ctg2 in self.contig_pairs:
                line_counts += 1
                try:
                    allele1 = allele_group_db[ctg1]
                    allele2 = allele_group_db[ctg2]
                except KeyError:
                    continue
                
                tmp_allele_pair1 = [ctg1, *allele1]
                tmp_allele_pair2 = [ctg2, *allele2]
                l1 = len(tmp_allele_pair1)
                l2 = len(tmp_allele_pair2)
        
                contig_edges = list(product(tmp_allele_pair1, tmp_allele_pair2))
                

                scores = [score_db.get(i, 0) for i in contig_edges]

                args.append([ctg1, ctg2, l1, l2, scores])
            
            res = Parallel(n_jobs=self.threads)(
                delayed(KPruner._remove_weak_with_allelic)(i, j, k, l, m)
                        for i, j, k, l, m in args)
    
        res = list(filter(lambda x: x is not None, res))
        weak_contacts = len(res)
        logger.info(f"Removed {weak_contacts} weak contacts with allelic.")
        self.weak_prune_list.extend(res)
    
    def save_prune_list(self, output, symmetric=False):
        self.allelic_prune_list = list(map(lambda x: tuple(
                                        map(lambda y: self.idx2contig[y], x)), 
                                        self.allelic_prune_list))
        allele_data = self.alleletable.data.set_index([1, 2])
        
        allelic_res = allele_data.reindex(self.allelic_prune_list) 
        allelic_res = allelic_res.reset_index()
        allelic_res.columns = self.alleletable.AlleleHeader2[2:]
        allelic_res['type'] = 0 

        if len(self.weak_prune_list) > 0:
            weak_res = pd.DataFrame(self.weak_prune_list)
            weak_res = weak_res.applymap(self.idx2contig.get)
            weak_res.columns = self.alleletable.AlleleHeader2[2:4]
            weak_res['type'] = 1 
            weak_res['mz1'] = 0
            weak_res['mz2'] = 0
            weak_res['mzShared'] = 0
            weak_res['similarity'] = 0
        
        if symmetric:
            allelic_res2 = allelic_res.rename(columns={'contig1': 'contig2',
                                                        'conrig2': 'conrig1'})
            if len(self.weak_prune_list) > 0:
                weak_res2 = weak_res.rename(columns={'contig1': 'contig2',
                                                            'conrig2': 'conrig1'})
                res = pd.concat([allelic_res, allelic_res2, 
                                    weak_res, weak_res2], axis=0)  
            else:
                res = pd.concat([allelic_res, allelic_res2], axis=0)
        else:
            if len(self.weak_prune_list) > 0:
                res = pd.concat([allelic_res, weak_res], axis=0)
            else:
                res = allelic_res
        res.to_csv(output, sep='\t', header=None, index=None)
        
      
    def run(self):
        self.remove_allelic()
        self.remove_weak_with_allelic()

class KPruneHyperGraph: 
    """
    generate a prune list by kmer similarity allele table
        rewrite of allhic prune to adapt the hypergraph.

    Params:
    --------
    alleletable: str
        alleletable file
    A: csr_matrix
        csr_matrix of hypergraph from clique expansion 
    contig2idx: dict
        contig to idx
        
    Examples:
    --------
    >>> kph = KPruneHyperGraph('allele.table', A, contig2idx)
    run allelic and weak with allelic contacts removing
    >>> kph.run()
    save the prune list
    >>> kph.save_prune_list('prune.contig.list')
    """
    def __init__(self, alleletable, A, contig2idx):
        self.alleletable = alleletable
        _alleletable = alleletable.data.copy()
        _alleletable[1] = _alleletable[1].map(contig2idx.get)
        _alleletable[2] = _alleletable[2].map(contig2idx.get)
        _alleletable = _alleletable.dropna()
        _alleletable[1] = _alleletable[1].astype(int)
        _alleletable[2] = _alleletable[2].astype(int)

        self.alleletable_idx = _alleletable

        self.A = A 
        self.A.setdiag(0)
        self.A = self.A.tocsr()
        
        self.contig2idx = contig2idx 

        self.allele_group_db = self.get_allele_group()
        print(self.allele_group_db)

        # self.pixeles = pd.Series.sparse.from_coo(self.A + self.A.T)
        self.contig_pairs = list(zip(triu(A).row, triu(A).col))

        self.contig_pairs = set(map(tuple, self.contig_pairs))


        self.allelic_prune_list = []
        self.weak_prune_list = []


    def get_allele_group(self):
        tmp_df = self.alleletable.data #.drop_duplicates(1)
        tmp_df = tmp_df.groupby(1, sort=False)[2].apply(lambda x: list(x)[0])
        tmp_df = tmp_df[tmp_df.index.isin(self.contig2idx)]
        tmp_df = tmp_df[tmp_df.isin(self.contig2idx)]
        
        tmp_df.index = tmp_df.index.map(self.contig2idx.get)
        tmp_df = tmp_df.map(self.contig2idx.get)
    

        return tmp_df.to_dict()

    def get_score_db(self):
        db = self.A.todok()
        db = defaultdict(lambda :0, db.items())

        return db 

    def remove_allelic(self):
        _allelic = self.alleletable_idx[self.alleletable_idx[1] < self.alleletable_idx[2]
                            ][[1, 2]]
        _allelic = _allelic.values.tolist()
        _allelic = list(map(tuple, _allelic))
        _allelic = list(set(_allelic).intersection(self.contig_pairs))
        logger.info(f"Removed {len(_allelic)} allelic contacts.")
        
        self.allelic_prune_list.extend(_allelic)

        _allelic_zipped = list(zip(*_allelic))
        self.A[_allelic_zipped[0], _allelic_zipped[1]] = 0
        
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
        
        flag = KPruneHyperGraph.is_strong_contact((ctg1, allele1), (ctg2, allele2), score_db)
        if not flag:
            return (ctg1, ctg2)
        else:
            return None

    def remove_weak_with_allelic(self):
        args = []
        res = []
        score_db = self.get_score_db()
        allele_group_db = self.allele_group_db

        for ctg1, ctg2 in self.contig_pairs:  

            
            res.append(KPruneHyperGraph._remove_weak_with_allelic(ctg1, ctg2, allele_group_db, score_db))

        res = list(filter(lambda x: x is not None, res))
        weak_contacts = len(res)
        logger.info(f"Removed {weak_contacts} weak contacts with allelic.")
        self.weak_prune_list.extend(res)
    
    def get_prune_table(self, symmetric=False):
        allele_data = self.alleletable_idx.set_index([1, 2])
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

        return res

    def run(self):
        self.remove_allelic()
        self.remove_weak_with_allelic()
