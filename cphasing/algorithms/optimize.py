#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the algorithms of travelling salesman problem
"""

import argparse
import logging
import os
import os.path as op
import sys

import cooler
import igraph as ig
import numpy as np
import pandas as pd

from collections import defaultdict
from joblib import Parallel, delayed
from itertools import combinations
from pytools import natsorted
from scipy.sparse import tril


from cphasing.core import CountRE, Clm

logger = logging.getLogger(__name__)

class OldOptimize0:
    
    def __init__(self, contigs, clm, threads=10):
        self.contigs = natsorted(contigs)
        self.contig_idx = dict(zip(self.contigs, range(len(self.contigs))))

        self.clm = clm

        self.data = self.parse()

    def parse(self):
        df = self.clm.dk_df 

        res = dict(zip(df.to_dict('split')['index'], 
                        df.to_dict('split')['data']))
        
        return res
    
    def filter(self, as_idx=False):
        """
        only retain the max score in different orientation pairs

        """
        score_res = {}
        orientation_res = {}
        for pair in self.data:
            scores = self.data[pair]
            max_index = np.argmax(scores)
            max_score = scores[max_index]
            if max_score == 0:
                continue

            score_res[pair] = max_score
            orientation_res[pair] = max_index
        
        graph_df = pd.DataFrame(score_res, index=['weight']).T
        graph_df = graph_df.reset_index()
        graph_df.columns = ['source', 'target', 'weight']
        
        if as_idx:
            graph_df['source'] = graph_df['source'].map(self.contig_idx.get)
            graph_df['target'] = graph_df['target'].map(self.contig_idx.get)
            graph_df = graph_df.dropna()
            graph_df['source'] = graph_df['source'].astype('int')
            graph_df['target'] = graph_df['target'].astype('int')

        return graph_df, orientation_res
    
    def save(self, output):
        """
        save the score data.

        Params:
        --------
        output: str 

        Examples:
        --------
        >>> so2.save("out.score.txt")

        """
        with open(output, "w") as out:
            for pair in self.data:
                for i, s in enumerate(self.data[pair]):
                    _pair = list(zip(pair, self.orientations[i]))
                    _pair = list(map(lambda x: "".join(x), _pair))
                    print(f"{_pair[0]} {_pair[1]}\t{s:f}", file=out)
        logger.info(f"Successful output the score of contig pairs into {output}")


    def graph(self):
        """
        construct a graph
        """
        graph_df, orientation_res = self.filter(as_idx=True)
        # graph_df['weight'] = 1 / graph_df['weight']
        G = ig.Graph.DataFrame(graph_df, directed=False)
        
        return G

    def _dfs(self, start):
        
        pass 

    def dfs(self):
        pass

    def minimum_spanning_tree(self):
        pass

class OldOptimize:
    orientations = ["++", "+-", "-+", "--"]
    
    def __init__(self, contigs, clm, threads=10):
        self.contigs = natsorted(contigs)
        self.contig_idx = dict(zip(self.contigs, range(len(self.contigs))))

        self.clm = clm

        self.data = self.parse()

    def parse(self):
        df = self.clm.dk_df 

        res = dict(zip(df.to_dict('split')['index'], 
                        df.to_dict('split')['data']))
        
        return res
    
    def filter(self, as_idx=False):
        """
        only retain the max score in different orientation pairs

        """
        score_res = {}
        orientation_res = {}
        for pair in self.data:
            scores = self.data[pair]
            max_index = np.argmax(scores)
            max_score = scores[max_index]
            if max_score == 0:
                continue

            score_res[pair] = max_score
            orientation_res[pair] = max_index
        
        graph_df = pd.DataFrame(score_res, index=['weight']).T
        graph_df = graph_df.reset_index()
        graph_df.columns = ['source', 'target', 'weight']
        
        if as_idx:
            graph_df['source'] = graph_df['source'].map(self.contig_idx.get)
            graph_df['target'] = graph_df['target'].map(self.contig_idx.get)
            graph_df = graph_df.dropna()
            graph_df['source'] = graph_df['source'].astype('int')
            graph_df['target'] = graph_df['target'].astype('int')

        return graph_df, orientation_res
    
    def save(self, output):
        """
        save the score data.

        Params:
        --------
        output: str 

        Examples:
        --------
        >>> so2.save("out.score.txt")

        """
        with open(output, "w") as out:
            for pair in self.data:
                for i, s in enumerate(self.data[pair]):
                    _pair = list(zip(pair, self.orientations[i]))
                    _pair = list(map(lambda x: "".join(x), _pair))
                    print(f"{_pair[0]} {_pair[1]}\t{s:f}", file=out)
        logger.info(f"Successful output the score of contig pairs into {output}")


    def graph(self):
        """
        construct a graph
        """
        graph_df, orientation_res = self.filter(as_idx=True)
        # graph_df['weight'] = 1 / graph_df['weight']
        G = ig.Graph.DataFrame(graph_df, directed=False)
        
        return G

    def _dfs(self, start):
        
        pass 

    def dfs(self):
        pass

    def minimum_spanning_tree(self):
        pass


class SimpleOptimize:
    """

    """
    orientations = ["++", "+-", "-+", "--"]
    
    def __init__(self, contigs, cool, threads=10):
        self.cool = cool
        self.contigs = sorted(contigs, key=lambda x: self.cool.chromnames.index(x))
        self.contig_idx = dict(zip(self.contigs, range(len(self.contigs))))

        self.matrix = cool.matrix(balance=False, sparse=True)
        
        self.threads = threads 
        
        self.data = self.parse()

    @staticmethod
    def _parse(matrix, pair):
        """
        single pair parser
        """
        res = []
        contig1, contig2 = pair 
        sub_matrix = matrix.fetch(contig1, contig2)
        if sub_matrix.getnnz() == 0:
            return
        sub_matrix = sub_matrix.tocsr()
        l1, l2 = sub_matrix.shape
        d1 = int(np.ceil(l1 / 2))
        d2 = int(np.ceil(l2 / 2))

        k = l2 - l1 if l1 < l2 else 0 

        ## contig+ contig+
        c = tril(sub_matrix[d1:, :d2], k).sum()
        s = c / (d1 * d2)
        res.append(s)
        # print(f"{contig1}+ {contig2}+\t{s}")
        
        ## contig+ contig-
        c = tril(sub_matrix[:, ::-1][d1:, :d2], k).sum()
        s = c / (d1 * d2)
        res.append(s)
        # print(f"{contig1}+ {contig2}-\t{s}")

        ## contig- contig+
        c = tril(sub_matrix[::-1, :][d1:, :d2], k).sum()
        s = c / (d1 * d2)
        res.append(s)
        # print(f"{contig1}- {contig2}+\t{s}")

        ##contig- contig-
        c = tril(sub_matrix[::-1, ::-1][d1:, :d2], k).sum()
        s = c / (d1 * d2)
        res.append(s)
        # print(f"{contig1}- {contig2}-\t{s}")

        return pair, res 

    def parse(self):
        """
        calculate the score between different contig pairs

        |------||--------|        
        """
        
        res = defaultdict(list)

        args = []
        for pair in combinations(self.contigs, 2):
            
            args.append((self.matrix, pair))
        
        res = Parallel(n_jobs=self.threads)(
                delayed(SimpleOptimize._parse)(i, j) for i, j in args )
        
        ## remove None value
        res = dict(filter(lambda x: x is not None, res))
           
        return res

    def filter(self, as_idx=False):
        """
        only retain the max score in different orientation pairs

        """
        score_res = {}
        orientation_res = {}
        for pair in self.data:
            scores = self.data[pair]
            max_index = np.argmax(scores)
            max_score = scores[max_index]
            if max_score == 0:
                continue

            score_res[pair] = max_score
            orientation_res[pair] = max_index
        
        graph_df = pd.DataFrame(score_res, index=['weight']).T
        graph_df = graph_df.reset_index()
        graph_df.columns = ['source', 'target', 'weight']
        
        if as_idx:
            graph_df['source'] = graph_df['source'].map(self.contig_idx.get)
            graph_df['target'] = graph_df['target'].map(self.contig_idx.get)
            

        return graph_df, orientation_res
    
    def save(self, output):
        """
        save the score data.

        Params:
        --------
        output: str 

        Examples:
        --------
        >>> so2.save("out.score.txt")

        """
        with open(output, "w") as out:
            for pair in self.data:
                for i, s in enumerate(self.data[pair]):
                    _pair = list(zip(pair, self.orientations[i]))
                    _pair = list(map(lambda x: "".join(x), _pair))
                    print(f"{_pair[0]} {_pair[1]}\t{s}", file=out)
        logger.info(f"Successful output the score of contig pairs into {output}")


    def graph(self):
        """
        construct a graph
        """
        graph_df, orientation_res = self.filter(as_idx=True)
        # graph_df['weight'] = 1 / graph_df['weight']
        G = ig.Graph.DataFrame(graph_df, directed=False)
        
        return G

    def _dfs(self, start):
        
        pass 

    def dfs(self):
        pass

    def minimum_spanning_tree(self):
        pass


class SimpleOptimize2:
    """

    """
    orientations = ["++", "+-", "-+", "--"]
    
    def __init__(self, contigs, cool, threads=10):
        self.cool = cool
        self.contigs = sorted(contigs, key=lambda x: self.cool.chromnames.index(x))
        self.contig_idx = dict(zip(self.contigs, range(len(self.contigs))))

        self.matrix = cool.matrix(balance=False, sparse=True)
        
        self.threads = threads 
        
        self.data = self.parse()

    @staticmethod
    def _parse(matrix, pair):
        """
        single pair parser
        """
        res = []
        contig1, contig2 = pair 
        sub_matrix = matrix.fetch(contig1, contig2)
        if sub_matrix.getnnz() == 0:
            return
        sub_matrix = sub_matrix.tocsr()
        l1, l2 = sub_matrix.shape
        d1 = int(np.ceil(l1 / 2))
        d2 = int(np.ceil(l2 / 2))

        d_matrix = np.zeros((l1, l2))
        for i in range(l1):
            for j in range(l2):
                d_matrix[i, j] = i + j + 1
        d_matrix = d_matrix[::-1]
        
    
        ## contig+ contig+
        s = (sub_matrix[d1:, :d2] / d_matrix[d1:, :d2]).sum()
        
        res.append(s)
        # print(f"{contig1}+ {contig2}+\t{s}")
        
        ## contig+ contig-
        s = (sub_matrix[:, ::-1][d1:, :d2] / d_matrix[d1:, :d2]).sum()
        res.append(s)
        # print(f"{contig1}+ {contig2}-\t{s}")

        ## contig- contig+
        s = (sub_matrix[::-1, :][d1:, :d2] / d_matrix[d1:, :d2]).sum()
        res.append(s)
        # print(f"{contig1}- {contig2}+\t{s}")

        ##contig- contig-
        s = (sub_matrix[::-1, ::-1][d1:, :d2] / d_matrix[d1:, :d2]).sum()
        res.append(s)
        # print(f"{contig1}- {contig2}-\t{s}")

        return pair, res 

    def parse(self):
        """
        calculate the score between different contig pairs

        |------||--------|        
        """
        
        res = defaultdict(list)

        args = []
        for pair in combinations(self.contigs, 2):
            
            args.append((self.matrix, pair))
        
        res = Parallel(n_jobs=self.threads)(
                delayed(SimpleOptimize2._parse)(i, j) for i, j in args )
        
        ## remove None value
        res = dict(filter(lambda x: x is not None, res))
           
        return res

    def filter(self, as_idx=False):
        """
        only retain the max score in different orientation pairs

        """
        score_res = {}
        orientation_res = {}
        for pair in self.data:
            scores = self.data[pair]
            max_index = np.argmax(scores)
            max_score = scores[max_index]
            if max_score == 0:
                continue

            score_res[pair] = max_score
            orientation_res[pair] = max_index
        
        graph_df = pd.DataFrame(score_res, index=['weight']).T
        graph_df = graph_df.reset_index()
        graph_df.columns = ['source', 'target', 'weight']
        
        if as_idx:
            graph_df['source'] = graph_df['source'].map(self.contig_idx.get)
            graph_df['target'] = graph_df['target'].map(self.contig_idx.get)
            

        return graph_df, orientation_res
    
    def save(self, output):
        """
        save the score data.

        Params:
        --------
        output: str 

        Examples:
        --------
        >>> so2.save("out.score.txt")

        """
        with open(output, "w") as out:
            for pair in self.data:
                for i, s in enumerate(self.data[pair]):
                    _pair = list(zip(pair, self.orientations[i]))
                    _pair = list(map(lambda x: "".join(x), _pair))
                    print(f"{_pair[0]} {_pair[1]}\t{s}", file=out)
        logger.info(f"Successful output the score of contig pairs into {output}")

    def graph(self):
        """
        construct a graph
        """
        graph_df, orientation_res = self.filter(as_idx=True)
        # graph_df['weight'] = 1 / graph_df['weight']
        G = ig.Graph.DataFrame(graph_df, directed=False)
        
        return G

    def _dfs(self, start):
        
        pass 

    def dfs(self):
        pass

    def minimum_spanning_tree(self):
        pass



def test(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('count_re', 
            help='contig group in countRE table')
    pReq.add_argument('cool', 
            help='Path to cool file')
    
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    contigs = CountRE(args.count_re).contigs
    cool = cooler.Cooler(args.cool)
    so =  SimpleOptimize(contigs, cool)
    so.parse()

if __name__ == "__main__":
    test(sys.argv[1:])