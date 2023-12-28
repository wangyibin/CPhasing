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
import tempfile

import cooler
import random
import igraph as ig
import gc
import numpy as np
import pandas as pd

from collections import defaultdict
from joblib import Parallel, delayed
from itertools import combinations, permutations
from pathlib import Path
from pytools import natsorted
from scipy.sparse import tril


from ..core import CountRE, ClusterTable, Clm
from ..utilities import choose_software_by_platform, run_cmd

logger = logging.getLogger(__name__)

class AllhicOptimize:

    def __init__(self, clustertable, count_re, clm, 
                    fasta=None, output="groups.agp", 
                    tmp_dir='scaffolding_tmp', threads=4):
        self.clustertable = ClusterTable(clustertable)
        self.count_re = CountRE(count_re, minRE=1)
        self.clm = pd.read_csv(clm, sep='\t', header=None, index_col=0)
        self.fasta = Path(fasta).absolute() if fasta else None
        self.output = output
        self.tmp_dir = tmp_dir 
        self.threads = threads 

        self.allhic_path = choose_software_by_platform("allhic")

    @staticmethod
    def extract_count_re(group, contigs, count_re):
        tmp_df = count_re.data.reindex(contigs)
        tmp_df.to_csv(f"{group}.txt", sep='\t', header=None)

        return f"{group}.txt"

    @staticmethod
    def extract_clm(group, contigs, clm):
        
        contig_pairs = list(permutations(contigs, 2))
        contig_with_orientation_pairs = []
        for pair in contig_pairs:
            for strand1, strand2 in [('+', '+'), ('+', '-'),
                                    ('-', '+'), ('-', '-')]:
                contig_with_orientation_pairs.append(f"{pair[0]}{strand1} {pair[1]}{strand2}")
        
        tmp_df = clm.reindex(contig_with_orientation_pairs).dropna().astype({1: int})
        tmp_df.to_csv(f"{group}.clm", sep='\t', header=None)

        return f"{group}.clm"

    @staticmethod
    def run_allhic_optimize(allhic_path, count_re, clm):
        cmd = [allhic_path, "optimize", count_re, clm]
        run_cmd(cmd, log=os.devnull, out2err=True)
        return count_re.replace(".txt", ".tour")

    @staticmethod
    def _run(allhic_path, count_re, clm, workdir):
        os.chdir(workdir)
        tmp_res = AllhicOptimize.run_allhic_optimize(allhic_path, count_re, clm)

        return tmp_res
    
    def run(self):
        from ..cli import build
        with tempfile.TemporaryDirectory(prefix=self.tmp_dir, dir='./') as tmpDir:
            logger.info('Working on temporary directory: {}'.format(tmpDir))
            os.chdir(tmpDir)
            workdir = os.getcwd()
            args = []
            for group in self.clustertable.data.keys():
                contigs = self.clustertable.data[group]
                tmp_clm = AllhicOptimize.extract_clm(group, contigs, self.clm)
                tmp_count_re = AllhicOptimize.extract_count_re(group, contigs, self.count_re)
                args.append((self.allhic_path, tmp_count_re, tmp_clm, workdir))
            
            del self.clm
            gc.collect()
            
            
            Parallel(n_jobs=min(len(args), self.threads))(delayed(
                        self._run)(i, j, k, l) for i, j, k, l in args)

            if not self.fasta:
                os.system(f"cp *tour ../")
            else:
                try:
                    build.main(args=[str(self.fasta), "--only-agp", "-oa", self.output], prog_name='build')
                except SystemExit as e:
                    exc_info = sys.exc_info()
                    exit_code = e.code
                    if exit_code is None:
                        exit_code = 0
                    
                    if exit_code != 0:
                        raise e
                    
                os.system(f"cp {self.output} ../")
            
            logger.info("Removed temporary directory.")
             
            os.chdir("../")

        logger.info("Done")


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
        mst_tree = self.G.spanning_tree()


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
        d1 = float(int(np.ceil(l1 / 2)))
        d2 = float(int(np.ceil(l2 / 2)))

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
    
    def __init__(self, contigs, cool, method="so", threads=10):
        self.cool = cool
        self.contigs = contigs #sorted(contigs, key=lambda x: self.cool.chromnames.index(x))
        self.contig_idx = dict(zip(self.contigs, range(len(self.contigs))))
        self.idx_to_contig = dict(zip(range(len(self.contigs)), self.contigs))

        self.matrix = cool.matrix(balance=False, sparse=True)
        
        self.method = method
        self.threads = threads 
        
        self.data = self.parse()
        self.score_df, self.orientation_res = self.filter(mode='score')


    @staticmethod
    def _parse_by_so(matrix, pair):
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

    @staticmethod
    def _parse_by_so2(matrix, pair):
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
        
        _parse = SimpleOptimize2._parse_by_so if self.method == "so" else SimpleOptimize2._parse_by_so2

        res = Parallel(n_jobs=self.threads)(
                delayed(_parse)(i, j) for i, j in args )
        
        ## remove None value
        res = dict(filter(lambda x: x is not None, res))
           
        return res

    def filter(self, mode='graph', as_idx=False):
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
        
        if mode == 'graph':
            graph_df = pd.DataFrame(score_res, index=['weight']).T
            graph_df = graph_df.reset_index()
            graph_df.columns = ['source', 'target', 'weight']
            
            if as_idx:
                graph_df['source'] = graph_df['source'].map(self.contig_idx.get)
                graph_df['target'] = graph_df['target'].map(self.contig_idx.get)
       
            return graph_df, orientation_res
        
        elif mode == 'score':
            score_df = pd.DataFrame(score_res, index=['score']).T
            score_df = score_df.reset_index()
            score_df.columns = ['contig1', 'contig2', 'score']
            score_df2 = score_df.rename(columns={"contig1": "contig2",
                                                    "contig2": "contig1"})
            score_df = pd.concat([score_df, score_df2], axis=0)

            score_df.set_index(["contig1", "contig2"], inplace=True)

            return score_df, orientation_res
        else:
            raise ValueError("mode must in {'graph', 'score'}")
        
        
    
    def save_score(self, output):
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

        graph_df, self.orientation_info = self.filter(as_idx=True)
    
        graph_df2 = graph_df[['target', 'source', 'weight']]
        graph_df2.columns = ['source', 'target', 'weight']

        graph_df = pd.concat([graph_df, graph_df2], axis=0)
        graph_df['weight'] = 1 / graph_df['weight']
        # graph_df['weight'] = max(graph_df['weight']) - graph_df['weight']
        # graph_df['weight'] = graph_df['weight'].astype(int)
        G = ig.Graph.DataFrame(graph_df, directed=False)
        

        return G, graph_df

    def get_global_score(self, order_list):
        order_pairs = [(i, j) for i, j in zip(order_list[:-1], order_list[1:])]
    
        scores = self.graph_df.reindex(order_pairs, fill_value=0)
        
        score = sum(scores['weight'])

        return score
    
    def bfs(self, start):
        res = self.G.bfs(start)[0]
        
        if len(res) == 1:
            return None
        # res = list(map(self.idx_to_contig.get, res[0]))
        return  self.get_global_score(res), res

    def dfs(self, start):
        
        res = self.G.dfs(start)[0]
        
        if len(res) == 1:
            return None
        # res = list(map(self.idx_to_contig.get, res[0]))
        return  self.get_global_score(res), res

    def order(self, algorithm='dfs'):
        
        args = []
        for idx in range(len(self.contigs)):
            args.append(idx)
        
        res = Parallel(n_jobs=self.threads)(
                delayed(getattr(self, algorithm))(i) for i in args)

        res = list(filter(lambda x: x is not None, res))
        res = sorted(res, key=lambda x: x[0])

        print(list(map(self.idx_to_contig.get, res[0][1])))
        return list(map(self.idx_to_contig.get, res[0][1]))

    def tsp_order(self):
        import networkx as nx
        tsp = nx.approximation.traveling_salesman_problem

        G = self.G.to_networkx()
        # SA_tsp = nx.approximation.simulated_annealing_tsp
        # method = lambda G, wt: SA_tsp(G, "greedy", weight=wt, move='1-0', temp=500)
        # path = tsp(G, cycle=False, method=method)
        
        path = tsp(G, cycle=False)

        return list(map(self.idx_to_contig.get, path))
    
    def lkh_order(self):
        import elkai
        import networkx as nx 
        G = self.G.to_networkx()
        M = nx.adjacency_matrix(G).todense()
        print(M)
        path = elkai.solve_int_matrix(M)

        return list(map(self.idx_to_contig.get, path))
    

    def nn_tsp(self, contigs, score_df, start=0):

        start_contig = "Chr5.ctg1" #contigs[start]
        candicate_contigs = set(set(contigs) - {start_contig})
        tour = [start_contig]
        print(start_contig)
        added_contigs = set() 
        added_contigs.add(start_contig)

        while candicate_contigs:
            next_contig = score_df.loc[tour[-1]].idxmax().values[0]
            tmp_df = score_df.loc[tour[-1]].sort_values(by=['score'], ascending=False)
            if next_contig in added_contigs:
                # print(score_df.loc[tour[-1]]['score'].nlargest(2, keep='all'))
                next_contig = score_df.loc[tour[-1]]['score'].nlargest(2, keep='all').index[1]
                i = 1
                while next_contig in added_contigs:
                    i += 1
                    next_contig = tmp_df.index[i]

            if next_contig not in contigs:
                continue
            
            print(next_contig)
            added_contigs.add(next_contig)
            tour.append(next_contig)

            candicate_contigs.remove(next_contig)

        return tour

    def orientation(self):
        
        for i, j in zip(self.ordering[:-1], self.ordering[1:]):
            try:
                print(i, j, self.orientations[self.orientation_info[(i, j)]])
            except:
                print(i, j)

    def save(self, output):
        with open(output, "w") as out:
            print("\n".join(self.ordering), file=out)

    def minimum_spanning_tree(self):
        pass


    
class SAOptimizer:
    """
    Simulated Annealing
    """
    def __init__(self) -> None:
        pass 

    def cooling_schedule(self, t):
        return 0.99 * t
    
    def acceptance_probability(self, energy, new_energy, temperature):
        if new_energy < energy:
            return 1.0
        return np.exp((energy - new_energy) / temperature)
    
    def simulated_annealing(self, graph, start, iterations=1000):
        """
        Simulated Annealing

        Params:
        -------
        graph: np.array
            adjacency matrix
        start: list 
            start node
        iterations: int
            number of iterations
        
        Returns:
        --------
        current_path: list
        current_cost: int

        Example:
        --------
        >>> graph = np.array([[0, 1, 2, 3],
                                [1, 0, 4, 5],
                                [2, 4, 0, 6],  
                                [3, 5, 6, 0]])
        >>> start = 0
        >>> iterations = 1000
        >>> sa = SA()
        >>> current_path, current_cost = sa.simulated_annealing(graph, start, iterations)
        >>> print(current_path, current_cost)
        [0, 1, 2, 3] 12
        """
        # initial state
        current_path = [start]
        current_cost = 0
        unvisited_nodes = set(range(len(graph))) - {start}
    
        # iteration
        for i in range(iterations):
            # cooling
            T = self.cooling_schedule(iterations / i)
    
            # random neighbour
            current_node = current_path[-1]
            next_node = random.sample(unvisited_nodes, 1)[0]
    
            # calculate current and neighbour path cost
            current_cost += graph[current_node][next_node]
            new_cost = current_cost - graph[current_node][next_node] + graph[next_node][current_node]
    
            # accept neighbour or not
            if self.acceptance_probability(current_cost, new_cost, T) > random.random():
                current_cost = new_cost
                current_path.append(next_node)
                unvisited_nodes.remove(next_node)
    
        return current_path, current_cost




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