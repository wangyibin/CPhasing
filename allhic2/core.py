#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
core functions of allhic2
"""


import argparse
import logging
import os
import os.path as op
import sys

import numpy as np
import pandas as pd 

from collections import OrderedDict
from pathlib import Path

logger = logging.getLogger(__name__)


class CountRE:
    """
    countRE table

    Params:
    --------

    Returns:
    --------
    Object of countRE

    Examples:
    --------
    >>> cr = CountRE("sample.bwa_mem_GATC.txt")
    """
    def __init__(self, countRE, minRE=3):
        self.filename = countRE
        self.minRE = minRE

        self.parse()

    def parse(self):
        self.data = pd.read_csv(self.filename, sep='\t', 
                                    header=0, index_col=0)
        self.data = self.data[self.data['RECounts'] >= self.minRE]
        
        logger.info(f'Load countRE `{self.filename}` (minRE={self.minRE}).')

    @property
    def header(self):
        """
        header of countRE

        Returns:
        list of header

        Examples:
        --------
        >>> cr.header
        ['#Contig', 'RECounts', 'Length']
        """
        return self.data.reset_index().columns.tolist()
    
    @property
    def contigs(self):
        """
        contigs of countRE

        Returns:
        --------
        list of contigs

        Examples:
        --------
        >>> cr.contigs
        ['utg001', 'utg002', ...]
        """
        return self.data.index.tolist()

    @property
    def ncontigs(self):
        """
        number of contigs in countRE

        Returns:
        --------
        `int`, number of contigs

        Examples:
        --------
        >>> cr.ncontigs
        100
        """
        return len(self.contigs)
    
    @property
    def length(self):
        """
        length of all contigs

        Returns:
        --------
        Length of all contigs

        Examples:
        >>> cr.length
        552000
        """
        return self.data['Length'].sum()


class AlleleLine:
    """
    Allele line

    Params:
    --------
    line: `str`, line of file

    Examples:
    --------
    >>> line = 'Chr01\tgene001\tctg1\tctg2'
    >>> al = AlleleLine(line)
    >>> al.chrom
    "Chr01"
    
    """

    def __init__(self, line):
        self.line = line.strip()
        line_list = line.strip().split()
        self.chrom = line_list[0]
        self.gene = line_list[1]
        self.contigs = line_list[2:]
        self.n = len(self.contigs)
    
    def __str__(self):
        return self.line 

class AlleleTable:
    """
    Allele table

    Params:
    --------
    infile: `str`, input file for allele table
    sort: `bool`, whether to sort contigs per row [True]
    Examples:
    --------
    >>> infile = 'Allele.ctg.table'
    >>> at = AlleleTable(infile)

    """
    def __init__(self, infile, sort=True):
        self.filename = infile
        self.sort = sort
        if not Path(self.filename).exists():
            logger.error(f'No such file of `{self.filename}`.')
            sys.exit()

        logger.info(f'Load AlleleTable: `{self.filename}`.')
        self.check()
        self.get_data()

    def check(self):
        """
        check file to resolve uneuqal fields in different lines

        Examples:
        --------
        >>> at.check()
        """
        with open(self.filename, 'r') as fp:
            n = 0
            for line in fp:
                if not line.strip():
                    continue
                al = AlleleLine(line)
                if al.n > n:
                    n = al.n

        self.columns = list(range(n + 2))

    def get_data(self, sort=True):
        def sort_row(row):
            row_filter = filter(lambda x: not pd.isna(x), row)
            return pd.Series(sorted(row_filter), name=row.name)

        df = pd.read_csv(self.filename, sep='\t', header=None, 
                            index_col=0, names=self.columns)
        df.index = df.index.astype('category')
        df = df.dropna(how='all', axis=1)
        df = df.drop(1, axis=1)
    
        df = df.apply(sort_row, axis=1) if sort else df
        df.columns = list(range(1, len(df.columns) + 1))

        self.data = df

    def filter(self, contigs):
        """
        remove contigs that not in input contigs list

        Params:
        --------
        contigs: list or list-like 
            contigs list
        
        Returns:
        --------
        Remove the contigs that not in input ocntigs

        Examples:
        --------
        >>> at.data
        Chr01 gene001 ctg001 ctg002 ctg003 ctg004
        >>> contigs = ['ctg001', 'ctg002', 'ctg003']
        >>> at.filter(contigs)
        >>> at.data 
        Chr01 gene001 ctg001 ctg002 ctg003
        """
        self.data = self.data[self.data.isin(contigs)]
        self.data.index.rename(0, inplace=True)

        if self.data.empty:
            logger.warn("After filter AlleleTable is empty")


    @property
    def groups(self):
        _groups = OrderedDict()
        _groupby = self.data.groupby(0)
        for chrom, item in _groupby:
            tmp = item.values.flatten()
            _groups[chrom] = tmp[~pd.isna(tmp)]

        return _groups

    @property
    def chromnames(self):
        """
        list of chromosomes
        """
        return self.groups.keys()

    @property
    def nchroms(self):
        """
        number of chromosomes
        """
        return len(self.chromnames)

    @property
    def contigs(self):
        """
        list of contigs
        """
        _contigs = set()
        for i in self.groups.values():
            _contigs.update(i)

        return list(_contigs)
    
    @property
    def ncontigs(self):
        """
        number of contigs
        """
        return len(self.contigs)
    
    def get_shared(self):
        
        return

    def save(self, output):
        self.data.to_csv(output, sep='\t', header=None, index=True)


class PairTable:
    """
    Object of allhic pairs file

    Params:
    --------
    infile: str
        input file of pairs table
    symmetrix: bool
        symmetric pairs [True]

    Examples:
    --------
    Constructing PairTable from a txt file

    >>> pt = PairTable('sample.pairs.txt')

    Detect the symmetric of PairTable
    >>> pt.is_symmetric()
    True

    """ 
    def __init__(self, infile, symmetric=True):
        self.filename = infile
        if not Path(self.filename).exists():
            logger.error(f'No such file of `{self.filename}`.')
            sys.exit()
            
        self.data = self.import_pairs()
        if symmetric is True:
            self.symmetric_pairs()


    def import_pairs(self):

        logger.info('Load pair table: `{}`'.format(self.filename))
        df = pd.read_csv(self.filename, header=0, index_col=None,
                        sep='\t')
        df = df.astype(
            {'Contig1': 'category', 'Contig2': 'category',
            '#X': 'int32', 'Y': 'int32', 'RE1': 'int64', 
            'RE2': 'int64', 'Label': 'category'}
        )

        return df
    
    @property
    def header(self):
        """
        return the header of pairtable

        Examples:
        >>> pt.header
        ['#X',
        'Y',
        'Contig1',
        'Contig2',
        'RE1',
        'RE2',
        'ObservedLinks',
        'ExpectedLinksIfAdjacent',
        'Label']
        """
        return self.data.columns.tolist()
    
    def symmetric_pairs(self):
        """
        symmetric pairs dataframe
        """
        _pairs_df = self.data.copy()
        _pairs_df['Contig1'], _pairs_df['Contig2'] = self.data['Contig2'], self.data['Contig1']
        _pairs_df['#X'], _pairs_df['Y'] = self.data['Y'], self.data['#X']
        try:
            _pairs_df['Group1'], _pairs_df['Group2'] = self.data['Group2'], self.data['Group1']
        except KeyError:
            pass
        _pairs_df['RE1'], _pairs_df['RE2'] = self.data['RE2'], self.data['RE1']
        
        self.data = pd.concat([self.data, _pairs_df])
        self.data = self.data.sort_values(by=['#X', 'Y'])
        self.data.reset_index(drop=True, inplace=True)
        
    def is_symmetric(self):
        """
        Detect the pairtable whether symmetry

        Returns:
        --------
        bool
            returns a boolean indicating whether the pairtable is symmetric

        Examples:
        --------
        >>> pt.is_symmetric()
        True
        """
        return (self.data['#X'] > self.data['Y']).any()


    def save(self, output):
        """
        save pairtable to file

        Params:
        --------
        output: `str` or `stdout`, output

        Examples:
        --------
        >>> pt.save('pairs.txt')
        """
        self.data.to_csv(output, sep='\t', header=True, index=False)
        logger.info(f'Output pairs table to `{output}`')
        

class ClusterGroup:
    def __init__(self, line):
        line_list = line.strip().split("\t")
        self.line = line.strip()
        self.group = line_list[0]
        self.nContigs = line_list[1]
        self.Contigs = line_list[2].split()
    
    def __repr__(self):
        return self.line 

    __str__ = __repr__
    
class ClusterTable:
    """
    Cluster table from partition results.

    Params:
    -------
    filename: str
        input file of cluster table

    Examples:
    --------
    >>> filename = "sample.clusters.txt"
    >>> ct = ClusterTable(filename)

    """
    def __init__(self, filename):
        self.filename = filename
        logger.info(f'Load ClusterTable: `{self.filename}`.')
        self.get_data()

    def parse(self):
        with open(self.filename, 'r') as fp:
            for line in fp:
                if line[0] == "#":
                    continue
                yield ClusterGroup(line)

    def get_data(self):
        self.data = OrderedDict()
        for i in self.parse():
            self.data[i.group] = i.Contigs

    @property
    def contig_groups(self):
        _res = OrderedDict()
        for group in self.data:
            for contig in self.data[group]:
                _res[contig] = group
        
        return _res
       
    @property
    def groups(self):
        _group = [i.group for i in self.parse()]
        return _group

    @property
    def contigs(self):
        _contigs = []
        for i in self.data:
            _contigs.extend(self.data[i])

        return _contigs
    @property
    def ncontigs(self):
        """
        Total number of contigs

        Returns:
        --------
        int
            number of contigs

        Examples:
        --------
        >>> ct.ncontigs
        1000
        """
        return len(self.contigs)
    
    @property
    def contig_groups_dict(self):
        return dict(self.contig_groups)
    
    @property
    def stat(self):
        """
        the statistics of clusterTable

        Returns:
        --------
        dictionary
            returns a dictionary of the number of contigs in 
            each groups
        
        Examples:
        --------
        >>> ct.stat
        OrderedDict([('Chr01g1', 184),
             ('Chr01g2', 180),
             ('Chr01g3', 189),
             ('Chr01g4', 179)])

        """
        _stat = OrderedDict()
        for i in self.data:
            _stat[i] = len(self.data[i])
        
        return _stat

    def from_frame(self, dataframe):
        """
        create clusterTable from DataFrame

        Params:
        --------
        dataframe: `pd.DataFrame`, dataframe of cluster results

        Examples:
        --------
        >> dataframe
            group1  group2
        0   NaN     utg001
        1   utg002  utg003
        >>> ct.from_frame(dataframe)
        """
        _data = dataframe.to_dict('list')
        for group in _data:
            _data[group] = list(set(filter(
                lambda x: pd.isna(x) is not True, _data[group])))
        
        self.data = _data

    def to_countRE(self, countRE):
        """
        convert clusterTable into different groupings countRE

        Params:
        --------
        countRE: `CountRE`, object of countRE

        Examples:
        --------
        >>> countRE = CountRE("sample.bwa_mem.GATC.txt")
        >>> ct.to_countRE(countRE)
        """
        self.CountRE = CountRE(countRE)
        prefix = Path(countRE).stem
        for group in self.groups:
            _contigs = self.data[group]
            tmp_df = self.CountRE.data.loc[_contigs]
            tmp_df.to_csv(f'{prefix}.{group}.txt', sep='\t',
                            header=True, index=True)
            logger.info(f'Output countRE file into `{prefix}.{group}.txt`')

    def save(self, output):
        """
        save cluster table to output file

        Params:
        output: `str` or `stdout`, output

        Examples:
        --------
        >>> ct.save('cluster.txt')
        """
        if isinstance(output, str):
            output = open(output, 'w')
        for group in self.data:
            print(group, f"{len(self.data[group])}", 
                    " ".join(self.data[group]), 
                    sep='\t', file=output)
        logger.info(f'Output cluster table to `{output}`')

    def get_max_group(self):
        res = 0
        for i in self.stat:
            if res < self.stat[i]:
                res = self.stat[i]
        
        return res 

    def __getitem__(self, key):
        return self.data[key]


class clmLine():
    def __init__(self, line):
        self.line = line.strip()
        self.line_list = self.line.split("\t")
        self.ctg1, self.ctg2 = self.line_list[0].split()
        self.strand1 = self.ctg1[-1]
        self.strand2 = self.ctg2[-1]
        self.ctg1 = self.ctg1[:-1]
        self.ctg2 = self.ctg2[:-1]
        self.count = int(self.line_list[1])
        self.distances = list(map(int, self.line_list[2].split()))
    
    @property
    def dk(self):
        
        reciprocal_values = list(map(lambda x: 1/(x + 1), self.distances))
        return sum(reciprocal_values)
    
    def __str__(self):
        return self.line
 
class clm(object):
    def __init__(self, clmfile, mem_cache='.'):
        
        self.file = clmfile
        self.parse()
        # self.mem_cache = mem_cache
        # self.memory = Memory(mem_cache, verbose=0)
        # self.dk_df = self.memory.cache(self._dk_df)

    def parse(self):
        self.data = {}
        with open(self.file, 'r') as fp:
            for line in fp:
                cl = clmLine(line)
                if (cl.ctg1, cl.ctg2) not in self.data:
                    self.data[(cl.ctg1, cl.ctg2)] = []
                
                self.data[(cl.ctg1, cl.ctg2)].append(cl.dk)
    
    @property
    def dk_df(self):
        res_dk_df = pd.DataFrame(self.data)
        res_dk_df = res_dk_df.T 
        res_dk_df.columns = self.strands
        
        return res_dk_df

    @property
    def strands(self):
        return [('+', '+'), ('+', '-'),
                ('-', '+'), ('-', '-')]
      