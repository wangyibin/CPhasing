#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
core functions of C-Phasing
"""

import gc
import logging
import os
import os.path as op
import re
import shutil
import sys

import dask.dataframe as dd
import numpy as np
import pandas as pd 
import pyranges as pr 
import tempfile

from collections import Counter, OrderedDict
from dask.distributed import Client
from joblib import Parallel, delayed
from itertools import combinations
from pathlib import Path
from subprocess import check_call

from ._config import PQ_ENGINE, PQ_VERSION
from .utilities import list_flatten, tail, xopen

logger = logging.getLogger(__name__)


class CountRE:
    """
    countRE table

    Params:
    --------
    countRE: str
        countRE file
    minRE: int
        minimum RE of countRE [3]

    Returns:
    --------
    object:
        object of countRE

    Examples:
    --------
    >>> cr = CountRE("sample.bwa_mem_GATC.txt")
    """
    def __init__(self, countRE, minRE=3):
        self.filename = countRE
        self._path = Path(self.filename)
        self.suffix = re.match(r'.*(counts_.*.txt)$', 
                                self.filename).groups()[0]
        self.minRE = minRE

        self.parse()

    def parse(self):
        self.data = pd.read_csv(self.filename, sep='\t', 
                                    header=0, index_col=0,
                                    dtype={'RECounts': int, 
                                            'Length': int})
        self.data = self.data[self.data['RECounts'] >= self.minRE]
        
        logger.info(f'Load count RE `{self.filename}` (minRE={self.minRE}).')

    @property
    def header(self):
        """
        header of countRE

        Returns:
        list:
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
        Length of all contigs

        Returns:
        --------
        int:
            Length of all contigs

        Examples:
        --------
        >>> cr.length
        552000
        """
        return self.data['Length'].sum()

    def get_group_length(self, contigs: list) -> int:
        """
        Get length of a group contigs

        Returns:
        --------
        int:
           Length of contigs

        Examples:
        --------
        >>> contigs = ['utg001', 'utg002']
        >>> cr.get_group_length(contigs)
        523330
        """
        return self.data.loc[contigs]['Length'].sum()

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

        logger.info(f'Load allele table: `{self.filename}`.')
        self.check()
        self.get_data(sort)

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
                            index_col=0, names=self.columns,
                            usecols=self.columns)
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
    

    def get_shared(self, symmetrix=True):
        """
        get contig pairs sharing informations

        Params:
        --------
        symmetrix: bool
            symmetric pairs [True]

        Returns:
        --------
        pd.DataFrame:
            shared allelic number between two contigs

        Examples:
        --------
        >>> at.get_shared()
                                                            0
        level_0                  level_1                     
        utg000008l_137000_173999 utg000517l_691000_1025999  1
        utg000008l_173999_267139 utg000123l_932000_1434999  1
                                utg000517l_691000_1025999  3
        """
        def _func(row):
            return list(combinations(row.dropna(), 2))
        
        share_table = self.data.apply(_func, axis=1).values.flatten()
        share_table = list_flatten(share_table)
        share_table = Counter(share_table)
        share_table = pd.Series(share_table)
        if symmetrix:
            share_table = share_table.reset_index()
            share_table2 = share_table.copy()
            share_table2.columns = ['level_1', 'level_0', 0]
            share_table = pd.concat([share_table, share_table2])
            share_table = share_table.sort_values(by=['level_0', 'level_1'], axis=0)
            share_table = share_table.set_index(['level_0', 'level_1'])

        return share_table
        

    def save(self, output):
        """
        save allele tabel to output

        Params:
        --------
        output: str
            output file

        Examples:
        --------
        >>> at.save('output.allele.table')
        """
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
    index_contig: bool
        index by contig [True]

    Examples:
    --------
    Constructing PairTable from a txt file

    >>> pt = PairTable('sample.pairs.txt')

    Detect the symmetric of PairTable
    >>> pt.is_symmetric()
    True

    """ 
    HEADER = ['#X',
        'Y',
        'Contig1',
        'Contig2',
        'RE1',
        'RE2',
        'ObservedLinks',
        'ExpectedLinksIfAdjacent',
        'Label']

    def __init__(self, infile, symmetric=True, index_contig=True):
        self.filename = infile
        if not Path(self.filename).exists():
            logger.error(f'No such file of `{self.filename}`.')
            sys.exit()
            
        self.data = self.import_pairs()
        if symmetric:
            self.symmetric_pairs()

        if index_contig:
            self.data = self.data.set_index(['Contig1', 'Contig2'])

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
    
    @property
    def Contig1(self):
        """
        contigs in Contig1

        Returns:
        --------
        np.array:
            contigs in Contig1
        """
        return pd.unique(self.data.index.get_level_values(0))

    @property
    def Contig2(self):
        """
        contigs in Contig2

        Returns:
        --------
        np.array:
            contigs in Contig2
        """
        return pd.unique(self.data.index.get_level_values(1))

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

    def unsymmetric_pairs(self):
        """
        unsymmetric pairs dataframe

        """
        _pairs_df = self.data.copy()
        _pairs_df = _pairs_df[_pairs_df['#X'] < _pairs_df['Y']]
        self.data = _pairs_df.sort_values(by=['#X', 'Y'])

    def is_symmetric(self):
        """
        Detect the pairtable whether symmetry

        Returns:
        --------
        bool
            returns a boolefffan indicating whether the pairtable is symmetric

        Examples:
        --------
        >>> pt.is_symmetric()
        True
        """
        return (self.data['#X'] > self.data['Y']).any()

    def get_contacts(self, contig1, contig2):
        """
        get contacts dataframe between two contig lists, 
        the pairs table must be symmetric and index by Contig1
        and Contig2

        Params:
        --------
        contig1: list or list-like
            list of contigs
        contig2: list or list-like
            list of contigs
        
        Returns:
        --------
        pd.DataFrame:
            dataframe of two contig lists

        Examples:
        --------
        >>> contig1 = ['utg001']
        >>> contig2 = ['utg002', 'utg003']
        >>> pt.get_contacts(contig1, contig2)
                            #X  Y  RE1   RE2  ObservedLinks  ExpectedLinksIfAdjacent     Label
            Contig1 Contig2                                                                                    
            utg001 utg002   3  20  1247  532  46             426.6                       ok
                   utg003   3  29  1247  558  28             430.0                       ok
        """
        assert self.is_symmetric, "The pair table must be symmetric"

        tmp_df = self.data.reindex((i, j) for i in contig1 for j in contig2)
        tmp_df = tmp_df.dropna(how='all', axis=0)
        
        return tmp_df
        
    def get_normalized_contact(self, contig1, contig2):
        """
        get normalized contact between two contig lists

        Params:
        --------
        contig1: list or list-like
            list of contigs
        contig2: list or list-like
            list of contigs
        
        Returns:
        --------
        float:
            normalized contact

        Examples:
        --------
        >>> contig1 = ['utg001']
        >>> contig2 = ['utg002', 'utg003']
        >>> pt.get_normalized_contact(contig1, contig2)
        0.05
        """
        tmp_df = self.get_contacts(contig1, contig2)
        
        if tmp_df.empty:
            return 0
        
        L = tmp_df['ObservedLinks'].sum()
        C = tmp_df['RE2'].sum()
        Cu = tmp_df['RE1'][0]
        S = L/(Cu + C)

        return S 
    
    def normalize(self):
        
        self.data = (self.data
                    .eval("NormalizedLinks = ObservedLinks / (RE1 + RE2)")
        )
    
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
        try:
            self.data[self.HEADER]
        except KeyError:
            self.data.reset_index(inplace=True)
            self.data = self.data[self.HEADER]
        
        if self.is_symmetric:
            self.unsymmetric_pairs()
        
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
        logger.info(f'Load cluster table: `{self.filename}`.')
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
        >>> countRE = CountRE()
        >>> ct.to_countRE(countRE)
        """
        self.CountRE = CountRE(countRE, minRE=1)
        prefix = Path(countRE).stem
        for group in self.groups:
            _contigs = self.data[group]
            tmp_df = self.CountRE.data.loc[_contigs]
            tmp_df.to_csv(f'{prefix}.{group}.txt', sep='\t',
                            header=True, index=True)
            logger.info(f'Output countRE file into `{prefix}.{group}.txt`')

    def to_assembly(self, fasta, output):
        from pyfaidx import Fasta
        fasta = Fasta(fasta)

        contig_idx_db = dict(zip(self.contigs, range(len(self.contigs))))
        for i, contig in enumerate(self.contigs, 1):
            print(f">{contig} {i} {len(fasta[contig])}", file=output)
        
        for group in self.groups:
            _contigs = self.data[group]
            _idx = list(map(lambda x: contig_idx_db[x], _contigs))
            _idx = list(map(lambda x: int(x)+1, _idx))
            print(" ".join(map(str, _idx)), file=output)

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


class TourSingle:
    """
    Single contig for tour

    Params:
    --------
    contig: str or tuple
        contig end with orientation

    Returns:
    --------
    str:
        conitg end without orientation

    Examples:
    --------
    from string
    >>> contig = 'utg0001+'
    >>> TourSingle(contig)
    'utg0001+'

    from tuple
    >>> contig = ('utg001', '+')
    >>> TourSingle(contig)
    'utg0001+'
    """
    def __init__(self, contig):
        if isinstance(contig, str):
            self._contig = contig
            self.contig = contig[:-1]
            self.orient = contig[-1]
        else:
            self._contig = ''.join(contig)
            self.contig = contig[0]
            self.orient = contig[1]
        
        assert self.orient in ['+', '-'], \
                'Contig must end in a orientation, and must in {"+", "-"}'

    @property
    def data(self):
        return self.contig, self.orient

    def __repr__(self):
        return self._contig

    __str__ = __repr__

class Tour:
    """
    Object of tour file.

    Params:
    --------
    infile: str
        input file of tour.
    
    Returns:
    --------
    object:
        object of tour

    Examples:
    --------
    >>> tour = Tour('Chr01.tour')
    """
    def __init__(self, infile):
        self.filename = infile
        self.group = Path(self.filename).stem

        self.line = tail(self.filename, 1)[0].strip()
        self.data = list(map(TourSingle, self.line.split(' ')))

    def __str__(self):
        return ' '.join(map(str, self.data))

    def __len__(self):
        return len(self.data)
    
    def __iter__(self):
        return iter(self.to_tuples())

    def __contains__(self, contig):
        return contig in self.contigs

    def __reversed__(self):
        """
        reverse tour data.

        Params:
        --------
        None

        Returns:
        --------
        None

        Examples:
        --------
        >>> tour.data
        ['utg001+', 'utg002-']
        >>> tour.reverse()
        >>> tour.data
        ['utg002+', 'utg001-']
        """
        _data = self.to_tuples()[::-1]
        _new_data = []
        for contig, orient in _data:
            if orient == '+':
                orient = '-'
            else:
                orient = '+'
            _new_data.append((contig, orient))
        
        self.data = list(map(lambda x: TourSingle(''.join(x)), _new_data))

    @property
    def contigs(self, rm_orient=True):
        """
        list of contigs

        Params:
        --------
        rm_orient: bool
            remove orientation in the end of contigs

        Returns:
        --------
        list:
            list of contigs

        Examples:
        --------
        >>> tour.contigs
        ['utg001', 'utg002]
  
        """
        if rm_orient:
            _contigs = list(map(lambda x: x.contig, self.data))
        else:
            _contigs = list(map(lambda x: x._contig, self.data))
        return _contigs
    
    @property
    def orients(self):
        """
        list of orient

        Params:
        --------
        None

        Returns:
        --------
        list:
            list of orient
        
        Examples:
        --------
        >>> tour.orients
        ['+', '-']
        """
        return list(map(lambda x: x.orient, self.data))

    @property
    def ncontigs(self):
        """
        number of contigs.

        Params:
        --------
        None

        Returns:
        --------
        int:
            number of contigs

        Examples:
        --------
        >>> tour.ncontigs
        2
        """
        return len(self.contigs)

    def to_tuples(self):
        """
        Contert tour data to tuples of contig and orient.

        Params:
        --------
        None
        
        Returns:
        --------
        list:
            A list of contig and orient.

        Examples:
        --------
        >>> tour.data
        ['utg001+', 'utg002-']
        >>> tour.to_tuples()
        ('utg001', '+'), ('utg002', '-')
        """
        return list(map(lambda x: x.data, self.data))

    @classmethod
    def from_tuples(self, tuples):
        """
        Create or update Tour from a tuples.

        Params:
        --------
        tuples: list-like
            A tuple list of contig and orient.
        
        Returns:
        --------
        None

        Examples:
        --------
        >>> tour = Tour
        >>> tour.from_tuples(('utg001', '+'), ('utg002', '-'))
        >>> tour.data
        ['utg001+', 'utg002-']
        """
        self.data = list(map(lambda x: TourSingle(''.join(x)), tuples))
    
    def to_agp(self, fasta, gap_length=100):
        """
        Convert tour to agp component.

        Params:
        --------
        fasta: pyfaidx.Fasta
            contig-level fasta
        
        Returns:
        --------
        pd.DataFrame:
            dataframe of agp component.
        
        Examples
        --------
        >>> tour.to_agp(fasta)
        """
        res = []
        start = 1
        i = 1
        old_end = 0
        for contig, orient in self:
            length = len(fasta[contig])
            
            start = old_end + 1
            end = start + length - 1
            
            res.append((self.group, start, end, i, 'W', contig, 
                    1, length, orient))
            
            i += 1
            if i <= len(self) * 2 - 1:
                res.append((self.group, end + 1, end + gap_length,
                        i, 'U', gap_length, 'contig', 'yes', 'map'))
                old_end = end + gap_length
                i += 1
            
        return pd.DataFrame(res)

    def get_fasta(self, fasta):
        """
        Get sequences from a contig-level fasta by tour contigs.

        Params:
        --------
        fasta: pyfaidx.Fasta
            contig-level fasta

        Returns:
        --------
        list:
            list of FastaRecord

        Examples:
        --------
        >>> fasta = Fasta('contig.fata')
        >>> tour.get_fasta(fasta)
        ['AAGCTT', 'ACGAAAGCATGG', 'AATTGGGAAGCA', 
         'GATACAGATAGA', 'TTTTAAAAATAG']
        """
        _seqs = []
        for contig, orient in self:
            if orient == '+':
                _seqs.append(str(fasta[contig]))
            else:
                _seqs.append(fasta[contig][::-1].complement.seq)
        
        return _seqs
    
    def reverse(self):
        """
        reverse tour data.

        Params:
        --------
        None

        Returns:
        --------
        None

        Examples:
        --------
        >>> tour.data
        ['utg001+', 'utg002-']
        >>> tour.reverse()
        >>> tour.data
        ['utg002+', 'utg001-']
        """
        self.__reversed__()

class PairHeader:
    """
    Object for 4DN pairs header.

    Params:
    --------
    header: list
        line list of pairs header.
    
    """
    def __init__(self, header):
        self.header = header
        self.pairs_format = self._get_pairs_format()
        self.shape = self._get_shape()
        self.chromsize = self._get_chromsize()
        self.columns = self._get_columns()
    
    @property
    def chromnames(self):
        return self.chromsize.keys()

    def _get_pairs_format(self):
        _res = list(filter(lambda x: x[:2] == "##", self.header))
        
        return _res[0].replace("## pairs format ", "") if _res else None

    def _get_shape(self):
        _res = list(filter(lambda x: x[:6] == "#shape", self.header))

        return _res[0].split(":")[1].strip() if _res else None

    def _get_chromsize(self):
        _res = list(filter(lambda x: x[:10] == "#chromsize", self.header))
        
        _res = [record.split(":")[1].split() for record in _res]
        db = OrderedDict()
        for chrom, size in _res:
            db[chrom] = int(size)

        return db

    def _get_columns(self):
        _res = list(filter(lambda x: x[:8] == "#columns", self.header))
        if not _res:
            return None
        else:
            _res = _res[0].split(":")[1].strip()
            _res = _res.split(" ")
            
            if "chrom1" not in _res or "chrom2" not in _res:
                try:
                    _res[_res.index("chr1")] = "chrom1"
                    _res[_res.index("chr2")] = "chrom2"
                except ValueError:
                    return None
            
            return _res

    def update(self):
        self.header = []
        self.header.append(f'## pairs format {self.pairs_format}')
        self.header.append(f'#shape: {self.shape}')
        for chrom, size in self.chromsize.items():
            self.header.append(f'#chromsize: {chrom} {size}')
        _columns = " ".join(self.columns)
        self.header.append(f'#columns: {_columns}')
    
    def from_file(self, filename):
        self.header = []
        with xopen(filename, 'r') as fp:
            for line in fp:
                if line.startswith('#'):
                    self.header.append(line.strip())
                else:
                    break
        self.pairs_format = self._get_pairs_format()
        self.shape = self._get_shape()
        self.chromsize = self._get_chromsize()
        self.columns = self._get_columns()

    def from_chromsize_file(self, chromsize, 
                            columns=None):
        self.pairs_format = "1.0"
        self.shape = "upper triangle"
        self.chromsize = dict(i.strip().split()[:2] 
                            for i in open(chromsize) if i.strip())
        self.columns = columns if columns else ["readID", "chrom1", 
                                                "pos1", "chrom2",
                                                "pos2", "strand1", 
                                                "strand2"]
        self.update()

    def save(self, filename):
        with xopen(filename, 'w') as out:
            out.write(str(self) + "\n")

    def __str__(self):
        return "\n".join(self.header)
    

class Pairs:
    """
    Object of 4DN pairs file.

    Params:
    --------
    pairs: str
        Path of 4DN pairs file.
    
    Returns:
    --------
    object

    Examples:
    --------
    >>> p = Pairs('Lib.pairs')
    """
    def __init__(self, pairs):
        self.file = Path(pairs)
        self.filename = self.file.name
        self.header = self._get_header()
        self.data = self.read_table()

    def _get_header(self):
        ph = PairHeader([])
        ph.from_file(self.file)

        return ph

    @property
    def chromsize(self):
        return self.header.chromsize
    
    @property
    def chromnames(self):
        return list(self.chromsize.keys())
    
    @property
    def columns(self):
        return self.header.columns

    @property
    def meta(self):
        chrom_dtype = pd.CategoricalDtype(self.chromnames, ordered=True)
        strand_dtype = pd.CategoricalDtype(["+", "-"], ordered=False)

        meta = {
            "readID": str,
            "chrom1": chrom_dtype,
            "chrom2": chrom_dtype,
            "strand1": strand_dtype,
            "strand2": strand_dtype
        }
        return meta 

    def read_table(self):
        import dask.dataframe as dd
        
        logger.info(f'Load pairs file of `{self.filename}`.')
        df = dd.read_csv(self.file, sep='\t', comment='#', 
                            names=self.header.columns, 
                            header=None,  dtype=self.meta)

        return df

    @staticmethod
    def _chrom2contig(df, n, source_gr, columns, meta=None):
        def _pos_to_range(df):
            """convert single pos to range."""
            df = (df.reset_index(drop=True)
                    .assign(
                        start1=lambda x: x['pos1'] - 1,
                        end1=lambda x: x['pos1'],
                        start2=lambda x: x['pos2'] - 1,
                        end2=lambda x: x['pos2']
                    )
            )
            return df 

        df = df.get_partition(n).compute()
        df = df.reset_index()
        df = _pos_to_range(df)

        query_df = df[['chrom1', 'start1', 'end1', 'index']]
        query_df.columns = ['Chromosome', 'Start', 'End', 'Index']
        query_gr = pr.PyRanges(query_df)
        res_df1 = query_gr.join(source_gr).df
        res_df1 = (res_df1.set_index('Index')
                    .assign(
                        Start=lambda x: x['Start'] - x['Start_b'] + 1,
                    )
        )

        query_df = df[['chrom2', 'start2', 'end2', 'index']]
        query_df.columns = ['Chromosome', 'Start', 'End', 'Index']
        query_gr = pr.PyRanges(query_df)
        res_df2 = query_gr.join(source_gr).df
        res_df2 = (res_df2.set_index('Index')
                    .assign(
                        Start=lambda x: x['Start'] - x['Start_b'] + 1,
                    )
        )

        res_df = res_df1.join(res_df2, lsuffix='_1', rsuffix='_2')
        df = df.assign(
            chrom1=res_df['Name_1'],
            pos1=res_df['Start_1'],
            chrom2=res_df['Name_2'],
            pos2=res_df['Start_2']
            )[columns]
        df = df.astype(meta)

        del res_df1, res_df2, query_df, res_df
        gc.collect()

        trans_lower_idx = (df['chrom1'].cat.codes.values > 
                                df['chrom2'].cat.codes.values)
        trans_lower_data = df.loc[trans_lower_idx]
        trans_lower_idx = trans_lower_data.index
        cis_lower_idx = df['pos1'] > df['pos2']
        cis_lower_data = df.loc[cis_lower_idx]
        cis_lower_idx = cis_lower_data.index 

        lower_idx = np.r_[trans_lower_idx, cis_lower_idx]
        lower_data = df.loc[lower_idx]
        
        df = df.astype(
            {'chrom1': str,
            'chrom2': str}
        )

        df.loc[lower_idx, 'chrom1'] = lower_data['chrom2']
        df.loc[lower_idx, 'pos1'] = lower_data['pos2']
        df.loc[lower_idx, 'chrom2'] = lower_data['chrom1']
        df.loc[lower_idx, 'pos2'] = lower_data['pos1']

        df = df.astype(meta)

        out = f'temp.part_{n}.corrected.pairs'
        df.to_csv(out, sep='\t', index=None, header=None)
        
        del df
        gc.collect()

        return  out

    def chrom2contig(self, source, output='corrected.pairs', threads=4):
        """
        Convert chromosome-level pairs to contig-level, 
            or convet draft assembly to corrected assembly.
        
        Params:
        --------
        source: str, pd.DataFrame, pr.PyRanges
            contig position information, four columns
                (chrom, start, end, name)

        output: str, default contig.pairs
            output of contig-level pairs.

        threads: int, default 4
            Number of threads.
        Returns:
        --------
        None

        Examples:
        --------
        >>> p.chrom2contig(corrected_df)
        """
        ## import corrected bed
        if isinstance(source, str):
            logger.info(f'Load source infomation of `{source}`.')
            source_df = pd.read_csv(source, sep='\t', 
                                        header=None, 
                                        index_col=None, 
                                        names=['Chromosome', 'Start', 
                                                'End', 'Name'])
            source_gr = pr.PyRanges(source_df)
        elif isinstance(source, pd.DataFrame):
            source_df = source
            source_gr = pr.PyRanges(source_df)
        elif isinstance(source, pr.PyRanges):
            source_df = source.df
            source_gr = source
        else:
            raise TypeError("source must be string/dataframe/pyranges.")
        
        ## create corrected header
        new_header = self.header
        chrom_dtype = pd.CategoricalDtype(source_df.Name)
        meta = {
            'chrom1': chrom_dtype,
            'chrom2': chrom_dtype
        }

        new_chroms = (
                        source_df['Name']
                        .astype(chrom_dtype)
                        .cat
                        .categories
                        .tolist()
        )
        new_length = (source_df['End'] - source_df['Start']).values.tolist()
        new_header.chromsize = OrderedDict(zip(new_chroms, new_length))
        new_header.update()
        new_header_filename = 'temp.corrected.new.header'
        new_header.save(new_header_filename)
        
        logger.info('Converting the position of pairs ...')
        columns = self.columns
        
        args = []
        for i in range(self.data.npartitions):  
            args.append((self.data, i, source_gr, columns, meta))
             
        res = Parallel(n_jobs=min(threads, len(args)))(
            delayed(self._chrom2contig)(i, j, k, l, m)
                for i, j, k, l, m in args)

        del source_gr
        gc.collect()

        pairs_files = ' '.join(res)
        command = f"""cat {new_header_filename} {pairs_files} > {output}"""
        check_call(command, shell=True)
        logger.info(f'Written corrected pairs file into `{output}`')

        ## remove temp file
        logger.debug('Removing temporary files ...')
        Parallel(n_jobs=min(threads, len(res)))(
            delayed(os.remove)(i) for i in res
        )
        os.remove(new_header_filename)

        return output


class PAFLine:
    def __init__(self, line):
        self.line = line.strip()
        data = line.strip().split()
        self.query = data[0]
        self.query_length = int(data[1])
        self.query_start = int(data[2])
        self.query_end = int(data[3])
        self.strand = data[4]
        self.target = data[5]
        self.target_length = int(data[6])
        self.target_start = int(data[7])
        self.target_end = int(data[8])
        self.num_match = int(data[9])
        self.aln_length = int(data[10])
        self.mapq = int(data[11])

    def __str__(self):
        return self.line


class PAFRecords:
    
    def __init__(self, paf_file):
        self.paf_file = paf_file
    
    def parse(self):
        with open(self.paf_file) as fp:
            for line in fp:
                yield PAFLine(line)

class PAFTable:
    """
    dataframe for paf file.

    Params:
    --------
    paf: str
        Path of paf file from minimap2.
    min_quality: int, default 1
        Minimum of mapping quality.
    min_length: int, default 50
        Minumum of fragment length.
    threads: int, default 4
        Number of threads.

    Returns:
    --------

    Examples:
    --------
    >>> 
    """
    
    PAF_HEADER = [
        'read_name', 'read_length', 
        'read_start', 'read_end', 
        'strand', 'chrom', 
        'chrom_length', 'start', 
        'end', 'matches', 
        'aln_length', 'mapping_quality',
        'pass_filter', 'filter_reason'
        ]

    META = {
        "read_name": "category",
        "strand": pd.CategoricalDtype(["+", "-"], ordered=False),
        "chrom": "category",
        "mapping_quality": "uint8",
        "filter_reason": "category",
        "identity": "float32"
    }
    def __init__(self, paf: str, 
                    threads: int = 4):
        self.file = Path(paf)
        self.filename = self.file.name
        self.threads = threads

        self.data = self.read_table()

        #self.client = Client(n_workers=self.threads) if self.threads > 1 else None
        self.tmpdir = tempfile.mkdtemp(prefix="tmp", dir='./')
    
    def read_table(self):
        logger.info('Starting ...')
        df = dd.read_csv(self.file, sep='\t', usecols=range(12), 
                        header=None, names=self.PAF_HEADER[:12], 
                        dtype=self.META, blocksize=1e8)

        df['read_idx'] = df['read_name'].cat.as_known().cat.codes
        # df = (df.set_index('read_idx')
        #         .repartition(df.npartitions)
        #         .reset_index())
        
        logger.info(f'Load paf file of `{self.filename}`.')

        return df
    
    @staticmethod
    def _filter_low_mq(df,
                    min_quality: int=1, 
                    min_identity: float=0.75,
                    min_length: int=50):

     
        
        df["pass_filter"] = df.pass_filter & df.pass_mq
        df["filter_reason"] = df.apply(
            lambda x: x.filter_reason if x.pass_mq else "low_mq", axis=1)
        df = df.drop('pass_mq', axis=1)

        return df  

    def producer(self, args):
        for i in range(len(self.data.npartitions)):
            df = self.data.partitions[i].compute()
            yield df

    @staticmethod
    def _filter(i, df, low_mq_condition, tmpdir):
        df = df.partitions[i].compute()

        ## init
        df = (df
            .assign(identity=lambda x: (x.matches/x.aln_length).round(2).astype(np.float32))
            .assign(fragment_length=lambda x: (x.read_end - x.read_start))
            .assign(pass_filter=True)
            .assign(filter_reason="pass")
        )

        ## filter low mapping quality and too short alignments
        min_quality, min_identity, min_length = low_mq_condition
        df = PAFTable._filter_low_mq(df, min_quality, min_identity, min_length)

        ## filter singleton
        singleton_idx = (df.query("pass_filter == True")
                            .groupby("read_idx")["read_idx"].count() == 1)
        singleton_idx = singleton_idx.index.to_numpy()[singleton_idx]
        df = df.reset_index().set_index("read_idx")
        candicate_idx = (df.loc[singleton_idx]
                            .query("pass_filter == True")["index"])
        
        df = df.reset_index().set_index("index")
        df.loc[candicate_idx, "pass_filter"] = False
        df.loc[candicate_idx, "filter_reason"] = "singleton"

        df.to_csv(f"{tmpdir}/tmp.part{i}.csv", sep='\t', index=False, header=True)

    def filter(self, min_quality: int=1, 
                    min_identity: float=0.75,
                    min_length: int=50):
        
        low_mq_condition = (min_quality, min_identity, min_length)
        filter_condition = (f'map_quality <= {min_quality}'
                                     f' | identity <= {min_identity}'
                                     f' | fragment_length <= {min_length}')
        logger.info(f"Filter low mapping quality by: {filter_condition}")
        logger.info("Filter singleton ...")
        args = []   
        for i in range(self.data.npartitions):
            args.append((i, self.data, low_mq_condition, self.tmpdir))
        
        Parallel(n_jobs=self.threads, verbose=2)(
                delayed(PAFTable._filter)(i, j, k, l)
                    for i, j, k, l in args)

        logger.info("Done.")

        self.data = dd.read_csv(f"{self.tmpdir}/tmp.*.csv", sep='\t', 
                                    header=0, dtype=self.META)

    def to_pairs(self, chromsize, output):
        def _to_pairs_per_concatemer(df):
            if len(df) <= 1:
                return None 
            
            rows = list(
                df.sort_values(['chrom', 'pos']).itertuples()
            )
            read_name = df.index.values[0]

            res = []
            for i, (align_1, align_2) in enumerate(combinations(rows, 2)):

                _res = [f'{read_name}_{i}', 
                        align_1.chrom, align_1.pos,
                        align_2.chrom, align_2.pos,
                        align_1.strand, align_2.strand]
                res.append(_res)
            
            res_df = pd.DataFrame(res, columns=['readID', 
                                'chrom1', 'pos1', 'chrom2', 
                                'pos2', 'strand1', 'strand2'])

            del res

            return res_df

        def _to_pairs(df):
            return df.groupby('read_idx').apply(_to_pairs_per_concatemer)
    
        logger.info('Converting Pore-C alignment to pairs.')
        self.data = self.data.assign(
            pos=lambda x: np.rint(x.start + (x.end - x.start)/2).astype(int))
        dtypes = self.data.head(1).dtypes
        meta = {
            "readID": str, 
            "chrom1": dtypes["chrom"],
            "pos1": dtypes["start"],
            "chrom2": dtypes["chrom"],
            "pos2": dtypes["start"],
            "strand1": dtypes["strand"],
            "strand2": dtypes["strand"]
        }

        (self.data.set_index('read_idx')[['chrom', 'pos', 'strand']]
                .map_partitions(_to_pairs, meta=meta)
                .to_csv(f'{output}.body', single_file=True, 
                sep='\t', header=None, index=False)
        )
        
        ph = PairHeader([])
        ph.from_chromsize_file(chromsize)
        ph.save(f'{output}.header')


        command = f"""cat {output}.header {output}.body > {output}"""
        check_call(command, shell=True)
        os.remove(f'{output}.body')
        os.remove(f'{output}.header')

        logger.info(f'Successful written pairs file into `{output}`.')

    def to_pore_c_table(self):
        data = self.data[["read_idx", "chrom", "start", 
                            "end", "strand", "read_name", 
                            "read_start", "read_end",
                            "mapping_quality", "identity",
                            "pass_filter", "filter_reason"]]

        pct = PoreCTable()
        pct.from_dask(data)

        return pct

    def clean_tempoary(self):
        shutil.rmtree(self.tmpdir)
class PoreCTable:
    HEADER = ["read_idx", "chrom", "start", "end", 
                "strand", "read_name", "read_start", 
                "read_end", "mapping_quality", "identity", 
                "pass_filter", "filter_reason"]
    
    def __init__(self, threads=4, lower_memory=False):
        self.threads = threads
        self.lower_memory = lower_memory

        if self.lower_memory:
            self.client = Client(n_workers=self.threads)

    def read_table(self, pore_c_table):
        if self.lower_memory:
            self.data = dd.read_parquet(pore_c_table, engine=PQ_ENGINE)
            self.data = self.data.set_index('read_idx').reset_index()
        else:
            self.data = pd.read_parquet(pore_c_table)

    def from_dask(self, data):
        self.data = data

    def to_pairs(self, chromsize, output):
        def _to_pairs_per_concatemer(df):
            rows = list(
                df.sort_values(['chrom', 'pos']).itertuples()
            )
            read_name = df.index.values[0]

            res = []
            for i, (align_1, align_2) in enumerate(combinations(rows, 2)):

                _res = [f'{read_name}_{i}', 
                        align_1.chrom, align_1.pos,
                        align_2.chrom, align_2.pos,
                        align_1.strand, align_2.strand]
                res.append(_res)
            
            res_df = pd.DataFrame(res, columns=['readID', 
                                'chrom1', 'pos1', 'chrom2', 
                                'pos2', 'strand1', 'strand2'])

            del res

            return res_df

        def _to_pairs(df):
            df = df.compute()
            return df.groupby('read_idx').apply(_to_pairs_per_concatemer)
    
        logger.info('Converting Pore-C alignment to pairs.')
        df = self.data.query('pass_filter == True').assign(
            pos=lambda x: np.rint(x.start + (x.end - x.start)/2).astype(int))
        dtypes = self.data.head(1).dtypes
        meta = {
            "readID": str, 
            "chrom1": dtypes["chrom"],
            "pos1": dtypes["start"],
            "chrom2": dtypes["chrom"],
            "pos2": dtypes["start"],
            "strand1": dtypes["strand"],
            "strand2": dtypes["strand"]
        }

        args = []
        for i in range(df.npartitions):
            args.append(df.partitions[i])
        
        res = Parallel(n_jobs=min(self.threads, len(args)))(
            delayed(_to_pairs)(i) for i in args
        )
        res_df = pd.concat(res, axis=0)
        res_df.to_csv(f"{output}.body", sep='\t', header=None, index=None)
        ph = PairHeader([])
        ph.from_chromsize_file(chromsize)
        ph.save(f'{output}.header')


        command = f"""cat {output}.header {output}.body > {output}"""
        check_call(command, shell=True)
        os.remove(f'{output}.body')
        os.remove(f'{output}.header')

        logger.info(f'Successful written pairs file into `{output}`.')

    @staticmethod
    def _read_and_alignment_stat(df, lower_memory=True):
        if lower_memory:
            df = df.compute()

        read_group = df.groupby('read_idx')
        
        alignment_per_reads = read_group['read_idx'].count()
        mapping_reads = len(alignment_per_reads)

        singleton_reads = len(alignment_per_reads[alignment_per_reads == 1])

        alignment_per_reads2 = alignment_per_reads[read_group['pass_filter'].any() == False]
        low_mq_reads = len(alignment_per_reads2[alignment_per_reads2 != 1])
        
        pass_reads = mapping_reads - low_mq_reads - singleton_reads

        read_statistics_df = pd.DataFrame([mapping_reads, pass_reads, 
                                        singleton_reads, low_mq_reads], 
                                        index=['mapping', 'pass', 
                                                'singleton', 'low_mq'])

        alignment_counts = len(df)
        filter_reason = df['filter_reason']
        index_name = dict(zip(range(len(filter_reason.cat.categories)), 
                                filter_reason.cat.categories))

        align_statistics_df = (filter_reason
                            .cat.codes
                            .value_counts()
                            .rename(index=index_name)
                            .rename_axis("filter_reason")
                            .to_frame()
                            .rename(columns={0: "count"}))

        #statistics_df.index = statistics_df.index.add_categories(['total alignments'])
        align_statistics_df.loc['total alignments'] = alignment_counts

        return read_statistics_df, align_statistics_df

    def read_and_alignment_stat_dask(self, read_counts=None):
        if self.lower_memory:
            args = []
            tmp_df = self.data[['read_idx', 'pass_filter', 'filter_reason']]
            for i in range(self.data.npartitions):
                args.append(tmp_df.partitions[i])
            
            res = Parallel(n_jobs=self.threads)(
                    delayed(PoreCTable._read_and_alignment_stat_dask)(i) for i in args)

            read_statistics_df = pd.concat(list(zip(*res))[0], axis=1).sum(axis=1).to_frame()
            align_statistics_df = pd.concat(list(zip(*res))[1], axis=1).sum(axis=1).to_frame()
        
        else:
            read_statistics_df, align_statistics_df = \
                PoreCTable._read_and_alignment_stat(self.data, lower_memory=False)
    
        read_statistics_df.columns = ['count']
        read_statistics_df = read_statistics_df.rename_axis("item")

        if read_counts:
            read_statistics_df.loc['total reads'] = read_counts
            read_statistics_df.loc['unmapping'] = read_counts - read_statistics_df.loc['mapping']
            read_statistics_df['percent'] = read_statistics_df['count'] / read_counts 

            read_statistics_df = read_statistics_df.loc[['total reads', 
                                                        'pass', 'singleton', 
                                                        'low_mq', 'unmapping']]

        return read_statistics_df, align_statistics_df

    @staticmethod
    def _alignment_stat(df):
        df = df.compute()
        alignment_counts = len(df)
        filter_reason = df['filter_reason']
        index_name = dict(zip(range(len(filter_reason.cat.categories)), 
                                filter_reason.cat.categories))

        align_statistics_df = (filter_reason
                            .cat.codes
                            .value_counts()
                            .rename(index_name))

        #statistics_df.index = statistics_df.index.add_categories(['total alignments'])
        align_statistics_df.loc['total alignments'] = alignment_counts
    
        return align_statistics_df

    def alignment_stat(self):
        args = []
        tmp_df = self.data[['read_idx', 'filter_reason']]
        for i in range(self.data.npartitions):
            args.append(tmp_df.partitions[i])
        
        res = Parallel(n_jobs=self.threads)(
                delayed(PoreCTable._alignment_stat)(i) for i in args)

        statistics_df = pd.concat(res, axis=1).sum(axis=1).to_frame()
        statistics_df.columns = ['count']

        return statistics_df

    def contact_stat(self):
        if self.lower_memory:
            res = (self.data.query("pass_filter == True")
                            .groupby('read_idx')['read_idx']
                            .count()
                            .compute()
                            )
        else:
            res = (self.data.query("pass_filter == True")
                            .groupby('read_idx')['read_idx']
                            .count()
                            )
        length_bins = pd.IntervalIndex.from_breaks([1, 2, 3, 4, 6, 11, 
                                                    21, 50, int(1e9)])
        length_bin_labels = {}
        for i in length_bins:
            if i.length == 1:
                label = str(i.right)
            elif i.length >= int(1e8):
                label = f"gt_{i.left}"
            else:
                label = f"{i.left + 1}-{i.right}"
            length_bin_labels[i] = label

        read_order_hist = (pd.cut(res, length_bins, labels=length_bin_labels)
                                .value_counts()
                                .rename("count")
                                .sort_index()
                                )
        
        read_order_hist = (read_order_hist
                            .to_frame()
                            .rename(index=length_bin_labels)
                            .rename_axis("concatemer_order")
                            )
        
        read_order_hist["perc"] = (read_order_hist
                                    .div(read_order_hist["count"].sum(), axis=1) * 100)

        return read_order_hist

    def save(self, output):
        self.data.to_parquet(output, write_metadata_file=True, 
                                write_index=False, engine=PQ_ENGINE, 
                                version=PQ_VERSION)
        logger.info(f"Successfully output pore_c_table to `{output}`")