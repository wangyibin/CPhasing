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

import math
import numpy as np
import pandas as pd 
import polars as pl
import pyranges as pr 
import tempfile

from collections import Counter, OrderedDict, defaultdict
from joblib import Parallel, delayed
from itertools import combinations
from multiprocessing import Process, Pool
from pathlib import Path
from pandarallel import pandarallel
from subprocess import check_call, Popen, PIPE

from ._config import *
from .expection import *
from .utilities import (
    binnify,
    decompress_cmd,
    list_flatten, 
    tail, 
    xopen, 
    gz_reader,
    read_chrom_sizes,
    parse_corrected_contig
)

# from line_profiler import profile
# from memory_profiler import profile as mprofile

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
    def __init__(self, countRE, minRE=3, has_header=True):
        self.filename = countRE
        self._path = Path(self.filename)
        # self.suffix = re.match(r'.*(counts_.*.txt)$', 
        #                         self.filename).groups()[0]
        self.minRE = minRE
        self.has_header = has_header

        self.parse()

    def parse(self):
        if self.has_header:
            self.data = pd.read_csv(self.filename, sep='\t', 
                                        header=0, index_col=0,
                                        dtype={'RECounts': int, 
                                                'Length': int})
        else:
            self.data = pd.read_csv(self.filename, sep='\t', 
                                        header=None, index_col=0,
                                        names=['#Contigs','RECounts', 'Length'],
                                        dtype={'RECounts': int, 
                                                'Length': int},
                                        comment="#")
            
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

    @property
    def length_db(self):
        """
        Returns:
        -------- 
        dict:
            dict of length 
    
        """
        return self.data['Length'].to_dict()

    @property 
    def midpoint_db(self):
        """
        Calculate the midpoint and return dictory 

        Returns:
        --------
        dict: 
            dict of midpoint
        """

        return self.data['Length'].map(lambda x: x//2).to_dict()


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

    def __init__(self, line, format='allele1'):
        self.line = line.strip()
        line_list = line.strip().split()
        self.chrom = line_list[0]
        self.gene = line_list[1]
        if format == "allele1":
            self.contigs = line_list[2:]
        elif format == "allele2":
            self.contigs = line_list[2:5]
            self.score = int(line_list[7])
            try:
                self.strand = int(line_list[8])
            except IndexError:
                self.strand = 0

        self.n = len(self.contigs)
    
    def __str__(self):
        return self.line 

class AlleleInfo:
    """
    object function of Allele table additional information
    """
    def __init__(self, line):
        self.line_list = line.strip("\n").strip("#").split()
        self.contig = self.line_list[0]
        self.length = self.line_list[1]
        self.mz = self.line_list[2]
        self.mzUnique = self.line_list[3]
    

class AlleleTable:
    """
    Allele table

    Params:
    --------
    infile: str
        input file for allele table
    sort: bool
        whether to sort contigs per row [True]
    fmt: str
        specified the format of allele table
            allele1: old format for table only contains alleles
                Chr01 12345 contig1 contig2 contig3 ...
                or 1 1 contig1 contig2
            allele2: new format for table contains the scores 
                1 1 contig1 contig2 mzConsidered1 mzConsidered2 mzShared kmerSimilarity

    Examples:
    --------
    >>> infile = 'Allele.ctg.table'
    >>> at = AlleleTable(infile)

    """
    AlleleHeader2 = ['idx1', 'idx2', 'contig1', 'contig2', 
                     'mz1', 'mz2', 'mzShared', 'similarity', 
                     'strand']
    def __init__(self, infile, sort=True, fmt="allele1"):
        self.filename = infile
        self.sort = sort
        self.fmt = fmt 
        assert self.fmt in ("allele1", "allele2"), \
                "format must in ['allele1', 'allele2']"
        if not Path(self.filename).exists():
            logger.error(f'No such file of `{self.filename}`.')
            sys.exit()

        logger.info(f'Load allele table: `{self.filename}`.')
        self.info = OrderedDict()
        self.check(self.fmt)
        self.get_data()

    def check(self, fmt="allele1"):
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
                if line[0] == "#":
                    if fmt == "allele1":
                        continue
                    else:
                        _info = AlleleInfo(line)
                        self.info[_info.contig] = _info 

                al = AlleleLine(line)
                if al.n > n:
                    n = al.n
        
        self.columns = list(range(n + 2))

    def get_data(self):
        def sort_row(row):
            row_filter = filter(lambda x: not pd.isna(x), row)
            return pd.Series(sorted(row_filter), name=row.name)
        
        df = pd.read_csv(self.filename, sep='\t', header=None, 
                            index_col=0, names=self.columns,
                            usecols=self.columns, comment="#")
        df.index = df.index.astype('category')
        df = df.dropna(how='all', axis=1)
        if df.empty:
            self.data = df 
            return 
    
        if self.fmt == "allele1":
            df = df.drop(1, axis=1)
            df = df.apply(sort_row, axis=1) if self.sort else df
            df.columns = list(range(1, len(df.columns) + 1))
            df = df.drop_duplicates(list(range(1, len(df.columns) + 1)))
        elif self.fmt == "allele2":
            df.drop(1, axis=1, inplace=True)
            df.columns = [1, 2] + self.AlleleHeader2[4:]

            df.drop_duplicates([1, 2], inplace=True)
            
        df = df.reset_index(drop=True)

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
        assert self.fmt == "allele1", "only support for format `allele1`"
        self.data = self.data[self.data.isin(contigs)]
        self.data.index.rename(0, inplace=True)
        
        if self.data.empty:
            logger.warning("After filter AlleleTable is empty")

    @property
    def contig_pairs(self):
        if self.fmt == 'allele1':
            _contig_pairs = set()
            groups = self.groups
            for group, items in groups.items():
                for pair in list(combinations(items, 2)):
                    _contig_pairs.add(pair)
        else:
            _contig_pairs = set(map(tuple, self.data[[1, 2]].values.tolist()))

        return _contig_pairs

    @property
    def groups(self):
        assert self.fmt == "allele1", "only support for format `allele1`"
        _groups = OrderedDict()
        data = self.data.reset_index()
        _groupby = data.groupby('index')
        for chrom, item in _groupby:
            tmp = item.values.flatten()[1:]
            _groups[chrom] = tmp[~pd.isna(tmp)]

        return _groups

    @property
    def n(self):
        """
        the number of allelic contigs per row
        """
        if self.fmt == 'allele1':
            return len(self.columns) - 2
        else:
            return len(self.columns) - 2 - 4

    @property
    def chromnames(self):
        """
        list of chromosomes
        """
        assert self.fmt == "allele1", "only support for format `allele1`"
        return self.groups.keys()

    @property
    def nchroms(self):
        """
        number of chromosomes
        """
        assert self.fmt == "allele1", "only support for format `allele1`"
        return len(self.chromnames)

    @property
    def contigs(self):
        """
        list of contigs
        """
        if self.fmt == "allele1":
            _contigs = self.data[range(1, self.n + 1)].values.flatten()
            _contigs = _contigs[~pd.isna(_contigs)]
        else:
            _contigs = self.data[[1, 2]].values.flatten()
            _contigs = _contigs[~pd.isna(_contigs)]

        return list(set(_contigs))
    
    @property
    def ncontigs(self):
        """
        number of contigs
        """
        return len(self.contigs)
    
    @property
    def similarity(self):
        assert self.fmt == "allele2", \
            "only support for allele2 format"
        
        df = self.data.set_index(list(range(1, self.n + 1)))
        
        return df["similarity"].to_dict()

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
    
    def allele1_to_allele2(self):
        assert self.fmt == 'allele1', "only support for format `allele1`"

        contig_pairs = self.contig_pairs 
    
        res = []
        i = 0
        for pair in contig_pairs: 
            i += 1
            res.append([i, i, pair[0], pair[1], 0, 0, 0, 1.0, 0])
            i += 2
            res.append([i, i, pair[1], pair[0], 0, 0, 0, 1.0, 0])
            
        res_df = pd.DataFrame(res)
        self.fmt = 'allele2'
        self.data = res_df 

    def allele2_to_allele1(self):
        assert self.fmt == 'allele2', 'only support for format `allele1`' 

        res = defaultdict(set)
        contig_pairs = self.contig_pairs
        
        for contig1, contig2 in contig_pairs:
            res[contig1].add(contig1)
            res[contig1].add(contig2)

        data = []
        for i, item in enumerate(set(map(tuple, map(sorted, res.values())))):
            data.append([i, i, *item])
        
        res_df = pd.DataFrame(data)
        self.data = res_df
   
        self.fmt = 'allele1'

    def get_contigsizes_from_header(self):
        assert self.fmt == 'allele2', "Only support for fmt=allele2"
        res = dict()
        with open(self.filename, 'r') as fp:
            for line in fp:
                if line.startswith("#"):
                    contig, length, minimizer, unique_minimizer = line[1:].strip().split()
                    res[contig] = int(length)

        return res

    def get_high_similar_pairs(self, minimum_similarity=0.995):
        assert self.fmt == 'allele2', "Only support for fmt=allele2"
        contigsizes = self.get_contigsizes_from_header()

        res = set()
        for idx, row in self.data[self.data['similarity'] > minimum_similarity].iterrows():
            contig1, contig2 = row[1], row[2]
            minimizer_count1, minimizer_count2 = row['mz1'], row['mz2']
            length1, length2 = contigsizes[contig1], contigsizes[contig2]
            
            if math.isclose(length1, length2, rel_tol=0.2, abs_tol=0) \
                or math.isclose(minimizer_count1, minimizer_count2, rel_tol=0.2, abs_tol=0):
                if contig1 > contig2:
                    contig1, contig2 = contig2, contig1

                res.add((contig1, contig2))

        res = list(res)

        return sorted(res)

    def save(self, output, index=True):
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
        if len(self.data.columns) == 7:
            data = self.data.reset_index().reset_index() 
        else:
            data = self.data
        
        data.to_csv(output, sep='\t', header=None, index=index)
    
class PruneTable:
    Header = ['contig1', 'contig2', 
              'mz1', 'mz2',
              'mzShared', 'similarity',
              'type']
    META = {
        'contig1': 'object',
        'contig2': 'object',
        'mz1': 'int64',
        'mz2': 'int64',
        'mzShared': 'int64',
        'similarity': 'float32',
        'type': 'category'
    }
    def __init__(self, prunetable, symmetric=True):
        self.file = Path(prunetable)
        self.filename = self.file.name
        self.data = self.read_table()

    def read_table(self):
        try:
            data = pd.read_csv(self.file, sep='\t', 
                            header=None, index_col=None,
                            names=self.Header, 
                            dtype=self.META)
        except pd.errors.ParserError:
            logger.error("Incorrect format of prune table.")
            sys.exit(-1)
        except ValueError:
            logger.error("Incorrect format of prune table.")
            sys.exit(-1)

        return data
    
    def symmetric_table(self):
        
        tmp_df2 = self.data.rename(columns={"contig1": "contig2",
                                            "contig2": "contig1"})
        
        data = pd.concat([self.data, tmp_df2], axis=0)
        self.data = data.reset_index(drop=True)
        
        return data
        

        

class PairTable:
    """
    Object of allhic pairs file

    Params:
    --------
    infile: str
        input file of pairs table
    symmetric: bool
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

    def __init__(self, infile, symmetric=True, index_contig=True, remove_self=True):
        self.filename = infile
        self.remove_self = remove_self
        
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
        if self.remove_self:
            df = df[df['Contig1'] != df['Contig2']]
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

    @property
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
        assert self.is_symmetric, "The pair table must be symmetric."

        tmp_df = self.data.reindex((i, j) for i in contig1 for j in contig2)
        tmp_df = tmp_df.dropna(how='all', axis=0)
        
        return tmp_df
    
    def get_contact(self, contig1, contig2, count_re):
        """
        get unnormalized contact between two contig lists

        Params:
        --------
        contig1: list or list-like
            list of contigs
        contig2: list or list-like
            list of contigs
        count_re: CountRE
            count RE table

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
        C1 = count_re.data.reindex(contig1)['RECounts'].sum()
        C2 = count_re.data.reindex(contig2)['RECounts'].sum()


        return L * 1000 / (C1 * C2)
    

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

    @property
    def hap_groups(self):
        db = OrderedDict()
        for group in self.data:
            try:
                hap, idx = group.rsplit("g", 1)
            except:
                if group not in db:
                    db[group] = []
                db[group].append(self.data[group])
            finally:
                if hap not in db:
                    db[hap] = []

                db[hap].append(self.data[group])
            
        return db

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
        
        groups = list(filter(lambda x: x in self.data, self.groups))
        for group in groups:
            _contigs = self.data[group]
            tmp_df = self.CountRE.data.loc[_contigs]
            tmp_df.to_csv(f'{group}.txt', sep='\t',
                            header=True, index=True)
            logger.info(f'Output countRE file into `{group}.txt`')

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
    
    def to_agp(self, contigsizes, output, orientation_db=None):
        """
        Convert cluster to a pseudo agp 
        """
        contigsizes = read_chrom_sizes(contigsizes)
        contigsizes_db = contigsizes.to_dict()['length']
        agp_res = []
        idx = 0
        added_contigs = set()
        for group in self.groups:
            _contigs = self.data[group]
            chrom_start = 1 
            chrom_end = 0
            for i, contig in enumerate(_contigs):
                contig_length = contigsizes_db[contig]
                chrom_start = chrom_end if i != 0 else 1
                chrom_end = chrom_start + contig_length 
                if orientation_db:
                    orientation = orientation_db[contig]
                else:
                    orientation = "+"
                tmp_res = [group, chrom_start, chrom_end, idx, 
                           "W", contig, 1, contig_length, orientation]
                agp_res.append(tmp_res)

                chrom_start = chrom_end + 1 
                chrom_end = chrom_start + 100 - 1
                idx += 1
                if i < len(_contigs):
                    tmp_res = [group, chrom_start, chrom_end, idx, 
                            "U", 100, "contig", "yes", "map"]
                    agp_res.append(tmp_res)
                    idx += 1 

                added_contigs.add(contig)

        for contig in contigsizes_db:
            if contig not in added_contigs:
                contig_length = contigsizes_db[contig]
                tmp_res = [contig, 1, contig_length, idx, 
                           "W", contig, 1, contig_length, "+"]
                
                agp_res.append(tmp_res)
                
                idx += 1

        agp_df = pd.DataFrame(agp_res)
   
        agp_df.to_csv(output, sep='\t', index=None, header=None)



    def to_fasta(self, fasta, trim_length=25000, outdir="./"):
        """
        extract fasta by cluster table 
        """
        from Bio import SeqIO 
        
        cluster_file = Path(self.filename).absolute()
        os.chdir(outdir)
        cmd = [
                "cphasing-rs",
                "extract-fasta",
                str(fasta),
                str(cluster_file),
                "-l",
                str(int(trim_length)),
                " 2>/dev/null"
                ]
      
        flag = os.system(" ".join(cmd))    
        assert flag == 0, "Error in extracting fasta"
        os.chdir("../")

        file_db = {}
        files = []
        for group in self.groups:
            # shutil.move(f"{group}.contigs.fasta", f"{outdir}/{group}.contigs.fasta")
            if len(self.data[group]) > 1:
                file_db[group] = f"{outdir}/{group}.contigs.fasta"
                files.append(f"{outdir}/{group}.contigs.fasta")
            else:
                if Path(f"{outdir}/{group}.contigs.fasta").exists():
                    os.remove(f"{outdir}/{group}.contigs.fasta")
        # db = self.contig_groups
        # fasta = SeqIO.parse(xopen(fasta), 'fasta')

        # threshold = trim_length * 3
        # for record in fasta:
        #     if record.id in db:
        #         seq = record.seq
        #         length = len(record.seq)
        #         if length > threshold:
        #             seq = seq[trim_length: length - trim_length + 1]
        #         if length > 1000_000:
        #             seq1 = seq[:500_000]
        #             seq2 = seq[500_000:]
        #             seq = seq1 + seq2
                
        #         file_db[db[record.id]].write(f'>{record.id}\n{seq}\n')
        
        # for group in file_db:
        #     file_db[group].close()

        return files

    def to_tour(self):
        """
        convert cluster table to several pseudo tour files.
        
        """
        for group in self.groups:
            with open(f"{group}.tour", "w") as out:
                _contigs = self.data[group]
                _contigs = list(map(lambda x: f"{x}+", _contigs))
                print(" ".join(_contigs), file=out)
            
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


class ClmLine():
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
 
class Clm(object):
    def __init__(self, clmfile, mem_cache='.'):
        
        self.file = clmfile
        self.parse()
        # self.mem_cache = mem_cache
        # self.memory = Memory(mem_cache, verbose=0)
        # self.dk_df = self.memory.cache(self._dk_df)

    def count_db(self, count_re=None):
        """
        generate a dictionary containing the contacts of contig pairs

        Params:
        --------
        count_re: None (default) or CountRE

        Returns:
        db: dict
            database of contacts between two contigs

        Examples:
        --------
        >>> clm.count_db()
        {"ctg1": 20, "ctg2": 30}
        >>> clm.count_db(cr) ## normalized by the distance between contigs
        {"ctg1": 0.5, "ctg2": 0.8}
        """
        db = {}
        with open(self.file, 'r') as fp:
            for line in fp:
                cl = ClmLine(line)
                if (cl.ctg1, cl.ctg2) in db:
                    continue
                if count_re:
                    try:
                        l1 = count_re.data.loc[cl.ctg1]['Length']
                        l2 = count_re.data.loc[cl.ctg2]['Length']
                    except KeyError:
                        continue
                    value = cl.count / ((l1 + l2) / 2)
                else:
                    value = cl.count
                db[(cl.ctg1, cl.ctg2)] = value
            
        return db

    def parse(self):
        self.data = {}
        with open(self.file, 'r') as fp:
            for line in fp:
                cl = ClmLine(line)
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
    def contigs(self):
        return list_flatten(list(self.data.keys()))
                            
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
        # logger.info(f"Load tour file `{self.filename}`")
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
        _new_data = [(contig, '-' if orient == '+' else '+') for contig, orient in _data]

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

    def to_tuples(self, strand_value=0) -> tuple:
        """
        Contert tour data to tuples of contig and orient.

        Params:
        --------
        strand_value: int
            0: represent output orientation of "+" and "-"
            1: represent output orientation of "1" and "-1"

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
         >>> tour.to_tuples(strand_value=1)
        ('utg001', 1), ('utg002', -1)
        """
        assert strand_value in {0, 1}, "strand_value must be 0 or 1"

        res = list(map(lambda x: x.data, self.data))
        if strand_value == 0:
            return res
        else:
            res = list(map(lambda x:(x[0], 1) if x[1] == "+" else (x[0], -1), res))
            return res

    def to_dict(self, strand_value=0) -> OrderedDict:
        """
        Convert tour data to dict of contig and orient 

        Params:
        --------
        strand_value: int
            0: represent output orientation of "+" and "-"
            1: represent output orientation of "1" and "-1"
        
        Returns:
        --------
        dict: 
            A dict of contig and orient
        
        Examples:
        --------
        >>> tour.data
        >>> tour.to_dict()
        OrderedDict([('utg001', '+'), ('utg002', '-')])
        >>> tour.to_dict(strand_value=1)
        OrderedDict([('utg001', 1), ('utg002', -1)])
        
        """

        data = self.to_tuples(strand_value)

        return OrderedDict(data)

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
    
    def to_agp(self, fasta, gap_length=100, corrected=False):
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
        corrected_contigs = []
        for contig, orient in self:
            corrected_contig = parse_corrected_contig(contig)
           
            if corrected and corrected_contig:
                contig, tig_start, tig_end = corrected_contig 
                corrected_contigs.append(corrected_contig)
                
            else:
                try:
                    length = fasta.faidx.index[contig].rlen
                    # length = len(fasta[contig])
                except KeyError:
                    logger.warning(f"Could not find `{contig}` in `{fasta.filename}`, skipped."
                                   "If corrected you should add `--corrected` parameter for `build`.")
                    
                    continue
                except AttributeError:
                    try:
                        length = len(fasta[contig])
                    except KeyError:
                        continue 
                
                tig_start = 1 
                tig_end = length 

            start = old_end + 1
            end = (tig_end - tig_start + 1)
            
            res.append([self.group, start, end, i, 'W', contig, 
                    tig_start, tig_end, orient])
            

            i += 1

        j = 1
        res2 = []
        start = 1
        old_end = 0

        for i, tmp_res in enumerate(res, 1):
            tmp_res[3] = j
            start = old_end + tmp_res[1]
            end = start + tmp_res[2] - 1
            tmp_res[1] = start 
            tmp_res[2] = end
            res2.append(tmp_res)
            j += 1
            if i != len(res):
                res2.append([self.group, end + 1, end + gap_length,
                        j, 'U', gap_length, 'contig', 'yes', 'map'])
                old_end = end + gap_length
                j += 1
            
        return corrected_contigs, pd.DataFrame(res2)

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
    
    def rm(self):
        os.remove(str(self.filename))

    def backup(self, suffix="bak"):
        try:
            shutil.move(f"{self.filename}", f"{self.filename}.{suffix}")
        except FileNotFoundError:
            logger.warning(f"No such file of `{self.filename}`, will not backup it.")

    def save(self, output):
        with open(output, 'w') as out:
            print(" ".join(list(map(lambda x: str(x), self.data))), file=out)
    
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
        
        _res = [record.split(":", 1)[1].split() for record in _res]
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
    HEADER = ["readID", "chrom1", "pos1",
                "chrom2", "pos2", "strand1",
                "strand2"]
    
    def __init__(self, pairs, threads=1):
        self.file = Path(pairs)
        self.filename = self.file.name
        self.header = self._get_header()
        self.data = self.read_table()
        self.threads = threads

    def _get_header(self):
        ph = PairHeader([])
        ph.from_file(self.file)

        return ph

    @property
    def chromsize(self):
        return self.header.chromsize
    
    @property
    def chromsizes(self):
        return self.header.chromsize

    @property
    def chromsizes_df(self):
        return pd.DataFrame(self.chromsize.items(), 
                            columns=["chrom", "length"])  
    
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
            "pos1": CHROM_COORD_DTYPE,
            "chrom2": chrom_dtype,
            "pos2": CHROM_COORD_DTYPE,
            "strand1": strand_dtype,
            "strand2": strand_dtype
        }
        return meta 

    def read_table(self):
        import dask.dataframe as dd
        
        logger.info(f'Load pairs file of `{self.filename}`.')
        if self.file.suffix == ".gz":

            df = dd.read_csv(self.file, sep='\t', comment='#', 
                                names=self.header.columns, 
                                header=None,  dtype=self.meta,
                                compression="gzip")
        else:
            df = dd.read_csv(self.file, sep='\t', comment='#', 
                            names=self.header.columns, 
                            header=None,  dtype=self.meta)

        return df
    
    def is_pairs_with_mapq(self):
        """
        Check if the pairs file contains mapq.
        """
        with xopen(self.file) as fp:
            for line in fp:
                if line.startswith("#"):
                    continue
                line_list = line.strip().split("\t")
                if len(line_list) <= 7:
                    return False 
                else:
                    if line_list[7].isdigit():
                        return True
                    else:
                        return False
                

    @staticmethod
    def _pos_to_range(df, idx=0):
        """convert single pos to range."""
        if idx == 0:
            df = (df.reset_index(drop=True)
                    .assign(
                        start1=lambda x: x['pos1'] - 1,
                        end1=lambda x: x['pos1'],
                        start2=lambda x: x['pos2'] - 1,
                        end2=lambda x: x['pos2']
                    )
            )
        elif idx == 1:
            df = (df.reset_index(drop=True)
                    .assign(
                        start1=lambda x: x['pos1'] - 1,
                        end1=lambda x: x['pos1'],
                    )
            )
        elif idx == 2:
            df = (df.reset_index(drop=True)
                    .assign(
                        start2=lambda x: x['pos2'] - 1,
                        end2=lambda x: x['pos2']
                    )
            )
        else: 
            raise ValueError("idx must in [0, 1, 2]")

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
            source_df.columns = ['Chromosome', 'Start', 
                                                'End', 'Name']
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
             
        res = Parallel(n_jobs=min(threads, len(args) + 1))(
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
    
    def intersection(self, bed, output, invert=False):
        """
        select pairs according to several regions
        """
        if output.rsplit(".", 1)[-1] == ".gz":
            output = output.rsplit(".", 1)[0]
            gzip_flag = True 
        else:
            gzip_flag = False

        if isinstance(bed, str):
            logger.info(f'Load source infomation of `{bed}`.')
            bed_df = pd.read_csv(bed, sep='\s+', 
                                        header=None, 
                                        index_col=None,
                                        usecols=[0, 1, 2],
                                        names=['Chromosome', 'Start', 
                                                'End'])
            bed_gr = pr.PyRanges(bed_df)
        elif isinstance(bed, pd.DataFrame):
            bed_df = bed
            bed_df.columns = ['Chromosome', 'Start', 
                                                'End']
            bed_gr = pr.PyRanges(bed_df)
        elif isinstance(bed, pr.PyRanges):
            bed_df = bed.df
            bed_gr = bed
        else:
            raise TypeError("source must be string/dataframe/pyranges.")

        df = self.data[['chrom1', 'pos1']].compute()
        df = df.reset_index()
        df = Pairs._pos_to_range(df, idx=1)
        query_df = df[['chrom1', 'start1', 'end1', 'index']]

        del df
        gc.collect()
    
        query_df.columns = ['Chromosome', 'Start', 'End', 'Index']
        query_gr = pr.PyRanges(query_df)
        if not invert:
            res_df1 = query_gr.join(bed_gr, nb_cpu=self.threads).df
        else: 
            res_df1 = query_gr.intersect(bed_gr, invert=False).df

        res_index1 = res_df1['Index']

        del res_df1
        gc.collect()

        df = self.data[['chrom2', 'pos2']].compute()
        df = df.reset_index()
        df = Pairs._pos_to_range(df, idx=2)
        query_df = df[['chrom2', 'start2', 'end2', 'index']]
        
        del df 
        gc.collect()

        query_df.columns = ['Chromosome', 'Start', 'End', 'Index']
        query_gr = pr.PyRanges(query_df)
        
        if not invert:
            res_df2 = query_gr.join(bed_gr, nb_cpu=self.threads).df
        else:
            res_df2 = query_gr.intersect(bed_gr, invert=True).df
        res_index2 = res_df2['Index']
        
        del query_df, query_gr, res_df2 
        gc.collect()

        res_index = np.intersect1d(res_index1, res_index2)

        df = self.data.compute().loc[res_index]

        self.header.save(f"temp.{output}.header")
        df.to_csv(f"temp.{output}.body", sep='\t', header=None, index=None)

        command = f"cat temp.{output}.header temp.{output}.body > {output}"
        check_call(command, shell=True)
        
        os.remove(f"temp.{output}.header")
        os.remove(f"temp.{output}.body")

        if gzip_flag:
            command = f"pigz -p 4 {output}"
            check_call(command, shell=True)
            logger.info(f"Successful output result in `{output}.gz`")
        else:
            logger.info(f"Successful output result in `{output}`")

    def to_mnd(self, output, threads=4):
        from dask.distributed import Client

        logger.info("Converting pairs to mnd file ...")
        with Client(n_workers=threads):
            (self.data[['strand1', 'chrom1', 'pos1', 
                        'strand2', 'chrom2', 'pos2']]
                .assign(
                    strand1=lambda x: x.strand1.replace('+', 0).replace('-', -1),
                    strand2=lambda x: x.strand2.replace('+', 0).replace('-', -1),
                    frag1=0,
                    frag2=1,
                    mapq1=1,
                    cigar1="-",
                    sequence1="-",
                    mapq2=1,
                    cigar2="-",
                    sequence2="-",
                    readname1="-",
                    readname2="-")
                .rename(
                    columns={'strand1': 'str1',
                                'chrom1': 'chr1',
                                'strand2': 'str2',
                                'chrom2': 'chr2'}
                        )
                    [['str1', 'chr1', 'pos1', 'frag1', 
                    'str2', 'chr2', 'pos2', 'frag2', 
                    'mapq1', 'cigar1', 'sequence1',
                    'mapq2', 'cigar2', 'sequence2',
                    'readname1', 'readname2']]
            ).to_csv(output, sep=' ', header=False, index=False, single_file=True)

        logger.info(f"Successful output mnd file intp `{output}`.")

class Pairs2:
    def __init__(self, pairs, min_mapq=0,
                 chunksize=1e8, threads=10):
        
        self.file = Path(pairs)
        self.filename = self.file.name
        
        self.chunksize = int(chunksize)
        self.threads = threads
        os.environ['POLARS_THREADS'] = str(self.threads)
        self.min_mapq = min_mapq
        
        if self.min_mapq > 0:
            logger.info(f"Only load pairs with minimum quality >= {self.min_mapq}.")

        if Path(self.file).is_dir():
            raise IsNotPairs("Input file must be a file, not a directory.") 
        else:
            try:
                self.header = self._get_header()
            except:
                raise IsNotPairs("Cannot read header of pairs file.")
            
            if not self.is_pairs():
                raise IsNotPairs("Input file is not a pairs file.")
    
    def _get_header(self):
        ph = PairHeader([])
        ph.from_file(self.file)

        return ph
    
    @property
    def chromsize(self):

        return pd.DataFrame(self.header.chromsize.items(), 
                            columns=['chrom', 'length']).astype({'length': int})
    

    @property
    def chromsizes(self):

        return pd.DataFrame(self.header.chromsize.items(), 
                            columns=['chrom', 'length']).astype({'length': int})

    @property
    def chromsizes_df(self):
        return self.chromsizes

    @property
    def chromnames(self):

        return self.chromsize['chrom']
    
    @property
    def columns(self):

        return self.header.columns
    
    @property
    def meta(self):
        max_value = self.chromsize['length'].max()
        if max_value < 2**8:
            pos_dtype = pl.UInt8
        elif max_value < 2**16:
            pos_dtype = pl.UInt16
        elif max_value < 2**32:
            pos_dtype = pl.UInt32
        else:
            pos_dtype = pl.UInt64 
        
        return {
                'chrom1': pl.Categorical, 
                'pos1': pos_dtype,
                'chrom2': pl.Categorical, 
                'pos2': pos_dtype,
                'mapq': pl.Int8
                }
    
    def is_pairs(self):
        if Path(self.file).is_dir():
            return False 
        else:
            if self.header.pairs_format:
                return True
            else:
                return False
            
    
    def is_pairs_with_mapq(self):
        """
        Check if the pairs file contains mapq.
        """
        with xopen(self.file) as fp:
            for line in fp:
                if line.startswith("#"):
                    continue
                line_list = line.strip().split("\t")
                if len(line_list) <= 7:
                    return False 
                else:
                    if line_list[7].isdigit():
                        return True
                    else:
                        return False

    def is_contain_mapq(self):
        df = pd.read_csv(self.file, sep='\t', header=None, index_col=None, comment="#", nrows=1)
        if df.shape[1] >= 8 and (df[7].dtype == np.int8 or df[7].dtype == np.int16 or df[7].dtype == np.int32 or df[7].dtype == np.int64):
            return True
        else:
            return False
        
    def chunks(self):
        dtype = self.meta 
        offset = 0
        is_contain_mapq = self.is_contain_mapq()
        columns = (
                    [1, 2, 3, 4, 7] 
                   if (self.min_mapq > 0 and is_contain_mapq) 
                   else [1, 2, 3, 4]
                   )
        new_columns = (
                        ['chrom1', 'pos1', 'chrom2', 'pos2', 'mapq'] 
                       if (self.min_mapq > 0 and is_contain_mapq) 
                       else ['chrom1', 'pos1', 'chrom2', 'pos2']
                       )

        # try:
        #     reader = pl.read_csv_batched(self.file, separator='\t', has_header=False,
        #                                     comment_prefix='#', columns=columns,
        #                                     new_columns=new_columns, 
        #                                     dtypes=dtype,
        #                                     n_threads=self.threads)
        #     chunks = reader.next_batches(1000)
        #     while chunks:
        #         chunk = pl.concat(chunks)

        #         if self.min_mapq > 0:
        #             chunk = chunk.filter(pl.col('mapq') >= self.min_mapq).drop(['mapq'])
                
        #         chunks = reader.next_batches(1000)
        #         yield chunk
    
        #     yield chunk
        # except pl.exceptions.ComputeError:
        # while True:
        #     try:
        #         chunk = (pl.read_csv(self.file, separator='\t', has_header=False, 
        #                     comment_prefix='#', columns=columns, 
        #                     new_columns=new_columns, 
        #                     dtypes=dtype, 
        #                     n_rows=self.chunksize, 
        #                     skip_rows=offset,
        #                     n_threads=self.threads)
        #                 )
        
        #     except pl.exceptions.NoDataError: 
        #         break

        #     if chunk.height == 0:
        #         break
        #     if chunk.height < self.chunksize:
        #         break

        #     offset += chunk.height
        #     if self.min_mapq > 0 and is_contain_mapq:
        #         chunk = chunk.filter(pl.col('mapq') >= self.min_mapq).drop(['mapq'])

        #     yield chunk
        
        # if self.min_mapq > 0 and is_contain_mapq:
        #     chunk = chunk.filter(pl.col('mapq') >= self.min_mapq).drop(['mapq'])

        # yield chunk
        if self.min_mapq > 0 and is_contain_mapq:
            input_file = self.file
            # cmd = f"cphasing-rs pairs-filter {self.file}  -q {self.min_mapq}"
            # process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            # stdout, stderr = process.communicate()
            # input_file = stdout
            
        else:
            input_file = self.file

        chunk = (pl.read_csv(input_file, separator='\t', has_header=False, 
                            comment_prefix='#', columns=columns, 
                            new_columns=new_columns, 
                            dtypes=dtype, 
                            n_threads=self.threads)
                        )
        if self.min_mapq > 0 and is_contain_mapq:
            chunk = chunk.filter(pl.col('mapq') >= self.min_mapq).drop(['mapq'])
        
        return chunk

    def process_chunk(self, chunk, bin_offset_db, binsize):
        bin1_id = (chunk['pos1'] // binsize) + chunk['chrom1'].map_elements(bin_offset_db.get)
        bin2_id = (chunk['pos2'] // binsize) + chunk['chrom2'].map_elements(bin_offset_db.get)
        chunk = chunk.with_columns([
            bin1_id.alias('bin1_id'),
            bin2_id.alias('bin2_id')
        ]).drop(['chrom1', 'chrom2', 'pos1', 'pos2'])
        # .filter(pl.col('bin1_id') <= pl.col('bin2_id'))

        # # swap bin1_id and bin2_id if bin1_id > bin2_id
        chunk = chunk.with_columns([
                pl.when(chunk['bin1_id'] > chunk['bin2_id'])
                .then(chunk['bin2_id'])
                .otherwise(chunk['bin1_id'])
                .alias('bin1_id'),
                pl.when(chunk['bin1_id'] > chunk['bin2_id'])
                .then(chunk['bin1_id'])
                .otherwise(chunk['bin2_id'])
                .alias('bin2_id')
            ])
        # .filter(pl.col('bin1_id') <= pl.col('bin2_id'))
        

    
        # chunk = chunk.with_columns([
        #     (pl.col('pos1') // binsize + pl.col('chrom1').map_elements(
        #                     bin_offset_db.get)).cast(self.bin_id_dtype).alias('bin1_id'),
        #     (pl.col('pos2') // binsize + pl.col('chrom2').map_elements(
        #                     bin_offset_db.get)).cast(self.bin_id_dtype).alias('bin2_id'),
        # ]).drop(['chrom1', 'chrom2', 'pos1', 'pos2'])

        chunk = chunk.group_by(['bin1_id', 'bin2_id'], maintain_order=True).agg(
            pl.count('bin1_id').alias('count'))
        
        chunk = chunk.sort(['bin1_id', 'bin2_id'])
     
        return chunk.to_pandas()
    
    def process_chunk_to_cool(self, bin_offset_db, binsize=10000):
        
        # with Pool(4) as pool:
        #     results = []
        #     for chunk in self.chunks():
            
        #         # return self.process_chunk(chunk, bin_offset_db, binsize)
        #         result = pool.apply_async(self.process_chunk, 
        #                                   args=(chunk, bin_offset_db, binsize))
        #         results.append(result)
            
        #     for result in results:
        #         yield result.get()
        # for chunk in self.chunks():
        #     yield self.process_chunk(chunk, bin_offset_db, binsize)
        return  self.process_chunk(self.chunks(), bin_offset_db, binsize)
    
    def to_cool(self, output, binsize=10000):
        import cooler
        from cooler.create import create, create_cooler

        if self.chromsize is None or self.chromsize.empty:
            raise ValueError("Chromsize in pairs file is required.")
        

        bins = cooler.binnify(
                        pd.DataFrame(self.chromsizes.set_index('chrom'))['length'], 
                        binsize=binsize)
        # bins = binnify(pd.DataFrame(self.chromsizes.set_index('chrom'))['length'], binsize=binsize)
        bins.to_csv('temp.bins', sep='\t', header=None, index=None)
        if len(bins) > 2**32:
            self.bin_id_dtype = pl.UInt64
        elif len(bins) > 2**16:
            self.bin_id_dtype = pl.UInt32
        elif len(bins) > 2**8:
            self.bin_id_dtype = pl.UInt16
        else:
            self.bin_id_dtype = pl.UInt8

        bin_offset = bins.groupby('chrom').size().cumsum()
        bin_offset = bin_offset.shift(1).fillna(0).astype(int)
        bin_offset_db = bin_offset.to_dict()
    
        iterator = self.process_chunk_to_cool(bin_offset_db, binsize)
        
        create(output, bins, iterator, 
                triucheck=False, dupcheck=False, boundscheck=False,)
        
        logger.info(f"Successful output cooler file into `{output}`.")

class MndTable:
    HEADER = ["str1", "chr1", "pos1", "frag1",
                "str2", "chr2", "pos2", "frag2",
                "score"]
    META = {"str1": "category", 
            "chr1": "category",
            "pos1": CHROM_COORD_DTYPE,
            "frag1": "category",
            "str2": "category",
            "chr2": "category",
            "pos2": CHROM_COORD_DTYPE,
            "frag2": "category",
            "score": "int32"
            }
    def __init__(self, mndtable):
        self.file = Path(mndtable)
        self.filename = self.file.name

        self.data = self.read_table()

    def read_table(self):
        df = pd.read_csv(self.file, sep='\s+', dtype=self.META,
                            usecols=range(9), names=self.HEADER,
                            header=None, index_col=None)
        
        return df 

    def to_pairs(self, chromsize, output):
        (
        self.data.assign(str1=lambda x: x.str1.map(lambda y: "+" if y =="0" else "-"),
                        str2=lambda x: x.str1.map(lambda y: "+" if y =="0" else "-"),
                        readID=lambda x: x.index,
                        pos1=lambda x: x.pos1 + 1,
                        pos2=lambda x: x.pos2 + 1
                        )
                        .rename(columns={"chr1": "chrom1",
                                            "chr2": "chrom2",
                                            "str1": "strand1",
                                            "str2": "strand2"})[Pairs.HEADER]
                        .to_csv(f"{output}.body", sep="\t", 
                                header=None, index=None)
        )

        ph = PairHeader([])
        ph.from_chromsize_file(chromsize)
        ph.save(f'{output}.header')

        command = f"""cat {output}.header {output}.body > {output}"""
        check_call(command, shell=True)
        os.remove(f'{output}.body')
        os.remove(f'{output}.header')

        logger.info(f'Successful written pairs file into `{output}`.')

class Contact:
    def __init__(self, coolfile, symmetric=True):
        import cooler
        self.coolfile = coolfile 
        self.cool = cooler.Cooler(self.coolfile)

        self.data = self.get_symmetric_pixels()
    
    @property
    def count_db(self):
        """
        get the symmetric contact counts database 
        """

        db = self.data.to_dict()['count']
        db = defaultdict(lambda :0, db)
        
        return db

    def get_symmetric_pixels(self):
        """
        get the symmetric pixels
        """

        pixels1 = self.cool.matrix(balance=False, join=True, 
                                       as_pixels=True)[:][['chrom1', 'chrom2', 'count']]
        
        pixels2 = pixels1.copy()
        pixels2.columns = ['chrom2', 'chrom1', 'count']
     
        pixels = pd.concat([pixels1, pixels2], axis=0)
        pixels.set_index(['chrom1', 'chrom2'], inplace=True)

        return pixels

    @property
    def is_symmetric(self):
        """
        Detect the contact data whether symmetry

        Returns:
        --------
        bool
            returns a boolefffan indicating whether the contact data is symmetric

        Examples:
        --------
        >>> contact.is_symmetric
        True
        """
       
        return (self.data.index.get_level_values(0) > 
                    self.data.index.get_level_values(1)).any()



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
        >>> contact.get_contacts(contig1, contig2)

        """
        assert self.is_symmetric, "The contact data must be symmetric."
        
        tmp_df = self.data.reindex(((i, j) for i in contig1  for j in contig2), fill_value=0)
        tmp_df = tmp_df.dropna(how='all', axis=0)

        return tmp_df

    def get_contact(self, contig1, contig2, count_re=None):
        """
        get unnormalized contact between two contig lists

        Params:
        --------
        contig1: list or list-like
            list of contigs
        contig2: list or list-like
            list of contigs
        count_re: CountRE
            count RE table

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

        L = tmp_df['count'].sum()
        if count_re:
            
            C1 = count_re.data.reindex(contig1)['RECounts'].sum()
            C2 = count_re.data.reindex(contig2)['RECounts'].sum()


            return L * 1000 / (C1 * C2)
        else:
            return L
        


    

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

def _process(args):
    n, df, low_mq_condition, output, tmpdir = args
    df = (df
        .assign(identity=lambda x: 
                (x.matches/x.aln_length).round(2).astype(np.float32))
        .assign(fragment_length=lambda x: (x.end - x.start))
        .assign(pass_filter=True)
        .assign(filter_reason="pass")
    )

    
    df = PAFTable._filter(n, df, low_mq_condition, tmpdir)
    df.to_parquet(f"{output}/part.{n}.parquet")

def _get_chunk(scans, read_idx, chunks, low_mq_condition, output, tmpdir):

    idx = 0
    for i, scan in enumerate(scans):
        try:
            df = chunks.get_chunk(scan)
        except StopIteration:
            return
        df['read_idx'] = read_idx[idx: idx + scan]
        idx += scan

        yield (i, df, low_mq_condition, output, tmpdir)


class Contacts:
    def __init__(self, filename):
        self.filename = filename
        self.file = Path(filename)
        if not self.file.exists():
            raise FileNotFoundError(f"File {filename} not found.")

        self.data = self.get_data()

    def get_data(self):
        df = pl.read_csv(self.file, separator='\t', has_header=False,
                            comment_prefix='#',
                            new_columns=['contig1', 'contig2', 'count']
        )

        return df.to_pandas()

    def __copy__(self):
        return Contacts(self.filename)
    
            

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
        'align_length', 'mapping_quality',
        'pass_filter', 'filter_reason'
        ]

    META = {
        "read_name": "category",
        "read_length": READ_COORD_DTYPE,
        "read_start": READ_COORD_DTYPE,
        "read_end": READ_COORD_DTYPE,
        "strand": pd.CategoricalDtype(["+", "-"], ordered=False),
        "chrom": "category",
        "chrom_length": CHROM_COORD_DTYPE,
        "start": CHROM_COORD_DTYPE,
        "end": CHROM_COORD_DTYPE,
        "matches": READ_COORD_DTYPE,
        "align_length": READ_COORD_DTYPE,
        "mapping_quality": MAPPING_QUALITY_DTYPE,
        "filter_reason": "category",
        "identity": PERCENTAGE_DTYPE
    }
    def __init__(self, paf: str, 
                    output: str = "out.pq",
                    min_quality: int = 1,
                    min_identity: float = 0.75,
                    min_length: int = 50,
                    no_read: bool = False,
                    use_dask: bool = False,
                    chunksize: int = 100000,
                    threads: int = 4):
        self.file = Path(paf)
        self.filename = self.file.name
        self.output = output
        self.min_quality = min_quality
        self.min_identity = min_identity
        self.min_length = min_length

        self.use_dask = use_dask
        self.chunksize = chunksize
        self.threads = threads
        self.tmpdir = tempfile.mkdtemp(prefix="tmp", dir='./')
        
        if self.file.suffix == ".gz":
            self.is_gz = True
        else:
            self.is_gz = False

        if self.use_dask:
            import dask.dataframe as dd
            from dask.distributed import Client
            self.client = Client(n_workers=self.threads) if self.threads > 1 else None
        
        if not no_read:
            self.read_table()
    
    
    def read_split_info(self):
        """
        import read_name to get the split information
        """
        read_names = pd.read_csv(self.file, sep='\t', usecols=[0],
                                 header=None)
        read_names = pd.Categorical(read_names[0], pd.unique(read_names[0]), ordered=True)
        
        read_idx = read_names.codes 
        del read_names
        gc.collect() 

        _, counts = np.unique(read_idx, return_counts=True)

        scan = list(range(0, len(read_idx), self.chunksize)) + [len(read_idx)]  

        a = np.searchsorted(np.cumsum(counts), scan[1:], side='right')
        b = np.searchsorted(read_idx, a, side='right')

        c = b - np.r_[0, b[:-1]]
      
        return read_idx, c

    def read_table(self):
        logger.info(f'Load paf file of `{self.filename}`.')
        if self.use_dask:
            df = dd.read_csv(self.file, sep='\t', usecols=range(12), 
                            header=None, names=self.PAF_HEADER[:12], 
                            dtype=self.META, blocksize=1e8)
            df['read_idx'] = df['read_name'].cat.as_known().cat.codes 
            # df = (df.set_index('read_idx')
            #         .repartition(df.npartitions)
            #         .reset_index())
            df = (df
                .assign(identity=lambda x: 
                        (x.matches/x.aln_length).round(2).astype(np.float32))
                .assign(fragment_length=lambda x: (x.read_end - x.read_start))
                .assign(pass_filter=True)
                .assign(filter_reason="pass")
            )
        else:
            
            Path(self.output).mkdir(exist_ok=True)
            
            read_idx, scans = self.read_split_info()
           
            chunks = pd.read_csv(self.file, sep='\t', usecols=range(12), 
                            header=None, names=self.PAF_HEADER[:12], 
                            dtype=self.META, iterator=True)
            

            low_mq_condition = (self.min_quality, self.min_identity, self.min_length)
            filter_condition = (f'map_quality < {self.min_quality}'
                                     f' | identity < {self.min_identity}'
                                     f' | fragment_length < {self.min_length}')
            logger.info(f"Filter low mapping quality by: {filter_condition}")
            

            new_chunk_iter = _get_chunk(scans, read_idx, chunks, low_mq_condition, self.output, self.tmpdir)
            
            pool = Pool(self.threads)
            pool.map(_process, new_chunk_iter)
            
            # for args in new_chunk_iter:
            #     _process(args)

            logger.info(f"Successfully output pore_c_table to `{self.output}`")

    
    @staticmethod
    def _filter_low_mq(df,
                    min_quality: int=1, 
                    min_identity: float=0.75,
                    min_length: int=10,
                    use_dask: bool=False):
        if use_dask:
            df["pass_mq"] = ((df.mapping_quality >= min_quality)
                            & (df.identity >= min_identity)
                            & (df.fragment_length >= min_length))
            df["pass_filter"] = df.pass_filter & df.pass_mq
            df["filter_reason"] = df.apply(
                lambda x: x.filter_reason if x.pass_mq else "low_mq", axis=1)

            df = df.drop('pass_mq', axis=1)

        else:
            low_mq_idx = ((df.mapping_quality < min_quality)
                            | (df.identity < min_identity)
                            | (df.fragment_length < min_length))

            df.loc[low_mq_idx, "pass_filter"] = False
            df.loc[low_mq_idx, "filter_reason"] = "low_mq"
        
        return df  

    @staticmethod
    def _filter(i, df, low_mq_condition, tmpdir, use_dask=False):
        if use_dask:
            df = df.partitions[i].compute()

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

        if use_dask:
            df.to_csv(f"{tmpdir}/tmp.part{i}.csv", sep='\t', 
                            index=False, header=True)
        
        return df

    def filter(self):
        
        low_mq_condition = (self.min_quality, self.min_identity, self.min_length)
        filter_condition = (f'map_quality <= {self.min_quality}'
                                     f' | identity <= {self.min_identity}'
                                     f' | fragment_length <= {self.min_length}')
        logger.info(f"Filter low mapping quality by: {filter_condition}")
        logger.info("Filter singleton ...")
        
        if self.use_dask:
            args = []
            for i in range(self.data.npartitions):
                args.append((i, self.data, low_mq_condition, self.tmpdir))
            
            Parallel(n_jobs=self.threads, verbose=2)(
                    delayed(PAFTable._filter)(i, j, k, l)
                        for i, j, k, l in args)

            self.data = dd.read_csv(f"{self.tmpdir}/tmp.*.csv", sep='\t', 
                                    header=0, dtype=self.META)
        else:
            self.data = PAFTable._filter(0, self.data, low_mq_condition, self.tmpdir)

        logger.info("Filter done.")
        
    def to_pairs(self, chromsize, output):
        """
        developing function.
        """
        logger.info('Converting Pore-C alignment to pairs.')
        from shutil import which
        from subprocess import PIPE, Popen
        if which('u4falign') is None:
            raise ValueError(f"`u4falign`: command not found.")
        cmd = ["u4falign", "paf23ddna", "-q", str(self.min_quality),
                "-p", str(self.min_identity), "-l", str(self.min_length),
                "/dev/stdin"]

        pipelines = []
        try:
            if self.is_gz:
                read_cmd = ["gzip", "-c", "-d", str(self.file)]
            else:
                read_cmd = ["cat", str(self.file)]
            pipelines.append(
                Popen(read_cmd, stdout=PIPE, 
                        stderr=open(f'{self.filename}.paf2pairs.log', 'w'),
                        bufsize=-1)
            )
            pipelines.append(
                Popen(cmd, stdin=pipelines[-1].stdout, 
                        stdout=open(f'{self.tmpdir}/temp.mnd.txt', 'w'), 
                        stderr=open(f'{self.filename}.paf2pairs.log', 'w'),
                        bufsize=-1)
            )
            pipelines[-1].wait()
        finally:
            for p in pipelines:
                if p.poll() is None:
                    p.terminate()
            else:
                assert pipelines != [], \
                    "Failed to execute command, please check log."
        
        mt = MndTable(f"{self.tmpdir}/temp.mnd.txt")
        mt.to_pairs(chromsize, output)

    def to_pore_c_table(self):
        data = self.data[["read_idx", "chrom", "start", 
                            "end", "strand", "read_name", 
                            "read_start", "read_end",
                            "mapping_quality", "identity",
                            "pass_filter", "filter_reason"]]

        pct = PoreCTable()
        if self.use_dask:
            pct.from_dask(data)
        else:
            pct.from_pandas(data)

        return pct

    def clean_tempoary(self):
        shutil.rmtree(self.tmpdir)

class PoreCTable:
    import cooler
    HEADER = ["read_idx", "read_length",
              "read_start", "read_end", 
              "strand", "chrom", "start", 
              "end", "mapping_quality", 
              "identity", "filter_reason"]
    
    META = {
        "chrom": "category",
        "start": CHROM_COORD_DTYPE,
        "end": CHROM_COORD_DTYPE,
        "read_start": READ_COORD_DTYPE,
        "read_end": READ_COORD_DTYPE,
        "strand": pd.CategoricalDtype(["+", "-"], ordered=False),
        "mapping_quality": MAPPING_QUALITY_DTYPE,
        "filter_reason": "category",
        # "identity": PERCENTAGE_DTYPE, 
    }
    def __init__(self, threads=4, use_dask=False):
        self.threads = threads
        self.use_dask = use_dask

        if self.use_dask:
            import dask.dataframe as dd
            from dask.distributed import Client
            self.client = Client(n_workers=self.threads)
        else:
            pandarallel.initialize(nb_workers=self.threads, verbose=0)

    def read_table(self, pore_c_table, columns=None):
        if self.use_dask:
            self.data = dd.read_parquet(pore_c_table, engine=PQ_ENGINE, columns=None, )
            self.data = self.data.set_index('read_idx').reset_index()
        else:
            self.data = pd.read_csv(pore_c_table, sep='\s+', comment='#', 
                                    names=self.HEADER, header=None, index_col=None)
            
            self.data = self.data.astype(self.META)

    def from_dask(self, data):
        self.use_dask = True
        self.data = data
    
    def from_pandas(self, data):
        self.use_dask = False
        self.data = data.astype(self.META)

    def to_pairs(self, chromsize, output):
        """
        Deprecated, Please use cphasing-rs porec2pairs.
        """
        def _to_pairs_per_concatemer(df):
            rows = list(
                df.sort_values(['chrom', 'pos']).itertuples()
            )
            read_name = df['read_name'][0]

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

        def _to_pairs(df, use_dask=False):
            if use_dask:
                df = df.compute()
            return df.groupby('read_idx').parallel_apply(_to_pairs_per_concatemer)
    
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
        if self.use_dask:
            args = []
            for i in range(df.npartitions):
                args.append(df.partitions[i])
            
            res = Parallel(n_jobs=min(self.threads, len(args)))(
                delayed(_to_pairs)(i) for i in args
            )
            res_df = pd.concat(res, axis=0)
        else:
            
            res_df = _to_pairs(df)
        
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
    def _read_and_alignment_stat(df, use_dask=True):
        if use_dask:
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
        if self.use_dask:
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
                PoreCTable._read_and_alignment_stat(self.data, use_dask=False)
    
        read_statistics_df.columns = ['count']
        read_statistics_df = read_statistics_df.rename_axis("item")

        if read_counts:
            read_statistics_df.loc['total reads'] = read_counts
            read_statistics_df.loc['unmapping'] = \
                read_counts - read_statistics_df.loc['mapping']
            read_statistics_df['percent'] = \
                read_statistics_df['count'] / read_counts 

            read_statistics_df = read_statistics_df.loc[['total reads', 
                                                        'pass', 'singleton', 
                                                        'low_mq', 'unmapping']]

        return read_statistics_df, align_statistics_df

    def contact_stat(self):
        if self.use_dask:
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

    @staticmethod
    def _chrom2contig(df, source_gr, columns, threads, meta=None):
        df = df.reset_index()
        query_gr = pr.PyRanges(df.rename(
                                    columns={
                                        'chrom': 'Chromosome',
                                        'start': 'Start', 
                                        'end': 'End',
                                        'index': 'align_idx'
                                    }
                                ))

        overlaps = query_gr.join(source_gr, nb_cpu=1).new_position('intersection')
        overlaps = overlaps.df.assign(
            chrom=lambda x: x['Name'],
            start=lambda x: x.eval('Start - Start_b'),
            end=lambda x: x.eval('End - Start_b'),
            ).astype(PoreCTable.META)

        del query_gr, source_gr, df
        gc.collect()

        tmp_columns = ['align_idx', 'read_start', 'read_end']
        overlaps_duplicates_first = \
            overlaps[tmp_columns][overlaps.align_idx.duplicated(keep='first')]
        overlaps_duplicates_last = \
            overlaps[tmp_columns][overlaps.align_idx.duplicated(keep='last')]
        overlaps_duplicates_first = (
            overlaps_duplicates_first.reset_index()
                                     .set_index('align_idx'))
        overlaps_duplicates_last = (
            overlaps_duplicates_last.reset_index()
                                    .set_index('align_idx'))

        overlaps_duplicates_last['length'] = overlaps_duplicates_last.eval('read_end - read_start')
        overlaps_duplicates_last['read_end'] = overlaps_duplicates_last.eval('read_start + length')
        overlaps_duplicates_first['read_start'] = \
            overlaps_duplicates_first['read_start'] + overlaps_duplicates_last['length']
        
        overlaps_duplicates_first = overlaps_duplicates_first.set_index('index')
        overlaps_duplicates_last = overlaps_duplicates_last.set_index('index')

        overlaps.loc[overlaps_duplicates_last.index, 
                    'read_end'] = overlaps_duplicates_last['read_end']
        overlaps.loc[overlaps_duplicates_first.index, 
                    'read_start'] = overlaps_duplicates_first['read_start']

        del overlaps_duplicates_first, overlaps_duplicates_last
        gc.collect()

        overlaps = (overlaps.sort_values(['read_idx', 'read_start'])
                    .assign(align_idx=lambda x: range(len(x))))
        overlaps = overlaps[columns]
        overlaps = overlaps.reset_index(drop=True)

        
        return overlaps[columns]

    def chrom2contig(self, source):
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
            source_df.columns = ['Chromosome', 'Start', 
                                  'End', 'Name']
            source_gr = pr.PyRanges(source_df)
        elif isinstance(source, pr.PyRanges):
            source_df = source.df
            source_gr = source
        else:
            raise TypeError("source must be string/dataframe/pyranges.")
        
        res_df = PoreCTable._chrom2contig(self.data, source_gr, 
                                            self.data.columns, self.threads)
        
        return res_df

    def intersection(self, bed, output):
        """
        According several regions to select contact 

        """
        assert self.use_dask == False, "Not support for dask dataframe"

        if isinstance(bed, str):
            logger.info(f'Load source infomation of `{bed}`.')
            bed_df = pd.read_csv(bed, sep='\s+', 
                                        header=None, 
                                        index_col=None,
                                        usecols=[0, 1, 2],
                                        names=['Chromosome', 'Start', 
                                                'End'])
            bed_gr = pr.PyRanges(bed_df)
        elif isinstance(bed, pd.DataFrame):
            bed_df = bed
            bed_df.columns = ['Chromosome', 'Start', 
                                                'End']
            bed_gr = pr.PyRanges(bed_df)
        elif isinstance(bed, pr.PyRanges):
            bed_df = bed.df
            bed_gr = bed
        else:
            raise TypeError("source must be string/dataframe/pyranges.")

        query_gr = pr.PyRanges(self.data.reset_index()[['chrom', 'start', 'end', 'index']]
                                .rename(columns={'chrom': 'Chromosome', 
                                                'start': 'Start', 
                                                'end': 'End',
                                                'index': 'Index'}))

        res_df = query_gr.join(bed_gr, nb_cpu=self.threads).df
        if res_df.empty:
            logger.warning('Could not found anything.')
            return
        res_index = res_df['Index'].astype(int)

        del query_gr, res_df
        gc.collect()

        self.data.loc[res_index].sort_index().to_csv(output, sep='\t', index=False, header=False)
        logger.info(f"Successful output result in `{output}`")


    def binnify(self, contigsizes, binsize=10000):
        """
        Divide contig into several bins

        Params:
        --------
        contigsizes: `pd.DataFrame`,
                dataframe of contigsizes, and chromosome is the index
        binsize: `int`, 
                the size of a bin
        """

        bins = cooler.util.binnify(contigsizes.iloc[:, 0], binsize)
        # self.data = (self.data.eval('pos=(start + end )/2')
        #                 .astype({'pos': 'int'})
        #                 .eval(f'start=pos - pos % {binsize}')
        #                 .eval(f'end=start + {binsize}'))

        bins.columns = ['Chromosome', 'Start', 'End']
        bins_gr = pr.PyRanges(bins)


        tmp_data = (self.data.reset_index()
                                .eval("pos=(start + end)/2")
                                .astype({'pos': 'int'})[['chrom', 'pos', 'index']])
        tmp_data = tmp_data.eval('start = pos -1')

        tmp_data.columns = ['Chromosome',  'End', 'Index', 'Start',]
        tmp_data = tmp_data[['Chromosome', 'Start', 'End', 'Index']]
        tmp_gr = pr.PyRanges(tmp_data)
        res_df = tmp_gr.join(bins_gr).df
        del tmp_data, tmp_gr 
        gc.collect()

        res_df = res_df[['Chromosome', 'Start_b', 'End_b', 'Index']]
        res_df.columns = ['chrom', 'start', 'end', 'index']
        res_df = res_df.set_index('index')
        self.data['start'] = res_df['start']
        self.data['end'] = res_df['end']

        self.data = self.data[self.HEADER]


    def divide_contig_into_nparts(self, contigsizes, n=5):
        """
        Divide contig into n parts

        Params:
        --------
        contigsizes: `pd.DataFrame`
                dataframe of contigsizes, and contig is the index

        n: `int`, 
            number of parts       
        """
        bins = contigsizes.copy()
        bins['start'] = contigsizes.iloc[:, 0].map(lambda x: np.linspace(0, x, 6).astype('int')[:-1])
        bins['end'] = contigsizes.iloc[:, 0].map(lambda x: np.linspace(0, x, 6).astype('int')[1:])
        bins.drop("length", axis=1, inplace=True)
        bins = bins.explode(['start', 'end'])
        
        bins = bins.reset_index()
        bins.columns = ['Chromosome', 'Start', 'End']
        bins_gr = pr.PyRanges(bins)


        tmp_data = (self.data.reset_index()
                                .eval("pos=(start + end)/2")
                                .astype({'pos': 'int'})[['chrom', 'pos', 'index']])
        tmp_data = tmp_data.eval('start = pos -1')

        tmp_data.columns = ['Chromosome',  'End', 'Index', 'Start',]
        tmp_data = tmp_data[['Chromosome', 'Start', 'End', 'Index']]
        tmp_gr = pr.PyRanges(tmp_data)
        res_df = tmp_gr.join(bins_gr).df
        del tmp_data, tmp_gr 
        gc.collect()

        res_df = res_df[['Chromosome', 'Start_b', 'End_b', 'Index']]
        res_df.columns = ['chrom', 'start', 'end', 'index']
        res_df = res_df.set_index('index')
        self.data['start'] = res_df['start']
        self.data['end'] = res_df['end']

        self.data = self.data[self.HEADER]
        

    def to_pao_csv(self, output="output.pao.csv"):
        """
        Convert pore-c table to csv, which can import into paohviz
        """
        tmp_data = self.data[["read_idx", "chrom", "start", "end"]]

        tmp_data['node_name'] = tmp_data.apply(lambda x: f"{x.chrom}:{x.start}-{x.end}", axis=1)
        tmp_data['time_slot'] = tmp_data['read_idx']
        tmp_data['edge_name_description'] = tmp_data['read_idx']
        tmp_data = tmp_data[['read_idx', 'node_name', 'time_slot', 'edge_name_description', 'chrom']]
        tmp_data.columns = ['edge_id', 'node_name', 'time_slot', 
                            'edge_name_description', 'group_name']
        
        tmp_data.to_csv(output, header=None, index=None)

    
    def save(self, output, tmpdir="/tmp"):
        if not self.use_dask:
            # self.data.to_csv(f"{tmpdir}/tmp.pore_c.csv", header=True, index=None)
            # df = dd.read_csv(f"{tmpdir}/tmp.pore_c.csv", dtype=self.META,
            #             header=0)
            self.data.to_parquet(output, 
                            engine=PQ_ENGINE, 
                            # write_metadata_file=True, 
                            # write_index=False, 
                            version=PQ_VERSION)
        else:
            self.data.to_parquet(output, 
                                engine=PQ_ENGINE, 
                                write_metadata_file=True, 
                                write_index=False, 
                                version=PQ_VERSION)
    
        logger.info(f"Successfully output pore_c_table to `{output}`")
