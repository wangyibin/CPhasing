#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
recluster partition results by allele table
"""

import logging
import os
import os.path as op
import sys

import numpy as np
import pandas as pd
import pprint

from collections import OrderedDict, defaultdict
# from joblib import Parallel, delayed
# from pytools import natsorted

from ..core import (
    AlleleTable,
    ClusterTable,
    PairTable,
    CountRE
)
from ..utilities import list_flatten

logger = logging.getLogger(__name__)


class reCluster(object):
    """
    recluster partition results by allele table

    Params:
    --------
    cluster_file: str
        cluster file from extract
    at_file: str
        allele table
    ploidy: int
        ploidy of assembly
    countRE_file: str
        countRE file from extract [None]
    pairs_file: str
        pairs table file from extract [None]
    threads: int
        threads of programs [4]
    
    Returns:
    -------
    object:
        object of recluster

    Examples:
    -------
    >>> recluster("cluster.txt", "Allele.ctg.table", 4)
    """
    def __init__(self, cluster_file: str, at_file: str, ploidy: int,
                countRE_file=None, pairs_file=None, threads=4):

        self.clusterFile = cluster_file
        self.atFile = at_file
        self.countREFile = countRE_file
        self.pairFile = pairs_file
        self.ploidy = ploidy
        self.threads = threads

        self.AlleleTable = AlleleTable(self.atFile)
        self.ClusterTable = ClusterTable(self.clusterFile)

        if self.countREFile:
            self.CountRE = CountRE(self.countREFile)
            self.AlleleTable.filter(self.CountRE.contigs)

        if self.pairFile:
            self.PairTable = PairTable(self.pairFile)
        self._contig_groups = self.ClusterTable.contig_groups_dict
        self.groups = self.ClusterTable.groups
      
        self.groupIdx = dict(zip(self.groups, range(1, len(self.groups) + 1)))

        self._contig_groups_idx = self.convert_group2idx()

        self.reClusterTable =  self.sort_by_group().reset_index()

    @property
    def passed(self):
        """
        list of passed contigs
        """

        tmp = self.reClusterTable[self.groups].values.flatten()
        tmp = pd.unique(tmp)
        _passed = sorted(tmp[~pd.isna(tmp)])

        return _passed

    @property
    def passed_dict(self):
        """
        dict of passed contigs and its group
        """
        _passed_set = list_flatten(map(lambda x: set(zip(
                            self.group(x), [x] * len(self.group(x)))), self.groups))
        
        _passed_dict = dict(_passed_set)

        return _passed_dict

    @property
    def uncluster_header(self):
        """
        header of uncluster
        """
        
        _length = self.reClusterTable['uncluster'].str.split(",", expand=True).values.shape[1]

        return [f'uncluster{i}' for i in range(1, _length + 1)]

    @property
    def uncluster_df(self):
        """
        uncluster dataframe
        
        Returns:
        --------
        pd.DataFrame:
            dataframe of uncluster
        
        Examples:
        --------
        >>> rc.uncluster_df
        uncluster1               uncluster2 uncluster3 uncluster4
        0                            NaN                      NaN        NaN        NaN
        1       utg004980l_348000_379999                      NaN        NaN        NaN
        2                            NaN                      NaN        NaN        NaN
        """
        _uncluster_df = self.reClusterTable['uncluster'].str.split(",", expand=True)
        _uncluster_df.columns = self.uncluster_header
        _uncluster_df = _uncluster_df.fillna(value=np.nan)

        return _uncluster_df

    @property
    def uncluster(self):
        """
        uncluster contigs
        """
        _uncluster = self.reClusterTable['uncluster'].str.split(",", expand=True)
        _uncluster = _uncluster.values.flatten()
        _uncluster = pd.unique(_uncluster)
        _uncluster = _uncluster[~pd.isna(_uncluster)]

        return sorted(_uncluster)

    @property
    def stat(self):
        """
        The statistics of reCluster

        Returns:
        --------
        Dictionary of statistics

        Examples:
        --------
        >>> rc.stat
        {'Total allelic': 1896,
            'Full cluster allelic': 897,
            'Uncluster allelic': 4,
            'Total contigs': 2215,
            'Total length': 416977289,
            'Anchored contigs': 2205,
            'Anchored length': 415973835,
            'Anchored rate': ' 99.76%',
            'groups': {'Chr01g1': 554, 
                        'Chr01g2': 570, 
                        'Chr01g3': 537, 
                        'Chr01g4': 544},
            'lengths': {'Chr01g1': 102008283,
            'Chr01g2': 105882123,
            'Chr01g3': 104950420,
            'Chr01g4': 103133009}}
        """
        return {'Total allelic': len(self.reClusterTable),
                'Full cluster allelic': len(self.reClusterTable[self.groups].dropna(how='any', axis=0)),
                'Uncluster allelic': len(self.reClusterTable['uncluster'].dropna()),
                'Total contigs': self.CountRE.ncontigs,
                'Total length': self.CountRE.length, 
                'Anchored contigs': len(self.passed),
                'Anchored length': self.CountRE.get_group_length(self.passed),
                'Anchored rate': f'{self.CountRE.get_group_length(self.passed)/self.CountRE.length:.2%}',
                'groups': dict(map(lambda x: (x, len(self.group(x))), self.groups)),
                'lengths': dict(map(lambda x: (x, self.CountRE.get_group_length(self.group(x) )), self.groups))
                }

    def convert_group2idx(self):
        _contig_groups_tuple = list(zip(*self._contig_groups.items()))
        _contig_groups_idx = list(map(lambda x: self.groupIdx[x],
                                        _contig_groups_tuple[1]))

        return dict(zip(_contig_groups_tuple[0], _contig_groups_idx))

    def sort_by_group(self):
        """
        sort allelic contigs by cluster table and move duplicated contigs doubt, 
            unpartiton contigs to uncluster
        
        Examples:
        --------
        >>> reClusterTable = rc.sort_by_group()
        """
        return self.AlleleTable.data.apply(self._sort_by_group, axis=1)

    def _sort_by_group(self, row):
        """
        pandas row apply function of sort_by_group

        Returns:
        -------
        pd.Series
            contain groups, uncluster.

        Examples:
        --------
        >>> rc.AlleleTable.data.apply(rc._sort_by_group, axis=1)
        """
        res = OrderedDict()
        res_list = []
        ## contig unpartition are move to uncluster
        res['uncluster'] = set()

        for i, item in row.items():
            try:
                idx = self._contig_groups[item]
                res_list.append((idx, item))

            except KeyError:
                if pd.isna(item) is True:
                    continue
                res['uncluster'].add(item)
                continue

        ## find duplicated partition groups
        res_df = pd.DataFrame(res_list)
        res_dict = dict(res_list)
        if res_df.empty is not True:
            res_df = res_df.set_index(0)
            dup_cond = res_df.index.duplicated(keep=False)
            unique_df = res_df[~dup_cond]
            dup_df = res_df[dup_cond]
            if dup_df.empty is not True:
                dup_ctg = dup_df[1]
                res['uncluster'].update(dup_ctg)

            res_df = unique_df
            res_dict = res_df.to_dict()[1]

        res['uncluster'] = ",".join(sorted(set(res['uncluster']))) \
                                        if res['uncluster'] else np.nan

        res.update(res_dict)

        return pd.Series(res)

    def rescue_by_allele(self):
        """
        rescue by allele table, which have four allelic contigs and only one uncluster

        Returns:
        --------
        number of rescued allelic contigs

        Examples:
        --------
        >>> rc.rescue_by_allele()
        520
        """
        n = 0
        self.groups = list(filter(lambda x: x in self.reClusterTable.columns[1:-1], self.groups))

        db = OrderedDict()
        for i, item in self.reClusterTable.iterrows():
            tmp_df = item[self.groups]

            if len(tmp_df.dropna()) == self.ploidy - 1:
                if pd.isna(item.uncluster) is True:
                    continue
                gap_group = tmp_df[tmp_df.isna()].index.tolist()[0]
                uncluster = item['uncluster']
                if uncluster not in db:
                    db[uncluster] = dict(zip(self.groups, [0]*len(self.groups)))
                db[uncluster][gap_group] += 1

        for uncluster in db:
            sorted_list = sorted(db[uncluster].items(), key=lambda x: x[1], reverse=True)
            if sorted_list[0][1] != 0 and sorted_list[0][1] == sorted_list[1][1]:
                continue
            gap_group = sorted_list[0][0]
            tmp_df = self.where(uncluster)

            def clean(x):
                if pd.isna(x) is True:
                    return True
                if len(x) <= 1:
                    return True
                else:
                    return False

            tmp_df = tmp_df[tmp_df['uncluster'].str.split(',').map(clean)]

            tmp_gap_df = tmp_df[gap_group].copy()
            tmp_df[gap_group] = tmp_df['uncluster']
            tmp_df['uncluster'] = tmp_gap_df

            self.reClusterTable.loc[tmp_df.index] = tmp_df
            n += len(tmp_df)

        return n

    def check(self):
        """
        check reClusterTable and correct it, such as a contig partition into 
            different groups.

        Examples:
        ---------
        >>> rc.check()
        """

        incorrcet_db = defaultdict(list)
        db = {}
        l = list(map(lambda x: list(x.items()),
                    self.reClusterTable[self.groups].dropna(how='all').to_dict('records')))

        for i in set([i[::-1] for item in l for i in item
                                if pd.isna(i[1]) is not True]):
            contig, group = i
            if contig in db or contig in incorrcet_db:
                try:
                    incorrcet_db[contig].append(group)
                    incorrcet_db[contig].append(db.pop(contig))
                except KeyError:
                    print(contig)
            else:
                db[contig] = group

        for contig in incorrcet_db:
            tmp_df = self.where(contig)
            self.reClusterTable.loc[tmp_df.index, self.groups] = \
                tmp_df[self.groups].replace(contig, np.nan)
            self.reClusterTable.loc[tmp_df.index, 'uncluster'] = \
                self.reClusterTable.loc[tmp_df.index, 'uncluster'].map(
                                        lambda x: self._uncluster_append(x, contig))

        ## check passed contigs in uncluster contigs
        _passed_dict = self.passed_dict
        _uncluster_df = self.uncluster_df.dropna(how='all')
        for i, items in _uncluster_df.iterrows():
            _uncluster_list = items.dropna().tolist()
            for contig in _uncluster_list.copy():
                try:
                    _group = _passed_dict[contig]
                except KeyError:
                    continue
                else:
                    _old = self.reClusterTable.loc[i, _group]
                    if pd.isna(_old) is True:
                        self.reClusterTable.loc[i, _group] = contig
                    _uncluster_list.remove(contig)

            self.reClusterTable.loc[i, 'uncluster'] = ",".join(_uncluster_list) \
                                                        if _uncluster_list else np.nan

    def rescue_by_contacts(self):
        """
        rescue contigs by Hi-C contacts

        Returns:
        --------
        int:
            number of rescued allelic

        Examples:
        --------
        >>> rc.rescue_by_contacts()
        10
        """
        _uncluster_counts = self.uncluster_df.dropna(how='all').count(axis=1)
        _uncluster_idx = _uncluster_counts.sort_values(ascending=False).index
        _uncluster_df = self.uncluster_df.iloc[_uncluster_idx]

        n = 0
        _passed_dict = self.passed_dict 
        for i, item in _uncluster_df.iterrows():
            tmp_uncluster = item.dropna().tolist()
            tmp_item = tmp_uncluster.copy()
            for contig in tmp_item:
                if contig in _passed_dict:
                    suggest_group = _passed_dict[contig]
    
                else:
                    res = list(
                        map(
                            lambda x: self.PairTable.get_normalized_contact(
                                [contig], self.group(x)), self.groups))
                    if any(res) is False:
                        continue
                    res_dict = dict(zip(self.groups, res))
                    suggest_groups = sorted(res_dict,
                                            key=lambda x: res_dict[x],
                                            reverse=True)

                    suggest_group = suggest_groups[0]

                    _old = self.reClusterTable.loc[i, suggest_group]
                    j = 1 
                    while pd.isna(_old) is not True and j < len(suggest_groups):
                        suggest_group = suggest_groups[j]
                        _old = self.reClusterTable.loc[i, suggest_group]
                        j += 1
                        if res_dict[suggest_group] == 0:
                            break
                
                _passed_dict.update({contig: suggest_group})
                self.reClusterTable.loc[i, suggest_group] = contig
                tmp_uncluster.remove(contig)
              
            self.reClusterTable.loc[i, 'uncluster'] = ','.join(tmp_uncluster) \
                                                            if tmp_uncluster else np.nan
            
        return n

    def rescue_by_conflict(self, method='greedy'):
        """
        rescue contigs by conflict between shared contigs

        Params:
        --------
        method: str
            method of rescue_by_conflict {"greedy", "strict"} ["greedy"]
                greedy: greedy to rescue uncluster contigs, 
                        which contig rescued in first gap
                strict: strict to rescue uncluster contigs, 
                        which contig must only rescued to one gap
        
        Returns:
        --------
        int:
            number of rescued allelic

        Examples:
        --------
        >>> rc.rescue_by_conflict(method='greedy')
        20
        """
        assert method in {"greedy", "strict"}, \
                        f"method {method} not in {{'greedy', 'strict'}}" 

        _uncluster_df = self.uncluster_df.dropna(how='all')
        n = 0
        for contig in self.uncluster:
            tmp_df = _uncluster_df.loc[_uncluster_df.where(
                        _uncluster_df == contig).dropna(how='all').index]

            tmp_recluster_df = self.reClusterTable.loc[tmp_df.index, self.groups]
            bool_df = pd.isna(tmp_recluster_df).all(axis=0)
            
            candidate_groups = bool_df.index[bool_df]
            if candidate_groups.empty:
                continue
            
            if method == 'greedy':
                suggest_group = candidate_groups[0]
               
            if method == 'strict':
                if len(candidate_groups) > 1:
                    continue
                suggest_group = candidate_groups[0]
            
            self.reClusterTable.loc[tmp_df.index, suggest_group] = contig 
            self.reClusterTable.loc[tmp_df.index, 'uncluster'] = tmp_df.replace(
                contig, np.nan).apply(lambda x: ",".join(
                    x.dropna()) if not x.dropna().empty else np.nan, axis=1)
            
            n += len(tmp_df)

        return n

    def _uncluster_append(self, row, contig):
        if pd.isna(row) is True:
            return contig
        else:
            tmp_list = row.split(",")
            if contig not in tmp_list:
                tmp_list.append(contig)

            return ",".join(tmp_list)

    def _uncluster_remove(self, row, contig):
        if pd.isna(row) is True:
            return row
        else:
            tmp_list = row.split(',')
            tmp_list.remove(contig)
            if tmp_list:
                return ",".join(tmp_list)
            else:
                return np.nan

    def group(self, _group):
        """
        get contigs by group

        Params:
        --------
        _group: str
            group of cluster table

        Returns:
        --------
        list:
            list of contigs in a group
        
        Examples:
        --------
        >>> rc.group('Chr01g1')
        ['utg001', 
        'utg002',
        ...]
        """
        assert _group in self.groups or _group == 'uncluster', \
                    f"No such of group `{_group}` in {str(self.groups)}"

        _res = pd.unique(self.reClusterTable[_group])
        _res = _res[~pd.isna(_res)]
        return _res[np.argsort(_res)]
        #return sorted(set(self.reClusterTable[_group].dropna().tolist()))

    def where(self, contig, columns=None):
        """
        where the contig in the reClusterTable

        Params:
        -------
        contig: str
            contig in reClusterTable
        columns: array-like
            columns of reClusterTable
        
        Returns:
        -------
        pd.DataFrame
            dataframe containing `contig`

        Examples:
        -------
        >>> rc.where('utg00001')
                  0   Chr01g1  ...   Chr01g4     uncluster
        1136  Chr01  utg00001  ...  utg00001           NaN
        1137  Chr01  utg00001  ...  utg00001           NaN
        """
        df = self.reClusterTable.loc[self.reClusterTable.where(
                        self.reClusterTable == contig).dropna(how='all').index]
        if columns:
            df  = df[columns]

        return df


    def save(self, output):
        """
        save recluster table to output file

        Params:
        --------
        output: `str` or `stdout`
            output

        Examples:
        --------
        >>> rc.save('recluster.txt')
        """
        self.reClusterTable.to_csv(output, sep='\t', header=True, index=False, na_rep='NA')
        logger.info(f'Successful output reClusterTable into `{output}`')

    def to_clusterTable(self):
        """
        updata the new cluster results into ClusterTable

        Examples:
        --------
        >>> rc.to_clusterTable()
        """
        self.ClusterTable.from_frame(self.reClusterTable[self.groups])

    def run(self, method='greedy', iter_round=5):
        """
        run pipeline of reCluster

        Params:
        --------
        method: str
            method of rescue
        iter_round: int
            round of rescue iteration [5]

        Examples:
        --------
        >>> rc.run()
        """
        round = 1
        logger.info(f"[Round {round}]")
        iter_num = 1
        n = self.rescue_by_allele()
        self.check()

        while n > 0 and iter_num <= iter_round:
            iter_num += 1
            n += self.rescue_by_allele()
            self.check()

        logger.info(f"\tTotal {n} allelic were rescued")

        if self.pairFile:
            round += 1
            logger.info(f"[Round {round}]")
            n = self.rescue_by_contacts()
            
            self.check()
            n += self.rescue_by_allele()
            self.check()
            logger.info(f"\tTotal {n} allelic were rescued")
        
        round += 1
        logger.info(f"[Round {round}]")
        n = self.rescue_by_conflict(method=method)
        self.check()
        n += self.rescue_by_allele()
        self.check()
        n += self.rescue_by_contacts()
        self.check()
        logger.info(f"\tTotal {n} allelic were rescued")

        logger.info("Done")
        logger.info(pprint.pformat(self.stat))

        self.to_clusterTable()
        self.ClusterTable.to_countRE(self.countREFile)