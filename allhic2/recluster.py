#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
recluster artition results by allele table
"""

import argparse
import logging
import os
import os.path as op
import sys

import numpy as np
import pandas as pd

from collections import Counter, OrderedDict, defaultdict
from joblib import Parallel, delayed
from pytools import natsorted

from allhic2.core import (
    AlleleTable,
    ClusterTable,
    PairTable,
    CountRE
)
from allhic2.utilities import list_flatten

logger = logging.getLogger(__name__)


class reCluster(object):
    """

    Params:
    -------
    cluster_file: `str`, cluster file from allhic
    at_file: `str`, allele table from allhic
    ploidy: `int`, ploidy of assembly
    threads: `int`, threads of programs [4]
    iter_round: `int`, round of rescue iteration [5]

    Returns:
    -------

    Examples:
    -------
    >>> recluster("cluster.txt", "Allele.ctg.table", 4)

    """
    def __init__(self, cluster_file, at_file, ploidy,
                countRE_file=None, pairs_file=None,
                threads=4, iter_round=5):

        self.clusterFile = cluster_file
        self.atFile = at_file
        self.countREFile = countRE_file
        self.pairFile = pairs_file
        self.ploidy = ploidy
        self.threads = threads
        self.iter_round = iter_round

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
        _passed = sorted(set(tmp[~pd.isna(tmp)]))

        return _passed

    @property
    def passed_dict(self):
        """
        dict of passed contigs and its group
        """
        l = list(map(lambda x: list(x.items()),
                    self.reClusterTable[self.groups].to_dict('records')))

        _passed_set = set([i[::-1] for item in l for i in item
                                if pd.isna(i[1]) is not True])
        _passed_dict = dict(_passed_set)

        return _passed_dict

    @property
    def uncluster(self):
        """
        uncluster contigs
        """
        return sorted(set(list_flatten(
            self.reClusterTable['uncluster'].str.split(",").dropna())))

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
        """
        return {'Total allelic': len(self.reClusterTable),
                'Full cluster allelic': len(self.reClusterTable[self.groups].dropna(how='any', axis=0)),
                'Uncluster allelic': len(self.reClusterTable['uncluster'].dropna()),
                'groups': dict(map(lambda x: (x, len(self.group(x))), self.groups))}

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

        for i, item in row.iteritems():
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

        self.check()
        logger.info(f'Successful rescue `{n}` allelic.')

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
                    self.reClusterTable[self.groups].to_dict('record')))
        for i in set([i[::-1] for item in l for i in item
                                if pd.isna(i[1]) is not True]):
            contig, group = i
            if contig in db or contig in incorrcet_db:
                incorrcet_db[contig].append(group)
                incorrcet_db[contig].append(db.pop(contig))
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
        for i, items in self.reClusterTable['uncluster'].dropna().iteritems():
            _uncluster_list = items.split(',')
            for contig in items.split(','):
                if contig in self.passed_dict:
                    _group = self.passed_dict[contig]
                    _old = self.reClusterTable.loc[i, _group]
                    if pd.isna(_old) is True:
                        self.reClusterTable.loc[i, _group] = contig
                    _uncluster_list.remove(contig)

            self.reClusterTable.loc[i, 'uncluster'] = ",".join(_uncluster_list) \
                                                        if _uncluster_list else np.nan

    def rescue_by_contacts(self):
        """
        rescue contigs by Hi-C contacts
        """
        _uncluster = self.reClusterTable.uncluster.dropna().str.split(",").map(
            len)
        _uncluster_idx = _uncluster.sort_values().index
        _uncluster_df = self.reClusterTable.loc[_uncluster_idx]

        for i, item in _uncluster_df.uncluster.iteritems():
            tmp_uncluster = item.split(",")
            for contig in item.split(","):
                if contig in self.passed_dict:
                    suggest_group = self.passed_dict[contig]
                else:
                    res = list(
                        map(
                            lambda x: self.PairTable.get_normalized_contact(
                                [contig], self.group(x)), self.groups))
                    if any(res) is not True:
                        continue
                    res_dict = dict(zip(self.groups, res))
                    suggest_groups = sorted(res_dict,
                                            key=lambda x: res_dict[x],
                                            reverse=True)

                    suggest_group = suggest_groups[0]

                    _old = self.reClusterTable.loc[i, suggest_group]
                    if pd.isna(_old) is not True:
                        suggset_group = suggest_groups[1]
                        if res_dict[suggest_group] == 0:
                            continue

                self.reClusterTable.loc[i, suggest_group] = contig
                tmp_uncluster.remove(contig)
            self.reClusterTable.loc[i, 'uncluster'] = ','.join(tmp_uncluster) \
                                                            if tmp_uncluster else np.nan

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
        list
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

        return sorted(set(self.reClusterTable[_group].dropna().tolist()))

    def where(self, contig):
        """
        where the contig in the reClusterTable

        Params:
        -------
        contig: str
            contig in reClusterTable

        Returns:
        -------
        pd.DataFrame
            dataframe containing `contig`

        Examples:
        -------
        >>> rc.where('utg00001')
        """
        df = self.reClusterTable.loc[self.reClusterTable.where(
                        self.reClusterTable == contig).dropna(how='all').index]
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
        self.reClusterTable.to_csv(output, sep='\t', header=True, index=False)

    def to_clusterTable(self):
        """
        updata the new cluster results into ClusterTable

        Examples:
        --------
        >>> rc.to_clusterTable()
        """
        self.ClusterTable.from_frame(self.reClusterTable[self.groups])

    def run(self):
        """
        run pipeline of reCluster

        Examples:
        --------
        >>> rc.run()
        """
        iter_num = 1
        logger.info(f'------------- Round {iter_num} ------------')
        n = self.rescue_by_allele()

        while n > 0 and iter_num <= self.iter_round:
            iter_num += 1
            logger.info(f'------------- Round {iter_num} ------------')
            n = self.rescue_by_allele()

        logger.info(f'Done {iter_num} round rescue')


        if self.pairFile:
            self.rescue_by_contacts()
            self.check()
            self.rescue_by_allele()
            self.check()

        self.to_clusterTable()
        self.ClusterTable.to_countRE(self.countREFile)
