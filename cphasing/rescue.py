#!/usr/bin/env python
# -*- coding:utf-8 -*-


import logging
import os
import os.path as op
import sys

from .core import (
    ClusterTable, 
    CountRE, 
    PairTable
    )

logger = logging.getLogger(__name__)

class Rescuer:
    def __init__(self, cluster_file, countRE_file, pairs_file,
                    exclude=[], min_score=0.01):
        self.clusterFile = cluster_file
        self.countREFile = countRE_file
        self.pairFile = pairs_file

        self.ClusterTable = ClusterTable(self.clusterFile)
        self.countRE = CountRE(self.countREFile, minRE=1)
        self.PairTable = PairTable(self.pairFile)

        self.exclude = exclude
        self.min_score = min_score

    @property
    def contigs(self):
        return self.countRE.contigs
    
    @property
    def anchored_contigs(self):
        return self.ClusterTable.contigs

    @property
    def unanchor_contigs(self):
        return list(set(self.contigs) - set(self.anchored_contigs))

    def rescue(self):
        logger.info('Starting rescue ...')
        pairs_df = self.PairTable.data.reset_index()

        pairs_df['Group1'] = pairs_df['Contig1'].map(
                                self.ClusterTable.contig_groups_dict)
        pairs_df['Group2'] = pairs_df['Contig2'].map(
                                self.ClusterTable.contig_groups_dict)

        pairs_df = pairs_df[pairs_df['Label'] == 'ok']
        pairs_df = pairs_df.set_index(['Contig1', 'Contig2'])
        
        
        for contig in self.unanchor_contigs:
            try:
                tmp_df = pairs_df.loc[contig]
            except KeyError:
                continue

            tmp_df = tmp_df.sort_values(by='ObservedLinks', ascending=False)
            tmp_df = tmp_df[~tmp_df['Group2'].isna()]
            if tmp_df.empty:
                continue
            sum_df = tmp_df.groupby('Group2')[['RE2', 'ObservedLinks']].sum()
            sum_df['score'] = sum_df['ObservedLinks'] / (sum_df['RE2'] + tmp_df.iloc[0]['RE1'])
            sum_df = sum_df.sort_values(by='score', ascending=False)

            for best_result, items in sum_df.iterrows():
                best_score = items['score']
                if self.exclude and best_result in self.exclude:
                    continue
                if best_score > self.min_score:
                    
                    self.ClusterTable.data[best_result].append(contig)
                    break
        
        logger.info('Done.')
    
    def save(self, output):
        self.ClusterTable.save(output)
        