#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
rewrite ALLHiC_prune
"""

import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

from itertools import combinations

from .core import (
    AlleleTable, 
    CountRE, 
    PairTable
    )


def Prune(at: AlleleTable, cr: CountRE, pt: PairTable, normalize=False):
    """
    Prune allelic signal by allele table.

    Params:
    --------
    at: AlleleTable
        AlleleTable object of allhic allele.table.
    cr: CountRE
        CountRE objcet of allhic count_RE.txt.
    pt: PairTable, unsymmetric
        PairTable object of allhic pairs.txt.
    normalize: bool, default False
        If normalize the Hi-C signal.

    Returns:
    --------
    None:
        But pairtable will be modify, where allelic signal will be removed.
    
    Example:
    --------
    >>> at = AlleleTable(alleletable)
    >>> cr = CountRE(count_re)
    >>> pt = PairTable(pairs, symmetric=False)
    >>> Prune(at, cr, pt)
    >>> pt.save(pt.filename.replace(".txt", ".prune.txt"))
    """
    contigs = cr.contigs
    contigs_index = dict(zip(contigs, range(len(contigs))))
    filter_func = lambda x: x in contigs

    allelic_pairs = []
    # for chrom, item in at.data.groupby(0):
    for _, alleles in at.data.iterrows():
        if len(alleles) == 1:
            continue
        alleles = list(filter(filter_func, alleles))
        alleles = sorted(alleles, key=lambda x: contigs_index[x])
        tmp_pairs = list(combinations(alleles, 2))
            
        allelic_pairs.extend(tmp_pairs)
    
    allelic_pairs = list(set(allelic_pairs))
    header = pt.HEADER
    pairs_df = pt.data

    valid_allelic_pairs = [pair for pair in allelic_pairs 
                                if pair in pairs_df.index]
    pt.data = (pairs_df.drop(valid_allelic_pairs, axis=0)
                .reset_index()[header]
    )
    pt.symmetric_pairs()
    pt.data = pt.data.set_index(['Contig1', 'Contig2'])

    if normalize:
        pt.normalize()

    at.data = at.data.drop_duplicates(at.data.columns)

    pairs_contigs = set(pt.Contig1)
    
    remove_db = set()
    retain_db = set()
    for _, item in at.data.iterrows():
        item = item.dropna()
        if len(item) == 1:
            continue
        
        item_list = item.values.tolist().copy()
        res = []
        for contig in item:
            if pd.isna(contig):
                continue
            if contig not in pairs_contigs:
                continue
            
            tmp_df = pt.data.loc[contig]
            if tmp_df.empty:
                item_list.remove(contig)
                continue

            if normalize:
                res.append(tmp_df['NormalizedLinks'])
            else:
                res.append(tmp_df['ObservedLinks'])

        
        if not res:
            continue
        if len(item_list) == 1:
            continue

        res_df = pd.concat(res, axis=1, keys=item_list)
        res_df = res_df[res_df.count(axis=1) > 1]

        if res_df.empty:
            continue
        
        all_db = (res_df.stack()
                    .index.values
                    .tolist())
        tmp_retain_db = list(res_df.idxmax(axis=1)
                    .reset_index()
                    .itertuples(index=False, name=None))
        retain_db.update(tmp_retain_db)

        tmp_remove_db = set(all_db) - set(tmp_retain_db)
        tmp_remove_db = tmp_remove_db - retain_db
        remove_db.update(tmp_remove_db)

    remove_db = set(remove_db)
    remove_db2 = set(map(lambda x: x[::-1], remove_db))
    remove_db = remove_db | remove_db2

    remove_db = [pair for pair in remove_db if pair in pt.data.index]

    pt.data = pt.data.drop(remove_db, axis=0)
    if normalize:
        pt.data = pt.data.drop(['NormalizedLinks'], axis=1)

    allelic_pairs2 = set(map(lambda x: x[::-1], allelic_pairs))
    allelic_pairs = set(allelic_pairs) | allelic_pairs2
    remove_db = set(remove_db) | set(allelic_pairs)

    with open("prune.contig.list", "w") as output:
        for pair in remove_db:
            if pair[0] > pair[1]:
                continue
            print("\t".join(pair), file=output)

    return remove_db