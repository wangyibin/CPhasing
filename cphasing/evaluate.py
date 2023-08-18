#!/usr/bin/env python
# -*- coding:utf-8 -*-


import argparse
import logging
import os
import os.path as op
import sys

from itertools import combinations

from .core import AlleleTable, ClusterTable
from .utilities import read_chrom_sizes, list_flatten

def allelic_error(cluster_table, alleletable, contigsizes):
    ct = ClusterTable(cluster_table)
    at = AlleleTable(alleletable, sort=False, fmt='allele2')
    df = read_chrom_sizes(contigsizes)

    allelic_pairs = set(map(tuple, at.data.query(f"similarity > 0.85 & mz1 > 10000")[[1, 2]].values.tolist()))
    
    for group in ct.groups:
        contigs = ct.data[group]
        if len(contigs) <= 1:
            continue
        pairs = set(combinations(contigs, 2))
        
        group_length = df.loc[contigs]['length'].sum()

        error_pairs = allelic_pairs & pairs
        error_contigs = list(set(list_flatten(list(error_pairs))))
    
        length = df.loc[error_contigs]['length'].sum()
        print(group, length/group_length)

