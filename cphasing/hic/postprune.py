#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import logging
import os
import os.path as op
import sys

import numpy as np
import pandas as pd

from collections import defaultdict, OrderedDict
from itertools import combinations, permutations
from joblib import Parallel, delayed

from cphasing.core import AlleleTable, ClusterTable, PairTable, CountRE
from cphasing.utilities import list_flatten

def post_prune(group, contigs, at, pt, cr):
    """
    
    Params:
    --------
    contigs: list

    at: AlleleTable

    pt: PairTable

    Returns:
    --------
    new_group: list

    Examples:
    --------
    >>> post_prune()

    """
    error_contigs = set() 
    new_contigs = set() 
    
    allelic_list = set(map(tuple, at.data[range(1, at.n + 1)].values))
    score_db = at.scores

    contig_length_db = cr.data['Length'].to_dict()

    contigs = sorted(contigs)
    contig_pairs = set(combinations(contigs, 2))
    
    
    conflict_contig_pairs = contig_pairs.intersection(allelic_list)
    conflict_contigs = list_flatten(conflict_contig_pairs)
    correct_contigs = set(contigs) - set(conflict_contigs)
    

    conflict_contig_pairs2 = []
    for i, pair in enumerate(conflict_contig_pairs):
        l1 = contig_length_db[pair[0]]
        l2 = contig_length_db[pair[1]]

        if l1 < l2:
            conflict_contig_pairs2.append(pair[::-1])
        else:
            conflict_contig_pairs2.append(pair)
        
    conflict_contig_pairs = sorted(conflict_contig_pairs2, 
                                        key=lambda x: (contig_length_db[x[0]], contig_length_db[x[1]]), 
                                        reverse=True)                                         

    res_db = defaultdict(lambda : 0)
    contact_db = {}
    for pair in conflict_contig_pairs:
        
        v1 = pt.get_contact([pair[0]], correct_contigs, cr)
        v2 = pt.get_contact([pair[1]], correct_contigs, cr)
       
        contact_db[pair] = (v1, v2)
    
        if v1 == 0 and v2 == 0:
            continue
        elif v1 > v2 and v2 != 0:
            res_db[pair[0]] += score_db[pair]
            res_db[pair[1]] -= score_db[pair]
            correct_contigs.add(pair[0])
        elif v1 < v2 and v1 != 0:
            res_db[pair[0]] -= score_db[pair]
            res_db[pair[1]] += score_db[pair]
            correct_contigs.add(pair[1])
        else:
            continue
    
    for contig in res_db:
        if contig == "utg000386l" or contig == "utg000097l":
            print(res_db[contig], file=sys.stderr)
        if res_db[contig] > 0:
            new_contigs.add(contig)
        elif res_db[contig] < 0:
            error_contigs.add(contig)
    

    new_contigs = set(contigs) - error_contigs
    
    return group, new_contigs, error_contigs 

def post_rescue(groups, error_contigs, at, cr, pt, mzTotal):
    """
    
    rescue error contigs after post prune 

    Params:
    --------
    groups: dict

    error_contigs: list

    """
    groups_array = np.array(list(groups.keys()))
    allele_db = at.data[range(1, at.n + 1)].set_index(1)
    score_db = at.scores
    
    new_cluster_db = OrderedDict()
    for contig in error_contigs:
        contacts = []
        allelic_score = []
        contig_mz = mzTotal[contig]
        score = []
        allelic_contigs = set(allele_db.loc[contig].values.flatten())

        for group in groups:
            group_contigs = groups[group]
            c = pt.get_contact([contig], group_contigs, cr)
            contacts.append(c)
            s = allelic_contigs.intersection(group_contigs)
            s = sum(score_db[(contig, i)] for i in s)
            allelic_score.append(s)
            score.append(c * (1 - s/contig_mz))
        
        contacts = np.array(contacts)
        allelic_score = np.array(allelic_score)
        score = np.array(score)
        
        zero_score_index = np.where(allelic_score == 0)
        if len(zero_score_index[0]) == 1:
            # print(contig, groups_array[zero_score_index][0], file=sys.stderr)
            new_cluster_db[contig] = groups_array[zero_score_index][0]
        else:
            if len(np.where(contacts == 0)[0]) != len(groups_array):
                if all(score <= 0):
                    continue
                
                new_cluster_db[contig] = groups_array[np.argmax(score)]
                # print(contig, groups_array[np.argmax(score)], file=sys.stderr)
            
                # print(contig, contacts, allelic_score, score, file=sys.stderr)

    for contig, group in new_cluster_db.items():
        groups[group].add(contig)
    
    return groups


def test(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('allele_table', 
            help='')
    pReq.add_argument('cluster_table')
    pReq.add_argument('pair_table')
    pReq.add_argument('count_re')
    pReq.add_argument('mz_total')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    at = AlleleTable(args.allele_table, sort=False, fmt="allele2")
    ct = ClusterTable(args.cluster_table)
    pt = PairTable(args.pair_table)
    cr = CountRE(args.count_re, minRE=1)
    
    mzTotal = dict((i.strip().split("\t")[0], int(i.strip().split("\t")[1]))
                            for i in open(args.mz_total) if i.strip())

    args = []
    for group in ct.data:
        args.append((group, ct.data[group], at, pt, cr))
    
    res = Parallel(n_jobs=10)(delayed(
                    post_prune)(i, j, k, l, m) for i, j, k, l, m in args)
    
    
    pruned_db = OrderedDict()
    error_contigs = []
    for (group, new_contigs, errors) in res: 
        pruned_db[group] = new_contigs
        error_contigs.extend(errors)
    
    # for (group, new_contigs, error_contigs) in res:
    #     print(group, len(new_contigs), " ".join(sorted(new_contigs)), 
    #                 sep='\t', file=sys.stdout)
    #     # print(group, len(error_contigs), " ".join(sorted(error_contigs)), 
        #             sep='\t', file=sys.stdout)

    new_groups = post_rescue(pruned_db, error_contigs, at, cr, pt, mzTotal)

    for group, new_contigs in new_groups.items():
        print(group, len(new_contigs), " ".join(sorted(new_contigs)), 
                     sep='\t', file=sys.stdout)
  


if __name__ == "__main__":
    test(sys.argv[1:])