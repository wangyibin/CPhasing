#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
evaluate the assembly by comparied with chromosome-level assembly

P = 
"""

import argparse
import logging
import os
import os.path as op
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


from collections import OrderedDict, Counter, defaultdict

from cphasing.core import ClusterTable

def assign_cluster():
    pass

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('cluster', 
            help='cluster table')
    pReq.add_argument('contigsizes',
            help='contig sizes')
    pReq.add_argument('real_list',
            help='real cluster list, two columns (chrom, contig)')

    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    real_db = OrderedDict()
    with open(args.real_list, 'r') as fp:
        for line in fp:
            if not line.strip():
                continue
            line_list = line.strip().split()[:2]
        
            chrom, contig = line_list
            if chrom not in real_db:
                real_db[chrom] = []
            real_db[chrom].append(contig)
    
    prefix = args.cluster.rsplit(".", 1)[0]

    contigsizes = dict(i.strip().split()[:2] for i in open(args.contigsizes) if i.strip())
    contigsizes = dict(map(lambda x: (x[0], int(x[1])), contigsizes.items()))
    
    chromsizes = defaultdict(lambda: 0)
    for contig in contigsizes:
        _chrom = contig.rsplit(".", 1)[0]
        chromsizes[_chrom] += contigsizes[contig]


    ct = ClusterTable(args.cluster)
    
    group_assign = OrderedDict()
    for group in ct.groups:
        _contigs = ct.data[group]
        _contigs_chrom = list(map(lambda x: x.split("_")[0], _contigs))
        _contigs_chrom_count = Counter(_contigs_chrom)
        _contigs_size = sum(contigsizes[i] for i in _contigs)
        main_chrom = max(_contigs_chrom_count, key=lambda x: _contigs_chrom_count.get(x))
        main_contigs = list(filter(lambda x: main_chrom in x, _contigs))
        main_chrom_size = sum(contigsizes[i] for i in main_contigs)
        
        if main_chrom_size < 0.5 * _contigs_size:
            continue
        
        _contigs_chrom_db = dict(zip(_contigs, _contigs_chrom))

        if main_chrom in group_assign:
            if len(_contigs) < len(group_assign[main_chrom][1]):
                continue
        group_assign[main_chrom] = (group, _contigs, _contigs_chrom, _contigs_chrom_count)

    
    res = OrderedDict()
    zero_chrom = []
    PR_output = open(f"{prefix}.PR", 'w')
    for real_chrom in real_db:
        try:
            group, _contigs, _contigs_chrom, _contigs_chrom_count = group_assign[real_chrom]
            main_contigs = list(filter(lambda x: real_chrom in x, _contigs))
            P = sum(list(map(lambda x: contigsizes.get(x), main_contigs))) \
                / sum(list(map(lambda x: contigsizes.get(x), _contigs)))
            R = sum(list(map(lambda x: contigsizes.get(x), main_contigs))) \
                / sum(list(map(lambda x: contigsizes.get(x), real_db[real_chrom])))
            
            F1 = (2 * P * R) / (P + R)

            print(real_chrom, P, R, F1, file=PR_output)
            res[real_chrom] = (P, R)
        except KeyError:
            print(real_chrom, 0, 0, 0, file=PR_output)
            zero_chrom.append(real_chrom)
            res[real_chrom] = (0, 0)
    
    plt.rcParams['font.family'] = 'Helvetica'
    fig, ax = plt.subplots(figsize=(5, 5))
    P_list, R_list = list(zip(*list(res.values())))
    ax.scatter(R_list, P_list, alpha=0.3, color='#00426d')
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["0", "20", "40", "60", "80", "100"], fontsize=16, fontname='Helvetica')
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels(["0", "20", "40", "60", "80", "100"], fontsize=16, fontname='Helvetica')
    ax.set_ylabel("Recall (%)", fontsize=20)
    ax.set_xlabel("Precision (%)", fontsize=20)
    ax.annotate(f"n = {len(res)}", (0.75, 0.03), fontsize=16, fontname='Helvetica')
    zero_chrom_text = '\n'.join(zero_chrom)
    ax.annotate(f"{zero_chrom_text}", (0.03, 0.03), fontsize=10)
    plt.savefig(f"{prefix}.precision_recall.jpg", dpi=600, bbox_inches='tight')
    plt.savefig(f"{prefix}.precision_recall.pdf", dpi=600, bbox_inches='tight')

if __name__ == "__main__":
    main(sys.argv[1:])