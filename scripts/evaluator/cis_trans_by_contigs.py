#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
calculate the cis/trans ratio according by contigs
"""

import argparse
import logging
import os
import os.path as op
import sys

import cooler 
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('cool', 
            help='Path to cool file')
    pReq.add_argument('contig_pairs',
            help="the pairs of contigs")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    cool = cooler.Cooler(args.cool)
    matrix = cool.matrix(balance=False, sparse=True)
    contig_paris = [i.strip().split()[:2] for i in open(args.contig_pairs) if i.strip()]
    
    data = []
    total_cis = 0
    total_trans = 0
    for pair in contig_paris:
        if pair[0] > pair[1]:
            continue
        cis1 = matrix.fetch(pair[0]).sum()
        cis2 = matrix.fetch(pair[1]).sum()
        trans = matrix.fetch(*pair).sum()
        total_cis += cis1 
        total_cis += cis2
        total_trans += trans
        if trans == 0:
            continue
        cis_trans_1 = cis1/trans
        cis_trans_2 = cis2/trans
        cis_trans = ((cis1 + cis2)/2)/trans
        print(pair[0], pair[1], cis_trans_1, cis_trans_2, cis_trans, sep='\t', file=sys.stdout)
        

        data.append(cis_trans_1)
        data.append(cis_trans_2)
    
    print(total_cis, total_trans, total_cis/total_trans, file=sys.stderr)
    fig, ax = plt.subplots(figsize=(5.5, 5))
    boxprops = dict(color='black', linewidth=1.5)
    medianprops=dict(color='black', linewidth=2.5)
    whiskerprops = dict(linestyle='--')

    bplot = ax.boxplot(data, 
                       showfliers=False, 
                       patch_artist=True, 
                       notch=True, 
                       widths=0.35,
                       medianprops=medianprops,
                       whiskerprops=whiskerprops,
                       boxprops=boxprops)
    for patch, color in zip(bplot['boxes'], ['#a83836', '#df8384', '#8dc0ed']):
        patch.set_facecolor(color)
    
    plt.savefig('boxplot.png', dpi=600, bbox_inches='tight')

        
if __name__ == "__main__":
    main(sys.argv[1:])