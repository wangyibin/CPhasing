#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot dotplot
"""
import argparse
import logging
import os
import os.path as op
import sys


import matplotlib.pyplot as plt
from collections import defaultdict
import scipy.stats
import sys
import numpy as np
from matplotlib.ticker import MaxNLocator


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('fai', 
            help='')
    pReq.add_argument('agp')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    ref = args.fai 
    qurry = args.agp
    winsize=10000

    ref_dic=defaultdict(list)
    qurry_dic=defaultdict(list)
    copy_num_dic=defaultdict(list)
    chromsizes=defaultdict(list)
    with open(ref) as ref_file:
        for line in ref_file:
            data=line.split('\t')[0].split('.')
            ref_dic[data[0]].append(data[1])
            copy_num_dic[data[0]].append([data[1],(int(line.split('\t')[1])//winsize)+1])
            chromsizes[data[0]].append(int(line.split('\t')[1]))
    for chrom in chromsizes:
        chromsizes[chrom]=sum(chromsizes[chrom])

    for key in copy_num_dic:
        copy_num_dic[key]=dict(copy_num_dic[key])

    with open(qurry) as qurry_file:
        for line in qurry_file:
            data=line.split('\t')[5].split('.')
            if len(data) == 2:
                qurry_dic[data[0]].append(data[1]+line.strip().split('\t')[8])

    copy_ref_dic=defaultdict(list)
    for key in ref_dic:
        for value in ref_dic[key]:
            for i in range(1,copy_num_dic[key][value]+1):
                copy_ref_dic[key].append(f"{value}_w{i}")

    copy_qurry_dic=defaultdict(list)
    for key in qurry_dic:
        for value in qurry_dic[key]:
            if value[-1] == '+':
                for i in range(1,copy_num_dic[key][value[:-1]]+1):
                    copy_qurry_dic[key].append(f"{value[:-1]}_w{i}")
            if value[-1] == '-':
                for i in range(copy_num_dic[key][value[:-1]],0,-1):
                    copy_qurry_dic[key].append(f"{value[:-1]}_w{i}")

    def dotplot(chr,rho, x_data, y_data,chrom_sizes):
        
        plt.rcParams['font.family'] = 'Arial'

        labels = copy_ref_dic[chr]

        # x_data = copy_ref_dic[chr]
        # y_data = copy_qurry_dic[chr]

        # x_indices = [labels.index(label) for label in x_data]
        # y_indices = [labels.index(label) for label in y_data]

        max_value = len(labels)
        num_ticks = 5
        ticks = [i * (max_value / (num_ticks - 1)) for i in range(num_ticks)]


        plt.figure(figsize=(5, 5))


        plt.scatter(x_data, y_data, s=10,color="#B3CDE3")


        # plt.xticks(ticks, fontsize=20)
        # plt.yticks(ticks, fontsize=20)
        
        ax = plt.gca()
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.tick_params(axis='x', length=5, width=2)
        #print(ax.get_xticks())
        xticklabels = ax.get_xticks()[1:] * winsize / 1e6
        
        ax.tick_params(axis='y', length=5, width=2)
        yticklabels = ax.get_yticks()[1:] * winsize / 1e6
        
        ax.set_xticklabels(xticklabels, fontsize=16)
        ax.set_yticklabels(yticklabels, fontsize=16)
        #print(yticklabels)
        ax.set_xlim(0)
        ax.set_ylim(0)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)

        if rho > 0:
            plt.text(0.05 * max_value, 1 * max_value, r'$\rho = {0:.2f}$'.format(rho), fontsize=16)
        else:
            plt.text(0.92 * max_value, 1 * max_value, r'$\rho = {0:.2f}$'.format(rho), fontsize=16, ha='right')

        plt.xlabel(f'Assembly: {(chrom_sizes/1e6):.2f} Mb',labelpad=15, fontsize=20)
        plt.ylabel(f'{chr}: {(chrom_sizes/1e6):.2f} Mb', labelpad=15, fontsize=20)

        plt.savefig(f'{chr}_dotplot.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{chr}_dotplot.pdf', dpi=300, bbox_inches='tight')
        plt.close()

    with open(f"{'.'.join(qurry.split('.')[:-1])}_rho.txt",'w') as output:

        for chr in qurry_dic:
            temp_ref = copy_ref_dic[chr]
            temp_query = copy_qurry_dic[chr]
            idx_db = dict(zip(temp_ref, range(len(temp_ref))))
            temp_query = list(map(idx_db.get, temp_query))
            rho=(scipy.stats.spearmanr(list(range(len(temp_ref))), temp_query)[0])
            output.write(f'{chr}\t{rho}\n')
            dotplot(chr,rho, list(range(len(temp_ref))), temp_query,chromsizes[chr])
            

if __name__ == "__main__":
    main(sys.argv[1:])