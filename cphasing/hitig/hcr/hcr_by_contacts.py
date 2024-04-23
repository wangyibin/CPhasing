#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Identify the high confidence regions of contigs by Pore-C/Hi-C contacts
"""

import argparse
import logging
import os
import os.path as op
import sys

import cooler 
import numpy as np 
import pandas as pd
import pyranges as pr 

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.ticker import MaxNLocator


logger = logging.getLogger(__name__)

def plot(data, lower_value=0.1, upper_value=1.75, output="output"):
    fig, ax = plt.subplots(figsize=(5.7, 5))
    plt.rcParams['font.family'] = 'Arial'
    data = data[data <= np.percentile(data, 98) * 1.5]
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    kdelines = sns.kdeplot(data, ax=ax, color='#253761', linewidth=2)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    formatter = plt.gca().get_yaxis().get_major_formatter()
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.gca().yaxis.get_offset_text().set_fontsize(14)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    formatter = plt.gca().get_xaxis().get_major_formatter()
    plt.gca().xaxis.set_major_formatter(formatter)
    plt.gca().xaxis.get_offset_text().set_fontsize(14)
    plt.xlim(0, np.percentile(data, 95) * 2)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("Contacts", fontsize=24)
    plt.ylabel("Density", fontsize=24)

    x = kdelines.lines[0].get_xdata()
    y = kdelines.lines[0].get_ydata()
    max_idx = np.argmax(y)
    
    
    ax.fill_between((x[max_idx] * lower_value, x[max_idx] * upper_value), 
                    0, ax.get_ylim()[1], alpha=0.5 , color='#bcbcbc')
    ax.axvline(x[max_idx] * lower_value, linestyle='--', color='k')
    ax.axvline(x[max_idx] * upper_value, linestyle='--', color='k')

    # plt.plot(x[max_idx], y[max_idx], ms=10, color='r')
    plt.savefig(f'{output}.kde.plot.png', dpi=600, bbox_inches='tight')
    plt.savefig(f'{output}.kde.plot.pdf', dpi=600, bbox_inches='tight')
    logger.info(f"Output kde plot of contacts distribution in `{output}.kde.plot.png`")

    return x[max_idx] * lower_value, x[max_idx] * upper_value

def hcr_by_contacts(cool_file, output, lower_value=0.1, upper_value=1.75 ):

    cool = cooler.Cooler(cool_file)
    binsize = cool.binsize
    bins = cool.bins()[:]
    
    contigsizes = cool.chromsizes
    matrix = cool.matrix(balance=False, sparse=True)[:]
    sum_values = np.array(matrix.sum(axis=1).T[0])
    logger.debug("Adjusting small bins value ...")
    small_bins = bins[bins['end'] - bins['start'] < cool.binsize]
    small_bins_sum_values = sum_values.T[small_bins.index]
    adjust_small_bins_sum_values = small_bins_sum_values.T / \
        ((small_bins['end'] - small_bins['start']) / binsize).values
    sum_values[:, small_bins.index] = adjust_small_bins_sum_values

    sum_values_nonzero = sum_values[sum_values > 0]

    min_value, max_value = plot(sum_values_nonzero)
    #median = np.median(sum_values)
    # max_value = np.percentile(sum_values_nonzero, percent)
    # logger.debug(f"Percent{percent} value is {max_value}")
    res = np.where((sum_values <= max_value) & (sum_values >= min_value))
    
    hcr_regions = bins.loc[res[1]]
    total_length = contigsizes.sum()
    num_hcr_regions = len(hcr_regions)
    logger.debug(f"Identify {num_hcr_regions} regions")
    hcr_length = sum(hcr_regions["end"] - hcr_regions["start"])
    logger.info(f"Identified {hcr_length/total_length:.2%} high-confidence regions")
    
    hcr_regions = hcr_regions[['chrom', 'start', 'end']]
    hcr_regions.columns = ['Chromosome', 'Start', 'End']
    hcr_regions_pr = pr.PyRanges(hcr_regions)
    hcr_regions_pr.merge().df.to_csv(output, sep='\t', index=None, header=None)
    logger.info(f"Successful output HCRs into `{output}`.")


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('cool_file', 
            help='')
    pOpt.add_argument('-p', '--percent', type=int, default=95,
            help='percentile' )
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    hcr_by_contacts(args.cool_file, args.output, args.percent)

if __name__ == "__main__":
    main(sys.argv[1:])