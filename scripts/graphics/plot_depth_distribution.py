#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot the distribution of depth
"""

import argparse
import logging
import os
import os.path as op
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, MaxNLocator
import seaborn as sns

import pandas as pd
import numpy as np

from scipy.signal import find_peaks
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bedgraph', 
            help='')
    pOpt.add_argument('-o', '--output', type=str,
                default=None, help='output file [default: suffix of input file]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    depth_df = pd.read_csv(args.bedgraph, sep='\t', header=None, index_col=None)
    depth_df.columns = ['chrom', 'start', 'end', 'count']
    depth_df['count'] = depth_df['count'].map(np.round)
    depthHist = depth_df.groupby(['count'])['chrom'].count()
    depthHist = depthHist[depthHist.index > 0]
    max_values = round(depthHist.argmax()  * 4)
    depthHist = depthHist[depthHist.index < max_values]
    depthHist = depthHist.to_dict()
    ## find wave vally
    xxxyyy = list(depthHist.items())
    xxxyyy.sort(key = lambda x:x[0])

    xxx = [x[0] for x in xxxyyy]
    yyy = [y[1] for y in xxxyyy]
    xxx = xxx[:max_values]
    yyy = yyy[:max_values]

    peak_ind = find_peaks(yyy, distance=10)

    peak_value = xxx[np.argmax(yyy)]

    fig, ax = plt.subplots(figsize=(5.7, 5))
    plt.rcParams['font.family'] = 'Arial'
    plt.plot(xxx, yyy, '-', label='Depth', color="#209093")
    plt.ylim(0)
    plt.xlim(0)
    plt.axvline(x=peak_value, linewidth=1, color='#b02418', linestyle='--')
    plt.text(peak_value, ax.get_ylim()[1] * 0.98 , f"{str(peak_value)}",
             fontsize=14, color='k', )
    plt.axvline(x=peak_value / 2, linewidth=1, color='#253761', linestyle='--')
    plt.text(peak_value / 2, ax.get_ylim()[1] * 0.98 , f"{peak_value / 2:.1f}",
             fontsize=14, color='k', )
    plt.axvline(x=peak_value * 2, linewidth=1, color='#253761', linestyle='--')
    plt.text(peak_value * 2, ax.get_ylim()[1] * 0.98 , f"{str(peak_value * 2)}",
             fontsize=14, color='k', )
    # data = np.column_stack((xxx, yyy))

    # gmm = GaussianMixture(n_components=3, covariance_type='full', tol=1e-6, 
    #                       max_iter=10000, random_state=0).fit(data)
    # labels = gmm.predict(data)

    # centers = gmm.means_

    # plt.scatter(xxx, yyy, c=labels, cmap='viridis', label='Data Points', s=3)

    # x = np.linspace(min(xxx), max(xxx), 100)
    # for i in range(gmm.n_components):
    #     mean = gmm.means_[i]
    #     cov = gmm.covariances_[i]
    #     y = np.exp(-0.5 * ((x - mean[0]) ** 2) / cov[0, 0]) / np.sqrt(2 * np.pi * cov[0, 0])
    #     plt.plot(x, y * max(yyy), label=f'Gaussian {i+1}')

    
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    plt.tick_params(which='both', width=1.5, length=5)
    sns.despine()

    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ylabel("Frequency", fontsize=20)
    plt.xlabel("Coverage", fontsize=20)

    output = args.output if args.output else args.bedgraph
    plt.savefig(f'{output}.hist.plot.png', dpi=600, bbox_inches='tight')
    plt.savefig(f'{output}.hist.plot.pdf', dpi=600, bbox_inches='tight')
    
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    plt.tick_params(which='both', width=1.5, length=5)
    sns.despine()
    


if __name__ == "__main__":
    main(sys.argv[1:])