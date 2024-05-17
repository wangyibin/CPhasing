#!/usr/bin/env python3


import argparse
import sys
import os
import math
import collections
from scipy import signal
import numpy as np
import pandas as pd
import pyranges as pr

"""
workflow:
(1) calculate depth distribution
(2) decision an cutoff to filter reads
"""

def read_depth(depthFile, Max):
    depthDic = {}
    Max = int(Max)
    with open(depthFile, 'r') as fin:
        for line in fin:
            ctg, bins, bine, depth = line.rstrip().split('\t')
            depth =Max if float(depth) >= float(Max) else float(depth)
            [bins, bine, depth] = list(map(int, [bins, bine, float(depth)]))
            if ctg not in depthDic:
                depthDic[ctg] = {}
            depthDic[ctg][(bins, bine)] = depth
    return depthDic


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

## return  peak_id
def get_wave_vally(depthDic, outPre):
    ## convert depth dict to depth hist
    depthHist = {}
    for ctg in depthDic:
        for bin in depthDic[ctg]:
            depth = int(depthDic[ctg][bin])
            if depth not in depthHist:
                depthHist[depth] = 0
            depthHist[depth] += 1
    ## find wave vally
    xxxyyy = list(depthHist.items())
    xxxyyy.sort(key = lambda x:x[0])
    xxx = [x[0] for x in xxxyyy]
    yyy = [y[1] for y in xxxyyy]
    xxx = xxx[:100]
    yyy = yyy[:100]
    negYYY = [-y for y in yyy]
    xxx = np.array(xxx)
    negYYY = np.array(negYYY)
    ###
    z1 = np.polyfit(xxx, negYYY, 7)
    p1 = np.poly1d(z1) 
    # print(p1)
    yvals=p1(xxx) 
    # calculate peak
    peak_ind = signal.find_peaks(yvals, distance=10)
    peak_ind = peak_ind[0]
    trough_ind = signal.find_peaks(-yvals, distance=10)
    
    trough_ind = trough_ind[0]
    trough = trough_ind[np.argmin(np.array(yvals)[trough_ind])]


    # draw picture
    width = 7
    height = 6
    filepath = outPre + "_dist.pdf"
    filepath_png = outPre + "_dist.png"
    plt.figure(num=None, figsize=(width, height))
    plt.rcParams['font.family'] = 'Arial'
    #xxx = xxx[:50]
    #yyy = yyy[:50]
    plt.plot(xxx, yyy, '*', label='Original', color="#209093")
    plt.plot(xxx, yvals, 'r', label='Polyfit', color='#cb6e7f')
    plt.legend()
    colors = ['#b02418', '#253761', 'c']
    ax = plt.gca()
    
    for peaki in [0, 1]:
        peak=peak_ind[peaki]
        plt.text(peak, ax.get_ylim()[1] / 4, str(peak), 
                 fontsize = 10, color = colors[peaki],
                )
        plt.axvline(x=peak, linewidth=1, color = colors[peaki], linestyle='--')
    
    plt.text(trough, ax.get_ylim()[0] / 4 * 3 , str(trough),
             fontsize=10, color='k', )
    plt.axvline(x=trough, linewidth=1, color='#bcbcbc', linestyle='--')

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel("Counts", fontsize=14)
    plt.xlabel("Depth", fontsize=14)
    plt.savefig(filepath, bbox_inches='tight', dpi=600)
    plt.savefig(filepath_png, bbox_inches='tight', dpi=600)

    plt.show()
    #plt.plot(x, hists[xm:xM], label = "l", color="blue")
    return peak_ind, trough

## return region
def filter_depth(depthDic, peak_ind):
    # print(peak_ind)
    filteredRegion = {}
    low_coverage_region = {}
    high_coverage_region = {}
    minpeak, maxpeak = min(peak_ind), max(peak_ind)
    for ctg in depthDic:
        for bin in depthDic[ctg]:
            depth = depthDic[ctg][bin]
            if depth >= minpeak and depth <= maxpeak:
                if ctg not in filteredRegion:
                    filteredRegion[ctg] = []
                filteredRegion[ctg].append(bin)
            elif depth < minpeak:
                if ctg not in low_coverage_region:
                    low_coverage_region[ctg] = []
                low_coverage_region[ctg].append((bin, depth))
            elif depth > maxpeak:
                if ctg not in high_coverage_region:
                    high_coverage_region[ctg] = []
                high_coverage_region[ctg].append((bin, depth))

    return low_coverage_region, filteredRegion, high_coverage_region

def output(filteredRegion, outPre):
    with open("{}.hcr_depth.bed".format(outPre), 'w') as fout:
        for ctg in filteredRegion:
            for bin in filteredRegion[ctg]:
                fout.write("{}\t{}\t{}\n".format(ctg, *bin))
    return "{}.hcr_depth.bed".format(outPre)

def output_coverage_region(coverage_region, outPre, low_or_high):
    with open("{}.{}.bed".format(outPre, low_or_high), 'w') as fout:
        for ctg in coverage_region:
            for bin, depth in coverage_region[ctg]:
                fout.write("{}\t{}\t{}\t{}\n".format(ctg, bin[0], bin[1], depth))

    return "{}.{}.bed".format(outPre, low_or_high)

def workflow(depthFile, win, Max, outPre, 
             contigsizes=None, junk_coverage=0.5, collapsed_coverage=0.1):
    depthDic = read_depth(depthFile, Max)
    peak_ind, trough = get_wave_vally(depthDic, outPre)
    low_coverage_region, filterRegion, high_coverage_region = filter_depth(depthDic, peak_ind)

    low_coverage_bed = output_coverage_region(low_coverage_region, outPre, low_or_high="low_coverage")
    high_coverage_bed = output_coverage_region(high_coverage_region, outPre, low_or_high="high_coverage")
    cmd = ['bedtools', 'merge', '-i', low_coverage_bed, '-c', '4', '-o', 'mean', ">", low_coverage_bed.replace(".bed", ".merge.bed")]
    os.system(" ".join(cmd))
    cmd = ['bedtools', 'merge', '-i', high_coverage_bed, '-c', '4', '-o', 'mean', ">", high_coverage_bed.replace(".bed", ".merge.bed")]
    os.system(" ".join(cmd))

    if contigsizes:
        low_coverage_df = pd.read_csv(low_coverage_bed.replace(".bed", ".merge.bed"), sep='\t', header=None, index_col=None, names=['Chromosome', 'Start', 'End', 'Depth'])
        low_coverage_df = low_coverage_df.eval('size=End-Start').drop(['Start', 'End'], axis=1)
        low_coverage_df = low_coverage_df.groupby('Chromosome', as_index=False).agg({"size": 'sum', 'Depth': 'mean'})
        high_coverage_df = pd.read_csv(high_coverage_bed.replace(".bed", ".merge.bed"), sep='\t', header=None, index_col=None, names=['Chromosome', 'Start', 'End', 'Depth'])
        high_coverage_df = high_coverage_df.eval('size=End-Start').drop(['Start', 'End'], axis=1)
        high_coverage_df = high_coverage_df.groupby('Chromosome', as_index=False).agg({"size": 'sum', 'Depth': 'mean'})

        low_coverage_df['length'] = low_coverage_df['Chromosome'].map(contigsizes.get)
        high_coverage_df['length'] = high_coverage_df['Chromosome'].map(contigsizes.get)
        low_coverage_df = low_coverage_df.eval("Coverage = size / length").drop(['size', 'length'], axis=1)
        high_coverage_df = high_coverage_df.eval("Coverage = size / length").drop(['size', 'length'], axis=1)

        low_coverage_df = low_coverage_df[low_coverage_df['Coverage'] > junk_coverage]
        if not low_coverage_df.empty:
            low_coverage_df['CN'] = (low_coverage_df['Depth'] / trough)
            low_coverage_df.to_csv(f'{outPre}.junk.list', sep='\t', header=None, index=None)
        high_coverage_df = high_coverage_df[high_coverage_df['Coverage'] > collapsed_coverage]
        if not high_coverage_df.empty:
            high_coverage_df['CN'] = (high_coverage_df['Depth'] / trough).astype(int)
            high_coverage_df.to_csv(f'{outPre}.collapsed.list', sep='\t', header=None, index=None)
    
    return output(filterRegion, outPre)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is the script for filter genome region.")
    parser.add_argument('-d', '--depthFile', default=None, required=True,
                        help='<filepath>  UL-ONT reads/Hifi depth file, 4colum, Chr start end depth.')
    parser.add_argument('-w', '--win', default=5000,
                        help='<int> window size when calculating depth.')
    parser.add_argument('-M', '--max', default=100,
                        help='<int> maximum depth.')
    parser.add_argument('-o', '--output', default="output",
                        help='<str> output file prefix, default is output')
    args = parser.parse_args()
    workflow(args.depthFile, args.win, args.max, args.output)
