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


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cphasing.utilities import read_chrom_sizes

"""
workflow:
(1) calculate depth distribution
(2) decision an cutoff to filter reads
"""

def read_depth(depthFile):
    # depthDic = {}
    # Max = int(Max)
    # with open(depthFile, 'r') as fin:
    #     for line in fin:
    #         ctg, bins, bine, depth = line.rstrip().split('\t')
    #         depth = Max if float(depth) >= float(Max) else float(depth)
    #         [bins, bine, depth] = list(map(int, [bins, bine, float(depth)]))
    #         if ctg not in depthDic:
    #             depthDic[ctg] = {}
    #         depthDic[ctg][(bins, bine)] = depth

    depth_df = pd.read_csv(depthFile, sep='\t', index_col=None, header=None, names=['chrom', 'start', 'end', 'count'])
   
    return depth_df


## return  peak_id
def get_wave_vally(depth_df, outPre, Max=None):
    ## convert depth dict to depth hist
    # depthHist = {}
    # for ctg in depthDic:
    #     for bin in depthDic[ctg]:
    #         depth = int(depthDic[ctg][bin])
    #         if depth not in depthHist:
    #             depthHist[depth] = 0
    #         depthHist[depth] += 1
    depth_df['count'] = depth_df['count'].astype(int)
    depthHist = depth_df.groupby(['count'])['chrom'].count()
    max_values = depthHist.argmax()  * 2.5 if not Max else Max

    depthHist = depthHist[depthHist.index < max_values]
    depthHist = depthHist.to_dict()
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
    z1 = np.polyfit(xxx, negYYY, 10)
    p1 = np.poly1d(z1) 
    # print(p1)
    yvals=p1(xxx) 
    # calculate peak
    peak_ind = signal.find_peaks(yvals, distance=10)
    peak_ind = peak_ind[0]
    trough_ind = signal.find_peaks(-yvals, distance=10)
    
    trough_ind = trough_ind[0]
    trough = trough_ind[np.argmin(np.array(yvals)[trough_ind])]
    peak_ind = [peak_ind[peak_ind < trough][0], peak_ind[peak_ind > trough][0]]
    
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

    for i, peaki in enumerate([0, 1]):
        peak=peak_ind[peaki]
        plt.text(peak, ax.get_ylim()[1] / 4, str(peak), 
                 fontsize = 10, color = colors[i],
                )
        plt.axvline(x=peak, linewidth=1, color = colors[i], linestyle='--')
    
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
    return [peak_ind[0], peak_ind[1]], trough

## return region
def filter_depth(depth_df, peak_ind):
    # print(peak_ind)
    filteredRegion = {}
    minpeak, maxpeak = min(peak_ind), max(peak_ind)
    
    filteredRegion = depth_df.query('count >= @minpeak & count <= @maxpeak')

    return filteredRegion

def output(filteredRegion, outPre):
    filteredRegion.to_csv("{}.hcr_depth.bed".format(outPre), header=None, sep='\t', index=None)

    return "{}.hcr_depth.bed".format(outPre)

def output_coverage_region(coverage_region, outPre, low_or_high):
    with open("{}.{}.bed".format(outPre, low_or_high), 'w') as fout:
        for ctg in coverage_region:
            for bin, depth in coverage_region[ctg]:
                fout.write("{}\t{}\t{}\t{}\n".format(ctg, bin[0], bin[1], depth))

    return "{}.{}.bed".format(outPre, low_or_high)

def workflow(depthFile, win,  outPre,  Max=None,
             contigsizes=None, junk_coverage=0.5, collapsed_coverage=0.1):
    depth_df = read_depth(depthFile)
    
    peak_ind, trough = get_wave_vally(depth_df, outPre, Max=Max)

    filterRegion= filter_depth(depth_df, peak_ind)
    minpeak, maxpeak = min(peak_ind), max(peak_ind)

    contig_depth = depth_df.groupby('chrom')['count'].mean()
    high_coverage_df = contig_depth[contig_depth > maxpeak].to_frame()
    high_coverage_df['CN'] = (high_coverage_df / trough)['count'].map(np.round)
    high_coverage_df.to_csv(f'{outPre}.collapsed.contigs', sep='\t', header=None, index=True)

    low_coverage_df = contig_depth[contig_depth < minpeak].to_frame()
    low_coverage_df['CN'] = (low_coverage_df / trough)['count'].map(np.round)
    low_coverage_df.to_csv(f'{outPre}.lowcoverage.contigs', sep='\t', header=None, index=True)


    # low_coverage_bed = output_coverage_region(low_coverage_region, outPre, low_or_high="low_coverage")
    # high_coverage_bed = output_coverage_region(high_coverage_region, outPre, low_or_high="high_coverage")
    # cmd = ['bedtools', 'merge', '-i', low_coverage_bed, '-c', '4', '-o', 'mean', ">", low_coverage_bed.replace(".bed", ".merge.bed")]
    # os.system(" ".join(cmd))
    # cmd = ['bedtools', 'merge', '-i', high_coverage_bed, '-c', '4', '-o', 'mean', ">", high_coverage_bed.replace(".bed", ".merge.bed")]
    # os.system(" ".join(cmd))
    

    # if contigsizes:
    #     low_coverage_df = pd.read_csv(low_coverage_bed.replace(".bed", ".merge.bed"), sep='\t', header=None, index_col=None, names=['Chromosome', 'Start', 'End', 'Depth'])
    #     low_coverage_df = low_coverage_df.eval('size=End-Start').drop(['Start', 'End'], axis=1)
    #     low_coverage_df = low_coverage_df.groupby('Chromosome', as_index=False).agg({"size": 'sum', 'Depth': 'mean'})
    #     high_coverage_df = pd.read_csv(high_coverage_bed.replace(".bed", ".merge.bed"), sep='\t', header=None, index_col=None, names=['Chromosome', 'Start', 'End', 'Depth'])
    #     high_coverage_df = high_coverage_df.eval('size=End-Start').drop(['Start', 'End'], axis=1)
    #     high_coverage_df = high_coverage_df.groupby('Chromosome', as_index=False).agg({"size": 'sum', 'Depth': 'mean'})

    #     low_coverage_df['length'] = low_coverage_df['Chromosome'].map(contigsizes.get)
    #     high_coverage_df['length'] = high_coverage_df['Chromosome'].map(contigsizes.get)
    #     low_coverage_df = low_coverage_df.eval("Coverage = size / length").drop(['size', 'length'], axis=1)
    #     high_coverage_df = high_coverage_df.eval("Coverage = size / length").drop(['size', 'length'], axis=1)
    #     print(high_coverage_df)
    #     low_coverage_df = low_coverage_df[low_coverage_df['Coverage'] > junk_coverage]
    #     if not low_coverage_df.empty:
    #         low_coverage_df['CN'] = (low_coverage_df['Depth'] / trough)
    #         low_coverage_df.to_csv(f'{outPre}.lowcoverage.contigs', sep='\t', header=None, index=None)
    #     high_coverage_df = high_coverage_df[high_coverage_df['Coverage'] > collapsed_coverage]
    #     if not high_coverage_df.empty:
    #         high_coverage_df['CN'] = (high_coverage_df['Depth'] / trough).astype(int)
    #         high_coverage_df.to_csv(f'{outPre}.collapsed.contigs', sep='\t', header=None, index=None)
    
    return output(filterRegion, outPre)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is the script for filter genome region.")
    parser.add_argument('-d', '--depthFile', default=None, required=True,
                        help='<filepath>  UL-ONT reads/Hifi depth file, 4colum, Chr start end depth.')
    parser.add_argument('-c', '--contigsizes', help='contigsizes', default=None)
    parser.add_argument('-w', '--win', default=5000,
                        help='<int> window size when calculating depth.')
    parser.add_argument('-M', '--max', default=None,
                        help='<int> maximum depth.')
    parser.add_argument('-o', '--output', default="output",
                        help='<str> output file prefix, default is output')
    args = parser.parse_args()
    contigsizes = read_chrom_sizes(args.contigsizes)
    contigsizes = contigsizes.to_dict()['length']

    workflow(args.depthFile, args.win, args.output, int(args.max) if args.max else None, contigsizes)
