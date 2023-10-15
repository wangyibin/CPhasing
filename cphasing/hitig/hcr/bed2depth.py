#!/usr/bin/env python3


import argparse
import sys
import os
import math
import collections
from scipy import signal
import numpy as np


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
    # draw picture
    width = 8
    height = 6
    filepath = outPre + "_dist.pdf"
    plt.figure(num=None, figsize=(width, height))
    #xxx = xxx[:50]
    #yyy = yyy[:50]
    plt.plot(xxx, yyy, '*',label='original values')
    plt.plot(xxx, yvals, 'r',label='polyfit values')
    colors = ['b', 'g', 'c']
    for peaki in [0, 1]:
        peak=peak_ind[peaki]
        plt.text(peak, 0, str(peak), fontsize = 5, color = colors[peaki])
        plt.axvline(x=peak, linewidth=1, color = colors[peaki])
    plt.savefig(filepath)
    plt.show()
    #plt.plot(x, hists[xm:xM], label = "l", color="blue")
    return peak_ind

## return region
def filter_depth(depthDic, peak_ind):
    # print(peak_ind)
    filteredRegion = {}
    minpeak, maxpeak = min(peak_ind), max(peak_ind)
    for ctg in depthDic:
        for bin in depthDic[ctg]:
            if depthDic[ctg][bin] >= minpeak and depthDic[ctg][bin] <= maxpeak:
                if ctg not in filteredRegion:
                    filteredRegion[ctg] = []
                filteredRegion[ctg].append(bin)
    return filteredRegion

def output(filteredRegion, outPre):
    with open("{}.hcr_depth.bed".format(outPre), 'w') as fout:
        for ctg in filteredRegion:
            for bin in filteredRegion[ctg]:
                fout.write("{}\t{}\t{}\n".format(ctg, *bin))
    return "{}.hcr_depth.bed".format(outPre)

def workflow(depthFile, win, Max, outPre):
    depthDic = read_depth(depthFile, Max)
    peak_id = get_wave_vally(depthDic, outPre)
    filterRegion = filter_depth(depthDic, peak_id)
    
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
