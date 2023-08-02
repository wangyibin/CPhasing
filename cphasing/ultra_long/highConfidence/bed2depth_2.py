#!/usr/bin/env python3
#coding=utf-8

import sys
import os
import math
import collections
from scipy import signal
import numpy as np
import argparse

"""
workflow:
(1) read LIS correced paf & read correct break point
(2) disable paf which walk though break-point
(3) calculate depth
(4) decision an cutoff
"""

"""
read LIS paf
pafDic[originQn][qi].append([qs,qe,s,tn,ts,te,ql,tl,al1,al2,AS])
"""
## return pafDic
def read_paf(paf):
    pafDic = {}
    with open(paf, 'r') as fin:
        for line in fin:
            tmpLst = line.rstrip().split('\t')
            qn = tmpLst[0]
            pafDic[qn] = tmpLst[1:]
    return pafDic

## return bpDIc
def read_bp(bpFile):
    bpDic = {}
    with open(bpFile,'r') as fin:
        for line in fin:
            ctg, bpStr = line.rstrip().split('\t')
            bpDic[ctg] = bpStr.split(',')
    return bpDic

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

## return fastaDic
def read_faSize(fastaFile):
    fastaDic = {}
    with open(fastaFile,'r') as IN:
        fastaName = IN.readline().strip()[1:]
        fa = ''
        for line in IN:
            if line.startswith('>'):
                fastaDic[fastaName] = len(fa)
                fastaName = line.strip()[1:]
                fa = ''
            else:
                fa += line.rstrip()
        fastaDic[fastaName] = len(fa)
    return fastaDic

## return winDic
def slid_win(fastaDic, win):
    winDic = {}
    for ctg in fastaDic:
        ctgSize = fastaDic[ctg]
        binLst = list(range(0, ctgSize, win))
        binLst = [[i, i+win] for i in binLst]
        binLst[-1] = (binLst[-1][0], ctgSize)
        if ctg not in winDic:
            winDic[ctg] = collections.OrderedDict()
        for bin in binLst:
            winDic[ctg][tuple(bin)] = []
    return winDic

## return winDic
def bin_reads(pafDic, winDic, win):
    """
    winDic: (0, 5000) : [(reads1, 0, 1000) , (reads2, 1, 1001)]
    pafDic[qn] : ([qs,qe,s,tn,ts,te,ql,tl,al1,al2,AS])
    """
    for qn in pafDic:
        qInfo = pafDic[qn]
        tn, ts, te = qInfo[3:6]
        idx = math.floor(int(ts)/win)
        blst = list(winDic[tn].keys())
        winDic[tn][blst[idx]].append((qn, int(ts), int(te)))
    return winDic

## return pafDic  
def disable_reads(bpDic, winDic, pafDic, win):
    for ctg in bpDic:
        bpLst = bpDic[ctg]
        binLst = winDic[ctg]
        for bp in bpLst:
            bpIdx = math.floor(int(bp)/win)
            qnLst = winDic[ctg][binLst[bpIdx]]
            qnLst.sort(key=lambda x:x[1])
            for qn, ts, te in qnLst:
                if ts >= bp:
                    continue
                elif ts < bp and bp < te:
                    pafDic.pop(qn)
                else:
                    break
    return pafDic

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

## return depthDic
def cal_depth(winDic, pafDic, win):
    """
    winDic: tn: (0, 5000) : [(reads1, 0, 1000) , (reads2, 1, 1001)]
    pafDic[qn] : ([qs,qe,s,tn,ts,te,ql,tl,al1,al2,AS])
    """
    depthDic = {}
    for ctg in winDic:
        depthDic[ctg] = {}
        for bin in winDic[ctg]:
            readsLst = winDic[ctg][bin]
            winSize = int(bin[1] - bin[0])
            overlapLst = list(map(lambda x :get_overlap(x, bin, pafDic), readsLst))
            depth = sum(overlapLst)/winSize
            depthDic[ctg][bin] = depth
    return depthDic
            
def get_overlap(readsItem, bin, pafDic):
    s,e = bin
    rn, rs, re = readsItem
    if rn in pafDic:
        overlap = min(re, e) - max(rs, s)
        return overlap
    else:
        return 0
    
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
    ## 
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
    print("this")
    z1 = np.polyfit(xxx, negYYY, 7) # 用7次多项式拟合
    p1 = np.poly1d(z1) #多项式系数
    print(p1) # 在屏幕上打印拟合多项式
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
    print(peak_ind)
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
    with open("{}.high.bed".format(outPre), 'w') as fout:
        for ctg in filteredRegion:
            for bin in filteredRegion[ctg]:
                fout.write("{}\t{}\t{}\n".format(ctg, *bin))
    
def workflow(depthFile, breakpointFile, fasta, win, Max, outPre):
    #pafDic = read_paf(paf)
    #fastaDic = read_faSize(fasta)
    #winDic = slid_win(fastaDic, win)
    #winDic = bin_reads(pafDic, winDic, win)
    depthDic = read_depth(depthFile, Max)
    #print("1")
    if bool(breakpointFile):
        bpDic = read_bp(breakpointFile)
        pafDic = disable_reads(winDic, pafDic, win)
    else:
        pass
    #depthDic = cal_depth(winDic, pafDic, win)
#    with open("{}.depth".format(outPre), 'w') as fout:
#        for ctg in depthDic:
#            for binItem in depthDic[ctg]:
#                fout.write("{}\t{}\t{}\t{}\n".format(ctg, binItem[0], binItem[1], depthDic[ctg][binItem]))
   # print("2")
    peak_id = get_wave_vally(depthDic, outPre)
  #  print("3")
    filterRegion = filter_depth(depthDic, peak_id)
  #  print("4")
    output(filterRegion, outPre)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is the script for filter genome region.")
    parser.add_argument('-d', '--depthFile', default=None,
                        help='<filepath>  UL-ONT reads/Hifi depth file, 4colum, Chr start end depth.')
    parser.add_argument('-b', '--breakpoint', default=None,
                        help='<filepath>  break-point of chimeric contigs.')
    parser.add_argument('-f', '--fasta', default=None,
                        help='<filepath> fasta file.')
    parser.add_argument('-w', '--win', default=5000,
                        help='<int> window size when calculating depth.')
    parser.add_argument('-M', '--max', default=100,
                        help='<int> maximum depth.')
    parser.add_argument('-o', '--output', default="output",
                        help='<str> output file prefix, default is output')
    args = parser.parse_args()
    # paf, breakpointFile, fasta, win, outPre
    workflow(args.depthFile, args.breakpoint, args.fasta, args.win, args.max, args.output)
