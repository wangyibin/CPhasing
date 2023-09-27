#!/usr/bin/env python3

import argparse
import logging
import os
import os.path as op
import sys

import math
import collections


def read_paf(paf):
    pafDic = {}
    with open(paf, 'r') as fin:
        for line in fin:
            tmpLst = line.rstrip().split('\t')
            qn = tmpLst[0]
            pafDic[qn] = tmpLst[1:]
    return pafDic

def read_fasta(faFile):
    print("Reading genome...", file=sys.stderr)
    fastaDic = {}
    with open(faFile,'r') as IN:
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

def build_windic(genomeSize, win):
    step = int(win/5)
    binDic = collections.OrderedDict()
    windic = collections.OrderedDict()
    for ctg in genomeSize:
        binDic[ctg] = {}
        windic[ctg] = {}
        ctgSize = genomeSize[ctg]
        binLst = list(range(0, ctgSize - win, step))
        binLst = [(i, i+win) for i in binLst]
        if binLst[-1][0] != (ctgSize - win):
            binLst.append((binLst[-1][0] + 1000, ctgSize))
        for bin in binLst:
            binDic[ctg][bin] = []
            windic[ctg][bin] = 0
    return binDic, windic

def bin_reads(pafDic, winDic, win):
    """
    winDic: (0, 5000) : [(reads1, 0, 1000) , (reads2, 1, 1001)]
    pafDic[qn] : ([qs,qe,s,tn,ts,te,ql,tl,al1,al2,AS]) # ql, qs, qe, s, tn, tl, ts, te, al1, al2, mapq, AS
    """
    for qn in pafDic:
        qInfo = pafDic[qn]
        tn, tl, ts, te = qInfo[4:8]
        ts, te = int(ts), int(te)
        if tn not in winDic:
            continue
        blst = list(winDic[tn].keys())
        maxBIn = len(blst)
        stratIdx = math.ceil(float(ts - 5000)/float(1000)) # 向上取整
        endIdx = math.floor(float(te)/float(1000)) + 1  # 向下取整
        stratIdx = stratIdx if stratIdx >= 0 else 0
        endIdx = endIdx if endIdx <= maxBIn else maxBIn
        #idx = math.floor(int(ts)/win)
        #blst = list(winDic[tn].keys())
        for bini in blst[stratIdx: endIdx]:
            overlap = min(int(bini[1]), te) - max(int(bini[0]), ts)
            winDic[tn][bini] += overlap
    return winDic

def cal_depth(winDic):
    """
    winDic: tn: (0, 5000) : [(reads1, 0, 1000) , (reads2, 1, 1001)]
    pafDic[qn] : ([qs,qe,s,tn,ts,te,ql,tl,al1,al2,AS])
    """
    depthDic = {}
    for ctg in winDic:
        depthDic[ctg] = {}
        for bin in winDic[ctg]:
            #readsLst = winDic[bin]
            winSize = int(bin[1]) - int(bin[0])
            overlapValue = winDic[ctg][bin]
            #overlapLst = list(map(lambda x :get_overlap(x, bin, pafDic), readsLst))
            #depth = sum(overlapLst)/winSize
            depth = overlapValue/winSize
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

def output(ouPre, depthDIc):
    with open(ouPre + ".depth", 'w') as fout:
        for ctg in depthDIc:
            for bini in depthDIc[ctg]:
                fout.write("{}\t{}\t{}\t{:.3f}\n".format(ctg, bini[0], bini[1], depthDIc[ctg][bini]))
    return ouPre + ".depth"

def workflow(paf, fastaFile, win, outPre):
    ## 
    win = int(win)
    pafDic = read_paf(paf)
    genomeSize = read_fasta(fastaFile)
    binDic, winDic = build_windic(genomeSize, win)
    winDic = bin_reads(pafDic, winDic, win)
    depthDic = cal_depth(winDic)

    return output(outPre, depthDic)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is the script for filter genome region.")
    parser.add_argument('-p', '--paf', required=True,
                        help='<filepath>  the corrected UL-ONT reads/Hifi reads mapping result, must be .paf format.')
    parser.add_argument('-f', '--fasta', required=True,
                        help='<filepath>  break-point of chimeric contigs.')
    parser.add_argument('-w', '--win', default=5000,
                        help='<int> window size when calculating depth.')
    parser.add_argument('-o', '--output', default="output",
                        help='<str> output file prefix, default is output')
    args = parser.parse_args()
    # paf, breakpointFile, fasta, win, outPre
    workflow(args.paf, args.fasta, args.win, args.output)
