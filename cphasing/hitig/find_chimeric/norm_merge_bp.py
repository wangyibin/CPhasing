#!/usr/bin/env python3

import logging
import math
import argparse
import sys

"""
workflow:
(1) read LIS correced paf & read correct break point
(2) disable paf which walk though break-point
(3) calculate depth
(4) decision an cutoff
"""


logger = logging.getLogger(__name__)

"""
read LIS paf
pafDic[originQn][qi].append([qs,qe,s,tn,ts,te,ql,tl,al1,al2,AS])
"""
## return pafDic

def read_mergeBP(mergeBP):
    bpDic = {}
    with open(mergeBP, 'r') as fin:
        for line in fin:
            tmpLst = line.rstrip().split('\t')
            ctg, bi1, bi2, count, reads = tmpLst
            if ctg not in bpDic:
                bpDic[ctg] = {}
            bini = (int(bi1), int(bi2))
            bpDic[ctg][bini] = int(count)
    return bpDic
          
def read_depth(depthFile):
    depthDic = {}
    with open(depthFile, 'r') as fin:
        for line in fin:
            tmpLst = line.rstrip().split('\t')
            ctg, bi1, bi2, depth = tmpLst
            if ctg not in depthDic:
                depthDic[ctg] = {}
            bini = (int(bi1), int(bi2))
            depthDic[ctg][bini] = float(depth)
    return depthDic

def find_region_depth(subDepthDic, subBpRegion, win, minDepth):
    """
    subDepthDic[(1,2)] = 19
    subBpRegion = (3,4)
    """
    regionDepth = {}
    [ts, te] = list(map(int,subBpRegion))
    blst = list(subDepthDic.keys())
    maxBIn = len(blst)
    step = int(float(win)/5)
    stratIdx = math.ceil(float(ts - win)/float(step)) # 向上取整
    endIdx = math.floor(float(te)/float(step)) + 1  # 向下取整
    stratIdx = stratIdx if stratIdx >= 0 else 0
    endIdx = endIdx if endIdx <= maxBIn else maxBIn
    for bini in blst[stratIdx: endIdx]:
        if minDepth <= subDepthDic[bini]:
            regionDepth[bini] = subDepthDic[bini]
    return regionDepth

## return depthDic
            

def normaliz_split_alignment(depthDic, bpDic, win, minDepth, cutoff):
    normBpDic = {}
    minDepth = float(minDepth)
    win = int(win)
    cutoff = float(cutoff)
    for ctg in bpDic:
   
        if ctg not in normBpDic:
            normBpDic[ctg] = {}
        subDepth = depthDic[ctg]
        for bini in bpDic[ctg]:
            regionDepth = find_region_depth(subDepth, bini, win, minDepth)
            # if not regionDepth:
            #     logger.info(f"Skip {ctg}, which not found depth.")
            #     continue
            bpCount = bpDic[ctg][bini]
            try: 
                avaDepth = float(sum(list(regionDepth.values())))/float(len(list(regionDepth.values())))
            except:
                # print(regionDepth)
                # print(ctg, bini)
                step = int(win/5)
                [ts, te] = list(map(int,bini))
                blst = list(subDepth.keys())
                
                maxBIn = len(blst)
                step = int(float(win)/5)
                stratIdx = math.ceil(float(ts - win)/float(step)) # 向上取整
                endIdx = math.floor(float(te)/float(step)) + 1  # 向下取整
                stratIdx = stratIdx if stratIdx >= 0 else 0
                endIdx = endIdx if endIdx <= maxBIn else maxBIn
                # print(stratIdx, endIdx,step, win, len(blst))
                # sys.exit()
            normBpRatio = float(bpCount)/avaDepth
            chimericType = "chimeric" if normBpRatio >= cutoff else "Non-chimeric"
            normBpDic[ctg][bini] = (chimericType, bpCount, avaDepth, normBpRatio)

    return normBpDic


def outputFile(normBpDic, outPre):
    chimericDic = {}
    ## output norm file
    with open(outPre + ".norm.result", 'w') as fout0:
        for ctg in normBpDic:
            for bini in normBpDic[ctg]:
                fout0.write("{}\t{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\n".format(ctg, bini[0], bini[1], *normBpDic[ctg][bini]))
    #with open()
            #if ctg not in chimericDic:
            #    chimericDic[ctg] = []
                if normBpDic[ctg][bini][0] == "chimeric":
                    if ctg not in chimericDic:
                        chimericDic[ctg] = []
                    bpPos = int(float(sum(bini))/float(2))
                    chimericDic[ctg].append(bpPos)
    with open(outPre + ".breakPos.txt", 'w') as fout1:
        for ctg in chimericDic:
            posLst = list(map(str, chimericDic[ctg]))
            fout1.write("{}\t{}\n".format(ctg, ",".join(posLst)))

    return outPre + ".breakPos.txt"


def workflow(mergeBP, depthFile, win, minDepth, cutoff, outPre):
    win = int(win)
    bpDic = read_mergeBP(mergeBP)
    depthDic = read_depth(depthFile)
    # for ctg in depthDic:
    #     print(len(list(depthDic[ctg].keys())))
 
    normBpDic = normaliz_split_alignment(depthDic, bpDic, win, minDepth, cutoff)
 
    #mergeBpDic, bpDic = check_breakpoint(normSaDic, cutoff, edges)
    return outputFile(normBpDic, outPre)

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is the script for filter genome region.")
    parser.add_argument('-b', '--breakpoint', required=True,
                        help='<filepath>  the merged breakpoint file.')
    parser.add_argument('-d', '--depth', required=True,
                        help='<filepath>  depth file.')
    parser.add_argument('-min', '--minDepth', default=3,
                        help='<int>  minimum depth of windows.')
    parser.add_argument('-c', '--cutoff', default=0.5,
                        help='<float>  cutoff of identifying chimeric contigs, which equal (count of splited alignment reads)/(avarage depth in chimeric region).')
    parser.add_argument('-w', '--win', default=5000,
                        help='<int> window size when calculating depth.')
    parser.add_argument('-o', '--output', default="output",
                        help='<str> output file prefix, default is output')
    args = parser.parse_args()
    # paf, breakpointFile, fasta, win, outPre
    workflow(args.breakpoint, args.depth, args.win, args.minDepth, args.cutoff, args.output)
