#/usr/bin/env python
#coding=utf-8
import sys
import os
import math
import argparse
from copy import deepcopy
import collections


def read_LISFile(LisFile, feature):
    LisDic = collections.OrderedDict()
    with open(LisFile, 'r') as fin:
        #tmpLine = fin.readline().rstrip().split('\t')
        for line in fin:
            tmpLine = line.rstrip().split('\t')
            lisType = tmpLine[2]
            if lisType == "LIS":
                qn = tmpLine[-1].split(":")[1]
                tn, tl, s, e, string = tmpLine[0], tmpLine[1], tmpLine[3], tmpLine[4], tmpLine[6]
                if qn not in LisDic:
                    LisDic[qn] = {}
                   # lisLst = []
                #lisLst.append((tn, s, e, string, tl))
                LisDic[qn][(tn, s, e, string, tl)] = []
            elif lisType == "alignment":
                qi = int(tmpLine[-1].split(":")[-1][:-1].split('_')[-1])
                LisDic[qn][(tn, s, e, string, tl)].append(qi)
            else:
                pass
    newLisDic = {}
    for qn in LisDic:
        if qn not in newLisDic:
            newLisDic[qn] = []
        for tn, s, e, string, tl in LisDic[qn]:
            qiRange = LisDic[qn][(tn, s, e, string, tl)]
            qiRange.sort()
            newLisDic[qn].append((tn, s, e, string, tl, qiRange[0], qiRange[-1]))
    return newLisDic
            
### return bpDic2
def check_break_point(LisDic, win, outPre, edge):
    """
    嵌合有几种模式(bp: break-point):
    1. SWITCH:
        ctg 1A|B:
            ctg 1A |<---| bp |---->| ctg 1B
        ctg 1B|A:
            ctg 1B |<---| bp |---->| ctg 1A
        LIS:
            1A|B + 1B|A, 1B|A + 1A|B
    2. INDEL:
        ctg 1|5A:
            ctg1 A |<---| bp |---->| ctg 5A
        LIS:
            ctg1|5A + ctg4A, ctg1|5A + ctg2A .....
    """
    win = int(win)
    edge = int(edge)
    ### 过滤掉只有1个窗口的LIS
    for qn in LisDic:
        newLisDic = []
        for lis in LisDic[qn]:
            #print(lis)
            #sys.exit()
            tn, s, e, string, tl, qi0, qi1 = lis
            s, e, tl =int(s), int(e), int(tl)
            if float(abs(s - e)) <= 1.25 * float(win):
                continue
            newLisDic.append(lis)
        LisDic[qn] = newLisDic
    ### find split alignment
    splitAlign = {}
    genomeSize = {}
    for qn in LisDic:
    #    if len(LisDic[qn]) == 2:
        if len(LisDic[qn]) >= 2:   # EDIT HERE
            for i in range(len(LisDic[qn]) - 1):
                ftn, fs, fe, fstr, ftl, fqi0, fqi1 = LisDic[qn][i]
                ctn, cs, ce, cstr, ctl, cqi0, cqi1 = LisDic[qn][i+1]
                if ftn == ctn:
                    continue
                [fs, fe, ftl, cs, ce, ctl] = list(map(int, [fs, fe, ftl, cs, ce, ctl]))
                if ftn not in genomeSize:
                    genomeSize[ftn] = ftl
                if ctn not in genomeSize:
                    genomeSize[ctn] = ctl
                """
                |--->||---->|
                |--->||<----|
                |<---||---->|
                |<---||<----|
                """
                #print(LisDic[qn][i])
                #print(LisDic[qn][i+1])
                winDis = int(1.25 * float(cqi0 - fqi1 + 1) * float(win))
                if fstr == "+" and cstr == "+":
                    #if (ftl - fe) <= win and cs <= win:
                    if (ftl - fe) + cs <= winDis:
                        continue
                    elif (ftl - fe) < edge and cs < edge:
                        continue
                    elif (ftl - fe) > win and cs > win:
                        fbp, cbp = fe, cs
                    elif (ftl - fe) > win and cs <= win:
                        fbp, cbp = fe, 0
                    else:
                        fbp, cbp = 0, cs

                elif fstr == "+" and cstr == "-":
                    if (ftl - fe) + (ctl - ce) <= winDis:
                    #if (ftl - fe) <= win and (ctl - ce) <= win:
                        continue
                    elif (ftl - fe) < edge and (ctl - ce) < edge:
                        continue
                    elif (ftl - fe) > win and (ctl - ce) > win:
                        fbp, cbp = fe, ce
                    elif (ftl - fe) > win and (ctl - ce) <= win:
                        fbp, cbp = fe, 0
                    else:
                        fbp, cbp = 0, ce
                
                elif fstr == "-" and cstr == "+":
                    if fs + cs <= winDis:
                    #if fs <= win and cs <= win:
                        continue
                    elif fs < edge and cs < edge:
                        continue
                    elif fs > win and cs > win:
                        fbp, cbp = fs, cs
                    elif fs > win and cs <= win:
                        fbp, cbp = fs, 0
                    else:
                        fbp, cbp = 0, cs

                else:
                    #if fs <= win and (ctl - ce) <= win:
                    if fs + (ctl - ce) <= winDis:
                        continue
                    elif fs < edge and (ctl - ce) < edge:
                        continue
                    elif fs > win and (ctl - ce) > win:
                        fbp, cbp = fs, ce
                    elif fs > win and (ctl - ce) <= win:
                        fbp, cbp = fs, 0
                    else:
                        fbp, cbp = 0, ce
                ## add into splitAlign
                if fbp == 0:
                    ctg1st, ctg2rd = ctn, ftn
                    bp1st, bp2rd = cbp, fbp
                elif cbp == 0:
                    ctg1st, ctg2rd = ftn, ctn
                    bp1st, bp2rd = fbp, cbp
                else:
                    ctg1st, ctg2rd = ftn, ctn
                    bp1st, bp2rd = fbp, cbp
                if ctg1st not in splitAlign:
                    splitAlign[ctg1st] = {}
                if ctg2rd not in splitAlign[ctg1st]:
                    splitAlign[ctg1st][ctg2rd] = []
                if ctg2rd not in splitAlign:
                    splitAlign[ctg2rd] = {}
                if ctg1st not in splitAlign[ctg2rd]:
                    splitAlign[ctg2rd][ctg1st] = []
                if bp2rd != 0:
                    splitAlign[ctg1st][ctg2rd].append((bp1st, bp2rd, qn))
                    splitAlign[ctg2rd][ctg1st].append((bp2rd, bp1st, qn))
                else:
                    splitAlign[ctg1st][ctg2rd].append((bp1st, bp2rd, qn))

    ### sliding window
    """
    为了确定断点,将contig拆分成5000bp,1000bp步长的窗口,分段比对数超过阈值的窗口即为嵌合的窗口
    splitAlign : [ftn][ctn] -> [(fbp, cbp), (fbp, cbp)]
    genomeSize[ctn] = ctl
    """
    win = 2*win
    step = int(win/5)
    binDic = collections.OrderedDict()
    for ctg in genomeSize:
        binDic[ctg] = {}
        ctgSize = genomeSize[ctg]
        binLst = list(range(0, ctgSize - win, step))
        binLst = [(i, i+win) for i in binLst]
        if binLst[-1][0] != (ctgSize - win):
            binLst.append((binLst[-1][0] + step, ctgSize))
        for bini in binLst:
            binDic[ctg][bini] = {}
    ### locate splitAlign
    """
    binLst = [(0, 5000), (1000, 6000), (2000,7000), (3000, 7650)]
    posLst = [423, 1234, 3456, 6500]
    x[0] = idx * step
    x[i] = idx * step + 5000
    pos >= x[0]
    pos <= x[1]
    """
    #print(splitAlign)
    #sys.exit()
    for ctg1 in splitAlign:
        for ctg2 in splitAlign[ctg1]:
            posLst = [x[0] for x in splitAlign[ctg1][ctg2]]
            # posLst.sort(key=lambda x:x[0])
            binLst = list(binDic[ctg1].keys())
            maxBIn = len(binLst)
            for posIdx in range(len(posLst)):
                pos = posLst[posIdx]
                stratIdx = math.ceil(float(pos - win)/float(step)) # 向上取整
                endIdx = math.floor(float(pos)/float(step)) + 1  # 向下取整
                stratIdx = stratIdx if stratIdx >= 0 else 0
                endIdx = endIdx if endIdx <= maxBIn else maxBIn
                #print(ctg1, stratIdx, endIdx)
                for bini in binLst[stratIdx: endIdx]:
                    pos1, pos2, qn = splitAlign[ctg1][ctg2][posIdx]
                    if ctg2 not in binDic[ctg1][bini]:
                        binDic[ctg1][bini][ctg2] = []
                    binDic[ctg1][bini][ctg2].append(qn)
                    
    ### find break point
    with open(outPre + ".splitAlign.txt",'w') as fout:
        for ctg in binDic:
            for bini in binDic[ctg]:
        #    with open("out.splitAlign.txt",'w') as fout:
                #qnLst = [x[-1] for x in binDic[ctg][bini]]
                for ctg2 in binDic[ctg][bini]:
                    qnLst = binDic[ctg][bini][ctg2]
                    fout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(ctg, bini[0], bini[1], len(qnLst), ctg2, ",".join(qnLst)))
                #fout.write("{}\t{}\t{}\t{}\t{}\n".format(ctg, bini[0], bini[1], len(binDic[ctg][bini]), ",".join(qnLst)))
    return binDic

def merge_bp(binDic, outPre, minSA):
    """
    bpDic[ctg1] = {(bini):[qnLst]}
    """
    minSA = int(minSA)
    bpDic = {}
    for ctg1 in binDic:
        for bini in binDic[ctg1]:
            for ctg2 in binDic[ctg1][bini]:
                qnLst = binDic[ctg1][bini][ctg2]
                if len(qnLst) >= minSA:
                    if ctg1 not in bpDic:
                        bpDic[ctg1] = {}
                    if bini not in bpDic[ctg1]:
                        bpDic[ctg1][bini] = []
                    bpDic[ctg1][bini].extend(qnLst)
    #print(bpDic["utg000598l"])
    #sys.exit()
    mergeBpDic = {}
    for ctg in bpDic:
        if ctg not in mergeBpDic:
            mergeBpDic[ctg] = {}
        biLst = list(bpDic[ctg].keys())
        biLst.sort(key=lambda x:x[0])
        newBi = [list(biLst[0])]
        newQnLst = bpDic[ctg][biLst[0]]
        for i in range(1,len(biLst)):
            if biLst[i][0] < newBi[-1][1]:
                newBi[-1][1] = biLst[i][1]
                newQnLst.extend(bpDic[ctg][biLst[i]])
            else:
                # remove qn duplication
                deDupQn = []
                for qn in newQnLst:
                    if qn not in deDupQn:
                        deDupQn.append(qn)
                mergeBpDic[ctg][tuple(newBi[-1])] = deDupQn
                newBi.append(list(biLst[i]))
                newQnLst = bpDic[ctg][biLst[i]]
        deDupQn = []
        for qn in newQnLst:
            if qn not in deDupQn:
                deDupQn.append(qn)
        mergeBpDic[ctg][tuple(newBi[-1])] = deDupQn
    with open(outPre + ".mergedSplitAlign.txt", 'w') as fout:
        for ctg in mergeBpDic:
            for bini in mergeBpDic[ctg]:
                qnLst = mergeBpDic[ctg][bini]
                qnStr = ",".join(qnLst)
                fout.write("{}\t{}\t{}\t{}\t{}\n".format(ctg, bini[0], bini[1], len(qnLst), qnStr))



def workflow(LISFIle, win, minSA, edge, outPre):
    allLISDIc = read_LISFile(LISFIle, "LIS")
    binDic = check_break_point(allLISDIc, win, outPre, edge)
    merge_bp(binDic, outPre, minSA)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is the script for maaping window UL-ONT reads into fasta, and identifying chimeric contig through split alignment.")
    parser.add_argument('-l', '--lis', default=None,
                        help='<filepath>  draft contig assembly fasta file.')
    parser.add_argument('-w', '--window', default=5000,
                        help='<int> window size, default is 5000.')
    parser.add_argument('-m', '--minSA', default=3,
                        help='<int> Number of minimum split alignment in windows, default is 3.')
    parser.add_argument('-e', '--edge', default=20000,
                        help='<int> Number of minimum split alignment in windows, default is 20000.')
    parser.add_argument('-o', '--output', default="output",
                        help='<str> output file prefix, default is output')
    args = parser.parse_args()
    # fastqFile, win, part, threads
    workflow(args.lis, args.window, args.minSA, args.edge, args.output)
