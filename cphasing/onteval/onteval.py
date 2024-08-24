#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
evaluate the asssembly by ultra-long reads
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)

def check_Mcount(typeLst, minCount, minMapqCount):
    mCount = typeLst.count("M")
    if len(typeLst) < minCount:
        return False
    ratio = float(mCount)/float(len(typeLst))
    if ratio >= minMapqCount:
        return True
    else:
        return False

def read_LIS(LisFile, minCount=3, minMapqCount=0):
    outDic = {}
    sequential_LIS = []
    with open(LisFile, 'r') as fin:
        data = fin.readline().rstrip().split('\t')
        ctg, s, e, String = data[0], int(data[3]), int(data[4]), data[6]
        typeLst = []
        pReadsid = data[8].split(';')[0].split("id:")[1]
        alignLst = [(ctg, s, e, String, pReadsid)]
        for line in fin:
            data = line.rstrip().split('\t')
            ctg, s, e, String, aliType = data[0], int(data[3]), int(data[4]), data[6], data[7]
            readsID = data[8].split(';')[0].split("id:")[1]
            if readsID != pReadsid:
                if check_Mcount(typeLst, minCount, minMapqCount) == True:
                    for i in alignLst:
                        sequential_LIS.append(i)
                typeLst = []
                alignLst = [(ctg, s, e, String, readsID)]
                pReadsid = readsID
            elif data[2] == "LIS":
                alignLst.append((ctg, s, e, String, readsID))
            else:
                typeLst.append(aliType)
        if check_Mcount(typeLst, minCount, minMapqCount) == True:
            for i in alignLst:
                sequential_LIS.append(i)
    ## reForm sequential_LIS
    for ctg, s, e, String, readsID in sequential_LIS:
        if ctg not in outDic:
            outDic[ctg] = []
        outDic[ctg].append([s,e,String, readsID])
    ## merge region
    for ctg in outDic:
        lisLst = outDic[ctg]
        lisLst.sort(key=lambda x:x[0])
        newLst = [lisLst[0][:2]]

        for i in range(1,len(lisLst)):
            if lisLst[i][0] <= newLst[-1][1] <= lisLst[i][1]:
                newLst[-1][-1] = max(lisLst[i][1], newLst[-1][1])
            else:
                newLst.append(lisLst[i][:2])
        outDic[ctg] = newLst
    return outDic





def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('lis', 
            help='')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    print(read_LIS(args.lis))

if __name__ == "__main__":
    main(sys.argv[1:])