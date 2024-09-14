#/usr/bin/env python

import logging
import sys
import os
import math
import time
from joblib import Parallel, delayed
import argparse
from copy import deepcopy
import collections
import pandas as pd
import polars as pl 

from cphasing.utilities import xopen, is_empty

logger = logging.getLogger(__name__)

def get_file_type(input_file):
    with xopen(input_file) as fp:
        for line in fp:
            if line.startswith(">"):
                return "fasta"
            elif line.startswith("@"):
                return "fastq"
            else:
                return

"""
0. 比对！
1. 筛选最佳AS的mapping,以 (reads_index, pos_in_genome, AS) 记录, 并根据reads_index 进行排序
2. 根据最佳AS的 pos_in_genome 迭代求解最长上升子串 (Best LMP)
3. 对于同一个contig的多个最长上升子串, 尝试填补中间的空缺。填补标准：
    (0) Located in same string
    (1) 最长子串之间的distance < 2 * window
    (2) score of second-alignment >= 75% * maximum alignment score
4. For Best LMP, we extent the edge of alignment. The extent criteria is:
    (1) Per base alignment accuracy (PBAA) = alignmet score / mapping length * 100
5. For Retained LMP, caculate the chimeric alignment
"""
### return outPre.paf
def minimap_mapping(fasta, reads, window, min_windows, threads, outPre, is_hifi=False):
    """
    check minimap2 version
    """
    checkMinimap2CMD = "minimap2 --version"
    minimapVer = "".join(os.popen(checkMinimap2CMD).readlines())
    minimapVer = float(minimapVer.split('-')[0].split('.')[1])
    if minimapVer < float(24):
        print("Warnning: The minimap2 version should be 2.24 or higher.")
        sys.exit()
    preset = "map-ont" if not is_hifi else "map-hifi"

    if ',' in reads:
        readsLst = reads.split(',')
    else:
        readsLst = [reads]

    file_type = get_file_type(readsLst[0])
    if not file_type:
        logger.error("Input error: please input fasta or fastq.")
        sys.exit()

    if readsLst[0][-3:] == ".gz":
        cat_cmd = "pigz -p 4 -dc"
    else:
        cat_cmd = "cat"

    min_length = window * min_windows

    if len(readsLst) > 1:
        minimap2CMD = "{} {} | cphasing-rs slidefq - -w {} -l {} -f {}\
                            | minimap2 -I 16g -t {} --qstrand -cx {} \
                            -p.3 {} - > {}.paf".format(
                                cat_cmd, ' '.join(readsLst), window, min_length, 
                                        file_type, threads, preset, fasta,  outPre)
    else:
        minimap2CMD = "cphasing-rs slidefq {} -w {} -l {} -f {} \
                            | minimap2 -I 16g -t {} --qstrand -cx {} \
                            -p.3 {} - > {}.paf".format(
                                readsLst[0], window, min_length, file_type, threads, preset, fasta,  outPre)
    
    logger.info("Running Command:")
    logger.info(f"\t\t{minimap2CMD}")
    flag = os.system(minimap2CMD)
    assert flag == 0, "Failed to run the `minimap2`"

    return outPre + ".paf"

### return pafDic
def read_paf(paf, minAS, nhap):
    pafDic = {}
    minAS = int(minAS)
    maxAlign = int(nhap)
    logger.info(f"Loading `{paf}` ...")

    with xopen(paf,'r') as fin:
        for line in fin:
            tmpLst = line.split("\t")
            if len(tmpLst) < 15:
                continue
            qn, ql, qs, qe, s, tn, tl, ts, te, al1, al2, mapq, NM, ms, AS = tmpLst[:15]
            if not qn:
                continue
            qi = int(qn.split("_")[-1])
      
            originQn = "_".join(qn.split("_")[:-1])
            NM, AS = list(map(lambda x: int(x.split(":")[-1]), [NM, AS]))
            if AS < minAS:
                continue
            if originQn not in pafDic:
                pafDic[originQn] = {}
            if qi not in pafDic[originQn]:
                pafDic[originQn][qi] = []
            pafDic[originQn][qi].append([int(qs),int(qe),s,tn,int(ts),int(te),int(ql),int(tl),int(al1),int(al2), int(mapq), int(AS)])
    # os.environ['POLARS_MAX_THREADS'] = '4'
    # df = pl.read_csv(paf, has_header=False, 
    #                  separator='\t', columns=range(15))

    # df = df.filter(df['column_1'].is_not_null())
    # df = (df.with_columns(pl.col("column_1")
    #                 .str.extract_groups(r"(.+)_(\d+)$")
    #                 .struct.rename_fields(['originQn', 'qi'])
    #                 .alias('fields')
    #                 ).unnest('fields')).drop("column_1")
    # df = df.with_columns([
    #     df['originQn'].cast(pl.Categorical),
    #     df['qi'].cast(pl.UInt32)
    # ])

    # df = df.with_columns(pl.col("column_15").str.splitn(":", 3).struct[2].cast(pl.Int32).alias('AS')).drop("column_15")
    # df = df.with_columns(pl.col("column_13").str.splitn(":", 3).struct[2].cast(pl.Int32).alias('NM')).drop("column_13")
    # df = df.filter(df['AS'] >= minAS)

    # print(df)
    
    # print(df.group_by(['originQn', 'qi']).apply(lambda x: x.head(maxAlign)))
    ## sort pafDic
    
    for qn in pafDic:
        qnDic = pafDic[qn]
        for qi in qnDic:
            totalReadsLst = qnDic[qi]
            totalReadsLst.sort(key=lambda x:x[-1], reverse=True)
            totalReadsLst = totalReadsLst[:maxAlign+1] if len(totalReadsLst) >= maxAlign else totalReadsLst
            qnDic[qi] = totalReadsLst
        pafDic[qn] = qnDic

    return pafDic


def select_mapq(pafDic, minMapq):
    """
    bestASPafDic = {"reads1":{reads_pos:[reads_index1,reads_index2]}}
    """
    minMapq = int(minMapq)
    bestASPafDic = {}
    mapqPafDic = {}
    for qn in pafDic:
        qnDic = pafDic[qn]
        for qi in qnDic:
            totalReadsLst = qnDic[qi]
            maxAS = totalReadsLst[0][-1]
            for ri in range(len(totalReadsLst)):
                r = totalReadsLst[ri]
                mapq = r[10]
                if mapq >= minMapq:
                    if qn not in mapqPafDic:
                        mapqPafDic[qn] = {}
                    if qi not in mapqPafDic[qn]:
                        mapqPafDic[qn][qi] = []
                    mapqPafDic[qn][qi].append(ri)
                AS = r[11]
                if maxAS == AS:                    
                    if qn not in bestASPafDic:
                        bestASPafDic[qn] = {}
                    if qi not in bestASPafDic[qn]:
                        bestASPafDic[qn][qi] = []
                    bestASPafDic[qn][qi].append(ri)

    return bestASPafDic, mapqPafDic


def reads2ctg(ASreadsDic, mapqPafDic, allreadsDic):
    ctg2AnchorDic = {}
    addDic = {}
    for qi in ASreadsDic:
        for ri in ASreadsDic[qi]:
            addDic[(qi, ri)] = 'P'
    for qi in mapqPafDic:
        for ri in mapqPafDic[qi]:
            addDic[(qi, ri)] = 'M'
    for qi, ri in addDic:
        alignment = allreadsDic[qi][ri]
        tn= alignment[3]
        if tn not in ctg2AnchorDic:
            ctg2AnchorDic[tn] = []
        ctg2AnchorDic[tn].append((qi, ri, addDic[(qi, ri)]))

    return ctg2AnchorDic

### return mergedlLIS : LIS in one reads
def calculate_LMP_pipeline(bestReadsDic, allreadsDic, mapqPafDic, win, qn):
    LISBAK = []
    lisWin = int(float(win) * 1.25)
    overlapWin = int(float(win)*0.25)
    winDis = 2 * win
    ctg2AnchorDic = reads2ctg(bestReadsDic, mapqPafDic, allreadsDic)
    for tn in ctg2AnchorDic:
        alignLst = ctg2AnchorDic[tn]
        alignLst.sort(key=lambda x:x[0])
        ## get LIS
        LIS_result = LIS(alignLst, allreadsDic, lisWin, overlapWin)
        LIS_result.sort(key = lambda x : x[0][0][0])
        ## fill LIS gap
        if len(LIS_result) == 0:
            continue
        filledGapLIS = fill_gap(allreadsDic, LIS_result, winDis, tn)
        filledGapLIS.sort(key = lambda x : x[0][0][0])
        ## get Tl & find LIS which loacated in 5' or 3' of this contig.
        finalLIS = extent_LIS(filledGapLIS, allreadsDic)
        LISBAK.extend(finalLIS)
    ### merge LIS and remove overlap
    """
    LIS_list : [([(qi, ri, 'P'), (qi, ri, 'P')], '+'), ]
                 [(0 [0,(0)]]
    """
    LISBAK.sort(key = lambda x : x[0][0][0])
    mergedlLIS = merge_LIS(LISBAK)
    return (qn, mergedlLIS)
        
def explainLIS(allreadsDic, LISBAK):
    """
    allreadsDic[qiPre][riPre] : [qs,qe,s,tn,ts,te,al1,al2,AS]
    This function is disabled.
    """
    LISlocation = []
    for LIS, string in LISBAK:
        if len(LIS) == 1:
            try:
                [qi ,ri] = LIS[0][:2]
            except:
                print(LIS[0][:2])
                sys.exit()
            ctn, beginPos, endPos = allreadsDic[qi][ri][3:6]
           # LISlocation.append((ctn, beginPos, endPos, string))
        elif len(LIS) >= 2 and string == "-":
            beginQI, beginRI = LIS[0][:2]
            endQI, endRI = LIS[-1][:2]
            ctn, beginPos, endPos = allreadsDic[beginQI][beginRI][3], allreadsDic[endQI][endRI][4], allreadsDic[beginQI][beginRI][5]
        elif len(LIS) >= 2 and string == "+":
            beginQI, beginRI = LIS[0][:2]
            endQI, endRI = LIS[-1][:2]
            ctn, beginPos, endPos = allreadsDic[beginQI][beginRI][3], allreadsDic[beginQI][beginRI][4], allreadsDic[endQI][endRI][5]
        LISlocation.append((ctn, beginPos, endPos, string))
    return LISlocation

### return LIS_list
def LIS(nums, allreadsDic, winDis, overlapWin):
    """
    nums: [(qi, ri)]
    """
    """
        rst : [(pos_in_reads, pos_in_ref)]
        align : [qs,qe,s,tn,ts,te,ql,tl,al1,al2,AS]
        processNum : [(qi, ri)]
        LIS_list : [([(qi, ri, 'P'), (qi, ri, 'P')], '+'), ]
        alignment type: 'P' primary
                        'S' secondary aligment
                        'E' extend aligment
    """
    LIS_list = []
    nums_new = [None] * len(nums)
    while nums:
        nums_index = 0
        rst = [nums[0]]
        fqi, fri, ftype = rst[0]
        falign = allreadsDic[fqi][fri]
        fstring = falign[2]
        n = len(nums)
        for i in range(1, n):
            qi, ri, type = nums[i]
            align = allreadsDic[qi][ri]
            
            if align[2] != fstring:
                continue
            pos1 = int(align[4])
            preQi, preRi, preType = rst[-1]
            preAlign = allreadsDic[preQi][preRi]
            pos2 = int(preAlign[4])
            posDis = pos1 - pos2
            
            if fstring == "-" and (posDis + winDis) >= 0 and (posDis + overlapWin) <= 0:
                rst.append(nums[i])
            elif fstring == "+" and (winDis - posDis) >= 0 and (-posDis + overlapWin) <= 0:
                rst.append(nums[i])
            else:
                nums_new[nums_index] = nums[i]
                nums_index += 1
        
    
        nums = nums_new[:nums_index]
        if len(list(filter(lambda x: x[2] == 'M', rst))) >= 2:
            LIS_list.append((rst, fstring))
    

    return LIS_list
 
### return newLIS : [([(qi, ri), (qi, ri)], '+'), ]
def fill_gap(allreadsDic, LIS_result, winDis, tn):
    """
    LIS_list : [([(qi, ri, 'P'), (qi, ri, 'P')], '+'), ]
    allreadsDic : []
    allreadsDic[qiPre][riPre] : [qs,qe,s,tn,ts,te,al1,al2,AS]
    """
    n = len(LIS_result)
    if n == 0:
        print("Can't find LIS of best alignment! Please check!")
        sys.exit()
    if n == 1:
        return LIS_result

    fillSecAli = {}
    for i in range(0,n-1):
        preLIS = LIS_result[i]
        curLIS = LIS_result[i + 1]
        preString = preLIS[-1]
        curString = curLIS[-1]
        if preString != curString:
            continue
        lAlignPre = preLIS[0][-1]
        fAlignCur = curLIS[0][0]
        qiPre, riPre, typePre = lAlignPre
        qiCur, riCur, typeCur = fAlignCur
        if preString == '+':
            lEndPos = allreadsDic[qiPre][riPre][5]
            fStartPos = allreadsDic[qiCur][riCur][4]
        elif preString == '-':
            lEndPos = allreadsDic[qiCur][riCur][5]
            fStartPos = allreadsDic[qiPre][riPre][4]
        else:
            print("Strange String {}".format(preString))
            sys.exit()
        ###
        #lIdxInPre, lPosInPre
        insertFlag = False
        posDis = fStartPos - lEndPos
        #print(fStartPos, lEndPos, posDis, winDis, posDis + winDis)
        try:
            if (winDis - posDis) >= 0 and posDis > 0:
                insertFlag = True
        except:
            print(posDis, winDis)
            sys.exit()
        #print(insertFlag)
        if insertFlag == False:
            continue
        #lIdxInPre, lPosInPre = preLIS[0][-1]
        #fIdxInCur, fPosInCur = curLIS[0][0]
        if (qiCur - qiPre) <= 1:
            continue
        """
        align : [qs,qe,s,tn,ts,te,al1,al2,AS]
        """
        tmpLst = []
        for qi in range(qiPre+1, qiCur):
            #print(qiPre+1, qiCur)
            if qi not in allreadsDic:
                continue
            else:
                allAlignLst = allreadsDic[qi]
            maxAS = allAlignLst[0][-1]
            cutoffAS = int(0.75 * maxAS)
            #print(maxAS, cutoffAS)
            for alignIdx in range(1, len(allAlignLst)):
                align = allAlignLst[alignIdx]
                ts, ctg, s, AS = align[4], align[3], align[2], align[-1]
                #print(ts, ctg, s, AS)
                #print(ts, max(fStartPos, lEndPos), min(fStartPos, lEndPos))
                if AS == maxAS or ctg != tn or s != preString:
                    continue
                if AS < cutoffAS:
                    break
                if ts <= fStartPos and ts >= lEndPos:
                    tmpLst.append((i, i+1, qi, alignIdx))
                    break
        if bool(tmpLst) == True:
            fillSecAli["{}_{}".format(i, i+1)] = tmpLst
    #print(list(fillSecAli.keys()))

    """
    fillSecAli : [(0, 1, qi, ri), (1, 2, qi, ri)]
    LIS_result : [([(qi, ri, 'P'), (qi, ri, 'P')], '+'), ]
    merge continue index
    """
    idxLst = list(fillSecAli.keys())
    idxLst.sort()
    if bool(idxLst) == False:
       return LIS_result
    try:
        newIdxLst = [idxLst[0].split('_')]
    except:
        #print(LIS_result)
        print(idxLst)
        sys.exit()
    newSecAli = [fillSecAli[idxLst[0]]]
    for i in range(1, len(idxLst)):
        idx1, idx2 = idxLst[i].split('_')
        if newIdxLst[-1][1] == idx1:
            newIdxLst[-1][1] = idx2
            newSecAli[-1].extend(fillSecAli[idxLst[i]])
        else:
            newIdxLst.append(idxLst[i].split('_'))
            newSecAli.append(fillSecAli[idxLst[i]])
        #curIdx = idxLst[i].split('_')
    newLIS = []
    fillgapIndex = []
    allIdx = list(range(n))
    for fillIdx in range(len(newIdxLst)):
        #preIdx, curIdx = list(map(int, newIdxLst[fillIdx]))
        preIdx, curIdx = list(map(int, newIdxLst[fillIdx]))
        fillgapIndex.append(preIdx)
        tmpFillLst = newSecAli[fillIdx]
        tmpNewLis, tmpLisString = LIS_result[preIdx]
        for lisIdx in range(preIdx + 1, curIdx + 1):
            fillgapIndex.append(lisIdx)
            ## fill gap
            for idx1, idx2, qi, ri in tmpFillLst:
                if idx1 == lisIdx - 1 and idx2 == lisIdx:
                    tmpNewLis.append((qi, ri, 'S'))
            ## junction Lis
            curLis = LIS_result[lisIdx][0]
            tmpNewLis.extend(curLis)
        newLIS.append((tmpNewLis, tmpLisString))
    ## add unfill gap
    for idx in allIdx:
        if idx not in fillgapIndex:
            newLIS.append(LIS_result[idx])
        
    return newLIS


def extent_LIS(filledGapLIS, allreadsDic):
    """
    LIS_list : [([(qi, ri, 'P'), (qi, ri, 'P')], '+'), ]
    
    """
    qiLst = list(allreadsDic.keys())
    qiLst.sort(reverse=True)
    maxQi = qiLst[0]
    #print(filledGapLIS[0][0])
    #sys.exit()
    fqi, fri, ftype = filledGapLIS[0][0][0]
    ## [qs,qe,s,tn,ts,te,ql,tl,al1,al2,AS]
    tn, tl = allreadsDic[fqi][fri][3], allreadsDic[fqi][fri][7]
    for lis in filledGapLIS:
        fLis, lLis, lisString = lis[0][0], lis[0][-1], lis[1]
        fqi, fri = fLis[:2]
        lqi, lri = lLis[:2]
        fts, fte = allreadsDic[fqi][fri][4:6]
        lts, lte = allreadsDic[lqi][lri][4:6]
        if lisString == '+':
            LIS_s, LIS_e = fts, lte
        elif lisString == '-':
            LIS_s, LIS_e = lts, fte
        else:
            pass
        ## extend left
        if lisString == '+' and LIS_s <= 5000:
            if fqi == 0:
                pass
            else:
                leftQi = fqi - 1
                #bestAPB = AS_per_base(allreadsDic, leftQi, lisString)
                extendAI = AS_per_base(allreadsDic, leftQi, tn, LIS_s, lisString, 5000, 'left')
                """
                bestAPB : apb, tn, ai
                """
                if bool(extendAI):
                    lis[0].insert(0, (leftQi, extendAI, 'E'))
        elif lisString =='-' and LIS_s <= 5000:
            if lqi == maxQi:
                pass
            else:
                leftQi = lqi + 1
                extendAI = AS_per_base(allreadsDic, leftQi, tn, LIS_s, lisString, 5000, 'left')
                #bestAPB = AS_per_base(allreadsDic, leftQi, lisString)
                if bool(extendAI):
                    lis[0].append((leftQi, extendAI, 'E'))
        else:
            pass
        ## extend right
        if lisString == '+' and (tl - LIS_e) <= 5000:
            if lqi == maxQi:
                pass
            else:
                rightQi = lqi + 1
                #bestAPB = AS_per_base(allreadsDic, leftQi, lisString)
                extendAI = AS_per_base(allreadsDic, rightQi, tn, LIS_e, lisString, 5000, 'right')
                """
                bestAPB : apb, tn, ai
                """
                if bool(extendAI):
                    lis[0].append((rightQi, extendAI, 'E'))
        elif lisString =='-' and (tl - LIS_e) <= 5000:
            if fqi == 0:
                pass
            else:
                rightQi = fqi - 1
                extendAI = AS_per_base(allreadsDic, rightQi, tn, LIS_e, lisString, 5000, 'right')
                #bestAPB = AS_per_base(allreadsDic, leftQi, lisString)
                if bool(extendAI):
                    lis[0].insert(0, (rightQi, extendAI, 'E'))
        else:
            pass
    return filledGapLIS
        

def AS_per_base(allreadsDic, leftQi, lisTn, lisEdge, lisString, winCutoff, extendDirction):
    """
    extend cretiria: 
    (1) some string ,same contig
    (2) distanse < 1.25 * win
    (3) maximum apb, AS > 1000
    """
    apbLst = []
    if leftQi not in allreadsDic:
        return False
    for ai in range(len(allreadsDic[leftQi])):
        align = allreadsDic[leftQi][ai] 
        s, tn, ts, al ,AS = align[2], align[3], align[4], align[-3], align[-1]
        if AS < 1000:
            continue
        try:
            apb = int(float(AS) / float(al) * float(1000))
        except:
            print(s, tn, ts, al ,AS)
            sys.exit()
        tmpAPB = (apb, tn, s, ts, ai)
        apbLst.append(tmpAPB)
    apbLst.sort(key=lambda x:x[0], reverse = True)
    #print(apbLst)
    bestAPB = apbLst[0]
    if bestAPB[1] != lisTn or bestAPB[2] != lisString:
        return False
    if extendDirction == 'left':
        dis = lisEdge - bestAPB[3]
    elif extendDirction == 'right':
        dis = bestAPB[3] - lisEdge
    if dis <= winCutoff:
        return bestAPB[-1]
    else:
        return False

def merge_LIS(finalLIS):
    """
    1. 计算每个LIS的覆盖度 覆盖度 > 2 的即为dup LIS
    2. 首先需要去掉包含在其他LIS中的dup LIS
    3. 采用迭代去掉最短的dup LIS来实现LIS的去冗余 直到LIS的平均覆盖度<2
    4. 此时可能存在一定的重叠 需要找到重叠的区域 再进行比较
    """
    ### caculate depth
   # print(finalLIS)
    depthDic, resultDic = calDepth(finalLIS)
    """
    return : [(lis_index, lis_length, lis_ava_depth, count_of_mapq10, count_of_bestAlignment)]
    """
   # print(resultDic)
    dupLIS = []
    for li in resultDic:
        avaDepth = resultDic[li][1]
        if avaDepth >= 2:
            dupLIS.append((li, resultDic[li][0], resultDic[li][1], resultDic[li][2], resultDic[li][3]))
            #windows_depth[qi] += 1
    ## dedup LIS
    while len(dupLIS) >= 1:
        #print(finalLIS)
        dupLIS.sort(key=lambda x: (x[1], -x[2], x[3], x[4]))
        theMostDup = dupLIS[0][0]
       # try:
        finalLIS.pop(theMostDup)
        depthDic, resultDic = calDepth(finalLIS)
        dupLIS = []
        for li in resultDic:
            avaDepth = resultDic[li][1]
            if avaDepth >= 2:
                dupLIS.append((li, resultDic[li][0], resultDic[li][1], resultDic[li][2], resultDic[li][3]))

    ## find overlaping
    """
    LIS_list : [([(qi, ri, 'P'), (qi, ri, 'P')], '+'), ]
    """
    finalLIS.sort(key=lambda x:x[0][0][0])
    depthDic, resultDic = calDepth(finalLIS)
    overlapDic = {}
    for qi in depthDic:
        if len(depthDic[qi]) >=2:
            overlapLiLst = depthDic[qi]
            """
            [3,4,5]
            {3:{4:[1,2,3]}; 4:}
            """
            if overlapLiLst[0] not in overlapDic:
                overlapDic[overlapLiLst[0]] = {}
            for tmpLi in overlapLiLst[1:]:
                if tmpLi not in overlapDic[overlapLiLst[0]]:
                    overlapDic[overlapLiLst[0]][tmpLi] = []
                overlapDic[overlapLiLst[0]][tmpLi].append(qi)
    if bool(overlapDic) == False:
        return finalLIS
    delLIS = ""
    delRegion = []
    """
    Note: 这里可能要加上一个不能删改contig边界mapq的代码 2023/06/18
    """
    for li1 in overlapDic:
        for li2 in overlapDic[li1]:
            lq, lr, lt = finalLIS[li1][0][-1]
            nq, nr, nt = finalLIS[li2][0][0]
            if lt == "E" and nt != "E":
                delLIS = li2
                delRegion = overlapDic[li1][li2]
            elif nt == "E" and lt != "E":
                delLIS = li1
                delRegion = overlapDic[li1][li2]
            else:
                pass
        if bool(delRegion) == False:
                delLIS = li2
                delRegion = overlapDic[li1][li2]
        finalLIS = delete_align(finalLIS, delLIS, delRegion)
        """
        如果存在重叠 有几种情况：
        1. 末端延伸 以 E 作为tag, 这种情况下 我们认为E的比对更好, 所以把非 E tag的比对删除
        2. 次优比对的gap 以 S 作为tag, 这种情况一般发生在contig中段, 会被先前的迭代过滤所过滤掉
        3. 同源交换的纯和区域
        4. contig末端的重复, 3/4 的情况下为了保证断点的一致性, 将重叠区都分给前一个LIS
        """
    return finalLIS
    

def delete_align(finalLIS, delLIS, delRegion):
    """
    delete overlap LIS
    """
    newTmpLIS = []
    tmpLIS = finalLIS[delLIS]
    tmpAlign, tmpS = tmpLIS
    for i in range(len(tmpAlign)):
        tmpq, tmpr, tmpt = tmpAlign[i]
        if tmpq in delRegion:
            continue
        newTmpLIS.append((tmpq, tmpr, tmpt))
    finalLIS[delLIS] = (newTmpLIS, tmpS)
    return finalLIS
        

def calDepth(finalLIS):
    """
    return : [(lis_index, lis_length, lis_ava_depth, count_of_mapq10, count_of_bestAlignment)]
    """
    resultDic = {}
    depthDic = {}
    for li in range(len(finalLIS)):
        lis, lisS = finalLIS[li]
        #print(lis)
        #sys.exit()
        for qi, ri, aType in lis:
            if qi not in depthDic:
                depthDic[qi] = []
            depthDic[qi].append(li)
        resultDic[li] = [len(lis)]
    for li in range(len(finalLIS)):
        lis, lisS = finalLIS[li]
        lisDepth = 0
        countMapq10 = 0
        countmaxAS = 0
        for qi, ri, aType in lis:
            lisDepth += len(depthDic[qi])
            if aType == 'M':
                countMapq10 += 1
            elif aType == 'P':
                countmaxAS += 1
            else:
                pass
        countmaxAS += countMapq10
        lisAvaDepth = float(lisDepth)/float(len(lis))
        resultDic[li].extend([lisAvaDepth, countMapq10, countmaxAS])
    return depthDic, resultDic

def checkLIS(lisBak, pafDic):
    corLis = {}
    for qn in lisBak:
        corLis[qn] = []
        for align, string in lisBak[qn]:
            qi0, ri0, tp0 = align[0]
            tmpLst = [(qi0, ri0, tp0)]
            if string == '-':
                preS = preE = 9999999999999
            else:
                preS, preE = -1, -1
            for i in range(1, len(align)):
                qi, ri, tp = align[i]
                itemInPaf = pafDic[qn][qi][ri]
                itemInPaf = list(map(str,itemInPaf))
                curS, curE = itemInPaf[4:6]
                curS = int(curS)
                curE = int(curE)
                if preS < curS and preE < curE and string == '+':
                    tmpLst.append((qi, ri, tp))
                    preS, preE = curS, curE
                elif preS > curS and preE > curE and string == '-':
                    tmpLst.append((qi, ri, tp))
                    preS, preE = curS, curE
                else:
                    pass
            corLis[qn].append((tmpLst, string))
    return corLis


## return LISFIle, newPAF 
def outputLIS(lisBak, pafDic, outpre):
    """
    LIS_list : [([(qi, ri, 'P'), (qi, ri, 'P')], '+'), ]
    pafDic[originQn][qi].append([qs,qe,s,tn,ts,te,ql,tl,al1,al2,mapq, AS])
    qn, ql, qs, qe, s, tn, tl, ts, te, al1, al2, mapq, NM, ms, AS = tmpLst[:15]
    """
    mapqLIS = []
    mapqPaf = []

    with open(outpre + '.LIS.gtf', 'w') as fout1, \
        open(outpre + ".corrected.paf", 'w') as fout2:
        for qn in lisBak:
            for align, string in lisBak[qn]:
                fq, fr = align[0][:2]
                lq, lr = align[-1][:2]
                fInfo = pafDic[qn][fq][fr]
                lInfo = pafDic[qn][lq][lr]
                if string == '-':
                    tl, tn, alignS, alignE = lInfo[7], lInfo[3], lInfo[4], fInfo[5]
                else:
                    tl, tn, alignS, alignE = lInfo[7], fInfo[3], fInfo[4], lInfo[5]
                """
                seq_id	seq_length	type	start	end	Alignment-score	strand	alignment-type	attributes
                ctg12	400000	LIS/alignment	25132483	25132543	8000	+	P/S/E/M	reads_id "reads1"; window_id "reads1_1"; 
                """
                atrriInfoLIS = "reads_id:" + qn
                tmpLIS = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(tn, tl, "LIS", alignS, alignE, ".", string, ".", atrriInfoLIS)
                fout1.write(tmpLIS)
                alignLisLst = [tmpLIS]
                alignPafLst = []
                addFlag = False
                for qi, ri ,tp in align: 
                    if tp == 'M':
                        addFlag = True
                    itemInPaf = pafDic[qn][qi][ri]
                    itemInPaf = list(map(str,itemInPaf))
                    """
                    [qs,qe,s,tn,ts,te,ql,tl,al1,al2,mapq, AS]
                    qn, ql, qs, qe, s, tn, tl, ts, te, al1, al2, mapq, AS
                    """
                    _qn, ql, [qs, qe, s, tn], tl, [ts, te], [al1, al2, mapq, AS] = "{}_{}".format(qn, qi), itemInPaf[6], itemInPaf[:4], itemInPaf[7], itemInPaf[4:6], itemInPaf[8:]
                    
                    mapq = "2" if mapq == "0" else mapq
                    # _qn2 = _qn.rsplit("_", 1)[0]
                    tmpAliPaf = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tAS:{}\n".format(_qn, ql, qs, qe, s, tn, tl, ts, te, al1, al2, mapq, AS)
                    fout2.write(tmpAliPaf)
                    tl, tn, ts, te, AS = itemInPaf[7], itemInPaf[3], itemInPaf[4], itemInPaf[5], itemInPaf[-1]
                    atrriInfo = "reads_id:{0};window_id:{0}_{1};".format(qn, qi)
                    tmpAliLis = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(tn, tl, "alignment", ts, te, AS, string, tp, atrriInfo)
                    fout1.write(tmpAliLis)
                    alignLisLst.append(tmpAliLis)
                    alignPafLst.append(tmpAliPaf)
                if addFlag == True:
                    mapqPaf.extend(alignPafLst)
                    mapqLIS.extend(alignLisLst)

    with open(outpre + '.mapq.LIS.gtf', 'w') as fout3:
        #for i in mapqLIS:
        fout3.write("".join(mapqLIS))
        logger.info(f"Successfull output {fout3.name}")
    with open(outpre + '.mapq.corrected.paf', 'w') as fout4:
        #for i in mapqLIS:
        fout4.write("".join(mapqPaf))
        logger.info(f"Successful output {fout4.name}")
                        


def workflow(fasta, reads, threads, outPre, win, min_windows, nhap, minAS, minMapq, hifi=False):
    """
    1. get LIS for every single reads
    2. output corrected alignment & LIS path
    3. determining chimeric break point
    4. break chimeric and output
    """
    pafFile = outPre + ".paf"
    if os.path.exists(pafFile) == False or is_empty(pafFile):
        pafFile = minimap_mapping(fasta, reads, win, min_windows, threads, outPre, hifi)
    else:
        logger.warning(f"Using existed mapping results: `{pafFile}`")

    threads, win = int(threads), int(win)
    #pafFile = minimap_mapping(fasta, reads, threads, outPre)
    pafDic = read_paf(pafFile, minAS, nhap)
    bestASpafDic, mapqPafDic = select_mapq(pafDic, minMapq)

    allLISDIc = {}
    ## get LIS for every single reads
    qnLst = []
    allAlignLst = []
    bestAlignLst = []
    mapqAlignLst = []

    for qn in pafDic:
        qnLst.append(qn)
        allAlignLst.append(pafDic[qn])
        if qn in mapqPafDic:
            mapqAlignLst.append(mapqPafDic[qn])
        else:
            mapqAlignLst.append({})
        bestAlignLst.append(bestASpafDic[qn])
    AllmergedLIS = Parallel(n_jobs=min(threads, 12))(
                    delayed(calculate_LMP_pipeline)(c, b, d, int(win), a) 
                        for a,b,c,d in zip(qnLst, allAlignLst, bestAlignLst, mapqAlignLst))
    
    for qn, mergedLIS in AllmergedLIS:
        allLISDIc[qn] = mergedLIS
       
    ## check LIS
    corLISDic = checkLIS(allLISDIc, pafDic)
    ## output
    
    outputLIS(corLISDic, pafDic, outPre)

    return pafFile



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is the script for maaping window UL-ONT reads into fasta, and identifying chimeric contig through split alignment.")
    parser.add_argument('-f', '--fasta', default=None,
                        help='<filepath>  draft contig assembly fasta file.')
    parser.add_argument('-i', '--fastq', default=None,
                        help='<filepath>  window fastq file.')
    parser.add_argument('-a', '--minAS', default=2000,
                        help='<int> minimun Alignment Score, default is 2000.')
    parser.add_argument('-m', '--minMapq', default=1,
                        help='<int> minimun Alignment Score, default is 1.')
    parser.add_argument('-n', '--nhap', default=4,
                        help='<int> maximum supplement aligment records.')
    parser.add_argument('-w', '--window', default=5000,
                        help='<int> window size, default is 5000.')
    parser.add_argument('-t', '--threads', default=10,
                        help='<int> number of threads used when running minimap2, default is 10')
    parser.add_argument('-o', '--output', default="output",
                        help='<str> output file prefix, default is output')
    args = parser.parse_args()
    # fastqFile, win, part, threads
    workflow(args.fasta, args.fastq, args.threads, args.output, args.window, 10, args.nhap, args.minAS, args.minMapq, hifi=False)

