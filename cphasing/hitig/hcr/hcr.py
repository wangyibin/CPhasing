#!/usr/bin/env python

import argparse
import logging
import os
import os.path as op
import sys

import numpy as np
import pandas as pd
import pyranges as pr 

from collections import OrderedDict

from ...utilities import read_chrom_sizes

"""
1. min length
2. M count
3. overlap with depth
"""

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

def read_LIS(LisFile, SAreads, minCount, minMapqCount):
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
        if readsID in SAreads:
            continue
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


def read_depth(depthFile):
    # depthDic = {}
    # with open(depthFile, 'r') as fin:
    #     for line in fin:
    #         tmpLst = line.rstrip().split('\t')
    #         ctg, s, e = tmpLst
    #         s, e = map(int, [s,e])
    #         #depth = float(depth)
    #         if ctg not in depthDic:
    #             depthDic[ctg] =[]
    #         depthDic[ctg].append((s,e))

    # return depthDic
    df = pd.read_csv(depthFile, sep='\t', index_col=None, header=None, 
                        names=['Chromosome', 'Start', 'End', 'Count'])
    
    return df


def read_SAreads(SAFile):
    SAreads = {}
    with open(SAFile, 'r') as fin:
        for line in fin:
            tmpLst = line.rstrip().split('\t')
            for reads in tmpLst[4].split(','):
                SAreads[reads] = ""

    return SAreads

## Deprecated
def overlap_with_depth(sequential_LIS, depth_region):
    overlapDic = {}
    for ctg in sequential_LIS:
        if ctg not in depth_region:
            continue
        lisLst = sequential_LIS[ctg]
        depthLst = depth_region[ctg]
        lisLst.sort(key=lambda x:x[0])
        depthLst.sort(key=lambda x:x[0])
        i, j = 0, 0
        res = []
        while i < len(lisLst) and j < len(depthLst):
            a1, a2 = lisLst[i][0], lisLst[i][1]
            b1, b2 = depthLst[i][0], depthLst[i][1]
            if b2 >= a1 and a2 >=b1:
                res.append([max(a1, b1), min(a2, b2)])
            if b2 < a2:
                j += 1
            else:
                i += 1
        overlapDic[ctg] = res
    return overlapDic

## Deprecated
def output(overlapDic, outPre):
    with open(outPre + ".hcr_all.bed", 'w') as fout:
        for ctg in overlapDic:
            for region in overlapDic[ctg]:
                fout.write("{}\t{}\t{}\n".format(ctg, *region))

def range_dict2frame(dic):
    res = []
    for ctg in dic:
        for s, e in dic[ctg]:
            res.append((ctg, s, e))
    
    df = pd.DataFrame(res)
    df.columns = ['Chromosome', 'Start', 'End']
    return df 

def correct_hcr_by_break_pos(hcrs, break_pos, contig_sizes, output):
    contig_sizes = read_chrom_sizes(contig_sizes)
    contig_sizes = contig_sizes.to_dict()['length']

    break_pos_db = OrderedDict()
    with open(break_pos, 'r') as fp:
        for line in fp:
            if line.strip():
                line_list = line.strip().split()
                contig, positions = line_list 
                positions = positions.split(",")
                positions = list(map(int, positions))
                if contig not in break_pos_db:
                    break_pos_db[contig] = []
                for pos in positions:
                    break_pos_db[contig].append(pos)
    
    res = []
    for contig in break_pos_db:
        contig_length = contig_sizes[contig]
        pos_list = break_pos_db[contig]

        for pos1, pos2 in zip(np.r_[0, pos_list], 
                                np.r_[pos_list, contig_length]):
            res.append((contig, pos1, pos2))

    hcrs_df = hcrs.df
    hcrs_gr = pr.PyRanges(hcrs_df.reset_index())
    df = pd.DataFrame(res)
    if not df.empty:
        df.columns = ['Chromosome', 'Start', 'End']
        gr = pr.PyRanges(df)
        overlapped = hcrs_gr.join(gr).new_position('intersection')
        overlapped = overlapped.df
    
        corrected_hcrs = []
        drop_idx = []
        for i, row in overlapped.iterrows():
            new_contig_id = f"{row.Chromosome}_{row.Start_b}_{row.End_b}"
            new_start = row.Start - row.Start_b 
            new_end = row.End - row.Start_b 
            corrected_hcrs.append((new_contig_id, new_start, new_end))
            drop_idx.append(row['index'])
        
        corrected_hcrs_df = pd.DataFrame(corrected_hcrs, columns=hcrs_df.columns)
        all_df = pd.concat([hcrs_df.drop(list(set(drop_idx)), axis=0), 
                        corrected_hcrs_df], axis=0)
    else:
        all_df = hcrs_df
    
    all_df.to_csv(output, sep='\t', index=None, header=None)


def workflow(LisFile, SA, depthFile, minCount, minMapqCount, 
                outPre, break_pos=None, contig_sizes=None):
    
    SAreads = read_SAreads(SA)
    logger.info(f"Load `{depthFile}`.")
    depth_df = read_depth(depthFile)
    depth_gr = pr.PyRanges(depth_df)
    depth_gr = depth_gr.merge()
    logger.info(f"Identified `{len(depth_gr)}` normal depth regions.")
    sequential_LIS = read_LIS(LisFile, SAreads, minCount, minMapqCount)
    logger.info("Identifing continuous regions ...")

    sequential_gr = pr.PyRanges(range_dict2frame(sequential_LIS))
    sequential_gr = sequential_gr.merge()
    logger.info(f"Identified `{len(sequential_gr)}` continous regions.")
    logger.info("Overlapping ...")
    overlap_gr = sequential_gr.overlap(depth_gr).merge()
    
    ## output
    output = outPre + ".hcr_all.bed"
    if break_pos and contig_sizes:
        logger.info("Correcting HCRs by chimeric position ...")
        correct_hcr_by_break_pos(overlap_gr, break_pos, contig_sizes, output)
        logger.info(f"Successful output chimeric-corrected `{len(overlap_gr)}` HCRs in `{output}`.")
    else:
        
        overlap_gr.df.to_csv(output, sep='\t', 
                                index=None, header=None)
        logger.info(f"Successful output `{len(overlap_gr)}` HCRs in `{output}`.")

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is the script for filter genome region.")
    parser.add_argument('-d', '--depthFile', required=True,
                        help='<filepath>  UL-ONT reads/Hifi depth file, 4colum, Chr start end depth.')
    parser.add_argument('-sa', '--splitAlign', required=True,
                        help='<filepath>  splitAlign file, 4 col: <chr> <start> <end> <count of split-align> <split-aligned reads>.')
    parser.add_argument('-l', '--lis', required=True,
                        help='<filepath> fasta file.')
    parser.add_argument('-m', '--minCount', default=5,
                        help='<int> minimum count of windows in LIS.')
    parser.add_argument('-M', '--minMapq', default=0.6,
                        help='<float> minmum ratio of (mapq10 window)/(window counts).')
    parser.add_argument('-o', '--output', default="output",
                        help='<str> output file prefix, default is output')
    args = parser.parse_args()
    # paf, breakpointFile, fasta, win, outPre
    workflow(args.lis, args.splitAlign, args.depthFile,  args.minCount, args.minMapq, args.output)