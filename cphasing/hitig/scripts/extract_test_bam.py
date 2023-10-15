#!/usr/bin/env python3
#coding=utf-8

import pysam
import sys
import os
from multiprocessing import Pool


def read_bed(bedFile):
    bedDic = {}
    with open(bedFile, 'r') as fin:
        for line in fin:
            tmpLst = line.rstrip().split('\t')
            chr, bed1, bed2 = tmpLst[:4]
            if chr not in bedDic:
                bedDic[chr] = []
            bedDic[chr].append((int(bed1),int(bed2)))
    return bedDic


def check_in_region(chrRead, pos, bedDic):
    #flag = False
    if chrRead not in bedDic:
        return False
    for highPos1, highPos2 in bedDic[chrRead]:
        if highPos1 <= pos & pos <= highPos2:
            return True
        else:
            pass
    return False    



def parse_bam(wrkDir, bamFile, bedDic, chromLst, index):
    bf = pysam.AlignmentFile(bamFile, "rb")
    outBam = wrkDir + "/" + str(index) + '.bam'
    obf1 = pysam.AlignmentFile(outBam, "wb", template=bf)
    #print(chromLst[index])
    for chrom in chromLst[index]:
        #print(chrom)
        for i in bf.fetch(chrom, multiple_iterators=True):
            r1Name = i.reference_name
            r1Pos = int(i.reference_start)
            r2Name = i.next_reference_name
            r2Pos = int(i.next_reference_start)
            flag1 = check_in_region(r1Name, r1Pos, bedDic)
            flag2 = check_in_region(r2Name, r2Pos, bedDic)
            if flag1 == flag2 == True:
                obf1.write(i)

def calcu_size(bedDic):
    sizeDic = {}
    for ctg in bedDic:
        size=0
        for s1, e1 in bedDic[ctg]:
            size += e1-s1
        sizeDic[ctg] = size
    totalSize = 0
    for ctg in sizeDic:
        totalSize += sizeDic[ctg]
    return totalSize, sizeDic

def split_ctg(sizeDic, totalSize, number):
    num_of_part = round(float(totalSize)/float(int(number)))
    count = 0
    chromLst = [[]]
    #tmpLst = []
    for chrom in sizeDic:
        if count <= num_of_part:
            chromLst[-1].append(chrom)
            count += sizeDic[chrom]
        else:
            count = sizeDic[chrom]
            chromLst.append([chrom])
    #print(chromLst)
    return chromLst

try:
    bedFile, bamFile, wrkDir, output, thread = sys.argv[1:6]
except:
    print("Usage: python {} <high_confidence_region|bed_file> <bam_file> <working_directory> <prefix_of_output> <threads>".format(sys.argv[0]))
    sys.exit()
bedDic = read_bed(bedFile)
totalSize, sizeDic = calcu_size(bedDic)
# split_ctg(sizeDic, totalSize, number)

## check 
if os.path.exists(wrkDir):
    print("{} is exists, please remove it and re-run.".format(wrkDir))
    sys.exit()
else:
    pass

## check sambamba
checkCMD="sambamba -h"
if os.system(checkCMD) !=0:
    print("Cannot find sambamba, please install.")
    sys.exit()
else:
    pass

header = list(bedDic.keys())
chromLst = split_ctg(sizeDic, totalSize, thread)
pool = Pool(processes=int(thread))
for chrIndex in range(len(chromLst)):
    pool.apply_async(parse_bam, (wrkDir, bamFile, bedDic, chromLst, chrIndex, ))
#for index in range(len(chromLst)):
#    parse_bam(bamFile, bedDic, chromLst, index)
pool.close()
pool.join()

## merge bam
mergeCMD = "sambamba merge -t {} {}.merged.bam {}/*bam".format(thread, output, wrkDir)
os.system(mergeCMD)
#pysam.merge
