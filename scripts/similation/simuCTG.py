#!/usr/bin/env python
import sys
import argparse
import random
import copy
import time


def GetOpts():
    group = argparse.ArgumentParser()
    group.add_argument('--min', help='minimum length of contig, default: 15k, you can use both number or string end with k,m', default="15k")
    group.add_argument('--max', help='minimum length of contig, default: 5m, you can use both number or string end with k,m', default="5m")
    group.add_argument('-n', '--n50', help='size of N50, default: 500k, you can use both number or string end with k,m', default="500k")
    group.add_argument('-i', '--input', help='origin fasta file of genome', required=True)
    group.add_argument('-o', '--output', help='filename of simulated data', required=True)

    return group.parse_args()


def ReadFasta(inFasta):
    fastaDB = {}
    with open(inFasta, 'r') as fIn:
        id = ''
        seq = ''
        for line in fIn:
            if line[0] == '>':
                if seq != '':
                    fastaDB[id] = seq
                id = line.strip().split(" ")[0][1:]
                seq = ''
            else:
                seq += line.strip()
        fastaDB[id] = seq
    return fastaDB


def GenCtgLen(fastaLenDB, cntLowerDB, cntHigherDB, n50, minLen, maxLen):
    ctgLenDB = {}
    for chrn in cntLowerDB:
        ctgLenDB[chrn] = []
        totalLower = 0
        totalHigher = 0
        totalLen = 0
        for i in range(0, cntLowerDB[chrn]):
            tmpLen = random.randint(minLen, n50)
            totalLower += tmpLen
            if totalLower > fastaLenDB[chrn]/2:
                break
            ctgLenDB[chrn].append(tmpLen)
            totalLen += tmpLen
        for i in range(0, cntHigherDB[chrn]):
            tmpLen = random.randint(n50, maxLen)
            totalHigher += tmpLen
            if totalHigher > fastaLenDB[chrn]/2:
                break
            ctgLenDB[chrn].append(tmpLen)
            totalLen += tmpLen
        cntN50 = int((fastaLenDB[chrn]-totalLen)/n50)
        for i in range(0, cntN50):
            tmpLen = random.randint(int(n50-n50*0.1), int(n50+n50*0.1))
            totalLen += tmpLen
            if totalLen > fastaLenDB[chrn]:
                totalLen -= tmpLen
                break
            ctgLenDB[chrn].append(tmpLen)
        ctgLenDB[chrn].append(fastaLenDB[chrn]-totalLen)
    return ctgLenDB


def GenCtgRegions(fastaLenDB, ctgLenDB):
    ctgRegionsDB = {}
    for chrn in fastaLenDB:
        ctgRegionsDB[chrn] = []
        totalCtgLen = 0
        for ctgLen in ctgLenDB[chrn]:
            totalCtgLen += ctgLen
        cntCtg = len(ctgLenDB[chrn])
        lastPos = 0
        ctgLenList = copy.deepcopy(ctgLenDB[chrn])
        for i in range(0, cntCtg):
            index = random.randint(0, cntCtg-1)
            ctgRegionsDB[chrn].append([lastPos, lastPos+ctgLenList[index]])
            lastPos += ctgLenList[index]
            del ctgLenList[index]
            cntCtg -= 1
    return ctgRegionsDB


def SimGenomeCtg(inFasta, outFasta, n50, minLen, maxLen):
    random.seed()
    print("\033[32m%s\033[0m Reading fasta"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
    fastaDB = ReadFasta(inFasta)
    fastaLenDB = {}
    cntLowerDB = {}
    cntHigherDB = {}
    for chrn in fastaDB:
        fastaLenDB[chrn] = len(fastaDB[chrn])
        cntLowerDB[chrn] = int(fastaLenDB[chrn]/(minLen+n50))
        cntHigherDB[chrn] = int(fastaLenDB[chrn]/(maxLen+n50))
    
    print("\033[32m%s\033[0m Generating contigs"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
    ctgLenDB = GenCtgLen(fastaLenDB, cntLowerDB, cntHigherDB, n50, minLen, maxLen)
    ctgRegionsDB = GenCtgRegions(fastaLenDB, ctgLenDB)

    print("\033[32m%s\033[0m Statistics"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
    for chrn in sorted(fastaDB):
        print("           Chromosome:\t%s"%(chrn))
        print("           Chromosome size:\t%d"%(fastaLenDB[chrn]))
        print("           Contig counts:\t%d"%(len(ctgLenDB[chrn])))
        tmpLen = 0
        n50Len = 0
        for ctgLen in sorted(ctgLenDB[chrn], reverse=True):
            tmpLen += ctgLen
            if tmpLen >= fastaLenDB[chrn]/2 and n50Len == 0:
                n50Len = ctgLen
        print("           Contig total size:\t%d"%(tmpLen))
        print("           N50 size:\t%d"%(n50Len))
    
    print("\033[32m%s\033[0m Writing contigs"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
    with open(outFasta, 'w') as fOut, open('new_genome.posi.bed', 'w') as out:
        for chrn in sorted(fastaDB):
            base = 1
            for region in ctgRegionsDB[chrn]:
                s = region[0]
                e = region[1]
                ctgName = f"{chrn}.ctg{base}"
                out.write(f"{chrn}\t{s}\t{e}\n")
                fOut.write(">%s\n%s\n"%(ctgName, fastaDB[chrn][s: e]))    
                base += 1
    
    print("\033[32m%s\033[0m Finished"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))


if __name__ == "__main__":
    opts = GetOpts()
    inFasta = opts.input
    outFasta = opts.output
    n50 = opts.n50
    n50 = n50.lower()
    n50 = n50.replace('m', '000000')
    n50 = n50.replace('k', '000')
    n50 = int(n50)
    minLen = opts.min
    minLen = minLen.lower()
    minLen = minLen.replace('m', '000000')
    minLen = minLen.replace('k', '000')
    minLen = int(minLen)
    maxLen = opts.max
    maxLen = maxLen.lower()
    maxLen = maxLen.replace('m', '000000')
    maxLen = maxLen.replace('k', '000')
    maxLen = int(maxLen)
    SimGenomeCtg(inFasta, outFasta, n50, minLen, maxLen)