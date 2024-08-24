#!/usr/bin/env python3
#coding=utf-8

import argparse
import sys

def read_fa(fastaFile):
    fastaDic = {}
    with open(fastaFile,'r') as IN:
        fastaName = IN.readline().strip()[1:]
        fa = ''
        for line in IN:
            if line.startswith('>'):
                fastaDic[fastaName] = fa
                fastaName = line.strip()[1:]
                fa = ''
            else:
                fa += line.rstrip()
        fastaDic[fastaName] = fa
    return fastaDic

def read_bp(bpFile):
    bpDic = {}
    with open(bpFile, 'r') as fin:
        for line in fin:
            ctg, bpLst = line.split('\t')
            bpLst = bpLst.split(',')
            bpLst = list(map(int, bpLst))
            bpDic[ctg] = bpLst
    return bpDic

## return None
def break_contig(bpDic2, fa, outPre):
    faDic = read_fa(fa)
    regionDic = {}
    for ctg in bpDic2:
        bpLst = bpDic2[ctg]
        if len(bpLst) == 1:
            bpRegion = [(1, bpLst[0]), (bpLst[0] + 1, len(faDic[ctg]))]
            regionDic[ctg] = bpRegion
            continue
        bpRegion = [(1, bpLst[0])]
        for bpi in range(0, len(bpLst)-1):
            bpRegion.append((bpLst[bpi], bpLst[bpi+1]))
        bpRegion.append((bpLst[-1], len(faDic[ctg])))
        regionDic[ctg] = bpRegion
    with open(outPre + ".corrected.fasta", 'w') as fout:
        for faName in faDic:
            if faName in regionDic:
                corrLst = regionDic[faName]
                for s, e in corrLst:
                   # try:
                    fa = faDic[faName]
                    subfa = fa[s:e]
                    fout.write(">{}:{}-{}\n{}\n".format(faName, s, e, subfa))
        #            except:
        #                print(faDic[faName], s, e)
        #                sys.exit()
            else:
                fout.write(">{}\n{}\n".format(faName, faDic[faName]))

    return outPre + ".corrected.fasta"

def workflow(fastaFile, bpFile, outPre):
    #fastaDic = read_fa(fastaFile)
    bpDic = read_bp(bpFile)
    breaked_fasta = break_contig(bpDic, fastaFile, outPre)
    
    return breaked_fasta


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is the script for chimeric contig correction.")
    parser.add_argument('-f', '--fasta', required=True,
                        help='<filepath>  the raw fasta file.')
    parser.add_argument('-bp', '--breakpos', required=True,
                        help='<filepath>  chimeric position of contigs.')
    parser.add_argument('-o', '--output', default="output",
                        help='<str> output file prefix, default is output')
    args = parser.parse_args()
    # paf, breakpointFile, fasta, win, outPre
    workflow(args.fasta, args.breakpos, args.output)
