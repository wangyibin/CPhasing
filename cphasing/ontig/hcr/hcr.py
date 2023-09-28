import argparse


"""
1. min length
2. M count
3. overlap with depth
"""

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
        outDic[ctg].append((s,e,String, readsID))
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
    depthDic = {}
    with open(depthFile, 'r') as fin:
        for line in fin:
            tmpLst = line.rstrip().split('\t')
            ctg, s, e = tmpLst
            s, e = map(int, [s,e])
            #depth = float(depth)
            if ctg not in depthDic:
                depthDic[ctg] =[]
            depthDic[ctg].append((s,e))

    return depthDic

def read_SAreads(SAFile):
    SAreads = {}
    with open(SAFile, 'r') as fin:
        for line in fin:
            tmpLst = line.rstrip().split('\t')
            for reads in tmpLst[4].split(','):
                SAreads[reads] = ""
    return SAreads

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
        
def output(overlapDic, outPre):
    with open(outPre + ".hcr_all.bed", 'w') as fout:
        for ctg in overlapDic:
            for region in overlapDic[ctg]:
                fout.write("{}\t{}\t{}\n".format(ctg, *region))

def workflow(LisFile, SA, depthFile, minCount, minMapqCount, outPre):
    SAreads = read_SAreads(SA)
    depthDic = read_depth(depthFile)
    sequential_LIS = read_LIS(LisFile, SAreads, minCount, minMapqCount)   
    overlapDic = overlap_with_depth(sequential_LIS, depthDic)
    output(overlapDic, outPre)

    
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