#!/usr/bin/env python3
#coding=utf-8

import sys
import collections

from pathlib import Path



def read_LISFile(LisFile, feature):
    LisDic = collections.OrderedDict()
    with open(LisFile, 'r') as fin:
        for line in fin:
            tmpLine = line.rstrip().split('\t')
            lisType = tmpLine[2]
            if lisType == feature:
                qn = tmpLine[-1].split(":")[1]
                tn, tl, s, e, string = tmpLine[0], tmpLine[1], tmpLine[3], tmpLine[4], tmpLine[6]             
                #tmpLine[:2], tmpLine[3:5], tmpLine[6]
                if qn not in LisDic:
                    LisDic[qn] = []
                LisDic[qn].append((tn, int(s), int(e), string, int(tl)))
    return LisDic



def check_break_point(LisDic, win):
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
    ### 过滤掉只有1个窗口的LIS
    for qn in LisDic:
        newLisDic = []
        for lis in LisDic[qn]:
            tn, s, e, string, tl = lis
            if float(abs(s - e)) < 1.25 * float(win):
                continue
            newLisDic.append(lis)
        LisDic[qn] = newLisDic
    ### find split alignment
    linkAlign = {}
    genomeSize = {}
    for qn in LisDic:
        if len(LisDic[qn]) >= 2:   # EDIT HERE
            for i in range(len(LisDic[qn]) - 1):
                ftn, fs, fe, fstr, ftl = LisDic[qn][i]
                ctn, cs, ce, cstr, ctl = LisDic[qn][i+1]
                if ftn == ctn:
                    continue
                [fs, fe, ftl, cs, ce, ctl] = list(map(int, [fs, fe, ftl, cs, ce, ctl]))
                if ftn not in genomeSize:
                    genomeSize[ftn] = ftl
                if ctn not in genomeSize:
                    genomeSize[ctn] = ctl
                for ctg in [ftn ,ctn]:
                    if ctg not in linkAlign:
                        linkAlign[ctg] = {"f1":{}, "f2":{}}
                """
                |--->||---->|
                |--->||<----|
                |<---||---->|
                |<---||<----|
                """
                linkLst = [0,0,0,0]
                if fstr == "+" and cstr == "+":
                    if (ftl - fe) <= 2*win and cs <= 2*win:
                        linkLst = ["f2", ctn, "f1", ftn]
                    else:
                        continue

                elif fstr == "+" and cstr == "-":
                    if (ftl - fe) <= 2*win and (ctl - ce) <= 2*win:
                        linkLst = ["f2", ctn, "f2", ftn]
                    else:
                        continue
                
                elif fstr == "-" and cstr == "+":
                    if fs <= win and cs <= win:
                        linkLst = ["f1", ctn, "f1", ftn]
                    else:
                        continue

                else:
                    if fs <= 2*win and (ctl - ce) <= 2*win:
                        linkLst = ["f1", ctn, "f2", ftn]
                    else:
                        continue
                ## add into splitAlign
                ## add into linkAlign
                ## linkAlign[ctg] = {"f1":{}, "f2":{}} linkLst = ["f1", ctn, "f2", ftn]
                if linkLst == [0,0,0,0]:
                    continue
                fDirec, fLinkCtg, cDirec, cLinkCtg = linkLst
                if (fLinkCtg, cDirec) not in linkAlign[ftn][fDirec]:
                    linkAlign[ftn][fDirec][(fLinkCtg, cDirec)] = 0
                if (cLinkCtg, fDirec) not in linkAlign[ctn][cDirec]:
                    linkAlign[ctn][cDirec][(cLinkCtg, fDirec)] = 0
                linkAlign[ftn][fDirec][(fLinkCtg, cDirec)]+=1
                linkAlign[ctn][cDirec][(cLinkCtg, fDirec)]+=1
    #print(linkAlign)
    #sys.exit()
    return linkAlign


def check_linkage(linkAlign, cutoff):
    """
    |--->||---->|
    |--->||<----|
    |<---||---->|
    |<---||<----|
    linkAlign[ftn][fDirec][(fLinkCtg, cDirec)] += 1
    linkAlign[ctn][cDirec][(cLinkCtg, fDirec)] += 1
    """
    ## calculate propotion
    newLinkAlign = {}
    passLinkAlign = {}
    cutoff = float(cutoff)
    for ctg in linkAlign:
        newLinkAlign[ctg] = {}
        #alleleCtg = allicPair[ctg]
        for direc in linkAlign[ctg]:
            newLinkAlign[ctg][direc] = {}
            tmpLst = list(linkAlign[ctg][direc].items())
            tmpLst.sort(key=lambda x:x[1])
            countLst = [x[1] for x in tmpLst]
            try:
                allLinkCount = sum(countLst)
            except:
                print(tmpLst)
                sys.exit()
            if allLinkCount <=3:
                continue
            for idx in range(len(tmpLst)):
                linkCtg, linkDirec = tmpLst[idx][0]
                linkCount = tmpLst[idx][1]
                #if linkCtg in alleleCtg:
                #    continue
                ratio = float(linkCount)/float(allLinkCount)
                newLinkAlign[ctg][direc][(linkCtg, linkDirec)] = (ratio, linkCount, allLinkCount)
                if ratio >= cutoff:
                    if ctg not in passLinkAlign:
                        passLinkAlign[ctg] = {}
                    if direc not in passLinkAlign[ctg]:
                        passLinkAlign[ctg][direc] = (linkCtg, linkDirec, ratio)
                    else:
                        print("WARNNING: {}-{} is more than 1 links higher than cutoff {}".format(ctg, direc, cutoff))
    with open("link.result", 'w') as fout:
        for ctg in newLinkAlign:
            for direc in newLinkAlign[ctg]:
                for linkCtg, linkDirec in newLinkAlign[ctg][direc]:
                    fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(ctg, direc, linkCtg, linkDirec, *newLinkAlign[ctg][direc][(linkCtg, linkDirec)]))
    return passLinkAlign

def reciprocal_links(passLinkAlign):
    recLinks = {}
    for ctg1_1 in passLinkAlign:
        for direc1_1 in passLinkAlign[ctg1_1]:
            ctg1_2, direc1_2, ratio1 = passLinkAlign[ctg1_1][direc1_1]
            if ctg1_2 not in passLinkAlign:
                continue
            if direc1_2 not in passLinkAlign[ctg1_2]:
                continue
            ctg2_1, direc2_1, ratio2 = passLinkAlign[ctg1_2][direc1_2]
            if ctg2_1 == ctg1_1 and direc2_1 == direc1_1:
                #ctgSet = frozenset((ctg1_1, ctg1_2))
                if (ctg1_1, ctg1_2) in recLinks or (ctg1_2, ctg1_1) in recLinks:
                    continue
                recLinks[(ctg1_1, ctg1_2)] = ((direc1_1, direc1_2), (ratio1, ratio2))
    with open("reci.links", 'w') as fout:
        for ctg1, ctg2 in recLinks:
            direcLst, ratioLst = recLinks[(ctg1, ctg2)]
            fout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(ctg1,ctg2,*direcLst,*ratioLst ))
    #print(recLinks)
    return recLinks

def scaffold_Links(recLinks):
    ## check order
    ### statics link number
    nodeDic = {}
    for ctg1, ctg2 in recLinks:
        #for ctg in [ctg1, ctg2]:
        if ctg1 not in nodeDic:
            nodeDic[ctg1] = [ctg2]
        else:
            nodeDic[ctg1].append(ctg2)
        if ctg2 not in nodeDic:
            nodeDic[ctg2] = [ctg1]
        else:
            nodeDic[ctg2].append(ctg1)
    statNodeDic = {}
    for ctg in nodeDic:
        edgeCount = len(nodeDic[ctg])
        if edgeCount not in statNodeDic:
            statNodeDic[edgeCount] = []
        statNodeDic[edgeCount].append(ctg)
    ### threads nodes
    #addedNodes = {}
    singleNodes = statNodeDic[1]
    pathDic = {}
    for node in singleNodes:
        addedNodes = {}
        tmpNodeLst = [node, nodeDic[node][0]]
       # nextNodeLst = nodeDic[tmpNodeLst[-1]]
        for i in [node, nodeDic[node][0]]:
            addedNodes[i] = ""
        stopFlag = False
        loopFlag = False
        while stopFlag == False:
            nextNodeLst = nodeDic[tmpNodeLst[-1]]
            if len(nextNodeLst) == 1:
                stopFlag = True
            lastNode = tmpNodeLst[-2]
            if len(nextNodeLst) ==2:
                nextNode = nextNodeLst[0] if nextNodeLst[1] == lastNode else nextNodeLst[1]
                if nextNode in addedNodes:
                    stopFlag = True
                    loopFlag = True
                else:
                    tmpNodeLst.append(nextNode)
        if loopFlag == True:
            print("{} if loop, please check".format(node))
            continue
        else:
            pathLen = len(tmpNodeLst)
            if pathLen not in pathDic:
                pathDic[pathLen] = []
            pathDic[pathLen].append(tmpNodeLst)
    ### merge same path
    simplePathDic = {}
    for pathLen in pathDic:
        for path in pathDic[pathLen]:
            reversePath = tuple(reversed(path))
            if reversePath not in simplePathDic:
                simplePathDic[tuple(path)] = ""
    ### check oritation
    #print(simplePathDic)
    tourPath = []
    for path in simplePathDic:
        tmpPath = []
        for i in range(len(path)-1):
            ctg1, ctg2 = path[i], path[i+1]
            if (ctg1, ctg2) in recLinks:
                direc1, direc2 = recLinks[(ctg1, ctg2)][0]
            elif (ctg2, ctg1) in recLinks:
                direc2, direc1 = recLinks[(ctg2, ctg1)][0]
            else:
                pass
            if direc1 == "f2" and direc2 == "f1":
                ctg1 = "{}{}".format(ctg1, "+")
                ctg2 = "{}{}".format(ctg2, "+")
                #tmpPath.append("{}{}".format(ctg1, "+"))
                #tmpPath.append("{}{}".format(ctg2, "+"))
            elif direc1 == "f2" and direc2 == "f2":
                ctg1 = "{}{}".format(ctg1, "+")
                ctg2 = "{}{}".format(ctg2, "-")
                #tmpPath.append("{}{}".format(ctg1, "+"))
                #tmpPath.append("{}{}".format(ctg2, "-"))
            elif direc1 == "f1" and direc2 == "f1":
                ctg1 = "{}{}".format(ctg1, "-")
                ctg2 = "{}{}".format(ctg2, "+")
                #tmpPath.append("{}{}".format(ctg1, "-"))
                #tmpPath.append("{}{}".format(ctg2, "+"))
            else:
                ctg1 = "{}{}".format(ctg1, "-")
                ctg2 = "{}{}".format(ctg2, "-")
                #tmpPath.append("{}{}".format(ctg1, "-"))
                #tmpPath.append("{}{}".format(ctg2, "-"))
            if bool(tmpPath) == True:
                tmpPath.append(ctg2)
            else:
                tmpPath.append(ctg1)
                tmpPath.append(ctg2)
        tourPath.append(tmpPath)
    ### output
#    for i in range(len(tourPath)):
#        with open("Scaffold{}.tour".format(i), 'w') as fout:
#            fout.write(" ".join(tourPath[i]))

    Path("scaffold_tours").mkdir(exist_ok=True)
    with open("scaffold.clusters", 'w') as fout:
        for i in range(0, len(tourPath)):
            j = i + 1
        #with open("Scaffold{}.tour".format(i), 'w') as fout:
            fout.write("{}\t{}\t{}\n".format("scaffold{}".format(j), len(tourPath[i]), " ".join(tourPath[i])))
            with open(f"scaffold_tours/scaffold{j}.tour", "w") as out_tour:
                out_tour.write(f"{' '.join(tourPath[i])}")
                
    

        
def workflow(LisFile, win, cutoff):
    win=int(win)
    LisDic = read_LISFile(LisFile, "LIS")
    linkDIc = check_break_point(LisDic, win)
    passLinkAlign = check_linkage(linkDIc, cutoff)
    #print(passLinkAlign)
    recLinks = reciprocal_links(passLinkAlign)
    scaffold_Links(recLinks)


if __name__=="__main__":
    LisFile, win, cutoff = sys.argv[1:]
    workflow(LisFile, win, cutoff)
