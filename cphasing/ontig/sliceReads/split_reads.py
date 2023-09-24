import os
import sys
import argparse
import multiprocessing


"""
0. fastp split-window qc
1. read fastq from gz
2. store chunk into cache
3. split fastq with window & write file
"""


def split_fastq(fastq, part, threads):
    split_cmd="seqkit split2 -p {} -j {} -O split_fastq {}".format(part, threads, fastq)
    print("{} will be split to {} parts, which stored in split_fastq directory.".format(fastq, part))
    os.system(split_cmd)

def splitFa(fa, qual, win):
    # 13450 5000
    # [0, 5000, 10000, 15000]
    totalLen = len(fa)
    binLst = list(range(0,totalLen,win))
    faLst = []
    qualLst = []
    for i in range(len(binLst) - 1):
        subFa = fa[binLst[i] : binLst[i+1]]
        subQual = qual[binLst[i] : binLst[i+1]]
        faLst.append(subFa)
        qualLst.append(subQual)
    faLst.append(fa[binLst[-1] :])
    qualLst.append(qual[binLst[-1] :])
    return faLst, qualLst


def trans_winSize(win):
    if int(win) < 1000:
        return win
    elif int(win)/1000 >= 1 and int(win)/1000 < 1000:
        return int(win)/1000 + "k"
    else:
        return int(win)/1000000 + "m"

def slide_reads_by_window(reads, outpre, win):
    if not os.path.exists(reads):
        raise Exception("no such file %s" % reads)
    affix = trans_winSize(win)
    with open("{}_{}.fastq".format(outpre, affix), 'w') as fout:
        #if reads.endwith("gz"):
        rf = os.popen("pigz -d -c {}".format(reads))
        lineCount = 0
        oneRead = []
        writeMode = False
        for line in rf:
            lineCount += 1
        #line.rstrip()
            oneRead.append(line.rstrip())
            if lineCount == 4 and writeMode == False:
                lineCount = 0
                writeMode = True
            else:
                continue
            if writeMode == True:
                faName, fa, faString, qual = oneRead
                faLst, qualLst = splitFa(fa, qual, win)
                #qualLst = splitFa(qual, win)
                faName = faName.split(' ')
                if len(faLst[1]) >= 5:
                    for faInd in range(len(faLst)):
                        fout.write("{}_{}\n{}\n{}\n{}\n".format(faName, faInd, faLst[faInd], faString, qualLst[faInd]))
                oneRead = []
                writeMode = False

def pipe(fastqFile, win, part, threads):
    split_fastq(fastqFile, part, threads)
    parDir = "{}/split_fastq".format(os.getcwd())
    fastqLst = os.listdir(parDir)
    fastqLst = list(map(lambda x: parDir + '/' + x, fastqLst))
    oupreLst = list(map(lambda x: x.replace(".fastq.gz", ""), fastqLst)) 
    pool = multiprocessing.Pool(processes=int(threads))
   # time_print("parsing %s" % id)
    for idx in range(len(fastqLst)):
        fq = fastqLst[idx]
        op = oupreLst[idx]
        #slide_reads_by_window(fq, op, win)
        #read_minimizer(fn_prefix, seq, ks, a.shape, sm.name)
        pool.apply_async(slide_reads_by_window, (fq, op, win,))
    pool.close()
    pool.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is the script for split fastq file into window reads. \
                                     For accelerate this process, we split fastq into N parts at first, \
                                     then sliding fastq file with window.")
    parser.add_argument('-i', '--fastq', default=None,
                        help='<filepath>  fastq file.')
    parser.add_argument('-w', '--window', default=5000,
                        help='<int> window size, default is 5000.')
    parser.add_argument('-p', '--npart', default=10,
                        help='<int> number of part, default is 10')
    parser.add_argument('-t', '--threads', default=10,
                        help='<int> number of threads, default is 10')
    args = parser.parse_args()
    # fastqFile, win, part, threads
    pipe(args.fastq, int(args.window), int(args.npart), int(args.threads))
