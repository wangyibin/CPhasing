#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot dotplot
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from natsort import natsort_keygen
from subprocess import Popen, PIPE
from pathlib import Path 
from pyfaidx import Fasta
from cphasing.utilities import run_cmd 


PAF_HADER = ["contig1", "length1", "start1", "end1", "strand",
                 "contig2", "length2", "start2", "end2", "mismatch",
                 "matches", "mapq", "identity"]


def align(ref, qry, threads=20, paf='wfmash.paf'):
    Path("logs").mkdir(parents=True, exist_ok=True)
    cmd = ["wfmash", ref, qry, 
                    "-m", 
                    "-4", 
                     "-t", str(threads)]
    pipelines = []
    try:
        pipelines.append(
            Popen(cmd, stdout=open(paf, "w"),
                    stderr=open(f"logs/ref.align.log", "w"),
                    bufsize=-1)
        )
        pipelines[-1].wait()
    except:
        raise Exception('Failed to execute command:' 
                            f'\t{cmd}.')
    finally:
        for p in pipelines:
            if p.poll() is None:
                p.terminate()
            else:
                assert pipelines != [], \
                    "Failed to execute command, please check log."
            if p.returncode != 0:
                raise Exception('Failed to execute command:' 
                                    f'\t{" ".join(cmd)}.')

    return paf


def align2(ref, qry, threads=20, paf='wfmash.paf2'):
    Path("logs").mkdir(parents=True, exist_ok=True)
    cmd = ["wfmash", ref, qry, 
                    "-m", "-s", "50k", "-l", "250k", 
                    "-4", "-N",
                     "-t", str(threads)]
    pipelines = []
    try:
        pipelines.append(
            Popen(cmd, stdout=open(paf, "w"),
                    stderr=open(f"logs/ref.align.log", "w"),
                    bufsize=-1)
        )
        pipelines[-1].wait()
    except:
        raise Exception('Failed to execute command:' 
                            f'\t{cmd}.')
    finally:
        for p in pipelines:
            if p.poll() is None:
                p.terminate()
            else:
                assert pipelines != [], \
                    "Failed to execute command, please check log."
            if p.returncode != 0:
                raise Exception('Failed to execute command:' 
                                    f'\t{" ".join(cmd)}.')

    return paf


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('ref', 
            help='ref')
    pReq.add_argument('qry', help='qry')
    pOpt.add_argument('--hap', help="haplotype resolved assembly containing the trio information")
    pOpt.add_argument('-o', '--output', type=str,
            default=None, help='output file [default: %(default)s]')
    
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    if not Path("wfmash.paf").exists():
        paf = align(args.ref, args.qry)
    else:
        paf = "wfmash.paf"

   

    
    df = pd.read_csv(paf, sep='\t', header=None, usecols=range(13),
                         names=PAF_HADER, index_col=None)
    
    df['identity'] = df['identity'].map(lambda x: x.replace("id:f:", "")).astype(float)
    df = df.loc[df['identity'] >= 0.90]
    df = df.sort_values(by=['contig2', 'start2'], key=natsort_keygen())

    df = df.loc[~df['contig1'].str.startswith('utg')]
    df = df.loc[~df['contig1'].str.startswith('utig')]
    df = df.loc[~df['contig1'].str.startswith('ctg')]
    df = df.loc[~df['contig1'].str.startswith('contig')]
    
    df['align_length2'] = df['end2'] - df['start2']
    df['align_length1'] = df['end1'] - df['start1']   
    
    ref_chromsizes = df.groupby('contig2', sort=False)['length2'].unique().map(lambda x: x[0])
    
    ref_cumsum = ref_chromsizes.copy()
    ref_cumsum = np.cumsum(np.r_[[0], ref_cumsum.values])
    
    ref_chromsizes_db = dict(zip(ref_chromsizes.index, ref_cumsum[:-1]))  

    df2 = df.sort_values(by=['contig1', 'start2'], key=natsort_keygen())
    
    qry_chromsizes = df2.groupby('contig1', sort=False)['length1'].unique().map(lambda x:x[0])
    qry_cumsum = qry_chromsizes.copy()
    qry_cumsum = np.cumsum(np.r_[[0], qry_cumsum.values])
    qry_chromsizes_db = dict(zip(qry_chromsizes.index, qry_cumsum[:-1]))

    if args.hap:
        with open("tmp.contigs.fasta", 'w') as out:
            fasta = Fasta(args.qry)
            for i, row in df.iterrows():
                print(f">{i}\n{fasta[row.contig1][row.start1:row.end1]}", file=out)

        paf2 = align2(args.hap, "tmp.contigs.fasta")
        paf2 = pd.read_csv("wfmash.paf2", sep='\t', header=None, usecols=range(13),
                         names=PAF_HADER, index_col=None)
        
        paf2['trio'] = paf2['contig2'].str.split("_").map(lambda x: x[1])
        paf2['trio'] = paf2['trio'].map(lambda x: 1 if x== 'PATERNAL' else 0)
        trio_info = paf2[['contig1', 'trio']]
        trio_info = trio_info.set_index('contig1').to_dict()['trio']
 
    
    def func(row):
        v = qry_chromsizes_db[row.contig1]
        row['start1'] += v
        row['end1'] += v

        v = ref_chromsizes_db[row.contig2]
        row['start2'] += v
        row['end2'] += v

        return row 

    df = df.apply(func, axis=1)
    
    line_list = []

    color_list = []
    for i, row in df.iterrows():
        if row.strand == "+":
            line_list.append([(row.start1, row.end1), (row.start2, row.end2)])
        else:
            line_list.append([(row.end1, row.start1), (row.end2, row.start2)])

        if args.hap:
            try:
                color_list.append(trio_info[row.name])
            except:
                color_list.append(2)
            

    plt.rcParams['font.family'] = 'Arial'

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)

    color_db = ["#FBB4AE", "#B3CDE3", "#bcbcbc"]
    for i, line in enumerate(line_list):
        try:
            color_idx = color_list[i]
        except:
            color_idx = 2
        plt.plot(line[0], line[1], color=color_db[color_idx],)
    

    ref_half_pos = ref_cumsum[:-1] + (ref_cumsum[1:] - ref_cumsum[:-1])//2
    qry_half_pos = qry_cumsum[:-1] + (qry_cumsum[1:] - qry_cumsum[:-1])//2
    plt.xticks(qry_half_pos, qry_chromsizes.index, fontsize=6, rotation=90)

    plt.yticks(ref_half_pos, ref_chromsizes.index, fontsize=6)

    plt.vlines(qry_cumsum, ax.get_xlim()[0], ax.get_xlim()[1],
                 color="#bcbcbc", linewidth=0.5, alpha=0.3)
    plt.hlines(ref_cumsum, ax.get_ylim()[0], ax.get_ylim()[1],
                 color="#bcbcbc", linewidth=0.5, alpha=0.3)
    plt.xlim(0, max(qry_cumsum))
    plt.ylim(0, max(ref_cumsum))
    plt.xlabel("Query", fontsize=12)
    plt.ylabel("Reference", fontsize=12)

    output = Path(paf).stem if not args.output else args.output
    plt.savefig(f"{output}.dotplot.png", dpi=600, bbox_inches='tight')
    plt.savefig(f"{output}.dotplot.pdf", dpi=300, bbox_inches='tight')
        

if __name__ == "__main__":
    main(sys.argv[1:])
