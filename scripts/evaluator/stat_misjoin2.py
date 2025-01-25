#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
stat misjoin points from paftools.js misjoin -e 
"""

import argparse
import logging
import os
import os.path as op
import sys
import io

import pandas as pd
import numpy as np

from pathlib import Path
from pyranges import PyRanges
from shutil import which
from subprocess import Popen, PIPE

from cphasing.utilities import cmd_exists, list_flatten, calculate_Nx_from_contigsizes
from cphasing.agp import import_agp

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('tsv', 
            help='tsv from paftools.js misjoin -e ')
    pReq.add_argument("agp", help="agp file")
    pOpt.add_argument('-m', '--min_count', default=2, type=int, 
            help='minimum support of read count [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    flank_length = 5000

    cmd = f"grep -v '^C' {args.tsv}"
    process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
 
    csv = io.StringIO(stdout.decode())

    try:
        df = pd.read_csv(csv, sep='\t', header=None, index_col=None,
                        comment="#", engine="python", usecols=list(range(12)))
    except pd.errors.EmptyDataError:
        print(f"No. misassemblies: {0}", file=sys.stderr)
        return 

    df = df[df[0] != "M"]
    df = df.sort_values([1, 3])
    candicate_res = []
    for read, tmp_df in df.groupby(1, sort=False):
        for i in range(0, len(tmp_df), 2):
            tmp_df2 = tmp_df.iloc[i:i+2]


            if (tmp_df2[5] == '-').all():
                tmp_df2 = tmp_df2.sort_values(3, ascending=False)
        
        
            if tmp_df2.iloc[0][6] > tmp_df2.iloc[1][6]:
                res = (
                    tmp_df2.iloc[1][6], tmp_df2.iloc[1][8], tmp_df2.iloc[1][8] + flank_length,
                    tmp_df2.iloc[0][6], tmp_df2.iloc[0][9] - flank_length, tmp_df2.iloc[0][9])
            else:
                res = (tmp_df2.iloc[0][6], tmp_df2.iloc[0][9] - flank_length, tmp_df2.iloc[0][9], 
                    tmp_df2.iloc[1][6], tmp_df2.iloc[1][8], tmp_df2.iloc[1][8] + flank_length)
            candicate_res.append(res)
            # candicate_res.append(res2)

    res_df = pd.DataFrame(candicate_res, columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'])                
    res_df.reset_index(inplace=True)
    
    res_df1 = res_df[['chrom1', 'start1', 'end1', 'index']]
    res_df1.columns = ["Chromosome", "Start", "End", "Index"]

    res_df2 = res_df[['chrom2', 'start2', 'end2', 'index']]
    res_df2.columns = ["Chromosome", "Start", "End", "Index"]
    
    
    # res_df = PyRanges(res_df).merge(count=True).df
    # res_df = res_df.sort_values(['Chromosome', 'Start', 'End'])
    # res_df = res_df[res_df['Count'] >= args.min_count]
    output = Path(args.tsv).stem
    # res_df.to_csv(f'.{output}.misassembly.txt', sep='\t', header=None, index=None)

    agp_df, gap_df = import_agp(args.agp)
    
    _agp_df = agp_df.reset_index()

    contig_size = _agp_df[['id', 'tig_end']].set_index('id')
 
    tig_df = _agp_df[_agp_df['chrom'] == _agp_df['id']]
    unanchored_contigs = set(tig_df['chrom'].values.tolist())
  
    chrom_sizes = agp_df.groupby(agp_df.index)['end'].max()
    chrom_sizes.to_csv(f".{output}.sizes", sep='\t', header=None, index=True)
    
    def flanking(row):
        row2 = row.copy()
        row['start'] = row['start'] - 1
        row['end'] = row['start'] + flank_length

        if row['orientation'] == "-":
            row['id'] = f"{row['id']}_R"
        else:

            row['id'] = f"{row['id']}_L"

        
        row2['start'] = row2['end'] - flank_length
        if row2['orientation'] == "-":
            row2['id'] = f"{row2['id']}_L"
        else:
            row2['id'] = f"{row2['id']}_R"
        
        row = row.rename(index = {"start": "start1", "end": "end1", "id": "id1"})
        row2 = row2.rename(index = {"start": "start2", "end": "end2", "id": "id2"})

        return pd.concat([row[['start1', 'end1', 'id1']], row2[['start2', 'end2', 'id2']]], axis=0) 
        
        
        
    target_df = agp_df.apply(flanking, axis=1)

    target_df = pd.concat([target_df[['start1', 'end1', 'id1']].rename(columns={"start1": "start",
                                                                                "end1": "end",
                                                                                "id1": "id"}), 
                           target_df[['start2', 'end2', 'id2']].rename(columns={"start2": "start",
                                                                                "end2": "end",
                                                                                "id2": "id"})], axis=0)
    target_df.reset_index(inplace=True)
    
    target_df.columns = ["Chromosome", "Start", "End", "ID"]
    target_gr = PyRanges(target_df)

    gr1 = PyRanges(res_df1)
    gr2 = PyRanges(res_df2)
    gr1 = gr1.join(target_gr, report_overlap=True).df[['Chromosome', 'Start_b', "End_b", "Index", "ID", "Overlap"]]
    gr2 = gr2.join(target_gr, report_overlap=True).df[['Chromosome', 'Start_b', "End_b", "Index", "ID", "Overlap"]]

    overlap = flank_length // 2

    
    gr1.columns = ['chrom1', 'start1', 'end1', 'index', 'id1', 'overlap']
    gr2.columns = ['chrom2', 'start2', 'end2', 'index', 'id2', 'overlap']
    gr1 = gr1[gr1['overlap'] > overlap]
    gr2 = gr2[gr2['overlap'] > overlap]


    gr1.set_index('index', inplace=True)
    gr2.set_index('index', inplace=True)

    res_df = pd.concat([gr1, gr2], axis=1).dropna().drop('overlap', axis=1)
    # res_df = res_df.drop_duplicates()

    res_df = res_df[res_df.duplicated()].drop_duplicates()

    unanchor_df = res_df[res_df['chrom1'].isin(unanchored_contigs) | res_df['chrom2'].isin(unanchored_contigs)]
    un_unanchor_df = unanchor_df[unanchor_df['chrom1'].isin(unanchored_contigs) & unanchor_df['chrom2'].isin(unanchored_contigs)]
    unanchor_df = unanchor_df[~(unanchor_df['chrom1'].isin(unanchored_contigs) & unanchor_df['chrom2'].isin(unanchored_contigs))]

    reported_unanchor_contigs = list_flatten([unanchor_df['chrom1'].values.tolist(), unanchor_df['chrom2'].values.tolist()])
    reported_unanchor_contigs = list(filter(lambda x: x in unanchored_contigs, reported_unanchor_contigs))

    reported_un_unanchor_contigs = list_flatten([un_unanchor_df['chrom1'].values.tolist(), un_unanchor_df['chrom2'].values.tolist()])
    reported_un_unanchor_contigs = list(filter(lambda x: x in unanchored_contigs, reported_un_unanchor_contigs))


    res_df = res_df[~(res_df['chrom1'].isin(unanchored_contigs) | res_df['chrom2'].isin(unanchored_contigs))]

    res_df = res_df.astype({"start1": np.int64, "end1": np.int64, "start2": np.int64, "end2": np.int64})
    res_df.to_csv(f"{output}.misassembly.txt", sep='\t', header=None, index=None)

    res_gap_df = pd.concat([res_df[['chrom1', 'start1', 'end1', 'id1']].rename(
                        columns={
                            "chrom1": "chrom",
                            "start1": "start",
                            "end1": "end",
                            "id1": "id"
                        }
                    ), 
                     res_df[['chrom2', 'start2', 'end2', 'id2']].rename(
                        columns={
                            "chrom2": "chrom",
                            "start2": "start",
                            "end2": "end",
                            "id2": "id"
                        }
                    )], axis=0).drop_duplicates()
    

    contig_size = contig_size.rename(columns={"tig_end": "length"})
    misassembly_contig_n50 = calculate_Nx_from_contigsizes(contig_size.reindex(res_gap_df['id'].str.split("_", n=1).map(lambda x: x[0])), 50)
   
    print(f"No. of misassembly pairs: {int(len(res_df))}", file=sys.stderr)
    print(f"No. of misassembly contig N50: {misassembly_contig_n50}", file=sys.stderr)
    print(f"No. of misassembly gaps: {int(len(res_gap_df))}", file=sys.stderr)
    print(f"No. of gaps: {int(len(gap_df))}", file=sys.stderr)
    print(f"Proportion of false gaps: {int(len(res_gap_df))/int(len(gap_df)):.4f}", file=sys.stderr)
    print(f"No. of supported single unanchor contigs: {int(len(reported_unanchor_contigs))}", file=sys.stderr)
    print(f"Size. of supported single unanchor contigs: {contig_size.reindex(reported_unanchor_contigs).sum().values[0]}", file=sys.stderr)
    print(f"No. of supported unanchor paired contigs: {int(len(reported_un_unanchor_contigs))}", file=sys.stderr)
    print(f"Size. of supported unanchor paired contigs: {contig_size.reindex(reported_un_unanchor_contigs).sum().values[0]} ", file=sys.stderr)
    print(f"No. of total supported unanchor contigs: {int(len(reported_un_unanchor_contigs)) + int(len(reported_unanchor_contigs))}", file=sys.stderr)
    print(f"Size. of total supported unanchor contigs: {contig_size.reindex(reported_unanchor_contigs).sum().values[0] + contig_size.reindex(reported_un_unanchor_contigs).sum().values[0]}", file=sys.stderr)
    print(f"No. of total unancored contigs: {int(len(unanchored_contigs))}", file=sys.stderr)
    print(f"Size. of total unancored contigs: {contig_size.reindex(unanchored_contigs).sum().values[0]}", file=sys.stderr)

if __name__ == "__main__":
    main(sys.argv[1:])