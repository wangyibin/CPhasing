#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
correct all alignments
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 

def is_edge(row, min_length=100):
    if (row[8] - row[10] < min_length) or (row[9] - 0 < min_length):
        return True 
    else:
        return False

def is_overlap(range1, range2, slope=10):
    start1, end1 = range1
    start2, end2 = range2

    return start1 - slope <= end2 and end1 + slope >= start2

def is_link(range1, range2, slope=50):
    start1, end1 = range1
    start2, end2 = range2

    return abs(start2 - end1) < slope 

def find_linked_read(df):
    df.sort_values(2, inplace=True)
    for i in range(len(df) - 1):
        row1 = df.iloc[i]
        range1 = (row1[2], row1[3])
        for j in range(i, len(df)):
            row2 = df.iloc[j]
            range2 = (row2[2], row2[3])
             

            if is_link(range1, range2): 
                linked_df = pd.concat([row1, row2], axis=1).T
                # linked_AS_per_base = (((linked_df[3] - linked_df[2]) / (linked_df[3] - linked_df[2]).sum()) * linked_df[17]).sum()
                linked_AS_per_base = linked_df[18].sum() / linked_df[10].sum()
             
                return row1.name, row2.name, linked_AS_per_base
            else:
                continue
            
            # if is_overlap(range1, range2):
            #     continue
    
    return None, None, None

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('paf', 
            help='')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
   
    df = pd.read_csv(args.paf, sep='\t', index_col=None, header=None, usecols=range(21))
    
    df[14] = df[14].str.split(":").map(lambda x:x[-1]).astype(int)
    df[17] = df[14] / df[10]
    df[18] = df[14]
    df[14] = df[14].map(lambda x: f"AS:i:{x}")
    
    for read, tmp_df in df.groupby(0, sort=False):
        if tmp_df[11].max() == 0:
            tmp_df.drop([17, 18], axis=1, inplace=True)
            tmp_df.to_csv(args.output, sep='\t', index=None, header=None)
            continue 

        alignments = []
        for read, row in tmp_df.iterrows():
            
            if row[16] == "tp:A:P":
                if alignments:
                    tmp_df2 = pd.concat(alignments, axis=1).T

                    if tmp_df2.apply(is_edge, axis=1).any():
                        i, j, linked_AS_per_base = find_linked_read(tmp_df)
                       
                        if tmp_df2[11].max() == 0:
                            tmp_df2.drop([17, 18], axis=1, inplace=True)
                        
                            tmp_df2.to_csv(args.output, sep='\t', index=None, header=None)
                            
                        else:
                            if tmp_df2[17].idxmax() == primary_idx:
                                tmp_df2.drop([17, 18], axis=1, inplace=True)
                                
                                tmp_df2.to_csv(args.output, sep='\t', index=None, header=None)
                            
                            else:
                                if linked_AS_per_base:
                             
                                    tmp_df2.loc[i, 17] = linked_AS_per_base
                                    tmp_df2.loc[j, 17] = linked_AS_per_base
                                    if linked_AS_per_base > tmp_df2.loc[primary_idx, 17]:
                                        tmp_df2.sort_values(17, ascending=False, inplace=True)
                                        tmp_df2.reset_index(inplace=True, drop=True)
                                        tmp_df2.iloc[0, 16] = "tp:A:P"
                                        tmp_df2.iloc[0, 11] = 2

                                        tmp_df2.loc[primary_idx, 16] = "tp:A:S"
                                        tmp_df2.loc[primary_idx, 11] = 0

                                tmp_df2.drop([17, 18], axis=1, inplace=True)
                                tmp_df2.to_csv(args.output, sep='\t', index=None, header=None)
                        
                    else:
                        tmp_df2.drop([17, 18], axis=1, inplace=True)
                        tmp_df2.to_csv(args.output, sep='\t', index=None, header=None)
                   

                primary_idx = row.name
                alignments = []
            
            alignments.append(row)
            
        else:
         
            if alignments:
                tmp_df2 = pd.concat(alignments, axis=1).T

                if tmp_df2.apply(is_edge, axis=1).any():
                    i, j, linked_AS_per_base = find_linked_read(tmp_df)
                 
                    if tmp_df2[11].max() == 0:
                        tmp_df2.drop([17, 18], axis=1, inplace=True)
                    
                        tmp_df2.to_csv(args.output, sep='\t', index=None, header=None)
                
                    else:
                    
                        if tmp_df2[17].idxmax() == primary_idx:
                            tmp_df2.drop([17, 18], axis=1, inplace=True)
                            
                            tmp_df2.to_csv(args.output, sep='\t', index=None, header=None)
                           
                        else:
                            if linked_AS_per_base:
                             
                                tmp_df2.loc[i, 17] = linked_AS_per_base
                                tmp_df2.loc[j, 17] = linked_AS_per_base

                                if linked_AS_per_base > tmp_df2.loc[primary_idx, 17]:  
                                    tmp_df2.sort_values(17, ascending=False, inplace=True)
                                    tmp_df2.iloc[0, 16] = "tp:A:P"
                                    tmp_df2.iloc[0, 11] = 2
                                    tmp_df2.loc[primary_idx, 16] = "tp:A:S"
                                    tmp_df2.loc[primary_idx, 11] = 0

                            tmp_df2.drop([17, 18], axis=1, inplace=True)
                            tmp_df2.to_csv(args.output, sep='\t', index=None, header=None)
        
                else:
                    tmp_df2.drop([17, 18], axis=1, inplace=True)
                    tmp_df2.to_csv(args.output, sep='\t', index=None, header=None)
                
                
 


if __name__ == "__main__":
    main(sys.argv[1:])