#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
seperate four contig to two groups
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 

import numpy as np
import igraph as ig 

import re 



def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('contigs', nargs="*",
            help='four contigs')
    pOpt.add_argument("-c", "--contacts", required=True,
                      help="contacts file")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    contacts = pd.read_csv(args.contacts, sep='\t', 
                           index_col=None, header=None,
                           names=["contig1", "contig2", "contact"])
    contacts2 = contacts.copy()
    contacts2.columns = ['contig2', 'contig1', 'contact']
    contacts = pd.concat([contacts, contacts2], axis=0)
    contacts_db = contacts.set_index(['contig1', 'contig2']).to_dict()['contact']

    
    if len(args.contigs) != 4:
        print("Need four contigs", file=sys.stderr)
        sys.exit(-1)
    
    group1 = args.contigs[:2]
    group2 = args.contigs[2:]
    scores = []
    contig_edges = []
    for contig1 in group1:
        for contig2 in group2:
            contig_edge = (contig1, contig2)
            score = contacts_db[contig_edge]
            scores.append(score)
            contig_edges.append(contig_edge)
    
    edges =  [(0, 2), (0, 3), (1, 2), (1, 3)]
    g = ig.Graph.Bipartite([0, 0, 1, 1], edges)
    g.es['weight'] = scores
    matching = g.maximum_bipartite_matching(weights='weight')

    contig_pairs = [] 
    for i, j in enumerate(matching.matching[:2]):
          
        contig_pairs.append((args.contigs[i], args.contigs[j]))
        
    print(contig_pairs)






if __name__ == "__main__":
    main(sys.argv[1:])