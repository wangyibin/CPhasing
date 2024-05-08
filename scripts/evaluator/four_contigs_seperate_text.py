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

from pathlib import Path
from rich.console import Console
console = Console()

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    # pReq.add_argument('contigs', 
    #         help='contigs text file')
    pOpt.add_argument("-c", "--contacts", required=True,
                      help="contacts file")
    pOpt.add_argument("-n", "--normalize", action='store_true',
                      default=False, 
                      help="Normalized the contacts")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    contacts = pd.read_csv(args.contacts, sep='\t', 
                           index_col=None, header=None,
                           names=["contig1", "contig2", "contact"])
    if args.normalize:
        intra_contacts = contacts[contacts['contig1'] == contacts['contig2']]
        intra_contacts_db = intra_contacts[['contig1', 'contact']].set_index('contig1').to_dict()['contact']
        contacts['contact1'] = contacts['contig1'].map(intra_contacts_db.get)
        contacts['contact2'] = contacts['contig2'].map(intra_contacts_db.get)
        contacts['contact'] = contacts['contact'] / np.sqrt(contacts['contact1'] * contacts['contact2'])
        contacts.drop(['contact1', 'contact2'], axis=1, inplace=True)

    contacts2 = contacts.copy()
    contacts2.columns = ['contig2', 'contig1', 'contact']
    contacts = pd.concat([contacts, contacts2], axis=0)

    contacts_db = contacts.set_index(['contig1', 'contig2']).to_dict()['contact']

    contigs = "start"
    while contigs != "exit" or contigs != "":
        contigs = input("Input contigs: ")
        if contigs == "exit":
            break 
        regex = re.compile('(utg\d+l)')
        contigs = regex.findall(contigs)
    

        if len(contigs) != 3 and len(contigs) != 4 and len(contigs) != 5:
            print("Need 3 or 4 or 5 contigs", file=sys.stderr)
            continue
        if len(contigs) == 5:
            contigs.pop(2)
        
        if len(contigs) == 4:
            group1 = contigs[:2]
            group2 = contigs[2:]
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
            
            if matching.match_of(0) == 2:

                console.print("Yes", style="magenta")
            else:
                console.print("No", style="red")
            contig_pairs = [] 
            for i, j in enumerate(matching.matching[:2]):
                
                contig_pairs.append((contigs[i], contigs[j]))
                
            console.print_json(data=contig_pairs)
        
        if len(contigs) == 3:
            contig1 = contigs[0]

            scores = []
            for contig2 in contigs[1:]:
                contig_edge = (contig1, contig2)
                score = contacts_db[contig_edge]
                scores.append(score)
            
            console.print(contigs[np.argmax(scores) + 1])
            console.print_json(data=dict(zip(contigs[1:], scores)))





if __name__ == "__main__":
    main(sys.argv[1:])