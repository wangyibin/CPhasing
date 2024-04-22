#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
estimate the memory of hypergraph
"""

import argparse
import logging
import os
import os.path as op
import sys
import msgspec
from cphasing.hyperpartition import HyperPartition
from cphasing.algorithms.hypergraph import HyperEdges
from cphasing.algorithms.hypergraph import HyperGraph

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('hypergraph', 
            help='')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    hypergraph = msgspec.msgpack.decode(open(args.hypergraph, 'rb').read(), type=HyperEdges)
    hypergraph.to_numpy()
    HG = HyperGraph(hypergraph)
    H = HG.incidence_matrix()


if __name__ == "__main__":
    main(sys.argv[1:])    
