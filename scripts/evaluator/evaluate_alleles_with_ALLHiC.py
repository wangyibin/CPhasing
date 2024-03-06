#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
comparison the `C-Phasing alleles` and `ALLHiC gmap allele`
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 
import numpy as np 
from pathlib import Path

from matplotlib_venn import venn2
import matplotlib.pyplot as plt


from cphasing.core import AlleleTable

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('allele1', 
            help='allele table from cphasing')
    pReq.add_argument('allele2',
            help="Allele table from ALLHiC gmap")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    prefix1 = Path(args.allele1.replace(".allele.table", "")).stem
    prefix2 = "ALLHiC"


    at1 = AlleleTable(args.allele1, fmt='allele2', sort=False)
    at1.data = at1.data[at1.data[1] < at1.data[2]]
    at2 = AlleleTable(args.allele2, fmt='allele1', sort=True)
    contig_pairs_1 = at1.contig_pairs
    contig_pairs_2 = at2.contig_pairs

    unique1 = contig_pairs_1 - contig_pairs_2
    unique2 = contig_pairs_2 - contig_pairs_1
    shared = contig_pairs_1 & contig_pairs_2

    venn2([contig_pairs_1, contig_pairs_2],
           ('C-Phasing', 'ALLHiC'),
           set_colors=("#b02418", "#253761"),
           subset_label_formatter=lambda x: f"{x:,}",
           alpha=0.3)
    
    plt.savefig(f"{prefix1}.{prefix2}.venn.png", dpi=600, bbox_inches='tight')
    plt.savefig(f"{prefix1}.{prefix2}.venn.pdf", dpi=600, bbox_inches='tight')

    total = sum((len(unique1), len(unique2), len(shared)))
    print(len(unique1) / total, len(shared) / total, len(unique2) / total, file=sys.stdout)
    
    with open(f"{prefix1}_unique.{prefix2}.stat", 'w') as out:
        print(len(unique1) / total, len(shared) / total, len(unique2) / total, file=out)

    with open(f"{prefix1}_unique.{prefix2}.tsv", 'w') as out:
        print("\n".join(map(lambda x: "\t".join(x), unique1)), file=out)

    with open(f"{prefix1}.{prefix2}_unique.tsv", 'w') as out:
        print("\n".join(map(lambda x: "\t".join(x), unique2)), file=out)

    with open(f"{prefix1}.{prefix2}.shared.tsv", 'w') as out:
        print("\n".join(map(lambda x: "\t".join(x), shared)), file=out)

if __name__ == "__main__":
    main(sys.argv[1:])