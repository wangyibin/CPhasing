import argparse
import logging
import os
import os.path as op
import sys

import re 

import gffpandas.gffpandas as gffpd
from collections import defaultdict
from pyfaidx import Fasta
from cphasing.utilities import read_fasta


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('gff', 
            help='')
    pReq.add_argument('fasta')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    gff_df = gffpd.read_gff3(args.gff)
    gff_df = gff_df.filter_feature_of_type(['gene']).df

    fasta = Fasta(args.fasta)

    gff_df = gff_df[gff_df['seq_id'].str.match(r'utg.*_d(\d+)')]
    dup_genes = gff_df['attributes'].apply(lambda x: re.search(r'ID=(.*?);', x).group(1))
 
    dup_raw_genes = list(map(lambda x: x.rsplit("_", 1)[0], dup_genes))
    dup_genes_db = defaultdict(list)
    for i, gene in enumerate(dup_genes):
        dup_genes_db[dup_raw_genes[i]].append(gene)
 
    for gene in fasta.keys():
        print(f'>{gene}', file=args.output)
        print(fasta[gene][:].seq, file=args.output)
    
    for genes in dup_genes_db:
        for gene in dup_genes_db[genes]:
            raw_gene = gene.rsplit("_", 1)[0]
            print(f'>{gene}', file=args.output)
            print(fasta[raw_gene][:].seq, file=args.output)


if __name__ == "__main__":
    main(sys.argv[1:])