#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 

from Bio import SeqIO
from io import StringIO
from pyfaidx import Fasta, FastaRecord
from pysam import VariantFile
from subprocess import PIPE, Popen, DEVNULL

from cphasing.utilities import get_contigs

def generate_consensus2(vcf_file, ranges, seq):
    vcf = VariantFile(vcf_file)
    print(ranges)
    for rec in vcf.fetch(*ranges):
        print(rec.pos)


def generate_consensus(vcf_file, seq_record):
    command = [
        "bcftools",
        "consensus",
        vcf_file
    ]


    process = Popen(command,
                  stdin=PIPE, 
                  stdout=PIPE,
                  stderr=DEVNULL,
                  text=True)
    process.stdin.write(seq_record)
    process.stdin.close()

    output_data = process.stdout.read()

    return_code = process.wait()

    if return_code == 0:
        return output_data
    else:
        return seq_record

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('fasta', 
            help='ref fasta')
    # pReq.add_argument('vcf',
    #         help="vcf file, compressed and indexed")
    pReq.add_argument("bed",
            help="porec alignments with read idx in four columns")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    fasta = args.fasta
    # vcf = args.vcf 
    bed = args.bed
    
    seq_db = Fasta(fasta)
  
    bed_df = pd.read_csv(bed, sep='\t', index_col=None, header=None)
    bed_df.columns = ['chrom', 'start', 'end', 'read_idx']
    assert len(bed_df.columns) == 4, "bed file should be four columns"

    prev_read_idx = None 
    fragments = []
    names = []
    for i, row in bed_df.iterrows():
        seq_id = f"{row.chrom}:{row.start}-{row.end}"
        try:
            seq = seq_db[row.chrom][row.start:row.end]
        except KeyError:
            continue
        # generate_consensus2(vcf, (row.chrom, row.start, row.end), seq)
        # output_data = generate_consensus(vcf, f'>{seq_id}\n{seq}')
        # seq_record = SeqIO.parse(StringIO(output_data), "fasta")
        if row.read_idx != prev_read_idx and prev_read_idx is not None:
            output_seq = "".join(fragments)
            output_names = " ".join(names)
            print(f">{prev_read_idx} {output_names}\n{output_seq}", file=sys.stdout)
            fragments = []
            names = []
        fragments.append(str(seq))
        names.append(seq_id)
        
        prev_read_idx = row.read_idx
        
        



if __name__ == "__main__":
    main(sys.argv[1:])