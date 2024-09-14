#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
simulate switch error from vcf, which vcf is a simulated vcf from a monoploid reference
"""

import argparse
import logging
import os
import os.path as op
import sys
import random
import pandas as pd
import numpy as np



from collections import OrderedDict

from cphasing.utilities import get_genome_size, run_cmd, cmd_exists

def bgzip(vcf):
    cmd1 = ["bgzip", "-f", vcf]

    cmd2 = ["tabix", "-p", "vcf", f"{vcf}.gz"]
    run_cmd(cmd1)
    run_cmd(cmd2)

def vcf2fasta(vcf, fasta, rename=None):
    output_prefix = vcf.replace(".vcf.gz", "")
    cmd = ["bcftools", "consensus", vcf, "-f", fasta, "-o", f"{output_prefix}.fa"]

    run_cmd(cmd)

    if rename:
        cmd = ["seqkit", "replace", "-p", "hap1", "-r", "hap2", 
               f"{output_prefix}.fasta",
               "-o", f".{output_prefix}.fa"]
        run_cmd(cmd)
        
        os.rename(f".{output_prefix}.fa", f"{output_prefix}.fa")
    

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('vcf', 
            help='')
    pReq.add_argument('ref_fasta')
    pOpt.add_argument('-p', '--rate', help='switch error', type=float, default=0.01)
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    if not cmd_exists('bgzip'):
        print("No such command of bgzip", file=sys.stderr)
        sys.exit(-1)
    if not cmd_exists('tabix'):
        print("No such command of tabix", file=sys.stderr)
        sys.exit(-1)
    
    if not cmd_exists('bcftools'):
        print("No such command of bcftools", file=sys.stderr)
        sys.exit(-1)


    header = []
    with open(args.vcf, 'r') as fp:
        for line in fp:
            if not line.strip():
                continue

            if line.startswith("#"):
                header.append(line.strip())


    vcf_df = pd.read_csv(args.vcf, sep='\t', header=None, index_col=None, comment="#")
    
    snp_counts = len(vcf_df)

    genomesize = get_genome_size(args.ref_fasta)
    
    sw = int(np.ceil(args.rate / 100 * genomesize ))
  
    assert sw < snp_counts, "switch rate too high"
    sample_df = vcf_df.sample(n=sw, random_state=12345)
    sample_df.sort_values([0, 1], inplace=True)
    # sample_df.columns = [0, 1, 2, 4, 3, 5, 6, 7]
    
    vcf1 = sample_df
    vcf2 = vcf_df.loc[~vcf_df.index.isin(sample_df.index)]
    
    with open("1.vcf", 'w') as out:
        print("\n".join(header), file=out)
        vcf1.to_csv(out, sep='\t', header=None, index=None)

    bgzip("1.vcf")
    vcf2fasta("1.vcf.gz", args.ref_fasta)

    with open("2.vcf", 'w') as out:
        print("\n".join(header), file=out)
        vcf2.to_csv(out, sep='\t', header=None, index=None)

    bgzip("2.vcf")
    vcf2fasta("2.vcf.gz", args.ref_fasta, rename=True)

    






if __name__ == "__main__":
    main(sys.argv[1:])