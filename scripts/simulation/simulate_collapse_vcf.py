#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
simulate collapsed contigs from vcf 
"""

import argparse
import logging
import os
import os.path as op
import sys

import random 
import gzip
import numpy as np 
import pandas as pd 

from collections import defaultdict
from ncls import NCLS
from pathlib import Path

from cphasing.utilities import run_cmd, read_fasta

def read_bed(bed):
    df = pd.read_csv(bed, sep='\t', header=None, index_col=None)
    
    return df 

def vcf2contigs(vcf, bed_df, suffix):
    bed_df = bed_df.copy()
    bed_df[0] = bed_df[0].astype(str)
    bed_df[3] = bed_df[3].str.split("ctg").map(lambda x: x[-1])
    bed_df[3] = bed_df[3].astype(int)
    ncls_db = {}
    for i, tmp_df in bed_df.groupby(0):
        if i not in ncls_db:
            ncls_db[i] = NCLS(tmp_df[1].values, tmp_df[2].values, tmp_df[3].values)
    
    output = f"{Path(vcf).stem.rsplit('.', 1)[0]}.contigs.vcf"
    with gzip.open(vcf, 'r') as fp:
        with open(output, 'w') as out:
            for line in fp:
                line = line.decode('utf-8')
                if line[0] == "#":
                    print(line.strip(), file=out)
                else:
                    line_list = line.strip().split()
                    chrom, pos = line_list[:2]
                    pos = int(pos)
                    it = ncls_db[chrom].find_overlap(pos-1, pos)
                    if not it:
                        continue 

                    i = it.__next__()
                    chrom = f"{chrom}{suffix}.ctg{i[-1]}"
                    pos = pos - i[0]

                    line_list[0] = chrom
                    line_list[1] = str(pos)
                    print("\t".join(line_list), file=out)

    return output

        

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pOpt.add_argument('-f', '--fasta', help="reference fasta", required=True)
    pOpt.add_argument('--vcf', nargs="+", required=True,
            help='')
    pOpt.add_argument('--suffix', nargs="+", required=True,
                      help="suffix of the chromosome in each sample, must equal to vcf")
    pOpt.add_argument('-n', '--n50', default="500k")
    pOpt.add_argument('--ratio', type=float, default=0.2,
            help="the ratio of collapsed contigs. [default: %(default)s]")
    pOpt.add_argument('--seed', type=int, default=1212,
            help="random seed of program. [default: %(default)s]")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    print("[Warning]: Deprected", file=sys.stderr)
    random.seed(args.seed)
    simu_ctg_cmd = ["simuCTG.py", "-n", args.n50, "-i", args.fasta, "-o", f"{args.n50}.contigs.fasta"]
    vcftools_cmd = "vcftools --gzvcf {} --exclude-bed {} --recode --stdout | bgzip -c > {}"
    tabix_cmd = "tabix -p vcf {}"
    bcftools_cmd = "seqkit replace -p '\\.' -r '{}.' {} | bcftools consensus -f - {} > {}"
    bedtools_cmd = "bedtools getfasta -fi {} -bed {} -nameOnly -fo {}"
    run_cmd(simu_ctg_cmd, log=os.devnull, out2err=True)
    bed_df = read_bed(f"{args.n50}.contigs.new_genome.posi.bed")
    contig_beds = []
    prefix_list = []
    contig_vcf = []
    for i, vcf in enumerate(args.vcf):
        prefix = Path(vcf).stem.rsplit(".", 1)[0]
        prefix_list.append(prefix)
        contig_vcf.append(vcf2contigs(vcf, bed_df, args.suffix[i]))
     
   
    total_n = len(bed_df) * len(prefix_list)
    collapsed_n = int(total_n * args.ratio)
    collapsed_contigs = []
    collapsed_contigs_db = defaultdict(lambda : defaultdict(set))
    while len(collapsed_contigs) < collapsed_n:
        idx1 = random.randint(0, len(prefix_list) - 1)
        idx2 = random.randint(0, len(prefix_list) - 1)
        if idx1 == idx2:
            continue
        
        idx = (idx1, idx2)

        contig_idx = random.randint(0, len(bed_df) - 1)

        tmp_collapsed_contigs = tuple(bed_df.loc[contig_idx].values.tolist())

        if tmp_collapsed_contigs not in collapsed_contigs_db:
            collapsed_contigs_db[tmp_collapsed_contigs][idx1].add(idx2)
            collapsed_contigs.append((tmp_collapsed_contigs, idx))
            continue

        if idx1 not in collapsed_contigs_db[tmp_collapsed_contigs]:
            
            continue 
        else:
            if len(collapsed_contigs_db[tmp_collapsed_contigs][idx1]) >= len(prefix_list) - 1:
                continue
            collapsed_contigs.append((tmp_collapsed_contigs, idx))
            collapsed_contigs_db[tmp_collapsed_contigs][idx1].add(idx2)
            

    bed_df_list = [[] for _ in range(len(prefix_list))]
    collapsed_contigs_origin_id_db = defaultdict(list)
    collapsed_contigs_collapsed_id_db = []
    for item in collapsed_contigs_db:
        for origin in collapsed_contigs_db[item]:
            contig = f"{args.suffix[origin]}.".join(item[3].split("."))
            collapsed_contigs_origin_id_db[contig] = (origin, collapsed_contigs_db[item][origin])
            bed_df_list[origin].append(item)
            for collapsed in collapsed_contigs_db[item][origin]:
                bed_df_list[collapsed].append(item)
                contig = f"{args.suffix[collapsed]}.".join(item[3].split("."))
                collapsed_contigs_collapsed_id_db.append((collapsed, contig))


   
    bed_df_list = list(map(pd.DataFrame, bed_df_list))
    contig_fasta_list = []
    for i, tmp_df in enumerate(bed_df_list):
        tmp_df[0] = tmp_df[3].str.split(".").map(lambda x: f"{args.suffix[i]}.".join(x))
        tmp_df.drop(3, axis=1, inplace=True)
        tmp_df[2] = tmp_df[2] - tmp_df[1]
        tmp_df[1] = tmp_df[1] - tmp_df[1]

        tmp_df.columns = ["chrom", "chromStart", "chromEnd"]
        tmp_df.to_csv(f"{prefix_list[i]}.collapsed.contigs.bed", sep='\t', header=True, index=None)
        os.system(vcftools_cmd.format(contig_vcf[i], 
                                       f"{prefix_list[i]}.collapsed.contigs.bed", 
                                       f"{prefix_list[i]}.collapsed.vcf.gz"))
        os.system(tabix_cmd.format(f"{prefix_list[i]}.collapsed.vcf.gz"))
        os.system(bcftools_cmd.format(args.suffix[i], f"{args.n50}.contigs.fasta", 
                                        f"{prefix_list[i]}.collapsed.vcf.gz", f"{prefix_list[i]}.contigs.fasta" ))
        contig_fasta_list.append(dict(read_fasta(f"{prefix_list[i]}.contigs.fasta")))
        # os.system(bedtools_cmd.format())


    
    for item in collapsed_contigs_collapsed_id_db:
        idx, contig = item 
        contig_fasta_list[idx].pop(contig, None)

    for item in collapsed_contigs_origin_id_db:
        origin, collapsed = collapsed_contigs_origin_id_db[item]
        
        seq = contig_fasta_list[origin][item]
        contig_fasta_list[origin].pop(item, None)
        collapsed_suffix = [item.replace(args.suffix[origin], args.suffix[i]) for i in collapsed]
        collapsed_suffix = "|".join(collapsed_suffix)
        new_contig = f"{item};collapsed|{collapsed_suffix}"
        contig_fasta_list[origin][new_contig] = seq 
    

    for i, db in enumerate(contig_fasta_list):
        with open(f"{prefix_list[i]}.collapsed.contigs.fasta", 'w') as out:
            for contig in db:
                print(f">{contig}", file=out)
                print(f"{db[contig]}", file=out)


if __name__ == "__main__":
    main(sys.argv[1:])