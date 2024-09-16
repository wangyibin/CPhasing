#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
evaluate the retain and remove of prune
"""


import argparse
import logging
import os
import os.path as op
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import shutil
import pandas as pd 
import tempfile

import pickle 
import pandas as pd


from pathlib import Path
from itertools import combinations, product, permutations
from joblib import Parallel, delayed
from cphasing.utilities import run_cmd

def get_remove_retain(prunelist, contacts, contigs):
    remove_list = set([tuple(sorted(i.strip().split()[:2])) for i in open(prunelist) if i.strip()])
    contigs_df = pd.read_csv(contigs, sep='\t', usecols=(0,), header=None, index_col=None)
    chrom_contigs = pd.DataFrame(contigs_df[0].str.split(".").values.tolist(), columns=["chrom", "contig"])
    chrom_contigs['contig'] = contigs_df[0]
    chrom_contigs['hap'] = chrom_contigs['chrom'].str[:-1]
    chrom_contigs['hap'] = chrom_contigs['hap'].str.split("_").map(lambda x: x[0])

    contact_df = pd.read_csv(contacts, sep='\t', index_col=None, 
                                    header=None)
    contact_df.columns = ['contig1', 'contig2', 'count']
    contact_df = contact_df.set_index(['contig1', 'contig2'])
    contact_db = contact_df.to_dict()['count']

    inter_pairs = []
    for i, df in chrom_contigs.groupby('hap'):
        tmp_list = []
        for j, df2 in df.groupby('chrom'):
            tmp_list.append(df2['contig'].tolist())

        tmp_list.sort() 
        
        for n, m in list(combinations(tmp_list, 2)):
            inter_pairs.extend(list(set(product(n, m))))
    
    intra_pairs = []
    for i, df in chrom_contigs.groupby('chrom'):
        tmp_list = df['contig'].values.tolist()
        tmp_list = sorted(tmp_list)
        intra_pairs.extend(list(combinations(tmp_list, 2)))
   
    intra_pairs = set(intra_pairs)

    intra_contacts = contact_df.reindex(intra_pairs).dropna()
    inter_contacts = contact_df.reindex(inter_pairs).dropna()
    total_intra_contacts = intra_contacts['count'].sum()
    total_inter_contacts = inter_contacts['count'].sum()
    retain_contacts =  total_intra_contacts -  intra_contacts.reindex(remove_list).dropna()['count'].sum()
    remove_contacts = inter_contacts.reindex(remove_list).dropna()['count'].sum()
    
    

    return Path(contacts).stem, remove_contacts / total_inter_contacts, retain_contacts/ total_intra_contacts


def downsample(pairs, percent):
    if percent == 1.0:
        os.link(pairs, f"{percent:.2f}.pairs.gz")
    else:
        cmd = ["cphasing-rs", "pairs-downsample", pairs,
            "-p", f"{percent}", "-o", f"{percent:.2f}.pairs.gz"]
        run_cmd(cmd)

    cmd = ["cphasing-rs", "pairs2contacts", f"{percent:.2f}.pairs.gz",
           "-o", f"{percent:.2f}.contacts"]
    run_cmd(cmd)

    cmd = ["cphasing-rs", "pairs2bam", f"{percent:.2f}.pairs.gz", 
           "-o", f"{percent:.2f}.bam"]
    run_cmd(cmd)

    cmd = ["samtools", "collate", "-@", "10", f"{percent:.2f}.bam", "-o", f"{percent:.2f}.shuffle.bam"]
    run_cmd(cmd)

    os.remove(f"{percent:.2f}.shuffle.bam", f"{percent:.2f}.bam")

    return f"{percent:.2f}.pairs.gz", f"{percent:.2f}.contacts", f"{percent:.2f}.bam"

def evaluate_cphasing(alleletable, contacts):
    name = Path(contacts).stem
    cmd = ["cphasing-rs", "kprune", f"{alleletable}", f"{contacts}", f"{name}.C-Phasing.prune.table"]

    run_cmd(cmd)

    return f"{name}.C-Phasing.prune.table", contacts

def evaluate_ALLHiC(allhic_alleletable, bam, fasta):

    tmpDir = tempfile.mkdtemp(dir='./')
    
    name = Path(bam).stem
    bam_absolute = Path(bam).absolute()
    fasta_absolute = Path(fasta).absolute()
    allhic_alleletable = Path(allhic_alleletable).absolute()
    
    os.chdir(tmpDir)
    os.link(fasta_absolute, Path(fasta).name)
    os.link(bam_absolute, bam)
    cmd = ["allhic", "extract", f"{bam}", f"{Path(fasta).name}",
           "--RE", "AAGCTT", "--minLinks", "1"]
    run_cmd(cmd)
    os.remove(f"{name}.clm")

    cmd = ["~/software/ALLHiC_components/Prune/ALLHiC_prune",
           "-i", f"{allhic_alleletable}", "-b", f"{bam}"]
    os.system(" ".join(cmd))

    os.rename("prunning.bam", f"{name}.prune.bam")

    cmd = ["allhic", "extract", f"{name}.prune.bam", f"{Path(fasta).name}",
           "--RE", "AAGCTT", "--minLinks", "1"]
    run_cmd(cmd)
                                                                                                                                                                                                             
    
    os.remove(f"{name}.prune.clm")

    raw_pairs_txt = f"{name}.pairs.txt"
    prune_pairs_txt = f"{name}.prune.pairs.txt"

    raw_pairs_df = pd.read_csv(raw_pairs_txt, sep='\t', index_col=None, 
                                header=0)
    prune_pairs_df = pd.read_csv(prune_pairs_txt, sep='\t', index_col=None, 
                                header=0)
    
    raw_pairs = set(map(lambda x: tuple(x), raw_pairs_df[['Contig1', 'Contig2']].values.tolist()))
    prune_pairs = set(map(lambda x: tuple(x), prune_pairs_df[['Contig1', 'Contig2']].values.tolist()))

    remove_pairs = raw_pairs - prune_pairs

    os.chdir("../")
    
    with open(f"{name}.ALLHiC.remove.table", 'w') as out:
        for pair in remove_pairs:
            out.write(f"{pair[0]}\t{pair[1]}\n")
            out.write(f"{pair[1]}\t{pair[0]}\n")
    shutil.rmtree(tmpDir)

    return f"{name}.ALLHiC.remove.table", name


def evaluate_HapHiC(fasta, bam, ploidy, n, contigs):
    tmpDir = tempfile.mkdtemp(dir='./')
    
    name = Path(bam).stem
    contigs = Path(contigs).absolute()
    bam_absolute = Path(bam).absolute()
    fasta_absolute = Path(fasta).absolute()
  
    os.chdir(tmpDir)
    os.link(fasta_absolute, Path(fasta).name)
    os.link(bam_absolute, bam)
    cmd = f"conda run -n haphic ~/software/HapHiC/HapHiC/haphic pipeline {fasta_absolute} {bam} {n} --remove_allelic {ploidy} --steps 1 --verbose"
    os.system(cmd)

    data = pickle.load(open("01.cluster/full_links.pkl", 'rb'))
    all_pairs = [i.strip().split()[0] for i in open(contigs) if i.strip()]
    all_pairs = set(permutations(all_pairs, 2))
    all_pairs = set(filter(lambda x: x[0] < x[1], all_pairs))
    full_pairs = set(list(data.keys()))
    
    pt_pairs = all_pairs - full_pairs
    
    os.chdir("..")
    shutil.rmtree(tmpDir)

    with open(f"{name}.HapHiC.remove.table", 'w') as out:
        for pair in pt_pairs:
            out.write(f"{pair[0]}\t{pair[1]}\n")
            out.write(f"{pair[1]}\t{pair[0]}\n")
    
    return f"{name}.HapHiC.remove.table", name


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('pairs', 
            help='total pairs file')
    pReq.add_argument("alleletable", help="alleletable")
    pReq.add_argument('contigs', 
            help='contig list')
    pReq.add_argument('allhic_alleletable')
    pReq.add_argument('fasta')
    pReq.add_argument('ploidy', type=int)
    pReq.add_argument('n', type=int)
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    s = []
    for i in range(50, 51, 1):

        s.append((args.pairs, i*0.02))
    
    res = Parallel(n_jobs=10)(delayed(downsample)(i, j) for i, j in s)


    s2 = []
    for _, contacts, _ in res:
        s2.append((args.alleletable, contacts))

    res2 = Parallel(n_jobs=10)(delayed(evaluate_cphasing)(i, j) for i, j in s2)

    s3 = []
    for prunetable, contacts in res2:
        s3.append((prunetable, contacts, args.contigs))
    res3 = Parallel(n_jobs=10)(delayed(get_remove_retain)(i, j, k) for i, j, k in s3)
    

    s4 = []
    for _, _, bam in res:
        s4.append((args.allhic_alleletable, bam, args.fasta))
    
    
    allhic_res = Parallel(n_jobs=10)(delayed(evaluate_ALLHiC)
                                             (i, j, k) for i,j,k in s4)
    allhic_remove_list = []
    for remove_list, name in allhic_res:
        contacts = f"{name}.contacts"
        allhic_remove_list.append((remove_list, contacts, args.contigs))
    
    res4 = Parallel(n_jobs=10)(delayed(get_remove_retain)(i, j, k) for i, j, k in allhic_remove_list)


    s5 = []
    for _, _, bam in res:
        s5.append((args.fasta, bam, args.ploidy, args.n, args.contigs))

    haphic_res = Parallel(n_jobs=10)(delayed(evaluate_HapHiC)(i, j, k, l, n) for i,j,k, l, n in s5)

    haphic_remove_list = []
    for remove_list, name in haphic_res:
        contacts = f"{name}.contacts"
        haphic_remove_list.append((remove_list, contacts, args.contigs))
    
    res5 = Parallel(n_jobs=10)(delayed(get_remove_retain)(i, j, k) for i, j, k in haphic_remove_list)


    
    res3_df = pd.DataFrame(res3)
    res4_df = pd.DataFrame(res4)
    res5_df = pd.DataFrame(res5)
    res3_df[3] = 'C-Phasing'
    res4_df[3] = 'ALLHiC'
    res5_df[3] = 'HapHiC'
    result_df = pd.concat([res3_df, res4_df, res5_df])
    
    result_df.to_csv("output.tsv", sep='\t', index=None, header=None)
    


if __name__ == "__main__":
    main(sys.argv[1:])