#!/usr/bin/env python
# -*- coding:utf-8 -*-


import logging
import os
import os.path as op
import sys

import pandas as pd

from collections import Counter
from natsort import natsort_keygen
from joblib import Parallel, delayed
from pathlib import Path

from .core import AlleleTable, CountRE, PairTable
from .utilities import xopen 

logger = logging.getLogger(__name__)

def split_countRE(cr: CountRE, chrom: str, item: list, out: str) -> tuple:
    df = cr.data.reindex(item)
    df = (df.sort_index(key=natsort_keygen())
            .dropna()
            .astype({'RECounts': int, 
                    'Length': int})
            )
    df = df.reset_index()
    df.to_csv(out, sep='\t', header=True, index=False)
    
    return (chrom, dict(zip(df['#Contig'], df.index)))
    
def split_pairTable(pt: PairTable, contig2idx: dict, item: list, out: str) -> None:
    df = pt.get_contacts(item, item)
    df = df.astype({'#X': 'int32', 'Y': 'int32', 
                    'RE1': 'int64', 'RE2': 'int64', 
                    'ObservedLinks': 'int32'})
    df = df.reset_index()[pt.HEADER]
    df = df[df['#X'] < df['Y']]
    df['#X'] = df['Contig1'].map(lambda x: contig2idx[x])
    df['Y'] = df['Contig2'].map(lambda x: contig2idx[x])
    df = df.sort_values(["#X", "Y"])
    df.to_csv(out, sep='\t', header=True, index=None)


def read_fasta(fasta: str) -> dict:
    from Bio import SeqIO
    
    logger.info(f'Load fasta file: `{fasta}`.')
    db = {}
    
    fasta = SeqIO.parse(xopen(fasta), 'fasta')
    for record in fasta:
        if record.id not in db:
            db[record.id] = record.seq
        
    return db

def extract_fasta(fasta_db: dict, item: list, out: str) -> None:
    with open(out, 'w') as o:
        for contig in item:
            o.write(f'>{contig}\n{fasta_db[contig]}\n')

def pregroup(at: AlleleTable, cr: CountRE, pt: PairTable, 
            fasta: str, outdir: str, threads: int) -> None:
    """
    pregroup allhic countRE and pairs by homologous groups.

    Params:
    --------
    at: AlleleTable
        AlleleTable object of allhic allele.table.
    cr: CountRE
        CountRE objcet of allhic count_RE.txt.
    pt: PairTable, unsymmetric
        PairTable object of allhic pairs.txt.
    fasta: str
        Fasta path.
    outdir: str
        Output directory of results.
    threads: int
        Number of threads.
    
    Returns:
    --------
    None

    Examples:
    --------
    >>> pregroup(at, cr, pt, "outdir", 4)
    """
    ## reassign greedy contigs, which may grouped into different chromosomes.
    res = []
    chroms = []
    for chrom, df in at.data.groupby(0):
        if chrom.lower().startswith('utg') \
            or chrom.lower().startswith('ctg') \
                or chrom.lower().startswith('tig') \
                    or chrom.lower().startswith('scaffold') \
                        or chrom.lower().startswith('hic_scaffold'):
            continue
        chroms.append(chrom)
        res.append(pd.Series(Counter(df.values.flatten())))

    res_df = pd.concat(res, axis=1, keys=chroms)
    res_df = res_df.loc[res_df.index.dropna()]
    count_df = res_df.count(axis=1)
    definded_df = (
                    res_df[count_df == 1]
                    .stack()
                    .reset_index()
                    .drop(0, axis=1)
                    .rename(columns={'level_0': 'contig', 
                                    'level_1': 'chrom'})
                    )
    doubt_df = (
                res_df[count_df > 1]
                .idxmax(axis=1)
                .reset_index()
                .rename(columns={'index': 'contig', 
                                0: 'chrom'})
                )
    
    res_df = pd.concat([definded_df, doubt_df])
    
    del definded_df, doubt_df, count_df 

    res_df = res_df.groupby('chrom').agg(list)
    

    logger.info('Splitting count RE ...')
    args1 = []
    for chrom, item in res_df.iterrows():
        chrom_outdir = Path(f'{outdir}/{chrom}')
        chrom_outdir.mkdir(parents=True, exist_ok=True)
        item = item.values[0]
        args1.append((cr, chrom, item, f'{chrom_outdir}/{chrom}.{cr.suffix}'))
        
    r = Parallel(n_jobs=threads)(delayed(split_countRE)(i, j, k, l) 
                            for i, j, k, l in args1)
    
    contig_db = dict(r)

    logger.info('Splitting pairs table ...')
    args2 = []
    for chrom, item in contig_db.items():
        chrom_outdir = Path(f'{outdir}/{chrom}')

        contig2idx = contig_db[chrom]
        item = item.keys()

        args2.append((pt, contig2idx, item, f'{chrom_outdir}/{chrom}.pairs.txt'))

    Parallel(n_jobs=threads)(delayed(split_pairTable)(i, j, k, l) 
                            for i, j, k, l in args2)

    if fasta:
        logger.info('Splitting fasta ...')
        fasta_db = read_fasta(fasta)
        
        args = []
        for chrom, item in contig_db.items():
            chrom_outdir = Path(f'{outdir}/{chrom}')
            item = item.keys()
            
            args.append((fasta_db, item, f'{chrom_outdir}/{chrom}.fa'))

        Parallel(n_jobs=threads)(delayed(extract_fasta)(i, j, k)
                                for i, j, k in args)
    
    logger.info(f'Done, results are in `{outdir}`.')