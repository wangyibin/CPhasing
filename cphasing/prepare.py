#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
prepare for subsequence analysis
"""

import logging
import os
import os.path as op
import sys

import pandas as pd
import polars as pl

from collections import OrderedDict
from pathlib import Path

from .utilities import (
    decompress_cmd,
    digest, 
    get_contig_length,
    read_fasta,
    run_cmd
)

logger = logging.getLogger(__name__)

def write_chrom_sizes(fasta, output):
    length_db = get_contig_length(fasta)
    with open(output, 'w') as out:
        for contig, length in length_db.items():
            out.write(f"{contig}\t{length}\n")

    logger.info(f"Successful output contigs size file in `{output}`.")

def count_re_in_genome(fasta, enzyme, output=None):
    """
    count the RE sites.

    Params:
    --------
    fasta: str
        Path of fasta file.
    enzyme: str

    """
    
    fasta_records = read_fasta(fasta)

    logger.info(f"Starting count {enzyme} sites in {fasta}...")
    site_df = digest(fasta_records, enzyme)
    res_df = site_df.groupby('chrom', sort=False)['start'].count()
    length_db = OrderedDict(zip(fasta_records.keys(), map(len, fasta_records.values())))
    length_df = pd.DataFrame(length_db, index=['length']).T

    res_df = pd.concat([res_df, length_df], axis=1)
    res_df = res_df.reset_index()
    res_df.columns = ['#Contig', "RECounts", "Length"]
    
    if output:
        res_df.to_csv(output, sep='\t', header=True, index=None)
        logger.info(f"Successful output Count RE file in `{output}`.")
    return res_df

def split_contacts_to_contacts(split_contacts, output):
    df = pl.read_csv(split_contacts, separator='\t', has_header=False,
                     new_columns=['contig1', 'contig2', 'count'],
                     dtypes={
                             'count': pl.UInt32})
    df = df.with_columns(pl.col('contig1').str.strip_suffix("_0").str.strip_suffix("_1"),
                         pl.col('contig2').str.strip_suffix("_0").str.strip_suffix("_1"))
    df.group_by(['contig1', 'contig2']).sum().write_csv(output, separator='\t', include_header=False)

def pipe(fasta, pairs, pattern="AAGCTT", min_mapq=0, min_contacts=3, 
            threads=4, low_memory=False,
            outprefix=None, skip_pairs2clm=False,
            skip_pairs2contacts=False, log_dir="logs"):

    if pattern is None:
        pattern = "AAGCTT"
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
   
    pairs = Path(pairs)
    prefix = pairs.with_suffix('')
    while prefix.suffix in {'.pairs', 'gz', 'pairs'}:
        prefix = prefix.with_suffix('')
    prefix = Path(prefix).name
    outprefix = prefix if outprefix is None else outprefix

    ## count re
    cmd = ["cphasing-rs", "count_re", fasta, "-o", 
            f"{outprefix}.counts_{pattern}.txt", "-p", pattern]
    flag = run_cmd(cmd, log=f'{log_dir}/prepare.count_re.log')
    assert flag == 0, "Failed to execute command, please check log."

    ## pairs2clm
    if not skip_pairs2clm:
        if str(pairs).endswith('.gz'):
            cmd0 = decompress_cmd(str(pairs), str(threads))
            logger.info(f"\nGenerating clm file from pairs file `{pairs}`")
            cmd = ["cphasing-rs", "pairs2clm", "-", "-c", str(min_contacts),
                    "-t", str(threads), "-o", f"{outprefix}.clm.gz", "-q", str(min_mapq)]

            # print(" ".join(cmd0) + f" 2>{log_dir}/prepare.decompress.log" + " | " + " ".join(cmd) + f" 2>{log_dir}/prepare.pairs2clm.log")
            flag = os.system(" ".join(cmd0) + f" 2>{log_dir}/prepare.decompress.log" + " | " + " ".join(cmd) + f" 2>{log_dir}/prepare.pairs2clm.log")
        else:
            cmd = ["cphasing-rs", "pairs2clm", str(pairs), "-c", str(min_contacts),
                    "-t", str(threads), "-o", f"{outprefix}.clm.gz", "-q", str(min_mapq)]
            
            flag = run_cmd(cmd, log=f'{log_dir}/prepare.pairs2clm.log')

        assert flag == 0, "Failed to execute command, please check log."

    # ## pairs2split_contacts
    # cmd = ["cphasing-rs", "pairs2contacts", str(pairs), "-c", str(min_contacts), 
    #         "-n", "2", "-o", f"{outprefix}.split.contacts"]
    # run_cmd(cmd, log=f'{log_dir}/prepare.pairs2contacts.log')

    if not skip_pairs2contacts:
        # cmd = ["cphasing-rs", "pairs2contacts", str(pairs), "-c", str(min_contacts),
        #         "-o", f"{outprefix}.contacts" ]
        # flag = run_cmd(cmd, log=f'{log_dir}/prepare.pairs2contacts.log')
        # assert flag == 0, "Failed to execute command, please check log."
        split_contacts_to_contacts(f"{outprefix}.split.contacts", f"{outprefix}.contacts"  )
        logger.info(f"Output contacts `{outprefix}.contacts`")
    

def pairs2depth(pq):
    """
    backend function for pairs2depth
    """
    schema = {
        "read_idx": pl.Utf8,
        "chrom1": pl.Utf8,
        "pos1": pl.UInt32,
        "chrom2": pl.Utf8,
        "pos2": pl.UInt32,
        "strand1": pl.Categorical,
        "strand2": pl.Categorical,
        "mapq": pl.UInt8
    }

    os.environ["POLARS_MAX_THREADS"] = str(10)
    df = pl.scan_parquet(f"{pq}", schema=schema)

    df = df.select([
        "chrom1", "pos1", "chrom2", "pos2", "mapq"
    ])
    df = (df.filter(pl.col("mapq") >= 0)
            .with_columns(
        ((pl.col("pos1") - 1) // 10000).alias("pos1"),
        ((pl.col("pos2") - 1) // 10000).alias("pos2"),
            )
    )

    df1 = df.select(["chrom1", "pos1", "mapq"]).group_by(["chrom1", "pos1"]).agg([
        pl.col("mapq").count().alias("count")
    ])

    df2 = df.select(["chrom2", "pos2", "mapq"]).group_by(["chrom2", "pos2"]).agg([
        pl.col("mapq").count().alias("count")
    ])
    

    ## merge the two dataframes
    df = df1.join(df2, left_on=["chrom1", "pos1"], right_on=["chrom2", "pos2"], how="full")

    df = df.drop_nulls(subset=["chrom1", "pos1", "chrom2", "pos2"])
    df = df.with_columns(
        (pl.col("count") + pl.col("count_right")).alias("count")
    )
    df = df.drop(["count_right"])
    df  = df.with_columns(
        ((pl.col("pos1") + 1) * 10000 - 10000).alias("start"),
        ((pl.col("pos1") + 1) * 10000).alias("end")
    ).drop(["pos1"])

    df = df.sort(["chrom1", "start"])
    df = df.select(["chrom1", "start", "end", "count"])

    df.collect().write_csv("output.csv", separator="\t", include_header=False)
    