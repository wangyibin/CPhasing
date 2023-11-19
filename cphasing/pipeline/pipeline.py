#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
pipeline of C-Phasing
"""

import argparse
import logging
import os
import os.path as op
import sys

from pathlib import Path
from ..utilities import run_cmd 

logger = logging.getLogger(__name__)


def run(fasta,
        porec_data,
        porec_table,
        pairs, 
        motif="AAGCTT",
        mode="phasing",
        steps=set([0, 1, 2, 3, 4, 5, 6, 7]),
        skip_steps=set(),
        n="",
        resolution1=1,
        resolution2=1,
        allelic_similarity=0.85,
        min_allelic_overlap=0.3,
        factor=50,
        threads=4):
    from ..cli import (mapper,
                       alleles, 
                       prepare,
                       kprune,
                       hypergraph,
                       hyperpartition,
                       scaffolding,
                       pairs2cool,
                       plot
    )


    fasta_prefix = fasta.rsplit(".", 1)[0]
    if porec_data:
        porec_prefix = porec_data.replace(".gz", "").rsplit(".", 1)[0]
        pairs_prefix = porec_data.replace(".gz", "").rsplit(".", 1)[0]
        pairs = f"{pairs_prefix}.pairs.gz"
        hg_input = f"{porec_prefix}.porec.gz"
        hg_flag = ""
        
    else:
        if "0" in steps:
            steps = steps - set("0")
            if "0" not in skip_steps:
                logger.warning("Mapping step will not be run, because the porec data is not specified")
        if porec_table:
            porec_prefix = porec_table.replace(".gz", "").rsplit(".", 1)[0]
            pairs_prefix = porec_table.replace(".gz", "").rsplit(".", 1)[0]
            hg_input = porec_table
            hg_flag = ""
            if not pairs:
                pairs = f"{porec_prefix}.pairs.gz"

        else:
            if pairs:
                pairs_prefix = pairs.replace(".gz", "").rsplit(".", 1)[0]
                hg_input = pairs 
                hg_flag = "--pairs"
            

    
    contigsizes = f"{fasta_prefix}.contigsizes"

    if "0" not in skip_steps and "0" in steps:
        try:
            mapper.main(args=[fasta, 
                              porec_data,
                              "-t",
                              str(threads)],
                            prog_name='alleles')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e

    if "1" not in skip_steps and "1" in steps:
        try:
            alleles.main(args=["-f",
                                fasta],
                            prog_name='alleles')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
        
    allele_table = f"{fasta_prefix}.allele.table"

    if "2" not in skip_steps and "2" in steps:
        try:
            prepare.main(args=[fasta,
                            pairs, 
                            "-m",
                            motif],
                            prog_name='prepare')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
    
    count_re = f"{pairs_prefix}.counts_{motif}.txt"
    contacts = f"{pairs_prefix}.contacts"
    clm = f"{pairs_prefix}.clm"

    if "3" not in skip_steps and "3" in steps:
        try:
            kprune.main(args=[allele_table,
                            contacts,
                            count_re,
                            "-t",
                            str(threads)],
                            prog_name='kprune')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
        
    prune_table = "prune.contig.table" if mode == "phasing" else None

    hg = hg_input.replace(".gz", "").rsplit(".", 1)[0] + ".hg"

    if "4" not in skip_steps and "4" in steps:
        try:
            if hg_flag:
                hypergraph.main(args=[
                                hg_input,
                                contigsizes,
                                hg,
                                hg_flag
                            ],
                            prog_name='hypergraph')
            else:
                hypergraph.main(args=[
                                    hg_input,
                                    contigsizes,
                                    hg,
                                ],
                                prog_name='hypergraph')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
    
    output_cluster = "output.clusters.txt"
  
    if "5" not in skip_steps and "5" in steps:
        try:
            hyperpartition.main(args=[
                                hg,
                                contigsizes,
                                output_cluster,
                                "--mode",
                                mode,
                                "-pt",
                                prune_table,
                                "-n",
                                n,
                                "-r1",
                                resolution1,
                                "-r2",
                                resolution2,
                                "-as",
                                allelic_similarity,
                                "-mao",
                                min_allelic_overlap,
                                "-t",
                                threads
                            ],
                            prog_name='hyperpartition')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
    
    
    out_agp = "groups.agp"
    if "6" not in skip_steps and "6" in steps:
        try:
            scaffolding.main(args=[
                                output_cluster,
                                count_re,
                                clm,
                                "-f",
                                fasta,
                                "-t",
                                threads,
                                "-o",
                                out_agp
                            ],
                            prog_name='scaffolding')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
    
    out_small_cool = f"{pairs_prefix}.10000.cool"


    if "7" not in skip_steps and "7" in steps:
        if Path(out_small_cool).exists():
            logger.warning(f"`{out_small_cool}` exists, skipped `pairs2cool`.")
        else:
            try:
                pairs2cool.main(args=[
                                    pairs,
                                    contigsizes,
                                    out_small_cool
                                ],
                                prog_name='pairs2cool')
            except SystemExit as e:
                exc_info = sys.exc_info()
                exit_code = e.code
                if exit_code is None:
                    exit_code = 0
                
                if exit_code != 0:
                    raise e
                
        try:
            plot.main(args=[
                            "-a",
                            out_agp,
                            "-m",
                            out_small_cool,
                            "-o",
                            "groups.wg.png",
                            "--factor",
                            factor
                            ],
                            prog_name='plot')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e

