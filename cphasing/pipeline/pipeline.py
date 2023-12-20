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
        pattern="AAGCTT",
        mode="phasing",
        hic=False,
        steps=set([0, 1, 2, 3, 4, 5, 6, 7, 8]),
        skip_steps=set(),
        n="",
        resolution1=1,
        resolution2=1,
        first_cluster=None,
        allelic_similarity=0.85,
        min_allelic_overlap=0.3,
        whitelist=None,
        factor=50,
        threads=4):
    from ..cli import (mapper as porec_mapper,
                       hcr,
                       alleles, 
                       prepare,
                       kprune,
                       hypergraph,
                       hyperpartition,
                       scaffolding,
                       pairs2cool,
                       plot
    )
    from ..hic.cli import mapper as hic_mapper
    from ..cli import alignments 

    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    
    if mode == 'basal':
        skip_steps.add("3")
        skip_steps.add("4")
    fasta_prefix = fasta.rsplit(".", 1)[0]
    if porec_data:
        porec_prefix = porec_data.replace(".gz", "").rsplit(".", 1)[0]
        pairs_prefix = porec_data.replace(".gz", "").rsplit(".", 1)[0]
        pairs = f"{pairs_prefix}.pairs.gz"
        hg_input = f"{porec_prefix}.porec.gz"
        hg_flag = ""
        porec_table = hg_input
        
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
            


    if "0" not in skip_steps and "0" in steps:
        if not hic:
            try:
                porec_mapper.main(args=[fasta, 
                                porec_data,
                                "-p",
                                pattern,
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
    
    contigsizes = f"{fasta_prefix}.contigsizes"
    if not Path(contigsizes).exists():
        cmd = ["cphasing-rs", "chromsizes", fasta, "-o", contigsizes]
        flag = run_cmd(cmd, log=os.devnull)
        assert flag == 0, "Failed to execute command, please check log."


    if "1" not in skip_steps and "1" in steps:
        if porec_table:
            hg_input = f"{porec_prefix}_hcr.porec.gz"
            prepare_input = f"{porec_prefix}_hcr.pairs.gz"
            try:
                hcr.main(
                    args=["-pct",
                            porec_table,
                            "-cs",
                            contigsizes,
                    ],
                    prog_name="hcr"
                )
            except SystemExit as e:
                exc_info = sys.exc_info()
                exit_code = e.code
                if exit_code is None:
                    exit_code = 0
                
                if exit_code != 0:
                    raise e
        
        else:
            hg_input = f"{pairs_prefix}_hcr.pairs.gz"
            prepare_input = f"{pairs_prefix}_hcr.pairs.gz"
            try:
                hcr.main(
                    args=["-prs",
                            pairs,
                            "-cs",
                            contigsizes,
                    ],
                    prog_name="hcr"
                )
            except SystemExit as e:
                exc_info = sys.exc_info()
                exit_code = e.code
                if exit_code is None:
                    exit_code = 0
                
                if exit_code != 0:
                    raise e
    else:
        prepare_input = pairs

    if not Path(prepare_input).exists():
        cmd = ["cphasing-rs", "porec2pairs", porec_table, contigsizes,
                    "-o", pairs]
        
        flag = run_cmd(cmd, log=f"logs/porec2pairs.log")
        assert flag == 0, "Failed to execute command, please check log."

    if "2" not in skip_steps and "2" in steps:
        try:
            prepare.main(args=[fasta,
                            prepare_input, 
                            "-m",
                            pattern],
                            prog_name='prepare')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
    
    prepare_prefix = prepare_input.replace(".gz", "").rsplit(".", 1)[0]
    count_re = f"{prepare_prefix}.counts_{pattern}.txt"
    contacts = f"{prepare_prefix}.contacts"
    clm = f"{prepare_prefix}.clm"


    if "3" not in skip_steps and "3" in steps:
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

    if "4" not in skip_steps and "4" in steps:
        try:
            kprune.main(args=[allele_table,
                            contacts,
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

    if "5" not in skip_steps and "5" in steps:
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
  
    if "6" not in skip_steps and "6" in steps:
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
                                "-fc",
                                first_cluster,
                                "-as",
                                allelic_similarity,
                                "-mao",
                                min_allelic_overlap,
                                "-wl",
                                whitelist,
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
    if "7" not in skip_steps and "7" in steps:
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


    if "8" not in skip_steps and "8" in steps:
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

