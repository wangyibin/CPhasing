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
        hcr_flag=False,
        mode="phasing",
        hic=False,
        steps=set([0, 1, 2, 3, 4, 5]),
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

    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    steps = set(steps)
    if mode == 'basal':
        skip_steps.add("1")
        allele_table = None

    fasta_prefix = Path(fasta).with_suffix("")
    while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
        fasta_prefix = fasta_prefix.with_suffix("")

    if porec_data:
        porec_prefix = str(Path(porec_data).name).replace(".gz", "").rsplit(".", 1)[0]
        pairs_prefix = str(Path(porec_data).name).replace(".gz", "").rsplit(".", 1)[0]
        pairs = f"{pairs_prefix}.pairs.gz"
        hg_input = f"{porec_prefix}.porec.gz"
        hg_flag = ""
        porec_table = hg_input
        input_param = "--porec"
        
        if not Path(pairs).exists() and not Path(porec_table).exists():
            steps.add("0")
        
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
            input_param = "--porec"
            if not pairs:
                pairs = f"{porec_prefix}.pairs.gz"

        else:
            if pairs:
                pairs_prefix = pairs.replace(".gz", "").rsplit(".", 1)[0]
                hg_input = pairs 
                hg_flag = "--pairs"
                input_param = hg_flag

    if "0" not in skip_steps and "0" in steps:
        if not hic:
            logger.info("""#----------------------------------#
#      Running step 0. mapper      #
#----------------------------------#""")
            try:
                porec_mapper.main(args=[fasta, 
                                porec_data,
                                # "-p",
                                # pattern,
                                "-t",
                                str(threads)],
                                prog_name='mapper')
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


    if hcr_flag:
        if porec_table:
            hg_input = f"{porec_prefix}_hcr.porec.gz"
            prepare_input = f"{porec_prefix}_hcr.pairs.gz"
            if not Path(hg_input).exists() or not Path(hg_input).exists():
                
                try:
                    hcr.main(
                        args=["-pct",
                                porec_table,
                                "-cs",
                                contigsizes,
                                "-p",
                                95
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
                logger.warn(f"Use exists hcr porec table of `{hg_input}`")
        
        else:
            hg_input = f"{pairs_prefix}_hcr.pairs.gz"
            prepare_input = f"{pairs_prefix}_hcr.pairs.gz"
            if not Path(hg_input).exists() or not Path(hg_input).exists():
                
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
                logger.warn(f"Use exists hcr porec table of `{hg_input}`")
        
            
    else:
        prepare_input = pairs

    if not Path(prepare_input).exists():
        logger.info("Generating pairs file ...")
        cmd = ["cphasing-rs", "porec2pairs", porec_table, contigsizes,
                    "-o", pairs, "-q", "1"]
        
        flag = run_cmd(cmd, log=f"logs/porec2pairs.log")
        assert flag == 0, "Failed to execute command, please check log."

    if "1" not in skip_steps and "1" in steps:
        logger.info("""#----------------------------------#
#      Running step 1. alleles     #
#----------------------------------#""")
        try:
            alleles.main(args=["-f",
                                fasta,
                                ],
                            prog_name='alleles')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
        
    allele_table = None if mode == "basal" else f"{fasta_prefix}.allele.table" 

    if "2" not in skip_steps and "2" in steps:
        logger.info("""#----------------------------------#
#       Running step 2. prepare    #
#----------------------------------#""")
        try:
            prepare.main(args=[fasta,
                            prepare_input, 
                            "-p",
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

    clm = f"{prepare_prefix}.clm"
    
    contacts = f"{prepare_prefix}.contacts"
    split_contacts = f"{prepare_prefix}.split.contacts"
    # if "3" not in skip_steps and "3" in steps:
    #     try:
    #         kprune.main(args=[allele_table,
    #                         contacts,
    #                         "-t",
    #                         str(threads)],
    #                         prog_name='kprune')
    #     except SystemExit as e:
    #         exc_info = sys.exc_info()
    #         exit_code = e.code
    #         if exit_code is None:
    #             exit_code = 0
            
    #         if exit_code != 0:
    #             raise e
        
    # prune_table = "prune.contig.table" if mode == "phasing" else None

    # hg = hg_input.replace(".gz", "").rsplit(".", 1)[0] + ".hg"

    # if "5" not in skip_steps and "5" in steps:
    #     try:
    #         if hg_flag:
    #             hypergraph.main(args=[
    #                             hg_input,
    #                             contigsizes,
    #                             hg,
    #                             hg_flag
    #                         ],
    #                         prog_name='hypergraph')
    #         else:
    #             hypergraph.main(args=[
    #                                 hg_input,
    #                                 contigsizes,
    #                                 hg,
    #                             ],
    #                             prog_name='hypergraph')
    #     except SystemExit as e:
    #         exc_info = sys.exc_info()
    #         exit_code = e.code
    #         if exit_code is None:
    #             exit_code = 0
            
    #         if exit_code != 0:
    #             raise e
    
    output_cluster = "output.clusters.txt"
    
    if hg_flag == "--pairs":
        hyperpartition_contacts = contacts
    else:
        hyperpartition_contacts = None

    if "3" not in skip_steps and "3" in steps:
        logger.info("""#----------------------------------#
#  Running step 3. hyperpartition  #
#----------------------------------#""")
        try:
            hyperpartition.main(args=[
                                hg_input,
                                contigsizes,
                                output_cluster,
                                input_param,
                                "--mode",
                                mode,
                                "-at",
                                allele_table,
                                "-c",
                                hyperpartition_contacts,
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
    if "4" not in skip_steps and "4" in steps:
        logger.info("""#----------------------------------#
#    Running step 4. scaffolding    #
#----------------------------------#""")
        try:
            scaffolding.main(args=[
                                output_cluster,
                                count_re,
                                clm,
                                "-at",
                                allele_table,
                                "-sc",
                                split_contacts,
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


    if "5" not in skip_steps and "5" in steps:
        logger.info("""#----------------------------------#
#      Running step 5. plot        #
#----------------------------------#""")
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

