#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
pipeline of C-Phasing
"""

import logging
import os
import os.path as op
import re
import sys
import time 
# import tracemalloc

from datetime import date
from pathlib import Path
from .. import __version__
from ..utilities import (
    run_cmd, 
    calculate_Nx_from_contigsizes,
    read_chrom_sizes
    )


logger = logging.getLogger(__name__)


def run(fasta,
        porec_data,
        porec_table,
        pairs, 
        hic1=None,
        hic2=None,
        pattern="AAGCTT",
        mapper_k=15,
        mapper_w=10,
        mapping_quality=0,
        hic_mapper_k=17,
        hic_mapper_w=7,
        chimeric_correct=False,
        chimeric_corrected=False,
        hcr_flag=False,
        hcr_lower=0.1,
        hcr_upper=1.75,
        collapsed_contig_ratio=0.9,
        hcr_bs=10000,
        hcr_bed=None,
        hcr_invert=False,
        mode="phasing",
        hic=False,
        steps=set([0, 1, 2, 3, 4, 5]),
        skip_steps=set(),
        alleles_kmer_size=19,
        alleles_window_size=19,
        alleles_minimum_similarity=0.2,
        alleles_diff_thres=0.1,
        scaffolding_method="haphic",
        n="",
        use_pairs=False,
        resolution1=1,
        resolution2=1,
        init_resolution1=1,
        init_resolution2=1,
        first_cluster=None,
        normalize=False,
        exclude_group_to_second=None,
        allelic_similarity=0.85,
        min_allelic_overlap=0.3,
        min_weight=0.1,
        min_quality1=1,
        min_quality2=2,
        min_contacts=5,
        min_length=10000,
        Nx=100,
        min_scaffold_length=5e6,
        disable_misassembly_remove=False,
        whitelist=None,
        blacklist=None,
        factor=50,
        low_memory=False,
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
    from ..chimeric import run as chimeric_run 
    # tracemalloc.start()
    start_time = time.time()
    logger.info(f"C-Phasing version: {__version__}")
    today = date.today().strftime("%Y-%m-%d")
    logger.info(f"Pipeline begins in {today}.")


    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    steps = set(steps)

    raw_fasta = fasta
    
    # if n:
    #     if len(re.split(":|x|\|", n)) <= 1:
    #         mode = 'basal'

    mode = 'basal' if mode == 'haploid' else mode
    if mode == 'basal':
        
        skip_steps.add("1")
        allele_table = None

    fasta_prefix = Path(fasta).with_suffix("")
    while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
        fasta_prefix = fasta_prefix.with_suffix("")

    if porec_data:
        porec_prefix = str(Path(porec_data[0]).name).replace(".gz", "").rsplit(".", 1)[0]
        pairs_prefix = str(Path(porec_data[0]).name).replace(".gz", "").rsplit(".", 1)[0]
        pairs = f"{pairs_prefix}.pairs.gz"
        hg_input = f"{porec_prefix}.porec.gz"
        hg_flag = ""
        porec_table = hg_input
        input_param = "--porec"
        
        if not Path(pairs).exists() and not Path(porec_table).exists() and "0" not in skip_steps:
            steps.add("0")
    
    elif hic1 and hic2:
        pairs_prefix = Path(Path(hic1).stem).with_suffix('')
        while pairs_prefix.suffix in {'.fastq', 'gz', 'fq', '.fq', '.gz', '_R1', '_1', '_2'}:
            pairs_prefix = pairs_prefix.with_suffix('')
        pairs_prefix = str(pairs_prefix).replace('_R1', '')
        pairs = f"{pairs_prefix}.pairs.gz"
        hg_input = f"{pairs_prefix}.pairs.gz"
        porec_table = None
        hg_flag = ""
        input_param = "--pairs"

        if not Path(pairs).exists() and "0" not in skip_steps:
            steps.add("0")
        
    else:
        if "0" in steps:
            steps = steps - set("0")
            if "0" not in skip_steps:
                logger.warning("Mapping step will not be run")
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
            if porec_data:
                try:
                    porec_mapper.main(args=[fasta, 
                                    *porec_data,
                                    # "-p",
                                    # pattern,
                                    "-k",
                                    mapper_k,
                                    "-w",
                                    mapper_w,
                                    "-q",
                                    mapping_quality,
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
            
            if hic1 and hic2:
                try:
                    hic_mapper.main(args=[
                                    "-r",
                                    fasta, 
                                    "-1" ,
                                    hic1,
                                    "-2",
                                    hic2,
                                    "-k",
                                    hic_mapper_k,
                                    "-w",
                                    hic_mapper_w,
                                    "-q",
                                    mapping_quality,
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
    
    corrected = False
    if chimeric_correct or chimeric_corrected:
        if porec_table:
            if chimeric_corrected:

                corrected_items = (f"{fasta_prefix}.chimeric.contigs.bed",
                                    f"{fasta_prefix}.corrected.fasta", 
                                   f"{pairs_prefix}.corrected.pairs.gz")
                if not all(map(lambda x: Path(x).exists, corrected_items)):
                    corrected_items = ()

                if corrected_items and not Path(f"{porec_prefix}.corrected.porec.gz").exists:
                    break_bed, fasta, pairs = corrected_items
                    corrected = True
                    cmd = ["cphasing-rs", "porec-break", porec_table, 
                            break_bed, "-o", f"{porec_prefix}.corrected.porec.gz"]
                    flag = run_cmd(cmd, log="logs/porec-break.log")
                    assert flag == 0, "Failed to execute command, please check log."
                    porec_table = f"{porec_prefix}.corrected.porec.gz"
                    hg_input = porec_table

                if corrected_items:
                    logger.info("Using exists corrected results.")
                    break_bed, fasta, pairs = corrected_items
                    corrected = True
                    porec_table = f"{porec_prefix}.corrected.porec.gz"
                    hg_input = porec_table
                    
            else:
                corrected_items = chimeric_run(fasta, pairs, break_pairs=True, 
                                            outprefix=fasta_prefix, 
                                            low_memory=low_memory, threads=threads)
                
                if corrected_items:
                    corrected = True
                    break_bed, fasta, pairs = corrected_items
                    hg_input = pairs
                    cmd = ["cphasing-rs", "porec-break", porec_table, 
                            break_bed, "-o", f"{porec_prefix}.corrected.porec.gz"]
                    flag = run_cmd(cmd, log="logs/porec-break.log")
                    assert flag == 0, "Failed to execute command, please check log."
                    porec_table = f"{porec_prefix}.corrected.porec.gz"


            fasta_prefix = Path(fasta).with_suffix("")
            while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
                fasta_prefix = fasta_prefix.with_suffix("")

            pairs_prefix = Path(Path(pairs).stem).with_suffix('')
            while pairs_prefix.suffix in {'.pairs', '.gz'}:
                pairs_prefix = pairs_prefix.with_suffix('')
            
            porec_prefix = Path(Path(porec_table).stem).with_suffix('')
            while porec_prefix.suffix in {'.porec', '.gz'}:
                porec_prefix = porec_prefix.with_suffix('')

        elif pairs:
            if chimeric_corrected:
                corrected_items = (f"{fasta_prefix}.chimeric.contigs.bed",
                                    f"{fasta_prefix}.corrected.fasta", 
                                   f"{pairs_prefix}.corrected.pairs.gz")
                if not all(map(lambda x: Path(x).exists, corrected_items)):
                    corrected_items = ()
                else:
                    logger.info("Using exists corrected results.")
                
            else:
                corrected_items = chimeric_run(fasta, pairs, break_pairs=True, 
                                                outprefix=fasta_prefix, 
                                                low_memory=low_memory, threads=threads)
          
            
            if corrected_items:
                break_bed, fasta, pairs = corrected_items
                corrected = True
                fasta_prefix = Path(Path(fasta).name).with_suffix("")
                while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
                    fasta_prefix = fasta_prefix.with_suffix("")

                pairs_prefix = Path(Path(pairs).stem).with_suffix('')
                while pairs_prefix.suffix in {'.pairs', '.gz'}:
                    pairs_prefix = pairs_prefix.with_suffix('')


    contigsizes = f"{fasta_prefix}.contigsizes"
    if not Path(contigsizes).exists():
        
        cmd = ["cphasing-rs", "chromsizes", fasta, "-o", contigsizes]
        flag = run_cmd(cmd, log=os.devnull)
        assert flag == 0, "Failed to execute command, please check log."

    if Nx < 100 and Nx > 0:
        contigsizes_df = read_chrom_sizes(contigsizes)
        min_length = calculate_Nx_from_contigsizes(contigsizes_df, Nx)
        retain_contigs = contigsizes_df[contigsizes_df['length'] > min_length]
        retain_contigs = retain_contigs.index.tolist()
        
        if whitelist:
            whitelist_contigs = set(i.strip().split()[0] for i in open(whitelist) if i.strip())
            whitelist_contigs = set(retain_contigs).intersection(whitelist_contigs)
        else:
            whitelist_contigs = retain_contigs

        whitelist = f"N{Nx}.contigs.list"
        with open(whitelist, 'w') as out:
            out.write("\n".join(whitelist_contigs))
        logger.info(f"Filter `{len(contigsizes_df) - len(retain_contigs)}` contig which length < {min_length} (N{Nx})")

    if hcr_flag or hcr_bed:
        hcr_invert_string = "-v" if hcr_invert else ""

        if porec_table and not use_pairs:
            hg_input = f"{porec_prefix}_hcr.porec.gz"
            prepare_input = f"{porec_prefix}.pairs.gz"
            if not Path(hg_input).exists() or not Path(prepare_input).exists():
                try:
                    if hcr_invert_string:
                        hcr.main(
                            args=["-pct",
                                    porec_table,
                                    "-cs",
                                    contigsizes,
                                    "-l",
                                    hcr_lower,
                                    "-u",
                                    hcr_upper,
                                    "-cr",
                                    collapsed_contig_ratio,
                                    "-bs",
                                    hcr_bs,
                                    "-b",
                                    hcr_bed,
                                    "-q",
                                    min_quality1,
                                    hcr_invert_string
                            ],
                            prog_name="hcr"
                        )
                    else:
                        hcr.main(
                            args=["-pct",
                                    porec_table,
                                    "-cs",
                                    contigsizes,
                                    "-l",
                                    hcr_lower,
                                    "-u",
                                    hcr_upper,
                                    "-cr",
                                    collapsed_contig_ratio,
                                    "-bs",
                                    hcr_bs,
                                    "-b",
                                    hcr_bed,
                                   
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
                logger.warn(f"Using existed hcr porec table of `{hg_input}`")
        
        else:
            hg_input = f"{pairs_prefix}_hcr.pairs.gz"
            prepare_input = f"{pairs_prefix}.pairs.gz"
            input_param = "--pairs"
            if not Path(hg_input).exists() or not Path(prepare_input).exists():
                
                try:
                    if hcr_invert_string:
                        hcr.main(
                            args=["-prs",
                                    pairs,
                                    "-cs",
                                    contigsizes,
                                    "-b",
                                    hcr_bed,
                                    "-u",
                                    hcr_upper,
                                    "-l",
                                    hcr_lower,
                                    "-q",
                                    min_quality1,
                                    hcr_invert_string

                            ],
                            prog_name="hcr"
                        )
                    else:
                        hcr.main(
                            args=["-prs",
                                    pairs,
                                    "-cs",
                                    contigsizes,
                                    "-b",
                                    hcr_bed,
                                    "-u",
                                    hcr_upper,
                                    "-l",
                                    hcr_lower,
                                    "-q",
                                    min_quality1,

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
        if not porec_table and min_quality1 > 0: 
            hg_input = f"{pairs_prefix}.q{min_quality1}.pairs.gz"
            if not Path(hg_input).exists():
                cmd = ["cphasing-rs", "pairs-filter", prepare_input, 
                    "-o", hg_input , "-q", str(min_quality1)]
                flag = run_cmd(cmd, log=f'{log_dir}/pairs_filter.log')
                assert flag == 0, "Failed to execute command, please check log."
            else:
                logger.warning(f"Using exists filtered pairs file of `{hg_input}`")


    if not Path(prepare_input).exists() and (hic1 is None ):
        logger.info("Generating pairs file ...")
        cmd = ["cphasing-rs", "porec2pairs", porec_table, contigsizes,
                    "-o", pairs, "-q", "0"]
        
        flag = run_cmd(cmd, log=f"logs/porec2pairs.log")
        assert flag == 0, "Failed to execute command, please check log."

    if "1" not in skip_steps and "1" in steps:
        logger.info("""#----------------------------------#
#      Running step 1. alleles     #
#----------------------------------#""")
        try:
            alleles.main(args=["-f",
                                fasta,
                                "-k",
                                alleles_kmer_size,
                                "-w",
                                alleles_window_size,
                                "-m",
                                alleles_minimum_similarity,
                                "-d",
                                alleles_diff_thres,
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

    prepare_prefix = prepare_input.replace(".gz", "").rsplit(".", 1)[0]
    if "2" not in skip_steps and "2" in steps:
        logger.info("""#----------------------------------#
#       Running step 2. prepare    #
#----------------------------------#""")
        try:
            if porec_table:
                args = [fasta,
                                prepare_input, 
                                "-p",
                                pattern, 
                                "-q",
                                mapping_quality,
                                "--skip-pairs2contacts"]
                if low_memory:
                    args.append("--low-memory")
                prepare.main(args=args,
                                prog_name='prepare')
            else:
                args = [fasta,
                        prepare_input, 
                        "-p",
                        pattern,
                        "-q",
                        mapping_quality,
                        "--skip-pairs2contacts"]
                if low_memory:
                    args.append("--low-memory")
                
                prepare.main(args=args,
                                prog_name='prepare')
                 
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
        
    if not porec_table:
        if mode == 'phasing' or mode == 'basal_withprune':
            if not Path(f"{prepare_prefix}.q{min_quality2}.contacts").exists():
                cmd = ["cphasing-rs", "pairs2contacts", str(hg_input), 
                        "-q", str(min_quality2), 
                    "-o", f"{prepare_prefix}.q{min_quality2}.contacts" ]
                flag = run_cmd(cmd, log=f'{log_dir}/prepare.pairs2contacts.log')
                assert flag == 0, "Failed to execute command, please check log."
            else:
                logger.warning(f"Use exists contacts of `{prepare_prefix}.q{min_quality2}.contacts`")

            contacts = f"{prepare_prefix}.q{min_quality2}.contacts"
        else:
            if not Path(f"{prepare_prefix}.q{min_quality1}.contacts").exists():
                cmd = ["cphasing-rs", "pairs2contacts", str(hg_input), 
                        "-q", str(min_quality1),
                    "-o", f"{prepare_prefix}.q{min_quality1}.contacts" ]
                flag = run_cmd(cmd, log=f'{log_dir}/prepare.pairs2contacts.log')
                assert flag == 0, "Failed to execute command, please check log."
            else:
                logger.warning(f"Use exists contacts of `{prepare_prefix}.q{min_quality1}.contacts`")
            contacts = f"{prepare_prefix}.q{min_quality1}.contacts"
   

    count_re = f"{prepare_prefix}.counts_{pattern}.txt"

    clm = f"{prepare_prefix}.clm"
    
    
    split_contacts = f"{prepare_prefix}.split.contacts"
    
    output_cluster = "output.clusters.txt"
    
    hyperpartition_normalize = "-norm" if normalize else None
    disable_misassembly_remove = "--disable-misassembly-remove" if disable_misassembly_remove else None 


    if hg_flag == "--pairs":
        hyperpartition_contacts = contacts
    else:
        hyperpartition_contacts = None

    if use_pairs and porec_table and not hcr_flag:
        hg_input = pairs 
        hg_flag == "--pairs"
        input_param = "--pairs"

    if "3" not in skip_steps and "3" in steps:
        logger.info("""#----------------------------------#
#  Running step 3. hyperpartition  #
#----------------------------------#""")
        
        hyperpartition_args = [
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
                                "-ir1",
                                init_resolution1,
                                "-r2",
                                resolution2,
                                "-ir2",
                                init_resolution2,
                                "-fc",
                                first_cluster,
                                "--exclude-group-to-second",
                                exclude_group_to_second,
                                "-as",
                                allelic_similarity,
                                "-mao",
                                min_allelic_overlap,
                                "-mw",
                                min_weight,
                                "-q1",
                                min_quality1,
                                "-q2",
                                min_quality2,
                                "-mc",
                                min_contacts,
                                "-ml",
                                min_length,
                                "-ms",
                                min_scaffold_length,
                                "-wl",
                                whitelist,
                                "-t",
                                threads
                            ]
        if disable_misassembly_remove:
            hyperpartition_args.append(disable_misassembly_remove)
        if hyperpartition_normalize:
            hyperpartition_args.append(hyperpartition_normalize)

        try:
            hyperpartition.main(args=hyperpartition_args,
                            prog_name='hyperpartition')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e

    out_agp = "groups.agp"
    if corrected:
        corrected_agp = out_agp.replace("agp", "corrected.agp")
    
    if "4" not in skip_steps and "4" in steps:
        logger.info("""#----------------------------------#
#    Running step 4. scaffolding   #
#----------------------------------#""")
        try:
            if corrected:
                scaffolding.main(args=[
                                    output_cluster,
                                    count_re,
                                    clm,
                                    "-at",
                                    allele_table,
                                    "-sc",
                                    split_contacts,
                                    "-f",
                                    raw_fasta,
                                    "-t",
                                    threads,
                                    "-o",
                                    out_agp,
                                    "-m",
                                    scaffolding_method,
                                    "--corrected"
                                ],
                                prog_name='scaffolding')
            else:
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
                                    out_agp,
                                    "-m",
                                    scaffolding_method
                                ],
                                prog_name='scaffolding')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
    
    out_small_cool = f"{pairs_prefix}.q{min_quality1}.10000.cool"


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
                                    out_small_cool,
                                    "-q", 
                                    min_quality1
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
            if corrected:
                if Path(corrected_agp).exists():
                    input_agp = corrected_agp
                else:
                    input_agp = out_agp
            else:
                input_agp = out_agp
            plot.main(args=[
                            "-a",
                            input_agp,
                            "-m",
                            out_small_cool,
                            "-o",
                            "groups.wg.png",
                            "--factor",
                            factor,
                            "-oc",
                            ],
                            prog_name='plot')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e

    today = date.today().strftime("%Y-%m-%d")
    end_time = time.time() - start_time
    # peak_memory = tracemalloc.get_traced_memory()[1] / 1024 / 1024 / 1024
    # peak_memory = 0
    logger.info(f"Pipeline finished in {today}. Elapsed time {end_time:.2f} s. ") #Peak memory: {peak_memory:2f} Gb")
    # tracemalloc.stop()