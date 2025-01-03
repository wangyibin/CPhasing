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
import shutil
import time 
# import tracemalloc


from datetime import date
from logging.handlers import RotatingFileHandler
from pathlib import Path
from .. import __version__
from ..utilities import (
    pretty_cmd,
    run_cmd, 
    calculate_Nx_from_contigsizes,
    to_humanized2,
    recommend_binsize_by_genomesize,
    read_chrom_sizes,
    is_empty,
    is_compressed_table_empty,
    is_file_changed,
    generate_to_hic_cmd,
    generate_plot_cmd
    )


logger = logging.getLogger(__name__)

# file_handler = RotatingFileHandler("app.log", maxBytes=10**6, backupCount=3)
# file_handler.setLevel(logging.INFO)
# file_formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
# file_handler.setFormatter(file_formatter)
# logger.addHandler(file_handler)

def run(fasta,
        ul_data,
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
        use_existed_hitig=False,
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
        disable_merge_in_first=False,
        exclude_group_to_second=None,
        exclude_group_from_first=None,
        allelic_similarity=0.85,
        min_allelic_overlap=0.3,
        min_weight=0.1,
        min_quality1=1,
        min_quality2=2,
        min_contacts=5,
        min_length=10000,
        Nx=100,
        min_scaffold_length=5e6,
        enable_misassembly_remove=False,
        whitelist=None,
        blacklist=None,
        binsize="500k",
        low_memory=False,
        outdir="cphasing_output",
        threads=4):
    from ..cli import (mapper as porec_mapper,
                       hcr,
                       alleles, 
                       prepare,
                       hyperpartition,
                       scaffolding,
                       pairs2cool,
                       plot
    )
    from ..hic.cli import mapper as hic_mapper
    
    
    # tracemalloc.start()
    start_time = time.time()
    logger.info(f"C-Phasing version: {__version__}")
    today = date.today().strftime("%Y-%m-%d")
    logger.info(f"Pipeline is started on {today}.")

    _n = re.split(":|x|\|", n) if n else None
    if _n is not None:
        if len(_n) == 2:
            try: 
                _n[1] = int(_n[1])
            except ValueError:
                if Path(_n[1]).exists():
                    _n[1] = str(Path(_n[1]).absolute())
                    n = ":".join(_n)

    if _n is not None:
        if len(_n) <= 1 and mode != "basal_withprune":
            mode = "haploid"
            logger.info("Mode is set to `haploid` because of the second `n` is not specified.")
        

    outdir = Path(outdir)
    if outdir.exists():
        logger.info(f"Working on existing directory: `{outdir}`")
    else:
        logger.info(f"Working on new directory: `{outdir}`")
        outdir.mkdir(parents=True, exist_ok=True)

    fasta_path = str(Path(fasta).absolute())
    fasta = str(Path(fasta).name)
    if ul_data:
        ul_data = str(Path(ul_data).absolute())
    
    if porec_data:
        porec_data = list(porec_data)
        for i, j in enumerate(porec_data):
            porec_data[i] = str(Path(j).absolute())
    
    if hic1 and hic2:
        hic1 = str(Path(hic1).absolute())
        hic2 = str(Path(hic2).absolute())
    
    
    if porec_table:
        porec_table_path = Path(porec_table).absolute()
        if not pairs:
            pairs_prefix = porec_table.replace(".gz", "").rsplit(".", 1)[0]
            _pairs = Path(f"{pairs_prefix}.pairs.gz")
            if _pairs.exists():
                pairs = str(_pairs)
        
    if pairs:
        pairs_path = Path(pairs).absolute()

    if hcr_bed:
        hcr_bed = str(Path(hcr_bed).absolute())

    if first_cluster:
        first_cluster = str(Path(first_cluster).absolute())
    
    if whitelist:
        whitelist = str(Path(whitelist).absolute())

    os.chdir(outdir)

    try:
        Path(fasta).symlink_to(fasta_path)
    except FileExistsError:
        pass

    if porec_table:
        target_path = Path(Path(porec_table).name)
        try:
            target_path.symlink_to(porec_table_path)
        except FileExistsError:
            pass
        porec_table = Path(porec_table).name
    
    if pairs:
        target_path = Path(Path(pairs).name)
        try:
            target_path.symlink_to(pairs_path)
        except FileExistsError:
            pass 

        pairs = Path(pairs).name


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
        logger.info("The mode is `haploid`, skip step 1. alleles.")

    filtered_pairs = None
   
    if ul_data:
        from ..hitig.pipeline import run  as hitig_run
        input_fasta = str(Path(fasta).absolute())
        input_ul_data = str(Path(ul_data).absolute())              
        if not use_existed_hitig:
            Path("0_1.hitig").mkdir(exist_ok=True)
            os.chdir("0_1.hitig")
            hitig_run(input_fasta, input_ul_data)
           
            os.chdir("..")

            hcr_bed = Path("0_1.hitig/output.hcr_all.bed").absolute()
            # clean_fasta = "hitig/output.cleaned.fasta"
        else:
            logger.warning("Use existed hitig results.")
            hcr_bed = Path("0_1.hitig/output.hcr_all.bed").absolute()
            # clean_fasta = "hitig/output.cleaned.fasta"

        if not is_empty("0_1.hitig/output.breakPos.txt"):
            fasta = "0_1.hitig/output.corrected.fasta"
        # if Path(hcr_bed).exists() and Path(clean_fasta).exists():
        #     fasta = clean_fasta
        # else:
        #     logger.warning("Skipped hitig.")

            

    fasta_prefix = Path(Path(fasta).name).with_suffix("")
    while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
        fasta_prefix = fasta_prefix.with_suffix("")


    if porec_data:
        porec_prefix = str(Path(porec_data[0]).name).replace(".gz", "").rsplit(".", 1)[0]
        pairs_prefix = str(Path(porec_data[0]).name).replace(".gz", "").rsplit(".", 1)[0]
        paf = f"{porec_prefix}.paf.gz"
        pairs = f"{pairs_prefix}.pairs.gz"
        hg_input = f"{porec_prefix}.porec.gz"
        hg_flag = ""
        porec_table = hg_input
        input_param = "--porec"
        
        
        if not Path(pairs).exists() or not Path(porec_table).exists() and "0" not in skip_steps:
            steps.add("0")

        if Path(paf).exists():
            if is_compressed_table_empty(paf):
                logger.info(f"The existing paf `{paf}` is empty, rerun step 0.mapper.")
                steps.add("0")
                
                os.remove(paf)

        if Path(porec_table).exists():
            if is_compressed_table_empty(porec_table):
                logger.info(f"The existing porec table `{porec_table}` is empty, rerun step 0.mapper.")
                steps.add("0")
                
                os.remove(porec_table)

        
    
    elif hic1 and hic2:
        pairs_prefix = Path(Path(hic1).stem).with_suffix('')
        hic1 = str(Path(hic1).absolute())
        hic2 = str(Path(hic2).absolute())
       
        while pairs_prefix.suffix in {'.fastq', 'gz', 'fq', '.fq', '.gz', '_R1', '_1', '_2'}:
            pairs_prefix = pairs_prefix.with_suffix('')

        pairs_prefix = str(pairs_prefix).replace('_R1', '').replace('_1', '')
        pairs = f"{pairs_prefix}.pairs.gz"
        hg_input = f"{pairs_prefix}.pairs.gz"
        porec_table = None
        hg_flag = ""
        input_param = "--pairs"

        if not Path(pairs).exists() and "0" not in skip_steps:
            steps.add("0")
        
        if Path(pairs).exists():
            if is_compressed_table_empty(pairs):
                logger.info(f"The existing pairs `{pairs}` is empty, rerun step 0.mapper.")
                steps.add("0")
                os.remove(pairs)

    else:
        if "0" in steps:
            steps = steps - set("0")
            if "0" not in skip_steps:
                logger.warning("Mapping step will not be run, because not specified porec data with `-pcd` or hic data with `-hic1` and `-hic2`")
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
    

    contigsizes = f"{fasta_prefix}.contigsizes"
    if not Path(contigsizes).exists() or is_empty(contigsizes):
        
        cmd = ["cphasing-rs", "chromsizes", fasta, "-o", contigsizes]
        flag = run_cmd(cmd, log=os.devnull)
        assert flag == 0, "Failed to execute command, please check log."
        
    if porec_table is not None and not Path(pairs).exists():
        logger.info("Generating pairs file ...")
        cmd = ["cphasing-rs", "porec2pairs", porec_table, contigsizes,
                    "-o", pairs, "-q", mapping_quality]
        
        flag = run_cmd(cmd, log=f"logs/porec2pairs.log")
        assert flag == 0, "Failed to execute command, please check log."


    corrected = False
    if chimeric_correct or chimeric_corrected:
        from ..chimeric import run as chimeric_run 
        correct_dir = Path("0_2.correct")
        correct_dir.mkdir(exist_ok=True)
        correct_dir = str(correct_dir)
        os.chdir(correct_dir)
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
                    cmd = ["cphasing-rs", "porec-break", f"../{porec_table}", 
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
                corrected_items = chimeric_run(f"../{fasta}", f"../{pairs}", break_pairs=True, 
                                            outprefix=fasta_prefix, 
                                            low_memory=low_memory, threads=threads)
                
                if corrected_items:
                    corrected = True
                    break_bed, fasta, pairs = corrected_items
                    cmd = ["cphasing-rs", "porec-break", f"../{porec_table}", 
                            break_bed, "-o", f"{porec_prefix}.corrected.porec.gz"]
                    flag = run_cmd(cmd, log="logs/porec-break.log")
                    assert flag == 0, "Failed to execute command, please check log."
                    porec_table = f"{porec_prefix}.corrected.porec.gz"
                    hg_input = porec_table

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
                corrected_items = chimeric_run(f"../{fasta}", f"../{pairs}", break_pairs=True, 
                                                outprefix=fasta_prefix, 
                                                low_memory=low_memory, threads=threads)
            
            if corrected_items:
                break_bed, fasta, pairs = corrected_items
                corrected = True
                hg_input = pairs
                fasta_prefix = Path(Path(fasta).name).with_suffix("")
                while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
                    fasta_prefix = fasta_prefix.with_suffix("")

                pairs_prefix = Path(Path(pairs).stem).with_suffix('')
                while pairs_prefix.suffix in {'.pairs', '.gz'}:
                    pairs_prefix = pairs_prefix.with_suffix('')
        
        
        os.chdir("..")
        
        if corrected:
            if porec_table:
                try:
                    Path(porec_table).symlink_to(f"{correct_dir}/{porec_table}")
                except FileExistsError:
                    pass
            try:
                Path(pairs).symlink_to(f"{correct_dir}/{pairs}")
                Path(fasta).symlink_to(f"{correct_dir}/{fasta}")
            except FileExistsError:
                pass
        else:
            pass
            # logger.info(f"Do not need correct, removed `{correct_dir}`")
            # shutil.rmtree(correct_dir)
    

    ## filtered contig by min length
    filtered_contig_by_min_length = False
    if filtered_contig_by_min_length:
        cmd = f"seqkit seq -m {min_length} {fasta} 2>/dev/null > ml{min_length}.{fasta}"
        flag = os.system(cmd)
        assert flag == 0, "Failed to execute command, please check log."

        fasta = f"ml{min_length}.{fasta}"
        fasta_prefix = Path(Path(fasta).name).with_suffix("")
        while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
            fasta_prefix = fasta_prefix.with_suffix("")
        
        contigsizes = f"{fasta_prefix}.contigsizes"
        if not Path(contigsizes).exists() or is_empty(contigsizes):
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

    
    if hcr_flag:
        if porec_table and not use_pairs:
            prepare_input = f"{porec_prefix}.pairs.gz"
        else:
            prepare_input = f"{pairs_prefix}.pairs.gz"
            input_param = "--pairs"

    else:
        prepare_input = pairs
        if not porec_table and min_quality1 > 0: 
            if low_memory or mapping_quality > 0:
                hg_input = f"{pairs_prefix}.q{mapping_quality}.pairs.gz"
                filtered_pairs = hg_input
                if not Path(hg_input).exists():
                    cmd = ["cphasing-rs", "pairs-filter", prepare_input, 
                        "-o", hg_input , "-q", str(mapping_quality)]
                    flag = run_cmd(cmd, log=f'{log_dir}/pairs_filter.log')
                    assert flag == 0, "Failed to execute command, please check log."
                else:
                    logger.warning(f"Using exists filtered pairs file of `{hg_input}`")
                prepare_input = hg_input

    if porec_table or porec_data:
        if (not Path(prepare_input).exists() and (hic1 is None )) or ((hic1 is None) and is_compressed_table_empty(prepare_input)):
            logger.info("Generating pairs file ...")
            cmd = ["cphasing-rs", "porec2pairs", porec_table, contigsizes,
                        "-o", pairs, "-q", "0"]
            
            flag = run_cmd(cmd, log=f"logs/porec2pairs.log")
            assert flag == 0, "Failed to execute command, please check log."


    alleles_dir = str("1.alleles")
    
    genomesize = read_chrom_sizes(contigsizes)['length'].sum()
    logger.info(f"Total size of contig-level assembly: {to_humanized2(genomesize)}")
    if "1" not in skip_steps and "1" in steps:
        logger.info("""
#----------------------------------#
#      Running step 1. alleles     #
#----------------------------------#""")
        Path(alleles_dir).mkdir(exist_ok=True)
        os.chdir(alleles_dir)

        args = ["-f",
                f"../{fasta}",
                "-k",
                alleles_kmer_size,
                "-w",
                alleles_window_size,
                "-m",
                alleles_minimum_similarity,
                "-d",
                alleles_diff_thres,
                "-t",
                threads]
        
        with open('alleles.cmd.sh', 'w') as _out_sh:
            _out_sh.write("cphasing alleles ")
            _out_sh.write(" ".join(map(str, args)))
            _out_sh.write("\n")
        
        try:
            alleles.main(args=args,
                            prog_name='alleles')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
        
        os.chdir("../")

    allele_table = None if mode == "basal" else f"{alleles_dir}/{fasta_prefix}.allele.table" 


    prepare_prefix = prepare_input.replace(".gz", "").rsplit(".", 1)[0]
    prepare_dir = str("2.prepare")
    
    if "2" not in skip_steps and "2" in steps:
        logger.info("""
#----------------------------------#
#       Running step 2. prepare    #
#----------------------------------#""")
        Path(prepare_dir).mkdir(exist_ok=True)
        os.chdir(prepare_dir)

        if porec_table:
            args = [f"../{fasta}",
                    f"../{prepare_input}", 
                    "-p",
                    pattern, 
                    "-q",
                    mapping_quality,
                    "-t",
                    threads,
                    "--skip-pairs2contacts"]
        else:
            args = [f"../{fasta}",
                    f"../{prepare_input}", 
                    "-p",
                    pattern,
                    "-q",
                    mapping_quality,
                    "-t",
                    threads,
                    "--skip-pairs2contacts"]
            
        with open("prepare.cmd.sh", 'w') as _out_sh:
            _out_sh.write("cphasing prepare ")
            _out_sh.write(" ".join(map(str, args)))
            _out_sh.write("\n")
                
        try:
            prepare.main(args=args,
                        prog_name='prepare')
            
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e

        os.chdir("..")
    
    # if not porec_table:
    #     if mode == 'phasing' or mode == 'basal_withprune':
    #         if not Path(f"{prepare_prefix}.q{min_quality2}.contacts").exists():
    #             cmd = ["cphasing-rs", "pairs2contacts", str(hg_input), 
    #                     "-q", str(min_quality2), 
    #                 "-o", f"{prepare_prefix}.q{min_quality2}.contacts" ]
    #             flag = run_cmd(cmd, log=f'{log_dir}/prepare.pairs2contacts.log')
    #             assert flag == 0, "Failed to execute command, please check log."
    #         else:
    #             logger.warning(f"Use exists contacts of `{prepare_prefix}.q{min_quality2}.contacts`")

    #         contacts = f"{prepare_prefix}.q{min_quality2}.contacts"
    #     else:
    #         if not Path(f"{prepare_prefix}.q{min_quality1}.contacts").exists():
    #             cmd = ["cphasing-rs", "pairs2contacts", str(hg_input), 
    #                     "-q", str(min_quality1),
    #                 "-o", f"{prepare_prefix}.q{min_quality1}.contacts" ]
    #             flag = run_cmd(cmd, log=f'{log_dir}/prepare.pairs2contacts.log')
    #             assert flag == 0, "Failed to execute command, please check log."
    #         else:
    #             logger.warning(f"Use exists contacts of `{prepare_prefix}.q{min_quality1}.contacts`")
    #         contacts = f"{prepare_prefix}.q{min_quality1}.contacts"
   
    count_re = f"{prepare_dir}/{prepare_prefix}.counts_{pattern}.txt"
    clm = f"{prepare_dir}/{prepare_prefix}.clm.gz"
    
    
    split_contacts = f"{prepare_dir}/{prepare_prefix}.split.contacts"

    output_cluster = "output.clusters.txt"
    
    hyperpartition_normalize = "-norm" if normalize else None
    enable_misassembly_remove = "--enable-misassembly-remove" if enable_misassembly_remove else None 


    if hg_flag == "--pairs":
        # hyperpartition_contacts = contacts
        hyperpartition_contacts = None
    else:
        hyperpartition_contacts = None

    if use_pairs and porec_table and not hcr_flag:
        hg_input = pairs 
        hg_flag == "--pairs"
        input_param = "--pairs"
    
    hyperpartition_dir = str("3.hyperpartition")
    
    if "3" not in skip_steps and "3" in steps:
        logger.info("""
#----------------------------------#
#  Running step 3. hyperpartition  #
#----------------------------------#""")
        Path(hyperpartition_dir).mkdir(exist_ok=True)
        os.chdir(hyperpartition_dir)

        if hcr_flag and "3" not in skip_steps and "3" in steps:
            hcr_invert_string = "-v" if hcr_invert else ""

            try:
                if hcr_invert_string:
                    args = ["-prs",
                                f"../{pairs}",
                                "-cs",
                                f"../{contigsizes}",
                                "-p", pattern,
                                "-f", f"../{fasta}",
                                "-b",
                                hcr_bed,
                                "-u",
                                hcr_upper,
                                "-l",
                                hcr_lower,
                                # "-q",
                                # min_quality1,
                                hcr_invert_string]
                    hcr.main(
                        args=args,
                        prog_name="hcr"
                    )
                else:
                    args = ["-prs",
                                f"../{pairs}",
                                "-cs",
                                f"../{contigsizes}",
                                "-p", pattern,
                                "-f", f"../{fasta}",
                                "-b",
                                hcr_bed,
                                "-u",
                                hcr_upper,
                                "-l",
                                hcr_lower,
                                # "-q",
                                # min_quality1,
                                ]
                    hcr.main(
                        args=args,
                        prog_name="hcr"
                    )

            except SystemExit as e:
                exc_info = sys.exc_info()
                exit_code = e.code
                if exit_code is None:
                    exit_code = 0
                
                if exit_code != 0:
                    raise e
                
            hcr_bed = f"{pairs_prefix}.{hcr_bs}.hcr.bed"


        hyperpartition_args = [
                                f"../{hg_input}",
                                f"../{contigsizes}",
                                output_cluster,
                                input_param,
                                "--mode",
                                mode,
                                "-at",
                                f"../{allele_table}" if allele_table else None,
                                "-c",
                                hyperpartition_contacts,
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
                                "--exclude-group-from-first",
                                exclude_group_from_first,
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
                                threads,
                                "-n",
                                n,
                            ]
        
        if disable_merge_in_first:
            hyperpartition_args.append("--disable-merge-in-first")
        if enable_misassembly_remove:
            hyperpartition_args.append(enable_misassembly_remove)
        if hyperpartition_normalize:
            hyperpartition_args.append(hyperpartition_normalize)

    
        if hcr_bed:
            hyperpartition_args.extend(["--hcr-bed", f"{hcr_bed}"])
            if hcr_invert:
                hyperpartition_args.append("--hcr-invert")

        with open("hyperpartition.cmd.sh", 'w') as _out_sh:
            _out_sh.write('cphasing hyperpartition \\\n    ')
            _hyperpartition_args = []
            remove_index = []
            for i in range(4, len(hyperpartition_args), 2):
                try:
                    if hyperpartition_args[i+1] is None:
                        remove_index.append(i)
                        remove_index.append(i+1)
                except IndexError:
                    continue
            
            for i, j in enumerate(hyperpartition_args):
                if i not in remove_index:
                    _hyperpartition_args.append(j)

            _out_sh.write(" ".join(pretty_cmd(map(str, _hyperpartition_args))))
            _out_sh.write("\n")

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

        os.chdir("../")


    out_agp = "groups.agp"
    if corrected:
        corrected_agp = out_agp.replace("agp", "corrected.agp")
    
    scaffolding_dir = str("4.scaffolding")

    if "4" not in skip_steps and "4" in steps:
        logger.info("""
#----------------------------------#
#    Running step 4. scaffolding   #
#----------------------------------#""")
        Path(scaffolding_dir).mkdir(exist_ok=True)
        os.chdir(scaffolding_dir)

        if corrected:
            args = [
                f"../{hyperpartition_dir}/{output_cluster}",
                f"../{count_re}",
                f"../{clm}",
                " ",
                "-at",
                f"../{allele_table}" if allele_table else None,
                "-sc",
                f"../{split_contacts}",
                "-f",
                f"../{raw_fasta}",
                "-t",
                threads,
                "-o",
                out_agp,
                "-m",
                scaffolding_method,
                "--corrected"]

        else:
            args=[
                f"../{hyperpartition_dir}/{output_cluster}",
                f"../{count_re}",
                f"../{clm}",
                " ",
                "-at",
                f"../{allele_table}" if allele_table else None,
                "-sc",
                f"../{split_contacts}",
                "-f",
                f"../{fasta}",
                "-t",
                threads,
                "-o",
                out_agp,
                "-m",
                scaffolding_method]
            
        with open("scaffolding.cmd.sh", "w") as _out_sh:
            _out_sh.write("cphasing scaffolding \\\n    ")
            _out_sh.write(" ".join(pretty_cmd(map(str, args), n=2)))
            _out_sh.write("\n")

        args = list(filter(lambda x: x != " ", args))
        try:
            scaffolding.main(args=args,
                            prog_name='scaffolding')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
        
        # shutil.copy('groups.agp', '..')
        generate_to_hic_cmd('groups.agp', f"../{pairs}", n=n )
        if corrected:
            # shutil.copy('groups.corrected.agp', '..')
            generate_to_hic_cmd('groups.corrected.agp', f"../{pairs}", n=n)
        
        os.chdir("..")

   

    plot_dir = str("5.plot")
    
    if "5" not in skip_steps and "5" in steps:
        logger.info("""
#----------------------------------#
#      Running step 5. plot        #
#----------------------------------#""")
        Path(plot_dir).mkdir(exist_ok=True)

        os.chdir(plot_dir)

        if binsize == "auto" or binsize is None:
            init_binsize, binsize = recommend_binsize_by_genomesize(genomesize)
            
            logger.info(f"Recommended cool's binsize: `{to_humanized2(init_binsize)}`, heatmap's binsizs: `{to_humanized2(binsize)}`")
        else:
            init_binsize = 10000

        out_small_cool = f"{pairs_prefix}.q{min_quality1}.{to_humanized2(init_binsize)}.cool"

        if not is_file_changed(pairs) and Path(out_small_cool).exists():
            logger.warning(f"`{out_small_cool}` exists, skipped `pairs2cool`.")
            args = [
                    f"../{pairs}",
                    f"../{contigsizes}",
                    out_small_cool,
                    "-q", 
                    min_quality1,
                    "-bs",
                    init_binsize,
                    ]

        else:
            if filtered_pairs:
                pairs = filtered_pairs
            
            args = [
                    f"../{pairs}",
                    f"../{contigsizes}",
                    out_small_cool,
                    "-q", 
                    min_quality1,
                    "-bs",
                    init_binsize,
                    ]

            try:
                
                pairs2cool.main(args=args,
                                prog_name='pairs2cool')
            except SystemExit as e:
                exc_info = sys.exc_info()
                exit_code = e.code
                if exit_code is None:
                    exit_code = 0
                
                if exit_code != 0:
                    raise e
                
        
        if corrected:
            if Path(f"../{scaffolding_dir}/{corrected_agp}").exists():
                input_agp = corrected_agp
            else:
                input_agp = out_agp
        else:
            input_agp = out_agp

        generate_plot_cmd(f"../{pairs}",
                        f"{pairs_prefix}",
                        f"../{contigsizes}",
                        f"../{scaffolding_dir}/{input_agp}", 
                        min_quality1, init_binsize,
                        binsize)
        args = [
                "-a",
                f"../{scaffolding_dir}/{input_agp}",
                "-m",
                out_small_cool,
                "-o",
                f"groups.q{min_quality1}.{to_humanized2(binsize)}.wg.png",
                "-bs",
                binsize,
                "-oc",
                ]
        
       
        try:
            plot.main(args=args,
                            prog_name='plot')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
            
        os.chdir("../")


    today = date.today().strftime("%Y-%m-%d")
    end_time = time.time() - start_time
    # peak_memory = tracemalloc.get_traced_memory()[1] / 1024 / 1024 / 1024
    # peak_memory = 0
    logger.info(f"Pipeline finished in {today}. Elapsed time {end_time:.2f} s. ") #Peak memory: {peak_memory:2f} Gb")
    logger.info(f"Results are store in `{outdir}`")
    # tracemalloc.stop()