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
import glob
import time 

import pandas as pd

from datetime import datetime, date
from pathlib import Path

from .. import __version__
from ..cli import (mapper as porec_mapper,
                    hcr,
                    alleles, 
                    prepare,
                    hyperpartition,
                    rescue as collapsed_rescue_cli,
                    from_gfa as collapsed_from_gfa,
                    pairs_dup,
                    agp_dup,
                    agp2fasta,
                    scaffolding,
                    pairs2cool,
                    pairs2pqs,
                    plot
)
from ..core import Pairs2
from .._config import *
from ..expection import *
from ..hic.cli import mapper as hic_mapper
from ..pqs import PQS
from ..utilities import (
    pretty_cmd,
    run_cmd, 
    calculate_Nx_from_contigsizes,
    to_humanized2,
    to_humanized3,
    recommend_binsize_by_genomesize,
    read_chrom_sizes,
    is_empty,
    is_compressed_table_empty,
    is_file_changed,
    generate_to_hic_cmd,
    generate_plot_cmd,
    generate_curation_cmd,
    MemoryMonitor
    )

monitor = MemoryMonitor()
logger = logging.getLogger(__name__)


def run(fasta,
        gfa,
        ul_data,
        porec_data,
        paf,
        porec_table,
        pairs, 
        hic1=None,
        hic2=None,
        pattern="AAGCTT",
        mapper_k=15,
        mapper_w=10,
        aligner="minimap2",
        mm2_params="-x map-ont",
        mapping_quality=0,
        hic_aligner="_chromap",
        hic_mapper_k=None,
        hic_mapper_w=None,
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
        alleles_trim_length=25000,
        kprune_norm_method="auto",
        scaffolding_method="precision",
        enable_haplotype_cluster=False,
        n="",
        use_pairs=False,
        edge_length=2e6,
        split_length='auto',
        resolution1=1,
        resolution2=1,
        init_resolution1=1,
        init_resolution2=1,
        first_cluster=None,
        normalize=False,
        allelic_factor=-1,
        disable_merge_in_first=False,
        disable_merge_use_allele=True,
        allelic_positive_factor=3.0,
        exclude_group_to_second=None,
        exclude_group_from_first=None,
        allelic_similarity=0.85,
        min_allelic_overlap=0.3,
        min_weight=0.1,
        min_cis_weight=1.0,
        min_quality1=1,
        min_quality2=2,
        min_contacts=5,
        min_length=10000,
        Nx=100,
        min_scaffold_length=5e6,
        cluster_method="louvain",
        enable_misassembly_remove=False,
        disable_recluster_by_linkage=False,
        refine=False,
        whitelist=None,
        blacklist=None,
        collapsed_rescue=False,
        collapsed_contigs=None,
        disable_gfa_collapsed=False,
        disable_conflict_check=False,
        cn_offset=0.0,
        binsize="500k",
        colormap="viridis",
        whitered=False,
        balance=False,
        low_memory=False,
        output_hg=False,
        outdir="cphasing_output",
        threads=4):
    
    is_pairs2pqs = False
    start_time = time.time()
    logger.info(f"C-Phasing version: {__version__}")
    today = date.today().strftime("%Y-%m-%d")
    logger.info(f"Pipeline is started on {today}.")
    
    outdir = Path(outdir)
    if outdir.exists():
        logger.info(f"Working on existing directory: `{outdir}`")
    else:
        logger.info(f"Working on new directory: `{outdir}`")
        outdir.mkdir(parents=True, exist_ok=True)
    
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
    

    fasta_path = str(Path(fasta).absolute())
    fasta = str(Path(fasta).name)
    if gfa:
        gfa = str(Path(gfa).absolute())

    if ul_data:
        ul_data = str(Path(ul_data).absolute())
    
    if porec_data:
        porec_data = list(porec_data)
        for i, j in enumerate(porec_data):
            porec_data[i] = str(Path(j).absolute())
    
    if hic1 and hic2:
        hic1 = list(hic1) 
        hic2 = list(hic2)
        for i, j in enumerate(hic1):
            hic1[i] = str(Path(j).absolute())
        for i, j in enumerate(hic2):
            hic2[i] = str(Path(j).absolute())

    if paf:
        paf_path = Path(paf).absolute()
    
    
    if porec_table:
        porec_table_path = Path(porec_table).absolute()
        if not pairs:
            pairs_prefix = porec_table.replace(".gz", "").rsplit(".", 1)[0]
            _pqs = Path(f"{pairs_prefix}.pairs.pqs")
            _pqs2 = Path(f"{outdir}/{_pqs}")
            if _pqs.exists():
                _p = PQS(str(_pqs))
                if _p.is_pqs(str(_pqs)):
                    is_pairs2pqs = True
                    pqs_file = _pqs
                    _pairs = _pqs
                else:
                    _pairs = Path(f"{pairs_prefix}.pairs.gz")
            elif _pqs2.exists():
                _p = PQS(str(_pqs2))
                if _p.is_pqs(str(_pqs2)):
                    is_pairs2pqs = True
                    pqs_file = _pqs
                    _pairs = _pqs
                else:
                    _pairs = Path(f"{pairs_prefix}.pairs.gz")
            else:
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

    if blacklist:
        blacklist = str(Path(blacklist).absolute())

    if collapsed_contigs:
        collapsed_contigs = str(Path(collapsed_contigs).absolute())
        collapsed_rescued_contigs = f"{str(Path.cwd())}/{outdir}/3.hyperpartition/collapsed.rescue.contigs.list"

    os.chdir(outdir)

    try:
        Path(fasta).symlink_to(fasta_path)
    except FileExistsError:
        pass

    def safe_remove(target):
        target_str = str(target)
        try:
            if not os.path.lexists(target_str):
                return
            
            if os.path.islink(target_str) or os.path.isfile(target_str):
                os.unlink(target_str)
            elif os.path.isdir(target_str):
                shutil.rmtree(target_str, ignore_errors=True)
        except Exception:
            pass
    
    mapper_dir = Path("0.mapper")
    mapper_dir.mkdir(exist_ok=True)
    Path("logs").mkdir(exist_ok=True)
    Path(f"{mapper_dir}/logs").mkdir(exist_ok=True, parents=True)

    
    if paf:
        m_paf = mapper_dir / paf_path.name
        if m_paf.resolve() != paf_path.resolve():
            if m_paf.exists() or m_paf.is_symlink(): 
                safe_remove(m_paf)
            try: 
                m_paf.symlink_to(paf_path)
            except Exception: 
                pass

        target = Path("input.paf.gz")
        if target.exists() or target.is_symlink(): 
            safe_remove(target)
        try: 
            target.symlink_to(m_paf.resolve())
        except Exception: 
            pass
        source_paf = paf_path.name
        paf = "input.paf.gz"

    if porec_table:
        m_porec = mapper_dir / Path("input.porec.gz")
        if m_porec.resolve() != porec_table_path.resolve():
            if m_porec.exists() or m_porec.is_symlink(): 
                safe_remove(m_porec)
            try: 
                m_porec.symlink_to(porec_table_path)
            except Exception: 
                pass

        target = Path("input.porec.gz")
        if target.exists() or target.is_symlink(): 
            safe_remove(target)
        try: 
            target.symlink_to(m_porec.resolve())
        except Exception: 
            pass
        porec_table = "input.porec.gz"
    
    if pairs:
        source_pairs_name = pairs_path.name
        m_pairs = mapper_dir / source_pairs_name
        
        if m_pairs.resolve() != pairs_path.resolve():
            if m_pairs.exists() or m_pairs.is_symlink(): 
                safe_remove(m_pairs)
            try: 
                m_pairs.symlink_to(pairs_path)
            except Exception: 
                pass
        
        target_name = "input.pairs.pqs" if str(pairs).endswith(".pqs") else "input.pairs.gz"
        target = Path(target_name)
        if target.exists() or target.is_symlink(): 
            safe_remove(target)
        try: 
            target.symlink_to(m_pairs.resolve())
        except Exception: 
            pass 
        pairs = target_name


    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    steps = set(steps)

    raw_fasta = fasta

    mode = 'basal' if mode == 'haploid' else mode
    if mode == 'basal':
        skip_steps.add("1")
        allele_table = None
        logger.info("The mode is `haploid`, skip step '1. alleles.'")
    elif mode == "phasing":
        if "1" in steps:
            skip_steps.add("1")
            allele_table = None
            logger.info("The mode is `phasing`, the step '1. alleles' will be intergated into '3. hyparpartion'.")
        else:
            allele_table = None
    
    if mode == "basal_withprune":
        split_length = 'none'
        logger.warning("The mode is `basal_withprune`, the `split_length` is set to `none`. If you want to split contigs, please using the parameters `-n 1:x --mode phasing`")
        
         
    filtered_pairs = None
   
    if ul_data:
        from ..hitig.pipeline import run  as hitig_run
        input_fasta = str(Path(fasta).absolute())
        input_ul_data = str(Path(ul_data).absolute())              
        if not use_existed_hitig:
            Path("0.1.hitig").mkdir(exist_ok=True)
            os.chdir("0.1.hitig")
            hitig_run(input_fasta, input_ul_data)
           
            os.chdir("..")

            hcr_bed = Path("0.1.hitig/output.hcr_all.bed").absolute()
            # clean_fasta = "hitig/output.cleaned.fasta"
        else:
            logger.warning("Use existed hitig results.")
            hcr_bed = Path("0.1.hitig/output.hcr_all.bed").absolute()
            # clean_fasta = "hitig/output.cleaned.fasta"

        if not is_empty("0.1.hitig/output.breakPos.txt"):
            fasta = "0.1.hitig/output.corrected.fasta"
        # if Path(hcr_bed).exists() and Path(clean_fasta).exists():
        #     fasta = clean_fasta
        # else:
        #     logger.warning("Skipped hitig.")

            

    fasta_prefix = Path(Path(fasta).name).with_suffix("")
    while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
        fasta_prefix = fasta_prefix.with_suffix("")

    contigsizes = f"{fasta_prefix}.contigsizes"
    if not Path(contigsizes).exists() or is_empty(contigsizes):
        
        cmd = ["cphasing-rs", "chromsizes", fasta, "-o", contigsizes]
        flag = run_cmd(cmd, log=os.devnull)
        assert flag == 0, "Failed to execute command, please check log."

    # if porec_data:
    #     porec_prefix = str(Path(porec_data[0]).name).replace(".gz", "").rsplit(".", 1)[0]
    #     pairs_prefix = str(Path(porec_data[0]).name).replace(".gz", "").rsplit(".", 1)[0]
    #     paf = f"{porec_prefix}.paf.gz"
    #     pairs = f"{pairs_prefix}.pairs.pqs"
    #     is_pairs2pqs = True
    #     pqs_file = f"{pairs_prefix}.pairs.pqs"
    #     hg_input = f"{porec_prefix}.porec.gz"
    #     hg_flag = ""
    #     porec_table = hg_input
    #     input_param = "--porec"
        
        
    #     if not Path(f"{mapper_dir}/{pairs}").exists() or not Path(f"{mapper_dir}/{porec_table}").exists() and "0" not in skip_steps:
    #         steps.add("0")

    #     if Path(paf).exists():
    #         if is_compressed_table_empty(paf):
    #             logger.info(f"The existing paf `{paf}` is empty, rerun step 0.mapper.")
    #             steps.add("0")
                
    #             os.remove(paf)

    #     if Path(porec_table).exists():
    #         if is_compressed_table_empty(porec_table):
    #             logger.info(f"The existing porec table `{porec_table}` is empty, rerun step 0.alignment.")
    #             steps.add("0")
                
    #             os.remove(porec_table)
    # elif paf:
    #     porec_prefix = paf.replace(".gz", "").rsplit(".", 1)[0]
    #     pairs_prefix = paf.replace(".gz", "").rsplit(".", 1)[0]

    #     expected_porec = f"{porec_prefix}.porec.gz"
    #     expected_pqs = f"{pairs_prefix}.pairs.pqs"
    #     original_porec_path = paf_path.parent / expected_porec
    #     original_pqs_path = paf_path.parent / expected_pqs

    #     if original_porec_path.exists():
    #         try:
    #             Path(expected_porec).symlink_to(original_porec_path)
    #             logger.info(f"Existing porec table found, linked: {expected_porec}")
    #         except FileExistsError:
    #             pass
        
    #     if original_pqs_path.exists():
    #         try:
    #             Path(expected_pqs).symlink_to(original_pqs_path)
    #             logger.info(f"Existing pairs pqs found, linked: {expected_pqs}")
    #         except FileExistsError:
    #             pass

    #     pqs_file = expected_pqs
    #     is_pairs2pqs = True
    #     hg_input = expected_porec
    #     hg_flag = ""
    #     porec_table = hg_input
    #     input_param = "--porec"
    #     skip_steps.add("0")
    #     pairs = pqs_file
    #     if not Path(porec_table).exists():
    #         logger.info(f"Generating porec table from existing paf: {paf}")
    #         cmd = ["cphasing-rs", "paf2porec", paf, "-o", porec_table]
    #         run_cmd(cmd, log=f"logs/paf2porec.log")
        
    
    # elif hic1 and hic2:
    #     pairs_prefix = Path(Path(hic1).stem).with_suffix('')
    #     hic1 = str(Path(hic1).absolute())
    #     hic2 = str(Path(hic2).absolute())
       
    #     while pairs_prefix.suffix in {'.fastq', 'gz', 'fq', '.fq', '.gz', '_R1', '_1', '_2'}:
    #         pairs_prefix = pairs_prefix.with_suffix('')

    #     pairs_prefix = str(pairs_prefix).replace('_R1', '').replace('_1', '')
       
    #     pairs = f"{pairs_prefix}.pairs.pqs"
    #     pqs_file = f"{pairs_prefix}.pairs.pqs"
    #     is_pairs2pqs = False
    #     hg_input = f"{pairs_prefix}.pairs.pqs"
    #     porec_table = None
    #     hg_flag = ""
    #     input_param = "--pairs"

    #     if not Path(f"{mapper_dir}/{pairs}").exists() and not Path(f"{mapper_dir}/{pqs_file}").exists() and "0" not in skip_steps:
    #         steps.add("0")
        
    #     if Path(pqs_file).exists() and not is_compressed_table_empty(pqs_file):
    #         pairs = pqs_file
    #         is_pairs2pqs = True
            

    #     if Path(pairs).exists():
    #         if is_compressed_table_empty(pairs):
    #             logger.info(f"The existing pairs `{pairs}` is empty, rerun step 0.mapper.")
    #             steps.add("0")
    #             shutil.rmtree(pairs)

    # else:
    #     if "0" in steps:
    #         if not paf:
    #             steps = steps - set("0")
    #             if "0" not in skip_steps:
    #                 logger.warning("Mapping step will not be run, because not specified porec data with `-pcd` or hic data with `-hic1` and `-hic2`")
    #     if porec_table:
    #         porec_prefix = porec_table.replace(".gz", "").rsplit(".", 1)[0]
    #         pairs_prefix = porec_table.replace(".gz", "").rsplit(".", 1)[0]
    #         hg_input = porec_table
    #         hg_flag = ""
    #         input_param = "--porec"
    #         if not pairs:
    #             _pqs = Path(f"{porec_prefix}.pairs.pqs")
    #             if _pqs.exists():
    #                 _pairs = _pqs
    #             elif Path(f"{porec_prefix}.pairs.gz").exists():
    #                 _pairs = Path(f"{porec_prefix}.pairs.gz")

    #             else:
    #                 _pairs = None

    #             pairs = _pairs
    #             if use_pairs:
    #                 hg_input = pairs
    #                 hg_flag = "--pairs"
    
    #     else:
    #         if pairs:
    #             pairs_prefix = Path(Path(pairs).name)
    #             while pairs_prefix.suffix in {'.pairs', '.gz', ".pqs"}:
    #                 pairs_prefix = pairs_prefix.with_suffix('')
    #             pairs_prefix = str(pairs_prefix)
    #             hg_input = pairs 
    #             hg_flag = "--pairs"
    #             input_param = hg_flag
    mapper_dir = Path("0.mapper")
    mapper_dir.mkdir(exist_ok=True)
    if "0" not in skip_steps and "0" in steps:
        if not hic:
            logger.info("""#----------------------------------#
#      Running step 0. mapper      #
#----------------------------------#""")
         
            os.chdir(mapper_dir)
            if porec_data:
                try:
                    porec_mapper.main(args=[f"../{fasta}", 
                                    *porec_data,
                                    "-k",
                                    mapper_k,
                                    "-w",
                                    mapper_w,
                                    "-q",
                                    mapping_quality,
                                    "--aligner",
                                    aligner,
                                    "--mm2-params",
                                    mm2_params,
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
                
                porec_name = Path(porec_data[0]).name
                clean_name = re.sub(r'(\.gz|\.bam|\.fastq|\.fq|\.fa|\.fasta)$', '', porec_name)
                porec_prefix = clean_name.rsplit(".", 1)[0] if "." in clean_name else clean_name
                pairs_prefix = porec_prefix
                source_paf = str(Path(f"{porec_prefix}.paf.gz").absolute())
                source_pairs_name = str(Path(f"{pairs_prefix}.pairs.pqs").absolute())
            
            if hic1 and hic2:
                hic1_args = []
                for h in hic1:
                    hic1_args.extend(["-1", h])
                hic2_args = []
                for h in hic2:
                    hic1_args.extend(["-2", h])
                try:
                    hic_mapper.main(args=[
                                    "-r",
                                    f"../{fasta}", 
                                    *hic1_args,
                                    *hic2_args,
                                    "-a",
                                    hic_aligner,
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
                    
                pairs_prefix = Path(Path(hic1[0]).stem).with_suffix('')
                while pairs_prefix.suffix in {'.fastq', 'gz', 'fq', '.fq', '.gz', '_R1', '_1', '_2'}:
                    pairs_prefix = pairs_prefix.with_suffix('')
                pairs_prefix = str(pairs_prefix).replace('_R1', '').replace('_1', '')
                source_pairs_name = str(Path(f"{pairs_prefix}.pairs.pqs").absolute())

                    
            os.chdir("..")

    porec_prefix = "input"
    pairs_prefix = "input"
    prepare_prefix = "input"

    if paf:

        if not Path(f"{mapper_dir}/input.porec.gz").exists() or Path(f"{mapper_dir}/input.porec.gz").stat().st_size == 0:
            logger.info("Generating porec table from PAF ...")
            os.chdir(mapper_dir)
            cmd = ["cphasing-rs", "paf2porec", f"{source_paf}", "-q", "0", "-o", "input.porec.gz"]
            flag = run_cmd(cmd, log="logs/paf2porec.log")
            assert flag == 0, "Failed to execute paf2porec, please check log."
            os.chdir("..")
        if not Path(f"{mapper_dir}/input.pairs.pqs").exists():
            logger.info("Generating pairs file from porec ...")
            os.chdir(mapper_dir)
            cmd = ["cphasing-rs", "porec2pairs", "input.porec.gz", f"../{contigsizes}",
                        "-o", "input.pairs.pqs", "-q", str(mapping_quality)]
            flag = run_cmd(cmd, log=f"logs/porec2pairs.log")
            assert flag == 0, "Failed to execute porec2pairs, please check log."
            os.chdir("..")

        if Path(paf).resolve() != Path("input.paf.gz").resolve():
            run_cmd(["ln", "-sfn", str(Path(f"../{source_paf}").absolute()), "input.paf.gz"])
            paf = "input.paf.gz"

        porec_table = "input.porec.gz"
        pqs = "input.pairs.pqs"
        source_pairs_name = pqs
        pairs = pqs

        hg_input = porec_table


    if porec_table is not None and pairs is None:
        logger.info("Generating pairs file ...")
        _pairs = "input.pairs.pqs"
        pqs_file = "input.pairs.pqs"
        source_pairs_name = _pairs
        is_pairs2pqs = True
        os.chdir(mapper_dir)
        if not Path(_pairs).exists():
            cmd = ["cphasing-rs", "porec2pairs", porec_table, f"../{contigsizes}",
                        "-o", str(_pairs), "-q", str(mapping_quality)]

            flag = run_cmd(cmd, log=f"logs/porec2pairs.log")
            assert flag == 0, "Failed to execute command, please check log."
        else:
            logger.warning(f"The pairs file `{_pairs}` already exists, skipped generating pairs file.")

        pairs = str(_pairs)
        os.chdir("..")
    
    if pairs:
        os.chdir(mapper_dir)
        try:
            _pairs = Pairs2(source_pairs_name)
            if _pairs.is_pairs():

                is_skip = True if skip_steps.intersection({"2", "3", "4", "5"}) else False
                if is_skip:
                    pairs_changed = False
                else:
                    pairs_changed = is_file_changed(str(source_pairs_name))

                pqs_path = "input.pairs.pqs"
                pqs_exists = Path(pqs_path).exists()
                if not pqs_exists or pairs_changed:
                    logger.info("Coverting pairs or pairs.gz to pairs.pqs ...")
                    args = [source_pairs_name, "-t", str(threads), "-o", pqs_path]
                    try:
                        pairs2pqs.main(args=args, prog_name='pairs2pqs')
                    except SystemExit as e:
                        exc_info = sys.exc_info()
                        exit_code = e.code
                        if exit_code is None:
                            exit_code = 0
                        
                        if exit_code != 0:
                            raise e
                    is_pairs2pqs = True
                    pqs_file = pairs = f"{pairs_prefix}.pairs.pqs"
                    if not porec_table or use_pairs:
                        hg_input = pqs_file
                        hg_flag = "--pairs"
                else:
                    _p = PQS(pqs_path)
                    if _p.is_pqs(pqs_path):
                        logger.info("Use exists pairs.pqs file.")
                        is_pairs2pqs = True
                        pqs_file = pairs = pqs_path
                        if not porec_table or use_pairs:
                            hg_input = pqs_file
                            hg_flag = "--pairs"
                    else:
                        logger.info("The pqs file is incorrect, use the pairs.gz file.")

        except IsNotPairs:
            pass
        os.chdir("..")



    def get_latest_file(pattern):
        files = list(Path("0.mapper").glob(pattern))
        if not files:
            return None
        files.sort(key=lambda x: x.stat().st_mtime, reverse=True)
        if len(files) > 1:
            logger.warning(f"Found multiple `{pattern}` files in `0.mapper/`: {[f.name for f in files]}. "
                           f"Auto-selected the newest one: `{files[0].name}`.")
        return files[0]
    
    _paf_latest = get_latest_file("*.paf.gz")
    _porec_latest = get_latest_file("*.porec.gz")
    _pqs_latest = get_latest_file("*.pairs.pqs")
    _gz_latest = get_latest_file("*.pairs.gz")

    if _paf_latest:
        run_cmd(["ln", "-sfn", f"0.mapper/{_paf_latest.name}", "input.paf.gz"])
        paf = "input.paf.gz"
        if porec_data:
            porec_table = None
            pairs = None

    if _porec_latest:
        run_cmd(["ln", "-sfn", f"0.mapper/{_porec_latest.name}", "input.porec.gz"])
        porec_table = "input.porec.gz"
        hg_input = porec_table

    if _pqs_latest:
        run_cmd(["ln", "-sfn", f"0.mapper/{_pqs_latest.name}", "input.pairs.pqs"])
        pairs = "input.pairs.pqs"
        pqs_file = pairs
        is_pairs2pqs = True
        hg_input = pairs

    elif _gz_latest:
        run_cmd(["ln", "-sfn", f"0.mapper/{_gz_latest.name}", "input.pairs.gz"])
        pairs = "input.pairs.gz"
        hg_input = pairs


    corrected = False
    if "0.2" in skip_steps:
        chimeric_correct = False
        chimeric_corrected = True 
    if chimeric_correct or chimeric_corrected:
        logger.info("""#----------------------------------#
#      Running step 0.2 correct    #
#----------------------------------#""")
        from ..chimeric import run as chimeric_run 
        correct_dir = Path("0.2.correct")
        correct_dir.mkdir(exist_ok=True)
        correct_dir = str(correct_dir)
        os.chdir(correct_dir)
        if porec_table:
            if chimeric_corrected:

                corrected_items = (f"{fasta_prefix}.chimeric.contigs.bed",
                                    f"{fasta_prefix}.corrected.fasta", 
                                   f"{pairs_prefix}.corrected.pairs.pqs")
                if not all(map(lambda x: Path(x).exists(), corrected_items)):
                    corrected_items = ()

                
                if corrected_items and not Path(f"{porec_prefix}.corrected.porec.gz").exists():
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
                    if is_pairs2pqs:
                        pqs_file = pairs
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
                    if is_pairs2pqs:
                        pqs_file = pairs
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
            while pairs_prefix.suffix in {'.pairs', '.gz', '.pqs'}:
                pairs_prefix = pairs_prefix.with_suffix('')
            
            porec_prefix = Path(Path(porec_table).stem).with_suffix('')
            while porec_prefix.suffix in {'.porec', '.gz'}:
                porec_prefix = porec_prefix.with_suffix('')

        elif pairs:
            if chimeric_corrected:
                if is_pairs2pqs:
                    corrected_items = (f"{fasta_prefix}.chimeric.contigs.bed",
                                        f"{fasta_prefix}.corrected.fasta", 
                                        f"{pairs_prefix}.corrected.pairs.pqs")
                else:
                    corrected_items = (f"{fasta_prefix}.chimeric.contigs.bed",
                                        f"{fasta_prefix}.corrected.fasta", 
                                        f"{pairs_prefix}.corrected.pairs.gz")
                if not all(map(lambda x: Path(x).exists(), corrected_items)):
                    corrected_items = ()
                    logger.warning("The corrected files are not exists, please specify `--chimeric-correct`.")
                else:
                    logger.info("Using exists corrected results.")
                
            else:
                corrected_items = chimeric_run(f"../{fasta}", f"../{pairs}", break_pairs=True, 
                                                outprefix=fasta_prefix, min_mapq=1,
                                                low_memory=low_memory, threads=threads)
            
            if corrected_items:
                break_bed, fasta, pairs = corrected_items
                corrected = True
                hg_input = pairs
                # pqs_file = f"{pairs_prefix}.corrected.pairs.pqs"
                fasta_prefix = Path(Path(fasta).name).with_suffix("")
                while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
                    fasta_prefix = fasta_prefix.with_suffix("")

                pairs_prefix = Path(Path(pairs).stem).with_suffix('')
                while pairs_prefix.suffix in {'.pairs', '.gz', '.pqs'}:
                    pairs_prefix = pairs_prefix.with_suffix('')
        
        
        os.chdir("..")
        
        if corrected:
            if porec_table:
                run_cmd(["ln", "-sfn", f"{correct_dir}/{porec_table}", "input.porec.gz"])
                porec_table = "input.porec.gz"
                hg_input = porec_table

            if pairs:
                target_pairs = "input.pairs.pqs" if pairs.endswith(".pqs") else "input.pairs.gz"
                run_cmd(["ln", "-sfn", f"{correct_dir}/{pairs}", target_pairs])
                pairs = target_pairs
                if is_pairs2pqs: 
                    pqs_file = pairs
                
            if fasta:
                run_cmd(["ln", "-sfn", f"{correct_dir}/{fasta}", "input.corrected.fasta"])
                fasta = "input.corrected.fasta"

            contigsizes = f"{fasta_prefix}.contigsizes"
            if not Path(contigsizes).exists() or is_empty(contigsizes):
                cmd = ["cphasing-rs", "chromsizes", fasta, "-o", contigsizes]
                flag = run_cmd(cmd, log=os.devnull)
                assert flag == 0, "Failed to execute command, please check log."
            
            fasta_path = Path(fasta).absolute()
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


    if whitelist or blacklist:
        filtered_fasta = f"{fasta_prefix}.filtered.fasta"
        if not Path(filtered_fasta).exists() or is_file_changed(fasta):
            logger.info(f"Generating filtered FASTA based on whitelist/blacklist: `{filtered_fasta}`")
            filter_cmd = ["seqkit", "grep"]
            if whitelist:
                filter_cmd.extend(["-f", str(whitelist)])
            if blacklist:
                filter_cmd.extend(["-v", "-f", str(blacklist)])
            
            filter_cmd.extend([fasta, "-o", filtered_fasta])
            
            flag = run_cmd(filter_cmd, log=os.devnull)
            assert flag == 0, "Failed to filter FASTA using seqkit."
        
        fasta = filtered_fasta
        fasta_path = Path(fasta).absolute()
        
        fasta_prefix = Path(Path(fasta).name).with_suffix("")
        while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
            fasta_prefix = fasta_prefix.with_suffix("")

        contigsizes = f"{fasta_prefix}.contigsizes"
        cmd = ["cphasing-rs", "chromsizes", fasta, "-o", contigsizes]
        flag = run_cmd(cmd, log=os.devnull)
        assert flag == 0, "Failed to update chromsizes for the filtered FASTA."


    if porec_table or porec_data:
        os.chdir(mapper_dir)
        target_pairs = "input.pairs.pqs"
        is_generate_pairs_from_porec = True
        if Path(target_pairs).exists() and not is_compressed_table_empty(target_pairs):
            is_generate_pairs_from_porec = False
        try:
            if Path(source_pairs_name).exists() and not is_compressed_table_empty(source_pairs_name):
                is_generate_pairs_from_porec = False
        except UnboundLocalError:
            pass

        if is_generate_pairs_from_porec and ((hic1 is None) or (len(hic1) == 1)) :
            logger.info("Generating pairs file ...")
            cmd = ["cphasing-rs", "porec2pairs", porec_table, f"../{contigsizes}",
                   "-o", str(target_pairs), "-q", "0"]
            flag = run_cmd(cmd, log=f"logs/porec2pairs.log")
            assert flag == 0, "Failed to execute command, please check log."
        
        if Path(target_pairs).exists():
            pairs = target_pairs
        
        os.chdir("../")

    hcr_wrkdir = "0.3.hcr"
    if hcr_flag and "0.3":
        steps.add("0.3")
    if collapsed_rescue or "0.3" in steps and "0.3" not in skip_steps:
        logger.info("""#----------------------------------#
#      Running step 0.3 hcr        #
#----------------------------------#""")
        Path(hcr_wrkdir).mkdir(exist_ok=True)
        os.chdir(hcr_wrkdir)
        Path("logs").mkdir(exist_ok=True)   
        with open("hcr.cmd.sh", 'w') as _out_sh:
            _out_sh.write(f"#!/bin/bash\n")
              
        os.chmod("hcr.cmd.sh", 0o755)

        if hcr_flag or collapsed_rescue or collapsed_contigs:
            run_hcr = hcr_flag or (not gfa)
            if run_hcr:
                hcr_invert_string = "-v" if hcr_invert else ""
                init_args = []
                current_pattern = pattern

                if paf and not corrected:
                    init_args = ["-paf", f"../{paf}"]
                    current_pattern = None 
                elif pairs and not corrected and Path(f"../2.prepare/input.{hcr_bs}.depth").exists():
                    init_args = ["-prs", f"../2.prepare/input.{hcr_bs}.depth"]
                else:
                    real_pairs = str(Path(f"../{pairs}").resolve())
                    run_cmd(["ln", "-sfn", real_pairs, pairs])
                    init_args = ["-prs", pairs]
                                
        
                hcr_args = init_args + [
                    "-cs", f"../{contigsizes}",
                    "-f", f"../{fasta}",
                    "-u", hcr_upper,
                    "-l", hcr_lower,
                    "-cr", collapsed_contig_ratio,
                    "-bs", hcr_bs
                ]

                
                if hcr_invert_string:
                    hcr_args.append(hcr_invert_string)
                if hcr_bed:
                    hcr_args.extend(["-b", f"{hcr_bed}"])

                hcr_args.extend(["-p", current_pattern])


                with open("hcr.cmd.sh", 'w') as _out_sh:
                    _out_sh.write("cphasing hcr \\")
                    _out_sh.write("\n")
                    _hcr_args = []
                    skip_next = False
                    for i in range(len(hcr_args)):
                        if skip_next:
                            skip_next = False
                            continue
                        
                        if i + 1 < len(hcr_args) and hcr_args[i+1] is None:
                            skip_next = True
                            continue
                        
                        if hcr_args[i] is not None:
                            _hcr_args.append(str(hcr_args[i]))
            
                    _out_sh.write(" ".join(pretty_cmd(map(str, _hcr_args))))
                    _out_sh.write("\n")

          
                try:
                    hcr.main(args=hcr_args, prog_name="hcr")

                except SystemExit as e:
                    exc_info = sys.exc_info()
                    exit_code = e.code
                    if exit_code is None:
                        exit_code = 0
                    
                    if exit_code != 0:
                        raise e
                    
                hcr_bed = str(Path(f"input.{hcr_bs}.hcr.bed").absolute())


            if gfa and not collapsed_contigs and not disable_gfa_collapsed:
                with open("hcr.cmd.sh", 'a') as _out_sh:
                    _out_sh.write("\n")
                    _out_sh.write("cphasing collapse from-gfa \\")
                    if corrected:
                        try:
                            _out_sh.write(f"--chimeric-bed ../0.2.correct/{break_bed} \\")
                        except UnboundLocalError:
                            pass
                    _out_sh.write(f"{gfa}\n")
                    _out_sh.write("\n")

                
                try:
                    collapsed_from_gfa.main(args=[
                         f"{gfa}" ],
                    prog_name="collapsed_from_gfa")
                except SystemExit as e:
                    exc_info = sys.exc_info()
                    exit_code = e.code
                    if exit_code is None:
                        exit_code = 0
                    
                    if exit_code != 0:
                        raise e
                

            if collapsed_rescue or collapsed_contigs:
                
                if collapsed_contigs is None:
                    collapsed_contigs = str(Path("contigs.collapsed.contig.list").absolute())
                    collapsed_dup_contigs = str(Path("contigs.collapsed.dup.contig.list").absolute())
                else:
                    collapsed_contigs = str(Path(collapsed_contigs).absolute())
                    collapsed_dup_contigs = str(Path("contigs.collapsed.dup.contig.list").absolute())
                    df = pd.read_csv(collapsed_contigs, header=None, sep="\t")
                    with open(collapsed_dup_contigs, 'w') as out:
                        for idx, row in df.iterrows():
                            contig = row[0]
                            cn = round(row[2])
                            if cn > 2:
                                for i in range(2, cn+1):
                                    out.write(f"{contig}\t{contig}_d{i}\n")


                if Path(collapsed_dup_contigs).exists() and Path(collapsed_dup_contigs).stat().st_size > 0:
                    if porec_table:
                        with open("hcr.cmd.sh", 'a') as _out_sh:
                            _out_sh.write(" ".join(map(str, ["cphasing-rs", "porec-dup", f"../{porec_table}", collapsed_dup_contigs, "-o", f"{porec_prefix}.collapsed.porec.gz"])))
                            _out_sh.write("\n")
                        new_porec = f"{porec_prefix}.dup.porec.gz"
                        is_porec_source_changed = is_file_changed(f"../{porec_table}")
                        is_dup_contigs_changed = is_file_changed(collapsed_dup_contigs)
                        if (not Path(new_porec).exists() or Path(new_porec).stat().st_size == 0 or
                            is_porec_source_changed or is_dup_contigs_changed):
                            cmd = ["cphasing-rs", "porec-dup", f"../{porec_table}", collapsed_dup_contigs, "-o", new_porec]
                            flag = run_cmd(cmd, log="logs/porec-dup.log")
                            assert flag == 0, "Failed to execute porec-dup, please check log."
                        else:
                            logger.warning(f"Using existing porec-dup output: {new_porec}")

                       
                        os.chdir("../")
                        safe_remove("input.porec.gz")
                        run_cmd(["ln", "-sfn", f"{hcr_wrkdir}/{new_porec}", "input.porec.gz"])
                        porec_table = "input.porec.gz"
                        hg_input = porec_table
                        os.chdir(hcr_wrkdir)
                    
            
                    with open("hcr.cmd.sh", 'a') as _out_sh:
                        _out_sh.write(" ".join(map(str, ["cphasing-rs", "pairs-dup", f"../{_pqs_latest}", collapsed_dup_contigs, "-o",  f"{pairs_prefix}.dup.pairs.pqs"])))
                        _out_sh.write("\n")

                    new_pairs = f"{pairs_prefix}.dup.pairs.pqs"
                    is_pairs_source_changed = is_file_changed(f"../{pairs}")
                    is_dup_contigs_changed_pairs = is_file_changed(collapsed_dup_contigs)
                    
                    if (not Path(new_pairs).exists() or Path(new_pairs).stat().st_size == 0 or
                        is_pairs_source_changed or is_dup_contigs_changed_pairs):
                        cmd = ["cphasing-rs", "pairs-dup", f"../{pairs}", collapsed_dup_contigs, "-o", new_pairs]
                        flag = run_cmd(cmd, log="logs/pairs-dup.log")
                        assert flag == 0, "Failed to execute pairs-dup, please check log."
                    else:
                        logger.warning(f"Using existing pairs-dup output: {new_pairs}")
                    
                   

                    os.chdir("../")
                    new_pairs = f"{pairs_prefix}.dup.pairs.pqs"
                    safe_remove("input.pairs.pqs")
                    run_cmd(["ln", "-sfn", f"{hcr_wrkdir}/{new_pairs}", "input.pairs.pqs"])
                    pairs = "input.pairs.pqs"
                    pqs_file = pairs
                    os.chdir(hcr_wrkdir)


        os.chdir("../")

    elif hcr_flag and "0.3" in skip_steps:
        hcr_bed_path = Path(f"{hcr_wrkdir}/input.{hcr_bs}.hcr.bed")
        if hcr_bed_path.exists():
            hcr_bed = str(hcr_bed_path.absolute())

    try:
        pqs_path = Path("input.pairs.pqs")
        
        if pqs_path.exists() and PQS.is_pqs(str(pqs_path)):
            logger.info(f"Using existing PQS file: {pqs_path}")
            is_pairs2pqs = True
            pqs_file = pairs = str(pqs_path)
            if not porec_table or use_pairs:
                hg_input = pqs_file
                hg_flag = "--pairs"
        else:
            if pairs and Path(pairs).exists():
                _pairs = Pairs2(pairs)
                if _pairs.is_pairs():
                    if not pqs_path.exists() or is_file_changed(str(pairs)):
                        logger.info("Converting pairs or pairs.gz to pairs.pqs ...")
                        args = [str(pairs), "-t", str(threads), "-o", str(pqs_path)]
                        try:
                            pairs2pqs.main(args=args, prog_name='pairs2pqs')
                        except SystemExit as e:
                            if e.code != 0:
                                raise e
                    
                    is_pairs2pqs = True
                    pqs_file = pairs = str(pqs_path)
                    
                    if not porec_table or use_pairs:
                        hg_input = pqs_file
                        hg_flag = "--pairs"
    except IsNotPairs:
        pass

    prepare_input = pairs or "input.pairs.pqs"

    alleles_dir = str("1.alleles")
    
    genomesize = read_chrom_sizes(contigsizes)['length'].sum()
    logger.info(f"Total size of contig-level assembly: {to_humanized3(genomesize)}")
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
                threads,
                "-tl",
                alleles_trim_length]
        
        with open('alleles.cmd.sh', 'w') as _out_sh:
            _out_sh.write("#!/bin/bash\n")
            _out_sh.write("cphasing alleles ")
            _out_sh.write(" ".join(map(str, args)))
            _out_sh.write("\n")
        os.chmod('alleles.cmd.sh', 0o755)
        
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

    allele_table_basename = f"{fasta_prefix}.allele.table"
    allele_table = None if mode in {"basal", "hapaware"} else f"{alleles_dir}/{allele_table_basename}" 

    prepare_prefix = Path(Path(prepare_input).name)
    while prepare_prefix.suffix in {".pairs", ".gz", ".pqs"}:
        prepare_prefix = prepare_prefix.with_suffix("")

    prepare_dir = str("2.prepare")
    
    if "2" not in skip_steps and "2" in steps:
        logger.info("""
#----------------------------------#
#       Running step 2. prepare    #
#----------------------------------#""")
        Path(prepare_dir).mkdir(exist_ok=True)
        os.chdir(prepare_dir)

        args = [f"../{fasta}",
                f"../{prepare_input}", 
                "-p",
                pattern, 
                "-bs",
                hcr_bs,
                "-q",
                mapping_quality,
                "-t",
                threads,
                "--skip-pairs2contacts"]
        
        with open("prepare.cmd.sh", 'w') as _out_sh:
            _out_sh.write("cphasing prepare ")
            _out_sh.write(" ".join(map(str, args)))
            _out_sh.write("\n")
        
        os.chmod("prepare.cmd.sh", 0o755)
                
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
    
    if pattern is None:
        pattern = DEFAULT_PATTERN
    
    pattern = pattern.upper().replace("^", "").replace(",", "_")
    count_re = f"{prepare_dir}/{prepare_prefix}.counts_{pattern}.txt"
    clm = f"{prepare_dir}/{prepare_prefix}.clm.gz"
    
    
    split_contacts = f"{prepare_dir}/{prepare_prefix}.split.contacts.gz"

    output_cluster = "final.clusters.txt"
    
    hyperpartition_normalize = "-norm" if normalize else None
    enable_misassembly_remove = "--enable-misassembly-remove" if enable_misassembly_remove else None 

        
    if use_pairs or not porec_table:
        hg_input = pairs
        input_param = "--pairs"
        hg_flag = "--pairs"
    else:
        hg_input = porec_table
        input_param = "--porec"
        hg_flag = ""

    hyperpartition_contacts = None
    

    hyperpartition_dir = str("3.hyperpartition")
    collapsed_rescued_contigs = None
    if "3" not in skip_steps and "3" in steps:
        logger.info("""
#----------------------------------#
#  Running step 3. hyperpartition  #
#----------------------------------#""")
        Path(hyperpartition_dir).mkdir(exist_ok=True)
        os.chdir(hyperpartition_dir)
        _output_cluster = "output.clusters.txt"
        _mode = mode
        _mode = "phasing" if _mode in {"hapaware", "phasing2"} else _mode

        hyperpartition_args = [
                                f"../{hg_input}",
                                f"../{contigsizes}",
                                _output_cluster,
                                input_param,
                                "--mode",
                                _mode,
                                "-e",
                                edge_length,
                                "-sl",
                                split_length,
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
                                "-mcw",
                                min_cis_weight,
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
        if kprune_norm_method != "auto":
            hyperpartition_args.extend(["-knm", kprune_norm_method])
            if kprune_norm_method == "re":
                logger.info(f"Using `re` kprune norm method.")
                hyperpartition_args.extend(["-cr", f"../{count_re}"])

        if allelic_factor != -1:
            hyperpartition_args.extend(["-af", allelic_factor])
        
        if disable_merge_in_first:
            hyperpartition_args.append("--disable-merge-in-first")
        
        if disable_merge_use_allele:
            hyperpartition_args.append("--disable-merge-use-allele")
        
        if allelic_positive_factor:
            hyperpartition_args.extend(["-apf", str(allelic_positive_factor)])
            
        if enable_misassembly_remove:
            hyperpartition_args.append(enable_misassembly_remove)
        if disable_recluster_by_linkage:
            hyperpartition_args.append("--disable-recluster-by-linkage")
        if refine:
            hyperpartition_args.append("--refine")
        if hyperpartition_normalize:
            hyperpartition_args.append(hyperpartition_normalize)

        if blacklist:
            hyperpartition_args.extend(["--blacklist", f"{blacklist}"])

        if hcr_bed:
            hyperpartition_args.extend(["--hcr-bed", f"{hcr_bed}"])
            if hcr_invert:
                hyperpartition_args.append("--hcr-invert")

        if mode == "phasing":
            hyperpartition_args.extend(["-f", fasta_path])
            hyperpartition_args.extend(["--alleles-k", alleles_kmer_size, 
                                        "--alleles-w", alleles_window_size,
                                        "--alleles-minimum-similarity", alleles_minimum_similarity,
                                        "--alleles-diff-thres", alleles_diff_thres,
                                        "--alleles-trim-length", alleles_trim_length])

        else:
            hyperpartition_args.extend(["-at",
                         f"../{allele_table}" if allele_table else None,])
        
        if cluster_method != CLUSTER_METHOD:
            hyperpartition_args.extend(["-cm", cluster_method])

        if output_hg:
            hyperpartition_args.extend(["--output-hg"])
            
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

        os.chmod("hyperpartition.cmd.sh", 0o755)

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


        if collapsed_contigs:
            _fp = Path(fasta_path).name
            while Path(_fp).suffix in {".fasta", ".gz", ".fa"}:
                _fp = Path(_fp).stem
            if mode == "phasing":
                actual_allele_table = f"{_fp}.allele.table" 
            elif mode == "base_withprune":
                actual_allele_table = f"{_fp}.allele.table" if allele_table else None
            else:
                actual_allele_table = None
            
            collapsed_rescue_args = [
                f"../{hg_input}",
                f"../{contigsizes}",
                _output_cluster,
                collapsed_contigs,
                "-n", str(_n[1] if _n[1] else 4) if _n and len(_n) > 1 else "2",
                "-at", actual_allele_table,
                "-mw", str(min_weight),
                "-q",  str(min_quality2),
                "-mcw", str(min_cis_weight),
                "-mc", str(min_contacts),
                "-as", str(allelic_similarity),
                "-t", str(threads),
                input_param,
            ]

            if disable_conflict_check:
                collapsed_rescue_args.append("--disable-conflict-check")

            if cn_offset:
                collapsed_rescue_args.extend(["--cn-offset", str(cn_offset)])
            logger.info("Rescuing collapsed contigs ...")

            with open("hyperpartition.cmd.sh", 'a') as _out_sh:
                _out_sh.write("\n")
                _out_sh.write("cphasing collapse rescue \\\n    ")
                _collapsed_rescue_args = []
                remove_index = []
                for i in range(4, len(collapsed_rescue_args), 2):
                    try:
                        if collapsed_rescue_args[i+1] is None:
                            remove_index.append(i)
                            remove_index.append(i+1)
                    except IndexError:
                        continue
                for i, j in enumerate(collapsed_rescue_args):
                    if i not in remove_index:
                        _collapsed_rescue_args.append(j)
                
                _out_sh.write(" ".join(pretty_cmd(map(str, _collapsed_rescue_args))))
                _out_sh.write("\n")

            try:
                collapsed_rescue_cli.main(args=collapsed_rescue_args,
                                      prog_name="rescue")
            except SystemExit as e:
                exc_info = sys.exc_info()
                exit_code = e.code
                if exit_code is None:
                    exit_code = 0
                
                if exit_code != 0:
                    raise e
            
            _output_cluster = "collapsed.rescue.clusters.txt"

            collapsed_rescued_contigs = str(Path("collapsed.rescue.contigs.list").absolute())


        if Path(_output_cluster).exists():
            logger.info(f"Renaming `{_output_cluster}` to `final.clusters.txt` ...")
            shutil.copy(_output_cluster, "final.clusters.txt")
            with open("hyperpartition.cmd.sh", 'a') as _out_sh:
                _out_sh.write("\n")
                _out_sh.write(f"cp {_output_cluster} final.clusters.txt\n")
            output_cluster = "final.clusters.txt"

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
        if mode == "phasing":
            _fp = Path(fasta_path).name
            while Path(_fp).suffix in {".fasta", ".gz", ".fa"}:
                _fp = Path(_fp).stem
            allele_table = f"{hyperpartition_dir}/{_fp}.allele.table"
        if corrected:
            args = [
                f"../{hyperpartition_dir}/{output_cluster}",
                f"../{count_re}",
                "--clm",
                f"../{clm}",
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
                "--clm",
                f"../{clm}",
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

        if allele_table:
            args.extend(["-at",
                        f"../{allele_table}"])
            
        if enable_haplotype_cluster:
            args.append("--enable-haplotype-cluster")
            
        with open("scaffolding.cmd.sh", "w") as _out_sh:
            _out_sh.write("cphasing scaffolding \\\n    ")
            _scaffolding_args = []
            remove_index = []
            for i in range(4, len(args), 2):
                try:
                    if args[i+1] is None:
                        remove_index.append(i)
                        remove_index.append(i+1)
                except IndexError:
                    continue
            for i, j in enumerate(args):
                if i not in remove_index:
                    _scaffolding_args.append(j)
            _out_sh.write(" ".join(pretty_cmd(map(str, _scaffolding_args))))
            _out_sh.write("\n")
        os.chmod("scaffolding.cmd.sh", 0o755)

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

        if gfa:
            from ..gfa_cli import reassemble
            input_agp = corrected_agp if corrected else out_agp
            logger.info("Reassembling scaffolds based on GFA ...")

            input_args = [input_agp, gfa, "-o", input_agp.replace(".agp", ".reassemble.agp")]
            if corrected:
                try:
                    input_args.extend(["--chimeric-bed", f"../0.2.correct/{break_bed}"])
                except UnboundLocalError:
                    pass

            with open("scaffolding.cmd.sh", 'a') as _out_sh:
                _out_sh.write("\n")
                _out_sh.write("cphasing gfa reassemble \\\n    ")

                _out_sh.write(" ".join(pretty_cmd(map(str, input_args), n=2)))
                _out_sh.write("\n")
            try:
                reassemble.main(args=input_args,
                                prog_name='scaffolding')
            except SystemExit as e:
                exc_info = sys.exc_info()
                exit_code = e.code
                if exit_code is None:
                    exit_code = 0
                
                if exit_code != 0:
                    raise e
            
            out_agp = input_agp.replace(".agp", ".reassemble.agp")


            with open("scaffolding.cmd.sh", 'a') as _out_sh:
                _out_sh.write("\n")
                _out_sh.write(f"cphasing agp2fasta {out_agp} ../{fasta} -o {out_agp.replace('.agp', '.asm.fasta')}\n")
            try:
                agp2fasta.main(args=[out_agp, f"../{fasta}", 
                                     "-o", f"{out_agp.replace('.agp', '.asm.fasta')}"],
                                prog_name='scaffolding')
            except SystemExit as e:
                exc_info = sys.exc_info()
                exit_code = e.code
                if exit_code is None:
                    exit_code = 0
                
                if exit_code != 0:
                    raise e

        if collapsed_contigs and not is_empty(collapsed_rescued_contigs):
            with open("scaffolding.cmd.sh", 'a') as _out_sh:
                _out_sh.write("\n")
                _out_sh.write("cphasing collapse agp-dup \\\n    ")
                _out_sh.write(" ".join(pretty_cmd(map(str, [
                    out_agp,
                    "-o",
                    out_agp.replace(".agp", ".dup.agp"),
                ]), n=2)))
                _out_sh.write("\n")
            try:
                agp_dup.main(args=[out_agp,
                                   "-o",
                                   out_agp.replace(".agp", ".dup.agp"),],
                                prog_name='scaffolding')
            except SystemExit as e:
                exc_info = sys.exc_info()
                exit_code = e.code
                if exit_code is None:
                    exit_code = 0
                
                if exit_code != 0:
                    raise e
                pass

            logger.info("Renamed rescued collapsed contigs to AGP ...")


        generate_to_hic_cmd(out_agp, f"../{fasta}", f"../{pairs}", n=n )
      
        
        os.chdir("..")

    if binsize == "auto" or binsize is None:
        init_binsize, binsize = recommend_binsize_by_genomesize(genomesize)
        
        logger.info(f"Recommended cool's binsize: `{to_humanized2(init_binsize)}`, heatmap's binsizs: `{to_humanized2(binsize)}`")
    else:
        init_binsize = 10000
        binsize = '500k'

    
    
    if corrected:
        if Path(f"{scaffolding_dir}/{corrected_agp}").exists():
            input_agp = corrected_agp
        else:
            input_agp = out_agp
    else:
        input_agp = out_agp

    if collapsed_contigs and Path(f"{scaffolding_dir}/{out_agp.replace('.agp', '.dup.agp')}").exists():
        input_agp = out_agp.replace('.agp', '.dup.agp')

    plot_mapq = min_quality2
    plot_dir = str("5.plot")
    if "5" not in skip_steps and "5" in steps:
        logger.info("""
#----------------------------------#
#      Running step 5. plot        #
#----------------------------------#""")
        Path(plot_dir).mkdir(exist_ok=True)

        os.chdir(plot_dir)
   

        out_small_cool = f"{pairs_prefix}.q{plot_mapq}.{to_humanized2(init_binsize)}.cool"

        if not is_file_changed(f"../{pairs}") and Path(out_small_cool).exists():
            logger.warning(f"`{out_small_cool}` exists, skipped `pairs2cool`.")
            args = [
                    f"../{pairs}",
                    f"../{contigsizes}",
                    out_small_cool,
                    "-q", 
                    plot_mapq,
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
                    plot_mapq,
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

        generate_plot_cmd(f"../{pairs}",
                        f"{pairs_prefix}",
                        f"../{contigsizes}",
                        f"../{scaffolding_dir}/{input_agp}", 
                        plot_mapq, init_binsize,
                        binsize, colormap, 
                        whitered, mode,)
        args = [
                "-a",
                f"../{scaffolding_dir}/{input_agp}",
                "-m",
                out_small_cool,
                "-o",
                f"groups.q{plot_mapq}.{to_humanized2(binsize)}.wg.png",
                "-o",
                f"groups.q{plot_mapq}.{to_humanized2(binsize)}.wg.pdf",
                "-bs",
                binsize,
                "-oc",
                ]
        
        if colormap != "redp1_r":
            args.extend(["--colormap", colormap])
        if whitered:
            args.append("--whitered")
        if balance:
            args.append("--balance")

        if mode == "phasing":
            args.extend(["--add-hap-border", "--no-lines"])

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

    if "6" not in skip_steps and "6" in steps:
        logger.info("""
#----------------------------------#
#    Running step 6. curation      #
#        (only generated cmd)      #
#----------------------------------#""")
        curation_dir = str("6.curation")
        Path(curation_dir).mkdir(exist_ok=True)
        os.chdir(curation_dir)
        
        cmd_output = generate_curation_cmd(
                        f"{input_agp}",
                        f"../{fasta}",
                        f"../{pairs}",
                        binsize,
                        init_binsize,
                        plot_mapq,
                        scaffolding_dir,
                        plot_dir,
                        n,
                    )
        
        logger.info("Output curation command script.")

        os.chdir("..")

    today = date.today().strftime("%Y-%m-%d")
    end_time = time.time() - start_time
    monitor.join()
    if monitor.memory_buffer:
        peak_memory = max(monitor.memory_buffer) / 1024 / 1024 / 1024
    else:
        peak_memory = 0.0

    columns = shutil.get_terminal_size((80, 20)).columns
    logger.info("\n")
    today = datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    logger.info(f"The pipeline finished in {today}")
    logger.info(f"    Elapsed time: {end_time:.2f} s")
    logger.info(f"    Peak memory: {peak_memory:.2f} Gb")
   
    logger.info("=" * (columns // 2))
    logger.info(f"Results are store in `{outdir}/`")
    
    if ("4" in steps or "5" in steps) and ("4" not in skip_steps or "5" not in skip_steps):
        logger.info("Please check these results:")
    if "4" in steps and "4" not in skip_steps:
        logger.info(f"    {outdir}/4.scaffolding/{out_agp}")
    if "5" in steps and "5" not in skip_steps:
        if corrected:
            logger.info(f"    {outdir}/5.plot/groups.q{plot_mapq}.{to_humanized2(binsize)}.wg.png")
        else:
            logger.info(f"    {outdir}/5.plot/groups.q{plot_mapq}.{to_humanized2(binsize)}.wg.png")


    if "6" in steps and "6" not in skip_steps:
        logger.info("\n")
        logger.info("-" * (columns // 2))
        logger.info("If you want to curate the assembly, please check:")
        logger.info(f"    {outdir}/6.curation/curation.cmd.sh")
    