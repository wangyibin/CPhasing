#!/usr/bin/env python



import logging
import sys 
import os
import os.path as op
import psutil 
import re

import numpy as np
import pandas as pd

import rich_click as click
from rich_click import RichCommand
from click_didyoumean import DYMGroup

from collections import defaultdict
from pathlib import Path
from pytools import natsorted
from shutil import which

from . import __version__, __epilog__
from .core import (
    ClusterTable, 
    CountRE, 
)
from .utilities import (
    run_cmd,
    humanized2numeric,
    to_humanized2,
    is_compressed_table_empty,
    )

logger = logging.getLogger("cphasing")

banner = """
   ____      ____  _               _             
  / ___|    |  _ \| |__   __ _ ___(_)_ __   __ _ 
 | |   _____| |_) | '_ \ / _` / __| | '_ \ / _` |
 | |__|_____|  __/| | | | (_| \__ \ | | | | (_| |
  \____|    |_|   |_| |_|\__,_|___/_|_| |_|\__, |
                                           |___/ 
"""

HIDDEN = True

# help_config = click.RichHelpConfiguration(
#     show_arguments=True,
    
# )
click.rich_click.USE_MARKDOWN = True
click.rich_click.STYLE_COMMANDS_TABLE_SHOW_LINES = False
click.rich_click.STYLE_COMMANDS_TABLE_PAD_EDGE = True
click.rich_click.STYLE_COMMANDS_TABLE_BOX = "SIMPLE"
click.rich_click.STYLE_COMMANDS_TABLE_BORDER_STYLE = "red"
click.rich_click.STYLE_USAGE_COMMAND = "bold red"
# click.rich_click.STYLE_COMMANDS_TABLE_ROW_STYLES = ["yellow", "green", "cyan"]
click.rich_click.MAX_WIDTH = 128

click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.STYLE_ERRORS_SUGGESTION = "magenta italic"
click.rich_click.ERRORS_SUGGESTION = "Try running the '--help' flag for more information."
click.rich_click.ERRORS_EPILOGUE = f"Version: {__version__} | To find out more, visit [https://wangyibin.github.io/CPhasing](https://wangyibin.github.io/CPhasing)"

click.rich_click.COMMAND_GROUPS = {
    "cphasing": [
        {
            "name": "Pipe (Recommended)",
            "commands": ["pipeline"],
            "table_styles": {
                "row_styles": ["yellow"],
            }
        },
        {
            "name": "Main Commands",
            "commands": ["hitig", "mapper", "chimeric", "prepare",
                         "alleles", "hyperpartition", "scaffolding",
                         "build", "rename", "pairs2cool", "plot"]
        },
        {
            "name": "Other Useful Commands",
            "commands": ["hic", "alignments", "hcr", "alleles2", 
                         "kprune", "hypergraph", "statagp", 
                         "agp2fasta", "utils", "evaluator"]
        }
    ],

    "cphasing hitig": [
        {
           "name": "Main Commands",
            "commands": ["pipeline"],
            "table_styles": {
                "row_styles": ["yellow"],
            } 
        }
    ],
    "hitig": [
        {
           "name": "Main Commands",
            "commands": ["pipeline"],
            "table_styles": {
                "row_styles": ["yellow"],
            } 
        }
    ]
}


click.rich_click.OPTION_GROUPS = {
    "cphasing pipeline": [
        {
            "name": "Options of Input Files",
            "options": ["--fasta", "--porec-data", "--porectable", 
                            "--pairs", "--hic1", "--hic2"]
        },
        {
            "name": "Advance Options of Pipeline",
            "options": ["--mode", "--steps", "--skip-steps"]
        },
        {
            "name": "Options of Hitig",
            "options": ["--ul-data", "--use-existed-hitig"]
        },
        {
            "name": "Options of Pore-C Mapper",
            "options": ["--mapper-k", "--mapper-w", "--mapping-quality"]
        },
         {
            "name": "Options of Hi-C Mapper",
            "options": ["--hic-mapper-k", "--hic-mapper-w", "--mapping-quality"]
        },
        {   "name": "Options of chimeric correction",
            "options": ["--chimeric-correct", "--chimeric-corrected"]

        },
        {
            "name": "Options of HCR",
            "options": ["--hcr", "--pattern", "--hcr-lower", "--hcr-upper",
                         "--collapsed-contigs-ratio", "--hcr-binsize",
                         "--hcr-bed", "--hcr-invert"]
        },
        {
            "name": "Options of Alleles",
            "options": ["--alleles-k", "--alleles-w", "--alleles-m", "--alleles-d"]
        },
    {
            "name": "Options of HyperPartition",
            "options": ["-n", "--use-pairs", "--min-contacts",
                        "--Nx", "--min-length", "--min-weight",
                        "--min-cis-weight", "--whitelist",
                        "--resolution1", "--resolution2",
                         "--init-resolution1",  "--init-resolution2",
                          "--first-cluster", "--normalize", 
                          "--min-quality1", "--min-quality2",
                           "--min-scaffold-length",
                          "--allelic-similarity", "--min-allelic-overlap",
                          "--disable-merge-in-first", "--exclude-group-from-first",
                          "--exclude-group-to-second", "--enable-misassembly-remove"
                          ]
        },
        {
            "name": "Options of Scaffolding",
            "options": ["--scaffolding-method"]
        },
        {
            "name": "Options of Plot",
            "options": ["--binsize"]
        },
        {
            "name": "Global Options",
            "options": ["--low-memory", "--threads", "--help"]
        }
       
    ],
    "cphasing pipe": [
        {
            "name": "Options of Input Files",
            "options": ["--fasta", "--porec-data", "--porectable", 
                            "--pairs", "--hic1", "--hic2"]
        },
        {
            "name": "Advance Options of Pipeline",
            "options": ["--mode", "--steps", "--skip-steps"]
        },     
        {
            "name": "Options of Hitig",
            "options": ["--ul-data", "--use-existed-hitig"]
        },
        {
            "name": "Options of Pore-C Mapper",
            "options": ["--mapper-k", "--mapper-w", "--mapping-quality"]
        },
         {
            "name": "Options of Hi-C Mapper",
            "options": ["--hic-mapper-k", "--hic-mapper-w", "--mapping-quality"]
        },
        {   "name": "Options of chimeric correction",
            "options": ["--chimeric-correct", "--chimeric-corrected"]

        },
        {
            "name": "Options of HCR",
            "options": ["--hcr", "--pattern", "--hcr-lower", "--hcr-upper", 
                        "--collapsed-contigs-ratio", "--hcr-binsize",
                         "--hcr-bed", "--hcr-invert"]
        },
        {
            "name": "Options of Alleles",
            "options": ["--alleles-k", "--alleles-w", "--alleles-m", "--alleles-d"]
        },
        {
            "name": "Options of HyperPartition",
            "options": ["-n", "--use-pairs", "--min-contacts", 
                        "--Nx", "--min-length", "--min-weight",
                        "--min-cis-weight", "--whitelist",
                        "--resolution1", "--resolution2",
                         "--init-resolution1",  "--init-resolution2",
                          "--first-cluster", "--normalize", 
                          "--min-quality1", "--min-quality2",
                           "--min-scaffold-length",
                          "--allelic-similarity", "--min-allelic-overlap",
                          "--disable-merge-in-first", "--exclude-group-from-first",
                          "--exclude-group-to-second", "--enable-misassembly-remove"
                          ]
        },
        {
            "name": "Options of Scaffolding",
            "options": ["--scaffolding-method"]
        },
        {
            "name": "Options of Plot",
            "options": ["--binsize"]
        },
        {
            "name": "Global Options",
            "options": ["--low-memory", "--threads", "--help"]
        }
       
    ],
    "cphasing plot": [
        {
            "name": "Options of Matrix Operation",
            "options": ["--matrix", "--binsize", 
                        "--only-coarsen",
                        "--balance", "--balanced"]
        },
        {
            "name": "Options of AGP Adjustment",
            "options": ["--agp", "--only-adjust"]
        },
        {
            "name": "Options of Heatmap",
            "options": [
                        "--chromosomes", 
                        "--disable-natural-sort",
                        "--per-chromosomes",
                        "--only-chr", "--chr-prefix",
                        "--chrom-per-row",
                        "--vmin", "--vmax", 
                        "--scale", 
                        "--triangle",
                        "--fontsize",
                        "--dpi", "--cmap", "--whitered",
                        "--no-lines", "--no-ticks", 
                        "--rotate-xticks", "--rotate-yticks"]
            
        },
        {
            "name": "Global Options",
            "options": ["--output", "--threads", "--help"]
        }
    ]
}

class CommandGroup(click.RichGroup, RichCommand):
    """
    List subcommand in the order there were added.
    """
    def list_commands(self, ctx):
        return list(self.commands)

    def get_command(self, ctx, cmd_name):
        """
        Alias command, https://stackoverflow.com/questions/46641928/python-click-multiple-command-names
        """
        try:
            cmd_name = ALIASES[cmd_name].name 
        except KeyError:
            pass 
        return super().get_command(ctx, cmd_name)



@click.version_option(__version__, "-V", "--version")
@click.group(context_settings={"help_option_names": ["-h", "--help", "-help"]},
             cls=CommandGroup,
             epilog=__epilog__)
@click.option(
    "-v", 
    "--verbose", 
    count=True, 
    help="Increase level of logging information, eg. -vvv"
)
@click.option(
    "--quiet", 
    is_flag=True, 
    default=False, 
    help="Turn off all logging.", 
    show_default=True
)
def cli(verbose, quiet):
    """
    **Phasing** and scaffolding polyploid genomes based on Pore-**C**, Ultra-long, or Hi-**C** data.    
    --------------------

    \f

    :param click.core.Context ctx: Click context.
    """
    if quiet:
        logger.setLevel(logging.CRITICAL)
    elif verbose > 0:
        LOG_LEVELS = [logging.CRITICAL, logging.ERROR, 
                        logging.WARNING, logging.INFO, 
                        logging.DEBUG]
        offset = 2 
        idx = min(len(LOG_LEVELS) - 1, offset + verbose)
        logger.setLevel(LOG_LEVELS[idx])
    else:
        logger.setLevel(logging.INFO)



@cli.command(cls=RichCommand, epilog=__epilog__)
@click.option(
    '-f',
    '--fasta',
    metavar="FASTA",
    help="Path to draft assembly",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '-ul',
    '--ul-data',
    '-ul-data',
    metavar="Ul Data",
    help="Provide UL or HiFi data to run hitig",
    default=None,
    show_default=True
)
@click.option(
    '-pcd',
    '--porec-data',
    'porec_data',
    metavar="PoreC Data",
    multiple=True,
    help="Pore-C data. Multiple files use multiple options, such as -pcd read1.fq.gz -pcd read2.fq.gz. "
    "However, the prefix of subsequence output files use the read1. "
    "If you want to run parallel of multi porec data, you can run `cphasing mapper` for each cell. "
    "And then use `cphasing-rs porec-merge` to merge porec table into a single file",
    type=click.Path(exists=True),
    default=None,
    show_default=True
)
@click.option(
    '-pct',
    '--porectable',
    metavar="PoreC-Table",
    help="PoreC Table file, which generated from `cphasing mapper`",
    type=click.Path(exists=True),
    default=None,
    show_default=True
)
@click.option(
    '-prs',
    '--pairs',
    metavar="Pairs",
    help="4DN pairs file, which generated from `cphasing mapper`, or `cphasing hic mapper`",
    type=click.Path(exists=True),
    default=None,
    show_default=True
)
@click.option(
    '-hic1',
    '--hic1',
    metavar="R1 Reads",
    help="Input hic read1 to run the pipeline by hic data, only support one file",
    type=click.Path(exists=True),
    default=None,
    show_default=True
)
@click.option(
    '-hic2',
    '--hic2',
    metavar="R2 Reads",
    help="Input hic read2 to run the pipeline by hic data, only support one file",
    type=click.Path(exists=True),
    default=None,
    show_default=True
)
@click.option(
    '-p',
    '--pattern',
    help="Pattern of restriction enzyme. Comma separated for multiple pattern."
    " Restriction enzyme used to normalization the hcr bins, if not use the `-hcr`, do not need specified.",
    metavar="STR",
    default=None,
    show_default=True,
    # hidden=True,
)
@click.option(
    '-mapper-k',
    "--mapper-k",
    "mapper_k",
    metavar="INT",
    help="kmer size for mapper",
    default=15,
    show_default=True
)
@click.option(
    '-mapper-w',
    "--mapper-w",
    "mapper_w",
    metavar="INT",
    help="window size for mapper",
    default=10,
    show_default=True
)
@click.option(
    '-q',
    '--mapping-quality',
    'mapping_quality',
    metavar='INT',
    help='minimum mapping quality of mapper or hicmapper, also control the minimum mapping quality for scaffolding.',
    default=0,
    show_default=True,
    type=click.IntRange(0, 60),
)
@click.option(
    '-hic-mapper-k',
    "--hic-mapper-k",
    "hic_mapper_k",
    metavar="INT",
    help="kmer size for mapper",
    default=17,
    show_default=True
)
@click.option(
    '-hic-mapper-w',
    "--hic-mapper-w",
    "hic_mapper_w",
    metavar="INT",
    help="window size for mapper",
    default=7,
    show_default=True
)
@click.option(
    "--chimeric-correct",
    is_flag=True,
    default=False,
    help="Correct the chimeric contigs by contacts.",
    show_default=True
)
@click.option(
    "--chimeric-corrected",
    is_flag=True,
    default=False,
    help="Use existed results of correction, overwrite of the `--chimeric-corrected`",
    show_default=True
)
@click.option(
    "-hcr",
    "--hcr",
    is_flag=True,
    default=False,
    help="Only retain high confidence regions to subsequence analysis by contacts.",
    show_default=True
)
@click.option(
    "-hcr-l",
    "--hcr-lower",
    "hcr_lower",
    help="Lower value of peak value.",
    type=float,
    default=.01,
    show_default=True,
)
@click.option(
    "-hcr-u",
    "--hcr-upper",
    "hcr_upper",
    help="Upper value of peak value.",
    type=float,
    default=1.5,
    show_default=True,
)
@click.option(
    "-cr",
    "--collapsed-contigs-ratio",
    "collapsed_contig_ratio",
    metavar="FLOAT",
    help="minimum regions of the collapsed contig, "
    " which mean the region of a contig higher than upper "
    "value will be regarded as collapsed and remove the whole contigs",
    type=click.FloatRange(0.0, 1.0),
    default=.0,
    show_default=True,
)
@click.option(
    "-hcr-bs",
    "--hcr-binsize",
    "hcr_bs",
    metavar="INT",
    help="Bin size for HCRs identification, only support for hcr by contacts",
    default=10000,
    type=int,
    show_default=True,
)
@click.option(
    "-hcr-bed",
    "--hcr-bed",
    "hcr_bed",
    metavar="STR",
    default=None,
    help="Only retain high confidence regions to subsequence analysis by input bed.",
    show_default=True,
    type=click.Path(exists=True)
)
@click.option(
    "-hcr-invert",
    "--hcr-invert",
    "hcr_invert",
    is_flag=True,
    default=False,
    help="Invert table by hcr bed.",
    show_default=True,
)
@click.option(
    '-m',
    '--mode',
    metavar="STR",
    help="mode of hyperpartition, the basal equal to haploid. "
    "`['basal', 'haploid', 'phasing', 'basal_withprune']`",
    default='phasing',
    show_default=True,
    type=click.Choice(['basal', 'haploid', 'phasing', 'phasing2', 'basal_withprune']),
)
@click.option(
    '-s',
    '--steps',
    metavar='STR',
    help="steps, comma seperate",
    default="1,2,3,4,5",
    show_default=True
)
@click.option(
    '-ss',
    '--skip-steps',
    'skip_steps',
    metavar="STR",
    help="skip following steps, comma seperate.",
    default=None,
    show_default=True
)
@click.option(
    '--use-existed-hitig',
    help="Use existed hitig results.",
    is_flag=True,
    default=False,
    show_default=True 
)
@click.option(
    '-alleles-k',
    '--alleles-k',
    '--alleles-kmer-size',
    'alleles_kmer_size',
    metavar="INT",
    help="kmer size for `alleles` similarity calculation.",
    default=19,
    show_default=True,
)
@click.option(
    '-alleles-w',
    '--alleles-w',
    '--alleles-window-size',
    'alleles_window_size',
    metavar="INT",
    help="minimizer window size for `alleles` similarity calculation.",
    default=19,
    show_default=True,
)
@click.option(
    '-alleles-m',
    '--alleles-m',
    '--alleles-minimum-similarity',
    'alleles_minimum_similarity',
    metavar="FLOAT",
    help="minimum k-mer similarity for `alleles` similarity calculation.",
    default=0.5,
    show_default=True,
)
@click.option(
    "-alleles-d",
    '--alleles-d',
    '--alleles-diff-thres',
    "alleles_diff_thres",
    help="minimum different threshold between contig pairs.",
    type=click.FloatRange(0, 1),
    default=.1,
    show_default=True
)
@click.option(
    "--use-pairs",
    "use_pairs",
    help="""
    Use pairs instead of porec table to construct the hypergraph/graph.
    """,
    default=False,
    is_flag=True
)
@click.option(
    "-n",
    help="""
    Number of groups. If set to 0 or None, the partition will be run automatically.
    If in phasing mode, you can set to n1:n2, 
    which meaning that generate n1 groups in the first round partition
    and generate n2 groups in the second round partition. 
    Also, you can set `-n 8:0` to set the second round partition to automatical mode.
     The n2 parameter also can be set a file, which containing two columns, 
      the firse column repersent the group index from `first.clusters.txt`, and the 
       second column repersent the n2 groups in each n1 group. [default: 0]
    """,
    default=None,
    show_default=True,
)
@click.option(
    "-r1",
    "--resolution1",
    metavar="FLOAT",
    help="Resolution of the first partition",
    type=click.FloatRange(-1.0, 10.0),
    default=-1.0,
    show_default=True
)
@click.option(
    "-r2",
    "--resolution2",
    metavar="FLOAT",
    help="Resolution of the second partition",
    type=click.FloatRange(-1.0, 10.0),
    default=1.0,
    show_default=True
)
@click.option(
    "-ir1",
    "--init-resolution1",
    "init_resolution1",
    metavar="FLOAT",
    help="Initial resolution of the automatic search resolution in first partition"
    ", only used when `--resolution1=-1`.",
    type=float,
    default=0.8,
    show_default=True
)
@click.option(
    "-ir2",
    "--init-resolution2",
    "init_resolution2",
    metavar="FLOAT",
    type=float,
    help="Initial resolution of the automatic search resolution in second partition"
    ", only used when `--resolution2=-1`.",
    default=0.8,
    show_default=True
)
@click.option(
    '-fc',
    '--first-cluster',
    'first_cluster',
    metavar="PATH",
    help='Use the existing first cluster results to second round cluster',
    default=None,
    show_default=True,
    type=click.Path(exists=True)
)
@click.option(
    '--disable-merge-in-first',
    '--disable-merge-first',
    help="""
    disable merge function in first round cluster results, if set the `n` that
     program select `n` group to second round cluster.
    """,
    default=False, 
    is_flag=True
)
@click.option(
    '--exclude-from-first',
    '--exclude-group-from-first',
    'exclude_group_from_first',
    metavar="STR",
    help='exclude several groups in first rount cluster, do not retain it to any group, comma seperate. 1-base.'
    '[default: None]',
    default=None,
    show_default=True
)
@click.option(
    '--exclude-to-second',
    '--exclude-group-to-second',
    'exclude_group_to_second',
    metavar="STR",
    help='exclude several groups that do not run in second round cluster, comma seperate. 1-base.'
    '[default: None]',
    default=None,
    show_default=True
)
@click.option(
    "-norm",
    "--normalize",
    help="Normalize expanded edges by contig length",
    is_flag=True,
    show_default=True,
    default=False,
)
@click.option(
    "-as",
    "--allelic-similarity",
    "allelic_similarity",
    metavar="FLOAT",
    help="The similarity of allelic, which used to merge clusters",
    default=0.85,
    type=click.FloatRange(0.0, 1.0),
    show_default=True
)
@click.option(
    '-mao',
    '--min-allelic-overlap',
    "min_allelic_overlap",
    metavar="FLOAT",
    help="Minimum overlap ratio bewteen two group when merging different groups",
    default=0.3,
    type=click.FloatRange(0.0, 1.0),
    show_default=True
)
@click.option(
    "-mw",
    "--min-weight",
    "min_weight",
    help="Minimum weight of graph",
    type=float,
    default=0.1,
    show_default=True
)
@click.option(
    "-mcw",
    "--min-cis-weight",
    "min_cis_weight",
    help="Minimum conitg's cis weight of graph",
    type=float,
    default=5.0,
    show_default=True
)
@click.option(
    "-q1",
    "--min-quality1",
    "min_quality1",
    help="Minimum quality of hyperedge in first cluster, also control the minimum mapping quality of plotting heatmap.",
    type=click.IntRange(0, 60),
    default=1,
    show_default=True, 
)
@click.option(
    "-q2",
    "--min-quality2",
    "min_quality2",
    help="Minimum quality of hyperedge in second cluster.",
    type=click.IntRange(0, 60),
    default=2,
    show_default=True,
)
@click.option(
    "-mc",
    "--min-contacts",
    "min_contacts",
    help="Minimum contacts of contigs",
    metavar="INT",
    type=int,
    default=25,
    show_default=True
)
@click.option(
    "-ml",
    "--min-length",
    "min_length",
    help="Minimum length of contigs",
    metavar="INT",
    type=int,
    default=10000,
    show_default=True
)
@click.option(
    "-Nx",
    "--Nx",
    metavar="INT",
    help="Only retain contigs which length longer than Nx.",
    type=click.IntRange(0, 100),
    default=100,
    show_default=True, 
)
@click.option(
    "-ms",
    "--min-scaffold-length",
    "min_scaffold_length",
    help="The minimum length of the output scaffolding.",
    type=float,
    default=5e5,
    show_default=True
)
@click.option(
    "--enable-misassembly-remove",
    "enable_misassembly_remove",
    help="enable misassembly after second round partition.",
    is_flag=True,
    default=False,
    show_default=True
)
@click.option(
    "-wl",
    "--whitelist",
    metavar="PATH",
    help="""
    Path to 1-column list file containing
    contigs to include in hypergraph to partition.  
    """,
    default=None,
    show_default=True,
    type=click.Path(exists=True)
)
@click.option(
    '-scaf-method',
    '--scaffolding-method',
    'scaffolding_method',
    metavar='STR',
    help="The method of scaffolding, `['precision', 'allhic', 'fast']`."
    "precision: haphic_fastsort + allhic, which will quicker than allhic only. "
    "It is worth noting that `allhic` in `C-Phasing` is parameterized to "
    "achieve better results than the previous version",
    default="precision",
    type=click.Choice(["precision", "allhic", "fast"]),
    show_default=True
)
# @click.option(
#     '--factor',
#     '-k',
#     metavar="INT",
#     help='Factor of plot matrix. '
#             'If you input 10k matrix and want to plot heatmap at 500k, '
#             'factor should be set with 50.',
#     type=int,
#     default=50,
#     show_default=True
# )
@click.option(
    '-bs',
    '--binsize',
    metavar='STR',
    help='Bin size of the heatmap you want to plot. Enabled suffix with [k, m].',
    default='auto',
    show_default=True
)
@click.option(
    '--low-memory',
    help="Reduce memory usage. Only used in chimeric correct and hyperpartition, depend on input data size.",
    is_flag=True,
    default=False,
)
@click.option(
    '-o',
    '--output',
    '--outdir',
    'outdir',
    hidden=True,
    help="Output directory of the pipeline",
    metavar="STR",
    default="cphasing_output",
    type=str,
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
def pipeline(fasta, 
             ul_data,
            porec_data, 
            porectable, 
            pairs, 
            hic1,
            hic2,
            pattern,
            mapper_k,
            mapper_w,
            mapping_quality,
            chimeric_correct,
            chimeric_corrected,
            hic_mapper_k,
            hic_mapper_w,
            hcr,   
            hcr_lower,
            hcr_upper,
            collapsed_contig_ratio,
            hcr_bs,
            hcr_bed,
            hcr_invert,
            mode, 
            steps,
            skip_steps,
            use_existed_hitig,
            alleles_kmer_size,
            alleles_window_size,
            alleles_minimum_similarity,
            alleles_diff_thres,
            n,
            use_pairs,
            resolution1,
            resolution2, 
            init_resolution1,
            init_resolution2,
            first_cluster,
            disable_merge_in_first,
            exclude_group_from_first,
            exclude_group_to_second,
            normalize,
            allelic_similarity,
            min_allelic_overlap,
            min_weight,
            min_cis_weight,
            min_quality1,
            min_quality2,
            min_contacts,
            min_length,
            nx,
            min_scaffold_length,
            enable_misassembly_remove,
            whitelist,
            scaffolding_method, 
            binsize,
            low_memory,
            outdir,
            threads):
    """
    A pipeline of diploid or polyploid phasing and scaffolding (a: pipe).\n
    Also support for haploid scaffolding.
    > **Steps:**\n
        0. mapper : mapping sequencing data to reference\n
        1. alleles : identify the allelic contigs by self comparison\n
        2. prepare : prepare some datas for subsequence analysis.\n
        3. hyperpartiton : partition contigs into several groups\n
        4. scaffolding : ordering and orientation contigs \n
        5. plot : plot the heatmap of assembly\n

    > **Usages:**\n
    > - Input pore-c data\n
    ```bash
    $ cphasing pipeline -f contigs.fasta -pcd sample.fastq.gz -t 10 -s "all"\n
    ```
    > - Input multiple pore-c datas. It will use the sample1 as the prefix of output\n
    ```bash
    $ cphasing pipeline -f contigs.fasta -pcd sample1.fastq.gz -pcd sample2.fastq.gz -t 10
    ```
    > - Input pore-c table\n
    ```bash
    $ cphasing pipeline -f contigs.fasta -pct sample.porec.gz -t 10\n
    ```
    > - Input hic reads \n
    ```bash
    cphasing pipeline -f contigs.fasta -hic1 sample_R1.fastq.gz -hic2 sample_R2.fastq.gz -t 10\n
    ```
    > - Input pairs\n
    ```bash
    $ cphasing pipeline -f contigs.fasta -prs sample.pairs.gz -t 10\n
    ```
    > - Skip step:\n
    ```bash
    $ cphasing pipeline -f contigs.fasta -pct sample.porec.gz -t 10 -ss 1\n
    ```
    > - Only run a step:\n
    ```bash
    $ cphasing pipeline -f contigs.fasta -pct sample.porec.gz -t 10 -s 1\n
    ```
    > - Haploid mode, will skip 1.alleles and only run a round of partition:\n
    ```bash
    $ cphasing pipeline -f contigs.fasta -pct sample.porec.gz -t 10 --mode haploid\n
    ```

        
    """
    from .pipeline.pipeline import run 
    
    assert any([(porec_data is not None), (porectable is not None), 
                (pairs is not None), ((hic1 is not None) and (hic2 is not None))]), \
        "PoreC data or PoreC table or Hi-C data or Pairs must specified"
    
    if all([(len(porec_data) > 0), ((hic1 is not None) and (hic2 is not None))]):
        logger.warning("Simulataneously process Pore-C and Hi-C is not yet supported, only use Pore-C data for subsequently steps.")
        
    
    if steps:
        if steps == "all":
            steps = set(["0", "1", "2", "3", "4", "5"])

        else:
            steps = steps.strip().split(",")
    else:
        steps = set(map(str, [0, 1, 2, 3, 4]))

    binsize = binsize.lower()

    if skip_steps:
        skip_steps = set(skip_steps.strip().split(","))
    else:
        skip_steps = set()
    
    run(fasta, 
        ul_data,
        porec_data,
        porectable, pairs, 
        hic1=hic1,
        hic2=hic2,
        pattern=pattern,
        mapper_k=mapper_k, 
        mapper_w=mapper_w,
        mapping_quality=mapping_quality,
        hic_mapper_k=hic_mapper_k,
        hic_mapper_w=hic_mapper_w,
        chimeric_correct=chimeric_correct,
        chimeric_corrected=chimeric_corrected,
        hcr_flag=hcr,
        hcr_lower=hcr_lower,
        hcr_upper=hcr_upper,
        collapsed_contig_ratio=collapsed_contig_ratio,
        hcr_bs=hcr_bs, 
        hcr_bed=hcr_bed,
        hcr_invert=hcr_invert,
        mode=mode, 
        steps=steps,
        skip_steps=skip_steps,
        use_existed_hitig=use_existed_hitig,
        alleles_kmer_size=alleles_kmer_size,
        alleles_window_size=alleles_window_size,
        alleles_minimum_similarity=alleles_minimum_similarity,
        alleles_diff_thres=alleles_diff_thres,
        scaffolding_method=scaffolding_method,
        n=n,
        use_pairs=use_pairs,
        resolution1=resolution1,
        resolution2=resolution2,
        init_resolution1=init_resolution1,
        init_resolution2=init_resolution2,
        first_cluster=first_cluster,
        normalize=normalize,
        disable_merge_in_first=disable_merge_in_first,
        exclude_group_to_second=exclude_group_to_second,
        exclude_group_from_first=exclude_group_from_first,
        whitelist=whitelist,
        allelic_similarity=allelic_similarity,
        min_allelic_overlap=min_allelic_overlap,
        min_weight=min_weight,
        min_cis_weight=min_cis_weight,
        min_quality1=min_quality1,
        min_quality2=min_quality2,
        min_contacts=min_contacts,
        min_length=min_length,
        Nx=nx,
        min_scaffold_length=min_scaffold_length,
        enable_misassembly_remove=enable_misassembly_remove,
        binsize=binsize,
        outdir=outdir,
        threads=threads,
        low_memory=low_memory)
    

## Subcommand of UL ONT pipeline
from .hitig.cli import hitig 

@cli.command(cls=RichCommand, epilog=__epilog__)
@click.argument(
    "reference",
    type=click.Path(exists=True)
)
@click.argument(
    "fastq",
    nargs=-1,
    type=click.Path(exists=True),
    required=True
)
@click.option(
    "-e",
    "--enzyme",
    metavar="STR",
    default="",
    help="Restrict site pattern, use comma to separate multiple patterns.",
    show_default=True,
)
@click.option(
    "-k",
    "kmer_size",
    help="kmer size for mapping.",
    metavar="INT",
    type=int,
    default=15,
    show_default=True
)
@click.option(
    "-w",
    "window_size",
    help="minimizer window size for mapping.",
    metavar="INT",
    type=int,
    default=10,
    show_default=True
)
@click.option(
    "--mm2-params",
    metavar="STR",
    help="additional parameters for minimap2",
    default="-x map-ont",
    show_default=True
)
@click.option(
    '-q',
    '--mapq',
    help='Minimum quality of mapping [0, 60].',
    metavar='INT',
    type=click.IntRange(0, 60, clamp=True),
    default=0,
    show_default=True
)
@click.option(
    '-p',
    '--min-identity',
    'min_identity',
    help='Minimum percentage identity of alignments [0, 1.0].',
    metavar='FLOAT',
    type=click.FloatRange(0.0, 1.0, clamp=True),
    default=.8,
    show_default=True
)
@click.option(
    '-l',
    '--min-length',
    'min_length',
    help='Minimum length of fragments.',
    metavar='INT',
    default=150,
    show_default=True
)
@click.option(
    '-me',
    '--max-edge',
    'max_edge',
    help='Maximum length of fragment located in the edge of contigs.',
    metavar='INT',
    default=2000,
    show_default=True
)
@click.option(
    '--realign',
    help="realign to rescue multiple alignments, only support for contig-level",
    is_flag=True,
    default=False,
    show_default=True,
    hidden=True,
)
@click.option(
    '-f',
    '--force',
    help='Force run all the command, ignore existing results.'
    ' The index file also will be removed.',
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    '-o',
    '--outprefix',
    help='output prefix, if none use the prefix of fastq',
    default=None,
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
def mapper(reference, fastq, enzyme, kmer_size, 
            window_size, mm2_params, mapq, 
            min_identity, min_length, max_edge,
            force, realign, outprefix, threads):
    """
    Mapper for pore-c reads.

        REFERENCE: Path of reference

        FASTQ: Path of pore-c reads, multiple file enabled, the prefix of output default only use sample 1.

    """
    from .mapper import PoreCMapper

    additional_arguments = mm2_params.strip().split()
 
    pcm = PoreCMapper(reference, fastq,
                        pattern=enzyme,
                        k = kmer_size,
                        w = window_size,
                        force=force,
                        realign=realign,
                        min_quality=mapq,
                        min_identity=min_identity,
                        min_length=min_length,
                        max_edge=max_edge,
                        additional_arguments=additional_arguments,
                        outprefix=outprefix,
                        threads=threads)
    pcm.run()


@cli.group(cls=CommandGroup, epilog=__epilog__, short_help='Process Pore-C alignments.')
@click.pass_context
def alignments(ctx):
    pass

@alignments.command(cls=RichCommand, hidden=True, deprecated=True)
@click.argument(
    "paf",
    metavar="PAF",
    type=click.Path(exists=True)
)
@click.argument(
    "chromsize",
    metavar="CHROMSIZE",
    type=click.Path(exists=True)
)
@click.argument(
    "output",
    metavar='OUTPUT_PATH',
)
@click.option(
    '-q',
    '--min-quality',
    'min_quality',
    help='Minimum quality of mapping [0, 255].',
    metavar='INT',
    type=click.IntRange(0, 255, clamp=True),
    default=1,
    show_default=True
)
@click.option(
    '-p',
    '--min-identity',
    'min_identity',
    help='Minimum percentage identity of alignments [0, 1.0].',
    metavar='FLOAT',
    type=click.FloatRange(0.0, 1.0, clamp=True),
    default=.75,
    show_default=True
)
@click.option(
    '-l',
    '--min-length',
    'min_length',
    help='Minimum length of fragments.',
    metavar='INT',
    default=10,
    show_default=True,
    type=int,
)
def paf2pairs(paf, chromsize, output, 
                min_quality, min_identity, min_length):
    """
    Convert Pore-C alignments to 4DN pairs file.

        PAF_PATH : Path of alignment file(paf format).

        CHROMSIZE : Path of chromosome sizes.

        OUTPUT : PATH of output pairs file.
    """
    from .core import PAFTable
    paf = PAFTable(paf, threads=1, no_read=True,
                        min_quality=min_quality,
                        min_identity=min_identity,
                        min_length=min_length)
    paf.to_pairs(chromsize, output)
    paf.clean_tempoary()


@alignments.command(cls=RichCommand, hidden=True, deprecated=True)
@click.argument(
    "paf",
    type=click.Path(exists=True)
)
@click.argument(
    "output"
)
@click.option(
    '-q',
    '--min_quality',
    help='Minimum quality of mapping [0, 255].',
    metavar='INT',
    type=click.IntRange(0, 255, clamp=True),
    default=1,
    show_default=True
)
@click.option(
    '-p',
    '--min_identity',
    help='Minimum percentage identity of alignments [0, 1.0].',
    metavar='FLOAT',
    type=click.FloatRange(0.0, 1.0, clamp=True),
    default=.75,
    show_default=True
)
@click.option(
    '-l',
    '--min_length',
    help='Minimum length of fragments.',
    metavar='INT',
    default=10,
    show_default=True,
    type=int,
)
@click.option(
    '-cs',
    '--chunksize',
    help='Chunk size of edges',
    type=float,
    default=1000000,
    metavar='Float',
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
def paf2porec(paf, output, 
                min_quality, min_identity, min_length, 
                chunksize, threads):
    """
    Convert paf to pore_c_table.
    """

    from .core import PAFTable

    chunksize = int(chunksize)
    paf = PAFTable(paf, 
                    output=output,
                    threads=threads, 
                    chunksize=chunksize,
                    min_quality=min_quality, 
                    min_identity=min_identity, 
                    min_length=min_length)
    # paf.filter()
    # pore_c_table = paf.to_pore_c_table()

    # pore_c_table.save(output, paf.tmpdir)
    paf.clean_tempoary()



@alignments.command(cls=RichCommand, hidden=True, deprecated=True)
@click.argument(
    "pore_c_table",
    metavar="Pore-C-Table",
    type=click.Path(exists=True)
)
@click.option(
    "-rn",
    "--read_number",
    metavar="INT",
    help="Number or total reads.",
    default=None,
    show_default=True,
    type=int
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
@click.option(
    '-ud',
    '--use_dask',
    help="Use dask to lower memory usage.",
    is_flag=True,
    default=False,
    show_default=True 
)
def summary(pore_c_table, read_number, threads, use_dask):
    """
    Summary the pore-c table.

    Pore-C-Table : Path to pore-c table.

    """
    from .core import PoreCTable

    prefix = Path(pore_c_table).stem
    pct = PoreCTable(threads=threads, use_dask=use_dask)
    pct.read_table(pore_c_table)
    read_stat, alignment_stat = pct.read_and_alignment_stat_dask(read_number)
    read_stat.to_csv(f"{prefix}.read.summary", sep='\t', 
                            header=True, index=True)
    alignment_stat.to_csv(f"{prefix}.alignment.summary", sep='\t', 
                            header=True, index=True)
    
    contact_df = pct.contact_stat()
    contact_df["perc"] = contact_df["perc"].round(2)
    contact_df.to_csv(f"{prefix}.concatemer.summary", sep='\t', 
                            header=True, index=True)
    
@alignments.command(cls=RichCommand, hidden=True)
@click.argument(
    'pore_c_table',
    type=click.Path(exists=True),
    metavar='Pore_C_table'
)
@click.argument(
    'contig_bed',
    type=click.Path(exists=True),
    metavar='Contig_bed'
)
@click.argument(
    'output',
    type=click.Path(exists=False),
    metavar='Output'
)
@click.option(
    '-t',
    '--threads',
    type=int,
    default=1,
    show_default=True,
    help='Number of threads.'
)
def porec_chrom2contig(
    pore_c_table, 
    contig_bed,
    output,
    threads):
    """
    Convert chromosome-level pore_c_table to contig-level.

    Pore_c_table : Path of pore-c table.

    Contig_bed : Path of contig bed, three columns or four columns with contig id.

    Output : Path of output.
    """
    from .core import PoreCTable
    # df = pd.read_parquet(pore_c_table)

    contig_df = pd.read_csv(contig_bed, 
                            sep='\t', 
                            header=None,
                            index_col=None, 
                            )
    names = ['chrom', 'start', 'end', 'id']
    dtypes = ['object', 'int64', 'int64', 'object']
    if len(contig_df.columns) == 3:
        contig_df.columns = names[:3]
        contig_df = contig_df.astype(
            dict(zip(names[:3], dtypes[:3]))
        )
        contig_df['id'] = (contig_df.groupby('chrom', as_index=False)
                            .cumcount() + 1).astype(str) 
        contig_df = contig_df.assign(
            id=lambda x: x['chrom'].astype(str) + '.ctg' + x['id']
        )
    elif len(contig_df.columns) >= 4:
        contig_df = contig_df.iloc[:, 0:4]
        contig_df.columns = names
        contig_df = contig_df.astype(
            dict(zip(names, dtypes))
        )
    else:
        logger.warning(f'`{contig_bed}` must in 3-columns or 4-columns.')
        sys.exit()
    
    
    logger.info('Starting convert chromosome-level to contig-level ...')
    pct = PoreCTable()
    pct.read_table(pore_c_table)
    res_df = pct.chrom2contig(contig_df)

    res_df.to_parquet(output)
    logger.info(f'Done, output new Pore-C record in `{output}`')

@alignments.command(cls=RichCommand, hidden=True)
@click.argument(
    "pairs",
    metavar="Pairs",
    type=click.Path(exists=True)
)
@click.argument(
    "contig_bed",
    metavar="ContigBed",
    type=click.Path(exists=True)
)
@click.argument(
    "output",
    metavar="Output"
)
def pairs_chrom2contig(pairs, contig_bed, output):
    """
    Convert chromosome-level pairs to contig-level.

        Pairs : Path of pairs file.

        ContigBed : Path of contig coordinate to chrom in bed format.

        Output : Path of output contig-level pairs.


    """
    from .core import Pairs
    
    contig_df = pd.read_csv(contig_bed, 
                            sep='\t', 
                            header=None,
                            index_col=None, 
                            )
    names = ['chrom', 'start', 'end', 'id']
    dtypes = ['object', 'int64', 'int64', 'object']
    if len(contig_df.columns) == 3:
        contig_df.columns = names[:3]
        contig_df = contig_df.astype(
            dict(zip(names[:3], dtypes[:3]))
        )
        contig_df['id'] = (contig_df.groupby('chrom', as_index=False)
                            .cumcount() + 1).astype(str) 
        contig_df = contig_df.assign(
            id=lambda x: x['chrom'].astype(str) + '.ctg' + x['id']
        )
    elif len(contig_df.columns) >= 4:
        contig_df = contig_df.iloc[:, 0:4]
        contig_df.columns = names
        contig_df = contig_df.astype(
            dict(zip(names, dtypes))
        )
    else:
        logger.warning(f'`{contig_bed}` must in 3-columns or 4-columns.')
        sys.exit()

    p = Pairs(pairs)
    p.chrom2contig(contig_df, output=output)


@alignments.command(cls=RichCommand)
@click.argument(
    "pairs",
    metavar="Pairs",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    metavar="OUTPUT",
    default="-",
    show_default=True
)
@click.option(
    "-q",
    "--min-mapq",
    "min_mapq",
    default=1,
    metavar="INT",
    help="Minimum mapping quality of alignments",
    type=click.IntRange(0, 60),
    show_default=True,
)
def pairs2mnd(pairs, output, min_mapq):
    """
    convert 4DN pairs to mnd file.

        Pairs : Path of pairs file.

        Output : Path of output mnd file.
    """
    # from .core import Pairs 
    # p = Pairs(pairs)
    # p.to_mnd(output, threads) 
    cmd = ["cphasing-rs", "pairs2mnd", pairs, "-o", output, "-q", str(min_mapq)]
    flag = run_cmd(cmd, log=os.devnull)
    assert flag == 0, "Failed to execute command, please check log."


@alignments.command(cls=RichCommand)
@click.argument(
    'pairs',
    metavar='Pairs',
    type=click.Path(exists=True)
)
@click.argument(
    'bed',
    metavar='bed',
    type=click.Path(exists=True)
)
@click.argument(
    'output',
    metavar='Output'
)
def pairs_intersect(pairs, bed, output):
    """
    According a bed file to intersection a pairs file.

    Pairs : Path of Pairs.

    Bed : Path of Bed.
    
    Output : Path of output.
    """
    # from .core import Pairs
    # p = Pairs(pairs)
    # p.intersection(bed, output)
    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["cphasing-rs", "pairs-intersect", pairs, bed, "-o", output]
    flag = run_cmd(cmd, log=f"logs/pairs-intersect.log")
    assert flag == 0, "Failed to execute command, please check log."

@alignments.command(cls=RichCommand)
@click.argument(
    'table',
    metavar='Pore-C-Table',
    type=click.Path(exists=True)
)
@click.argument(
    'bed',
    metavar='Bed',
    type=click.Path(exists=True)
)
@click.argument(
    'output',
    metavar='Output'
)
# @click.option(
#     '-t',
#     '--threads',
#     help='Number of threads.',
#     type=int,
#     default=1,
#     metavar='INT',
#     show_default=True,
#     hidden=True,
# )
def porec_intersect(table, bed, output, threads):
    """
    According several regions to select contact 

    Pore-C-Table : Path to pore-c table.

    Bed : Path to bed file.

    Output : Path to output pore-c table.
    """
    # from .core import PoreCTable
    # pct = PoreCTable(threads)
    # pct.read_table(table)
    # pct.intersection(bed, output)
    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["cphasing-rs", "porec-intersect", table, "-o", output]
    flag = run_cmd(cmd, log=f"logs/porec-intersect.log")
    assert flag == 0, "Failed to execute command, please check log."

@alignments.command(cls=RichCommand)
@click.argument(
    'table',
    metavar='Pore-C-Table',
    type=click.Path(exists=True)
)
@click.argument(
    "contigsizes",
    metavar="ContigSizes",
    type=click.Path(exists=True),
)
@click.option(
    "--method",
    type=click.Choice(["binnify", "nparts"]),
    default="binnify",
    show_default=True
)
@click.option(
    "-n",
    "--nparts",
    type=int,
    default=5,
    show_default=True,
)
@click.option(
    "-bs",
    "--binsize",
    help="Bin size in bp.",
    type=int,
    default=10000,
    show_default=True
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout
)
def porec2csv(table, contigsizes, method, nparts, binsize, output):
    """
    Convert pore-c table to csv, which can import to paohvis
    """
    from .core import PoreCTable
    from .utilities import read_chrom_sizes

    contigsizes = read_chrom_sizes(contigsizes)

    pct = PoreCTable()
    pct.read_table(table)
    if method == "binnify":
        pct.binnify(contigsizes, binsize,)
    else:
        pct.divide_contig_into_nparts(contigsizes, nparts)
        
    pct.to_pao_csv(output)

@cli.command(cls=RichCommand, epilog=__epilog__)
@click.option(
    '-f',
    '--fasta',
    metavar="FASTA",
    help="Path to draft assembly",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '-prs',
    '--pairs',
    metavar="Pairs",
    help="4DN pairs file, which generated from `cphasing mapper`, or `cphasing hic mapper`",
    type=click.Path(exists=True),
    default=None,
    show_default=True
)
@click.option(
    '-d',
    '--depth',
    metavar="Depth",
    help="Cis long-range counts.",
    type=click.Path(exists=True),
    default=None,
    show_default=True,
    hidden=True
)
@click.option(
    '-w',
    '--window-size',
    'window_size',
    metavar='INT',
    default=500,
    show_default=True,
)
@click.option(
    '-telo-m',
    '--telo-motif',
    metavar="STR",
    help="5'-end of telomere motif. Filter out telomere fragment, avoid trim the telomere from a contig",
    default="CCCTAA",
    show_default=True 
)
@click.option(
    '-bp',
    '--break-pairs',
    help="Convert raw pairs file into chimeric corrected pairs file.",
    is_flag=True,
    default=True, 

)
@click.option(
    '-od',
    '--output-depth',
    help="Output cis depth file.",
    is_flag=True,
    default=False,
)
@click.option(
    '-o',
    '--outprefix',
    help='output prefix, if none use the prefix of fasta',
    default=None,
    show_default=True
)
@click.option(
    '--low-memory',
    help="Reduce memory usage.",
    is_flag=True,
    default=False,
)
@click.option(
    '-t',
    '--threads',
    help="Number of threads. ",
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
def chimeric(fasta, pairs, depth, window_size, 
             break_pairs, output_depth,
             telo_motif, outprefix, 
             low_memory, threads):
    """
    Correct chimeric contigs by contacts.
    """
    from .chimeric import run, correct
    assert any([pairs, depth]), "Pairs or Depth file must be input."

    if outprefix is None:
        outprefix = Path(Path(fasta).absolute().stem).name
    

    if pairs:
        run(fasta, pairs, window_size, 
                break_pairs=break_pairs, outprefix=outprefix,
                output_depth=output_depth, 
                telo_motif=telo_motif,
                low_memory=low_memory,
                threads=threads)
        
    elif depth:
        df = pd.read_csv(depth, sep='\t', header=None, index_col=None,
                        names=['contig', 'start', 'end', 'depth'])
        df = df.sort_values(['contig', 'start'])
    
        depth_dict = df.groupby('contig')['depth'].apply(np.array).to_dict()
        break_point_res = correct(depth_dict, window_size=window_size, threads=threads)

        if len(break_point_res) != 0:
            break_point_res = break_point_res.to_csv("output.breakPos.txt", 
                                            sep='\t', index=None, header=None)


@cli.command(cls=RichCommand, epilog=__epilog__, 
                short_help='Only retain the HCRs from Pore-C data')
@click.option(
    "-f",
    "--fasta",
    metavar="FILE",
    help="polyploid contig-level fasta.",
    type=click.Path(exists=True),
)
@click.option(
    '-pct',
    '--porectable',
    metavar="PoreC-Table",
    help="PoreC Table file, which generated from `cphasing mapper`",
    type=click.Path(exists=True),
    default=None,
    show_default=True
)
@click.option(
    '-prs',
    '--pairs',
    metavar="Pairs",
    help="4DN pairs file, which generated from `cphasing mapper`, or `cphasing hic mapper`",
    type=click.Path(exists=True),
    default=None,
    show_default=True
)
@click.option(
    "-cs",
    "--contigsize",
    help="Two columns file of contig name and length.",
    metavar="CONTIG_SIZE",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '-p',
    '--pattern',
    help="Pattern of restriction enzyme. Comma separated for multiple pattern."
        "  If specified, depth will normalization by RE counts, "
        "which will increase the performance of `HindIII` or other distribution biased RE.",
    metavar="STR",
    default=None,
    show_default=True,
)
@click.option(
    "-bs",
    "--binsize",
    help="Bin size in bp.",
    type=int,
    default=10000,
    show_default=True
)
@click.option(
    "-l",
    "--lower",
    help="Lower value of peak value.",
    type=float,
    default=.01,
    show_default=True,
)
@click.option(
    "-u",
    "--upper",
    help="Upper value of peak value.",
    type=float,
    default=1.5,
    show_default=True,
)
@click.option(
    "-cr",
    "--collapsed-contigs-ratio",
    "collapsed_contig_ratio",
    metavar="FLOAT",
    help="minimum regions of the collapsed contig, "
    " which mean the region of a contig higher than upper "
    "value will be regarded as collapsed and remove the whole contigs",
    type=click.FloatRange(0.0, 1.0),
    default=.0,
    show_default=True,
)
@click.option(
    "-b",
    "--bed",
    help="HCRs bed file",
    default=None,
    show_default=True
)
@click.option(
    "-v",
    "--invert",
    help="Invert regions",
    is_flag=True,
    default=False,
    show_default=True
)
@click.option(
    '-q',
    '--min_quality',
    help='Minimum quality of mapping `[0, 60]`.',
    metavar='INT',
    type=click.IntRange(0, 60, clamp=True),
    default=0,
    show_default=True
)
@click.option(
    '--fofn',
    help="""
    If this flag is set then the SRC_ASSIGN_TABLES is a file of
    filenames corresponding to the contact tables you want to merge.
    This is workaround for when the command line gets too long.
    """,
    is_flag=True,
    default=False,
    show_default=True
)
@click.option(
    "-o",
    "--output",
    metavar="STR",
    help="Output high confidence contacts, '.gz' supported.",
    default=None,
    show_default=True
)
def hcr(fasta, porectable, pairs, contigsize, 
        pattern, binsize, 
        lower, upper, collapsed_contig_ratio,
          bed, invert, min_quality, fofn, output):
    """
    Only retain high confidence regions to subsequence analysis.
        High confidence regions are identified from contacts.

    """
    from .hitig.hcr.hcr_by_contacts import hcr_by_contacts
    from .utilities import is_file_changed

    assert porectable or pairs, "The Pore-C table or 4DN pairs file should be provided at least one."

    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    if porectable:
        if fofn:
            pore_c_tables = [i.strip() for i in open(porectable) if i.strip()]
        else:
            pore_c_tables = [porectable]

        if not pairs:
            pairs_files = []
            for pore_c_table in pore_c_tables:
                path = Path(pore_c_table).parent
         
                prefix = Path(Path(pore_c_table).stem).with_suffix("")
                
                while prefix.suffix in {'.gz', 'gz', 'porec', ".porec"}:
                    prefix = prefix.with_suffix('')
                
                pairs = f"{path}/{prefix}.pairs.gz"
                pairs_files.append(pairs)
                if not Path(pairs).exists() or is_compressed_table_empty(pairs):
                    cmd = ["cphasing-rs", "porec2pairs", porectable, contigsize,
                    "-o", pairs, "-q", "0"]
                
                    flag = run_cmd(cmd, log=f"logs/porec2pairs.log")
                    assert flag == 0, "Failed to execute command, please check log."
                else:
                    logger.warning(f"Use exists 4DN pairs file of `{pairs}`.")
        else:
            prefix = Path(pairs).with_suffix("")
            while prefix.suffix in {'.gz', 'gz', '.pairs'}:
                prefix = prefix.with_suffix('')
        
    else:
        if fofn:
            pairs_files = [i.strip() for i in open(pairs) if i.strip()]
        else:
            pairs_files = [pairs]
        
        for pairs in pairs_files:
            prefix = Path(Path(pairs).stem).with_suffix("")
            while prefix.suffix in {'.gz', 'gz', '.pairs'}:
                prefix = prefix.with_suffix('')

    depth_file = f"{prefix}.{binsize}.depth"

    if bed:
        if is_file_changed(str(pairs)) or is_file_changed(bed) or not Path(depth_file).exists():
            is_file_changed(bed)
            logger.info("Calculating the depth of pairs ...")
            cmd = ['cphasing-rs', 'pairs-intersect', 
                   '-q', str(min_quality),
                    str(pairs), str(bed),
                    f"2>logs/pairs2intersect.1.log", "|",
                   'cphasing-rs', 'pairs2depth', "-", "-o", depth_file, 
                    f"2>logs/pairs2depth.log",]

            flag = os.system(" ".join(cmd))
            assert flag == 0, "Failed to execute command, please check log."

        else:
            logger.warning(f"Use exists depth file of `{depth_file}`.")
    
    else:
        
        if is_file_changed(str(pairs)) or not Path(depth_file).exists():
            cmd = ['cphasing-rs', 'pairs2depth', 
                   '-q', str(min_quality),
                   str(pairs), "-o", depth_file]
            flag = run_cmd(cmd, log=f"logs/pairs2depth.log")
            assert flag == 0, "Failed to execute command, please check log."

        else:
            logger.warning(f"Use exists depth file of `{depth_file}`.")

    if pattern and fasta:
        logger.info("Normalizing each bin by RE counts ...")
        cmd = ["bedtools", "getfasta", "-fi", str(fasta), 
               "-bed", str(depth_file), "2>/dev/null",
               ">", "tmp.depth.fasta"]
        os.system(" ".join(cmd))
        cmd = ["cphasing-rs", "count_re", "--pattern", str(pattern), "tmp.depth.fasta", 
               "2>/dev/null", "-o", "tmp.depth.countre.txt"]
        os.system(" ".join(cmd))

        df1 = pd.read_csv(depth_file, sep='\t', header=None, index_col=None, 
                            names=['chrom', 'start', 'end', 'depth'])
        df2 = pd.read_csv("tmp.depth.countre.txt", sep='\s+', header=0, index_col=None, 
                            usecols=[0, 1],
                            names=['item', "count"])

        if Path("tmp.depth.fasta").exists():
            os.remove("tmp.depth.fasta")
        if Path("tmp.depth.countre.txt").exists():
            os.remove("tmp.depth.countre.txt")

        df1['item'] = df1['chrom'] + ':' + df1['start'].astype(str) + '-' + df1['end'].astype(str)
        df1.set_index('item', inplace=True)
        df2.set_index('item', inplace=True)
        
        df = pd.concat([df1, df2], axis=1).dropna().reset_index().drop('item', axis=1)

        df['depth'] = df['depth'] / df['count']
        df.fillna(0, inplace=True)
        df.drop('count', axis=1, inplace=True)
        depth_file = depth_file.replace(".depth", ".norm.depth")
        df.to_csv(depth_file, header=None, index=None, sep='\t')
        
       
    hcr_by_contacts(depth_file, f"{prefix}.{binsize}.hcr.bed" , 
                        lower, upper, collapsed_contig_ratio)
    bed =  f"{prefix}.{binsize}.hcr.bed"

    if invert:
        logger.info("Invert regions.")
        import pyranges as pr
        contigsize = pd.read_csv(contigsize, sep='\t', header=None, index_col=None)
        contigsize.columns = ['Chromosome', 'End']
        contigsize['Start'] = 0

        bed_df = pd.read_csv(bed, sep='\t', header=None, index_col=None, 
                                names=['Chromosome', 'Start', 'End'])

        bed_gr = pr.PyRanges(bed_df)
        contig_gr = pr.PyRanges(contigsize)
        bed_gr = contig_gr.subtract(bed_gr)
        bed_df = bed_gr.df

        if bed_df.empty:
            logger.warning("No inverted regions.")
            
        else:
            total_length = (bed_df['End'] - bed_df['Start']).sum()
            total_length_of_contigs = contigsize['End'].sum()
            logger.info(f"Total length of inverted regions is {total_length / total_length_of_contigs:.02%}.")
        bed_df.to_csv(bed, sep='\t', header=None, index=None)
        


@cli.command(cls=RichCommand, epilog=__epilog__, short_help='Prepare data for subsequence analysis.')
@click.argument(
    "fasta",
    metavar="INPUT_FASTA_PATH",
    type=click.Path(exists=True)
)
@click.argument(
    'pairs',
    metavar='INPUT_PAIRS_PATH',
    type=click.Path(exists=True)
)
@click.option(
    '-p',
    '--pattern',
    help='Pattern of restrict enzyme. Commam separate for multiple pattern.'
    ' Suggest using the default parameter, different restrict enzymes do '
    'not effect any results.',
    metavar="STR",
    default="AAGCTT",
    show_default=True,
)
@click.option(
    "-q",
    "--min-mapq",
    "min_mapq",
    default=0,
    metavar="INT",
    help="Minimum mapping quality of alignments",
    type=click.IntRange(0, 60),
    show_default=True,
)
@click.option(
    '-mc',
    '--min-contacts',
    'min_contacts',
    help='Minimum contacts between contig pairs',
    metavar='INT',
    default=1,
    show_default=True,
)
@click.option(
    '-spc',
    '--skip-pairs2contacts',
    'skip_pairs2contacts',
    help='skip convert pairs to contacts',
    default=False,
    is_flag=True,
)
@click.option(
    '--skip-pairs2clm',
    'skip_pairs2clm',
    help='skip convert pairs to clm',
    default=False,
    is_flag=True,
)
@click.option(
    '--low-memory',
    help="Reduce memory usage.",
    is_flag=True,
    default=True,
    hidden=True
)
@click.option(
    '-t',
    '--threads',
    help="Number of threads. ",
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
@click.option(
    '-o',
    '--outprefix',
    help='output prefix, if none use the prefix of pairs',
    default=None,
    show_default=True
)
def prepare(fasta, pairs, min_mapq,
            min_contacts, pattern,
            skip_pairs2contacts, skip_pairs2clm,
            low_memory, threads, outprefix):
    """

    """
    from .prepare import pipe 

    pipe(fasta, pairs, pattern, 
         min_mapq,
         min_contacts, skip_pairs2clm=skip_pairs2clm,
         skip_pairs2contacts=skip_pairs2contacts,
         low_memory=low_memory,
         threads=threads, outprefix=outprefix)


@cli.command(cls=RichCommand, epilog=__epilog__)
@click.option(
    "-f",
    "--fasta",
    metavar="FILE",
    help="polyploid contig-level fasta.",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    "-o",
    "--output",
    metavar="PATH",
    help="path of output allele table [default: fasta_prefix.allele.table]",
    default=None
)
@click.option(
    "-k",
    "kmer_size",
    help="kmer size for similarity calculation.",
    metavar="INT",
    type=int,
    default=19,
    show_default=True
)
@click.option(
    "-w",
    "window_size",
    help="minimizer window size for similarity calculation.",
    metavar="INT",
    type=int,
    default=19,
    show_default=True
)
@click.option(
    "-c",
    "max_occurance",
    help="max occurance of allelic pairs",
    metavar="INT",
    type=int,
    default=60,
    show_default=True
)
@click.option(
    "-n",
    "min_cnt",
    help="minimum chain end",
    metavar="INT",
    type=int,
    default=5,
    show_default=True
)
@click.option(
    "-m",
    "minimum_similarity",
    help="minimum k-mer similarity for similarity calculation.",
    metavar="FLOAT",
    type=click.FloatRange(0, 1),
    default=.5,
    show_default=True
)
@click.option(
    "-d",
    "diff_thres",
    help="minimum different threshold between contig pairs.",
    type=click.FloatRange(0, 1),
    default=.1,
    show_default=True
)
@click.option(
    "--min-length",
    "min_length",
    help="minimum the sum of valid kmer length",
    type=int,
    default=100,
    show_default=True,
    hidden=True
)
@click.option(
    "-wl",
    "--whitelist",
    metavar="PATH",
    help="""
    Path to 1-column list file containing
    contigs to include in hypergraph to partition.  
    """,
    default=None,
    show_default=True,
    type=click.Path(exists=True),
    hidden=True
)
@click.option(
    '-fc',
    '--first-cluster',
    'first_cluster',
    metavar="PATH",
    help='Use the first cluster results to execute alleles',
    default=None,
    show_default=True,
    hidden=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
def alleles(fasta, output, 
                kmer_size, window_size, 
                max_occurance, min_cnt,
                minimum_similarity, diff_thres,
                whitelist, first_cluster,
                min_length, threads):
    """
    Build allele table by kmer similarity.

    """
    
    from .alleles import PartigAllele
    from .core import AlleleTable, ClusterTable
    from .utilities import is_file_changed
    from joblib import Parallel, delayed 

    if not output:
        fasta_prefix = Path(Path(fasta).name).with_suffix("")
        while fasta_prefix.suffix in {".fasta", "gz", "fa", '.fa', '.gz'}:
            fasta_prefix = fasta_prefix.with_suffix("")
    
        output = f"{fasta_prefix}.allele.table"
    
    # if whitelist:
    #     whitelist = set([i.strip() for i in open(whitelist) if i.strip()])

    if not first_cluster:
        pa = PartigAllele(fasta, kmer_size, window_size, max_occurance, 
                            min_cnt, minimum_similarity, diff_thres,
                            filter_value=min_length, output=output, threads=threads)
        pa.run()
    else:
        ct = ClusterTable(first_cluster)
        fasta = Path(fasta).absolute()
        Path("alleles_workdir/").mkdir(exist_ok=True)
        fastas = ct.to_fasta(fasta, outdir="alleles_workdir")
        fastas = list(map(lambda x: Path(x).name, fastas))
        os.chdir("alleles_workdir")
        
        args = []
        allele_tables = []
        for fa in fastas:
            fasta_prefix2 = Path(Path(fa).name).with_suffix("")
            while fasta_prefix2.suffix in {".fasta", "gz", "fa", '.fa', '.gz'}:
                fasta_prefix2 = fasta_prefix2.with_suffix("")
            output2 = f"{fasta_prefix2}.allele.table"
            allele_tables.append(f"alleles_workdir/{output2}")
            if Path(output2).exists():
                is_file_changed(output2)
            if not is_file_changed(fa) and Path(output2).exists() and not is_file_changed(output2):
                continue
            
            args.append((fa, kmer_size, window_size, max_occurance, 
                            min_cnt, minimum_similarity, diff_thres,
                            min_length, output2, 4))
            
        def func(fa, kmer_size, window_size, max_occurance, 
                    min_cnt, minimum_similarity, diff_thres,
                    min_length, output,threads):
            
            pa = PartigAllele(fa, kmer_size, window_size, max_occurance, 
                            min_cnt, minimum_similarity, diff_thres,
                            min_length, output, threads)
            pa.run()

        logger.info("Identifing allelic contig pairs in each homologous group ...")
        if len(args) > 0:
            Parallel(n_jobs=min(threads, len(args)))(delayed(func)(*a) for a in args)
        else:
            logger.info("Load existing allele table to generate total allele table.")
           
        for fa in fastas:
            if Path(fa).exists():
                os.remove(fa)
            if Path(f"{fa}.fai").exists():
                os.remove(f"{fa}.fai")

        os.chdir('..')
        ## merge allele table 
        logger.setLevel(logging.WARNING)
        allele_df = pd.concat([AlleleTable(at, sort=False, fmt="allele2").data for at in allele_tables], axis=0)
        allele_df = allele_df.reset_index(drop=True).reset_index().reset_index()
        logger.setLevel(logging.INFO)

        with open(output, 'w') as out:
            for at in allele_tables:
                with open(at) as fp:
                    for line in fp:
                        if line.startswith("#"):
                            print(line.strip(), file=out)

            allele_df.to_csv(out, sep='\t', header=None, index=None)

        logger.info(f"Successful output allele table in `{output}`")
        # for at in allele_tables:
        #     if Path(at).exists():
        #         os.remove(at)
        #     if Path(f"{at.replace('.allele.table', '')}.similarity.res").exists():
        #         os.remove(f"{at.replace('.allele.table', '')}.similarity.res")


@cli.command(cls=RichCommand, short_help="Build allele table by self mapping, lower memory usage, but higher errors.")
@click.option(
    "-f",
    "--fasta",
    metavar="FILE",
    help="polyploid contig-level fasta.",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    "-n",
    "--ploidy",
    help="ploidy level of genome.",
    metavar="INT",
    type=int,
    required=True
)
@click.option(
    "-k",
    "--kmer",
    metavar="INT",
    type=int,
    help="kmer size for wfmash",
    default=19,
    show_default=True
)
@click.option(
    "-w",
    "--sketch-size",
    metavar="INT",
    type=int,
    help="window size for wfmash",
    default=19,
    show_default=True
)
@click.option(
    "-s",
    "--segment-length",
    "segment_length",
    help="segment length for wfmash",
    metavar="STR",
    default="5k",
    show_default=True
)
@click.option(
    "-l",
    "--block-length",
    "block_length",
    help="keep merged mappings supported by homologies of this total length for wfmash",
    metavar="STR",
    default="25k",
    show_default=True,
)
@click.option(
    "-p",
    "--map-pct-id",
    "identity",
    help="Percent identity in the wfmash",
    metavar="INT",
    type=click.IntRange(50, 100),
    default=95,
    show_default=True,
)
@click.option(
    "--kmer-threshold",
    "kmer_threshold",
    help="Ignore the top % most-frequent kmers for wfmash",
    metavar="FLOAT",
    type=click.FloatRange(0.0, 100.0),
    default=0.001,
    show_default=True,
)
@click.option(
    "-r",
    "--output-raw-table",
    "output_raw_table",
    help="Ouput the raw format allele table, which decripted in ALLHiC.",
    default=False,
    is_flag=True,
    show_default=True
)
@click.option(
    "-o",
    "--output",
    metavar="PATH",
    help="path of output allele table [default: fasta_prefix.allele.table]",
    default=None
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=8,
    metavar='INT',
    show_default=True,
)
def alleles2(fasta, ploidy, kmer, sketch_size, segment_length, block_length,
                identity, kmer_threshold, output_raw_table, output, threads):
    """
    Identify the allelic contigs by self-mapping, lower memory usage than alleles but more errors.
        
        And, `-mao` of `hyperpartition` should be set to `0.8` and `-as` set to `0.95`.

    """
    from .alleles import AlignmentAlleles
    no_raw_table = False if output_raw_table else True
    aa = AlignmentAlleles(fasta, ploidy, 
                            k=kmer, 
                            w=sketch_size,
                            s=segment_length,
                            l=block_length,
                            H=kmer_threshold,
                            p=identity,
                            no_raw_table=no_raw_table,
                            threads=threads, 
                            output=output)
    aa.run()

@cli.command(cls=RichCommand, hidden=True,
             short_help="Build allele table by self mapping, lower memory usage, but higher errors.")
@click.option(
    "-f",
    "--fasta",
    metavar="FILE",
    help="polyploid contig-level fasta.",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    "-n",
    "--ploidy",
    help="ploidy level of genome.",
    metavar="INT",
    type=int,
    default=10,
)
@click.option(
    "-p",
    "--map-pct-id",
    "identity",
    help="Percent identity in the wfmash",
    metavar="INT",
    type=click.IntRange(50, 100),
    default=95,
    show_default=True,
)
@click.option(
    "-o",
    "--output",
    metavar="PATH",
    help="path of output allele table [default: fasta_prefix.allele.table]",
    default=None
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=8,
    metavar='INT',
    show_default=True,
)
def alleles3(fasta, ploidy, 
                identity, output, threads):
    """
    Identify the allelic contigs by self-mapping, lower memory usage than alleles but more errors.
        
        And, `-mao` of `hyperpartition` should be set to `0.8` and `-as` set to `0.95`.

    """
    from .alleles import ANIAlleles

    aa = ANIAlleles(fasta, ploidy, 
                            percent=identity,
                            threads=threads, 
                            output=output)
    
    aa.run()


@cli.command(cls=RichCommand, epilog=__epilog__, short_help="Generate the allelic contig and cross-allelic contig table.")
@click.argument(
    'alleletable',
    metavar='AlleleTable',
    type=click.Path(exists=True)
)
@click.argument(
    'contacts',
    metavar='INPUT_CONTACTS_PATH',
    type=click.Path(exists=True)
)
# @click.argument(
#     "countRE",
#     metavar="INPUT_COUNT_RE_PATH",
#     type=click.Path(exists=True)
# )
# @click.option(
#     '-c',
#     '--countRE',
#     metavar="CountRE",
#     help="input count RE table to normalize the contacts",
#     type=click.Path(exists=True),
#     default=None,
#     show_default=True
# )
# @click.option(
#     "-wl",
#     "--whitelist",
#     metavar="PATH",
#     help="""
#     Path to 1-column list file containing
#     contigs to include in hypergraph to partition.  
#     """,
#     default=None,
#     show_default=True,
#     type=click.Path(exists=True),
#     hidden=True
# )
# @click.option(
#     '-ns',
#     '--no-sort',
#     'no_sort',
#     help='Do not sort the prune table by similarity.',
#     metavar="BOOL",
#     is_flag=True,
#     default=True,
#     show_default=True,
#     hidden=True
# )
@click.option(
    '-o',
    '--output',
    help='output prune table',
    default='prune.contig.table',
    show_default=True,
)
# @click.option(
#     '-cs',
#     '--chunksize',
#     help='chunksize of contacts data',
#     default=100000,
#     show_default=True,
# )
@click.option(
    '-m',
    '--method',
    help='Method of cross-allelic, fast or greedy.',
    default='precise',
    type=click.Choice(["fast", "precise", "greedy"]),
    show_default=True
)
@click.option(
    '-fc',
    '--first-cluster',
    'first_cluster',
    metavar="PATH",
    help='Use the first cluster results to execute kprune',
    default=None,
    show_default=True
)
@click.option(
    '-n',
    '--norm-method',
    'norm_method',
    metavar="STR",
    help="Normalization method of contacts for kprune: none|cis|cis_unique",
    default="none",
    show_default=True ,
    type=click.Choice(["none", "cis", "cis_unique", "auto"])
)
@click.option(
    '-t',
    '--threads',
    help="Number of threads. ",
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
# @click.option(
#     '--symmetric',
#     help='output the symmetric contig pairs',
#     is_flag=True,
#     default=False,
#     show_default=True
# )
def kprune(alleletable, contacts, 
            output, method, first_cluster,
            norm_method, threads):
    """
    Generate the allelic contig and cross-allelic contig pairs by sequences similarity.

        Note: it is already intergated in `hyperpartition`.

        AlleleTable : allele table from cphasing allele in allele2 format.

        INPUT_CONTACTS_PATH : path of whole contigs contacts from `prepare`.
    

    """
    from .kprune import KPrunerRust
    
    # if whitelist:
    #     whitelist = [i.strip() for i in open(whitelist) if i.strip()]
    # sort_by_similarity = False if no_sort else True
    
    kp = KPrunerRust(alleletable, contacts, 
                     output, method=method,
                     first_cluster=first_cluster, 
                     norm_method=norm_method,
                      threads=threads)
    kp.run()
    # kp.save_prune_list(output, symmetric)


@cli.command(cls=RichCommand, epilog=__epilog__)
@click.argument(
    "contacts",
    metavar="Contacts",
    type=click.Path(exists=True)
)
@click.argument(
    "contigsize",
    metavar="Contig_sizes",
    type=click.Path(exists=True),
)
@click.argument(
    "output",
    metavar="OUTPUT"
)
@click.option(
    "-min",
    "--min-order",
    "min_order",
    help="Minimum contig order of pore-c reads",
    metavar="INT",
    type=int,
    default=2,
    show_default=True,
    hidden=True
)
@click.option(
    "-max",
    "--max-order",
    "max_order",
    help="Maximum contig order of pore-c reads",
    metavar="INT",
    type=int,
    default=50,
    show_default=True,
    hidden=True
)
@click.option(
    "-ma",
    "--min-alignments",
    "min_alignments",
    help="Minimum length of pore-c alignments",
    metavar="INT",
    type=int,
    default=30,
    show_default=True,
    hidden=True
)
@click.option(
    '-q',
    '--min_quality',
    help='Minimum quality of mapping [0, 255].',
    metavar='INT',
    type=click.IntRange(0, 255, clamp=True),
    default=1,
    show_default=True
)
@click.option(
    '-prs',
    '--pairs',
    '--hic',
    help="""
    construct common graph from 4DN pairs file. 
    Phasing by Hi-C data, should add this parameter.
    """,
    is_flag=True,
    default=False,
    show_default=True
)
@click.option(
    '--fofn',
    help="""
    If this flag is set then the SRC_ASSIGN_TABLES is a file of
    filenames corresponding to the contact tables you want to merge.
    This is workaround for when the command line gets too long.
    """,
    is_flag=True,
    default=False,
    show_default=True
)
@click.option(
    '-s',
    '--split',
    metavar="INT",
    help="""
    Split the contig into specify number bins.
    """,
    type=int,
    default=None,
    show_default=True,
    hidden=True
)
@click.option(
    "-wl",
    "--whitelist",
    metavar="PATH",
    help="""
    Path to 1-column list file containing
    contigs to include in hypergraph to partition.  
    """,
    default=None,
    show_default=True,
    type=click.Path(exists=True)
)
@click.option(
    "-hcr-bed",
    "--hcr-bed",
    "hcr_bed",
    metavar="STR",
    default=None,
    help="Only retain high confidence regions to subsequence analysis by input bed.",
    show_default=True,
    type=click.Path(exists=True)
)
@click.option(
    "-hcr-invert",
    "--hcr-invert",
    "hcr_invert",
    is_flag=True,
    default=False,
    help="Invert table by hcr bed.",
    show_default=True,
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads. Only when multiple files are input.',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
def hypergraph(contacts,
            contigsize,
            output,
            min_order, 
            max_order, 
            min_alignments, 
            min_quality,
            pairs, 
            fofn,
            split,
            whitelist,
            hcr_bed,
            hcr_invert,
            threads):
    """
    Construct hypergraph from 4DN pairs or porec table (a: hg).

    The hypergraph or graph of pore-c or hic enabled extract from pore-c table or 4DN pairs. 

        Contacts : Path of Pore-C table or 4DN pairs.\n
        
        Contig_sizes : Path of contig sizes.\n

        OUTPUT : Path of output hypergraph.\n

    """
    from .hypergraph import (
        HyperExtractor, 
        HyperExtractorSplit,
        Extractor,
        ExtractorSplit
        )
    from .utilities import read_chrom_sizes

    contigsizes = read_chrom_sizes(contigsize)
    contigs = contigsizes.index.values.tolist()
    # contigs = natsorted(contigsizes.index.values.tolist())

    if whitelist:
        whitelist = set([i.strip() for i in open(whitelist) if i.strip()])
        contigs = list(filter(lambda x: x in whitelist, contigs))

    if not hcr_bed and hcr_invert:
        logger.warning("The `--hcr-invert` parameter is only valid when `--hcr-bed` is set.")

    contig_idx = defaultdict(None, dict(zip(contigs, range(len(contigs)))))

    if not pairs:
        if fofn:
            pore_c_tables = [i.strip() for i in open(contacts) if i.strip()]
        else:
            pore_c_tables = contacts

        if not split:
            he = HyperExtractor(pore_c_tables, contig_idx, contigsizes.to_dict()['length'], 
                                min_order, max_order, min_alignments, 
                                hcr_bed=hcr_bed, hcr_invert=hcr_invert,
                                min_quality=min_quality, threads=threads)
            he.save(output)
        else:
            he = HyperExtractorSplit(pore_c_tables, contig_idx, contigsizes.to_dict()['length'], 
                                     split, min_quality, threads)
            he.save(output)
    
    else:
        if fofn:
            pairs_files = [i.strip() for i in open(contacts) if i.strip()]

        else:
            pairs_files = contacts
        
        if not split:
            
            e = Extractor(pairs_files, contig_idx, contigsizes.to_dict()['length'], 
                            hcr_bed= hcr_bed, hcr_invert=hcr_invert,
                            min_quality=min_quality, threads=threads)
            e.save(output)
        else:
            e = ExtractorSplit(pairs_files, contig_idx, contigsizes.to_dict()['length'], split, threads=threads)
            e.save(output)


@cli.command(cls=RichCommand, epilog=__epilog__)
@click.argument(
    "hypergraph",
    metavar="HyperGraph",
    type=click.Path(exists=True)
)
@click.argument(
    "contigsizes",
    metavar="ContigSizes",
    type=click.Path(exists=True),
)
@click.argument(
    "output",
    metavar="Output",
)
@click.option(
    '-prs',
    '--pairs',
    '--hic',
    help="""
    construct common graph from 4DN pairs file. 
    Phasing by Hi-C data, should add this parameter.
    """,
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    '-pct',
    '--porec',
    help="""
    construct hypergraph from porec table. 
    Input porec table, should add this parameter.
    """,
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    '-c',
    '--contacts',
    help="""
    contacts file for kprune, generate from prepare
    """,
    default=None, 
    show_default=True,
    type=click.Path(exists=True), 
)
@click.option(
    "-ul",
    "--ultra-long",
    "ultra_long",
    metavar="HyperGraph",
    help="""
    Add another linking information which generated from ultra-long reads. 
    Input file must the hypergraph format.
    """,
    type=click.Path(exists=True),
    default=None,
    show_default=True,
    hidden=True,
)
@click.option(
    "-ul-w",
    "--ul-weight",
    metavar="Float",
    help="""
    The weight of ultra-long hypergraph.
    """,
    default=1,
    type=int,
    show_default=True,
    hidden=True,
)
@click.option(
    '--mode',
    help="mode of hyperpartition, conflict of `-inc`. The basal equal to phasing. [default: None]",
    type=click.Choice(["basal", "haploid", "phasing", "basal_withprune", "None"]),
    default=None,
    show_default=True
)
@click.option(
    "-n",
    help="""
    Number of groups. If set to 0 or None, the partition will be run automatically.
    If in phasing mode, you can set to n1:n2, 
    which meaning that generate n1 groups in the first round partition
    and generate n2 groups in the second round partition. 
    Also, you can set `-n 8:0` to set the second round partition to automatical mode.
     The n2 parameter also can be set a file, which containing two columns, 
      the firse column repersent the group index from `first.clusters.txt`, and the 
       second column repersent the n2 groups in each n1 group. [default: 0]
    """,
    default=None,
    show_default=True,
)
@click.option(
    "-f",
    "--fasta",
    metavar="FILE",
    help="polyploid contig-level fasta. if input fasta, allele table will generate in this step when alleletable is None and prune table is None.",
    type=click.Path(exists=True),
)
@click.option(
    '-alleles-k',
    '--alleles-k',
    '--alleles-kmer-size',
    'alleles_kmer_size',
    metavar="INT",
    help="kmer size for `alleles` similarity calculation.",
    default=19,
    show_default=True,
)
@click.option(
    '-alleles-w',
    '--alleles-w',
    '--alleles-window-size',
    'alleles_window_size',
    metavar="INT",
    help="minimizer window size for `alleles` similarity calculation.",
    default=19,
    show_default=True,
)
@click.option(
    '-alleles-m',
    '--alleles-m',
    '--alleles-minimum-similarity',
    'alleles_minimum_similarity',
    metavar="FLOAT",
    help="minimum k-mer similarity for `alleles` similarity calculation.",
    default=0.5,
    show_default=True,
)
@click.option(
    "-alleles-d",
    '--alleles-d',
    '--alleles-diff-thres',
    "alleles_diff_thres",
    help="minimum different threshold between contig pairs.",
    type=click.FloatRange(0, 1),
    default=.1,
    show_default=True
)
@click.option(
    "-at",
    "--alleletable",
    metavar="PATH",
    help="""
    Path to allele table that is generated from `alleles`. 
        Skip when prunetable parameter added. [default: None]
    """,
    default=None,
    show_default=True,
    type=click.Path(exists=True),
)
@click.option(
    "-pt",
    "--prunetable",
    metavar="PATH",
    help="Path to prune table that is generated from `kprune`. [defualt: None]",
    default=None,
    show_default=True,
    type=click.Path(exists=True)
)
@click.option(
    "-norm",
    "--normalize",
    help="Normalize expanded edges by contig length",
    is_flag=True,
    show_default=True,
    default=False,
)
@click.option(
    "-af",
    "--allelic-factor",
    "allelic_factor",
    metavar="INT",
    help="Factor of allelic weight.",
    default=-1,
    type=float,
    show_default=True
)
@click.option(
    "-caf",
    "--cross-allelic-factor",
    "cross_allelic_factor",
    metavar="INT",
    help="Factor of cross-allelic weight.",
    default=0,
    type=float,
    show_default=True
)
@click.option(
    "-as",
    "--allelic-similarity",
    "allelic_similarity",
    metavar="FLOAT",
    help="The similarity of allelic, which used to merge clusters",
    default=0.85,
    type=click.FloatRange(0.0, 1.0),
    show_default=True
)
@click.option(
    '-mao',
    '--min-allelic-overlap',
    "min_allelic_overlap",
    metavar="FLOAT",
    help="Minimum overlap ratio bewteen two group when merging different groups",
    default=0.3,
    type=click.FloatRange(0.0, 1.0),
    show_default=True
)
@click.option(
    '-inc',
    '--incremental',
    help='Use incremental partition algorithm',
    is_flag=True,
    default=False,
    show_default=True,
    hidden=True
)
@click.option(
    '-knm',
    '--kprune-norm-method',
    'kprune_norm_method',
    metavar="STR",
    help="Normalization method of contacts for kprune",
    default="auto",
    show_default=True ,
    type=click.Choice(["none", "cis", "cis_unique", "auto"])
)
# @click.option(
#     '-uc',
#     '--ultra-complex',
#     'ultra_complex',
#     metavar='FLOAT',
#     help='ultra-complex polyploid',
#     default=0,
#     type=click.FloatRange(0.0, 10.0),
#     show_default=True,
# )
@click.option(
    "--merge-cluster",
    "merge_cluster",
    metavar="PATH",
    help="Don't run hyperpartition, only merge this cluster into k group. "
    "[default: None]",
    default=None,
    show_default=True,
)
@click.option(
    '-fc',
    '--first-cluster',
    'first_cluster',
    metavar="PATH",
    help='Use the existing first cluster results to second round cluster. '
    '[default: None]',
    default=None,
    show_default=True
)
@click.option(
    '--disable-merge-in-first',
    '--disable-merge-first',
    help="""
    disable merge function in first round cluster results, if set the `n` that
     program select `n` group to second round cluster.
    """,
    default=False, 
    is_flag=True
)
@click.option(
    '--exclude-group-to-second',
    '--exclude-to-second',
    'exclude_group_to_second',
    metavar="STR",
    help='exclude several groups that do not run in second round cluster, comma seperate. 1-base.'
    '[default: None]',
    default=None,
    show_default=True
)
@click.option(
    '--exclude-from-first',
    '--exclude-group-from-first',
    'exclude_group_from_first',
    metavar="STR",
    help='exclude several groups in first rount cluster, do not retain it to any group, comma seperate. 1-base.'
    '[default: None]',
    default=None,
    show_default=True
)
@click.option(
    "-wl",
    "--whitelist",
    metavar="PATH",
    help="""
    Path to 1-column list file containing
    contigs to include in hypergraph to partition.  
    [default: None]
    """,
    default=None,
    show_default=True,
    type=click.Path(exists=True)
)
@click.option(
    "-bl",
    "--blacklist",
    metavar="PATH",
    help="""
    Path to 1-column list file containing
    contigs to exclude from hypergraph to partition.  
    [default: None]
    """,
    default=None,
    show_default=True,
    type=click.Path(exists=True)
)
@click.option(
    "-hcr-bed",
    "--hcr-bed",
    "hcr_bed",
    metavar="STR",
    default=None,
    help="Only retain high confidence regions to subsequence analysis by input bed.",
    show_default=True,
    type=click.Path(exists=True)
)
@click.option(
    "-hcr-invert",
    "--hcr-invert",
    "hcr_invert",
    is_flag=True,
    default=False,
    help="Invert table by hcr bed.",
    show_default=True,
)
@click.option(
    "-mc",
    "--min-contacts",
    "min_contacts",
    help="Minimum contacts of contigs",
    metavar="INT",
    type=int,
    default=25,
    show_default=True
)
@click.option(
    "-ml",
    "--min-length",
    "min_length",
    help="Minimum length of contigs",
    metavar="INT",
    type=int,
    default=10000,
    show_default=True
)
@click.option(
    "-r1",
    "--resolution1",
    metavar="FLOAT",
    help="Resolution of the first partition",
    type=click.FloatRange(-1.0, 10.0),
    default=-1.0,
    show_default=True
)
@click.option(
    "-r2",
    "--resolution2",
    metavar="FLOAT",
    help="Resolution of the second partition",
    type=click.FloatRange(-1.0, 10.0),
    default=1.0,
    show_default=True
)
@click.option(
    "-ir1",
    "--init-resolution1",
    "init_resolution1",
    metavar="FLOAT",
    help="Initial resolution of the automatic search resolution in first partition"
    ", only used when `--resolution1=-1`.",
    type=float,
    default=0.8,
    show_default=True
)
@click.option(
    "-ir2",
    "--init-resolution2",
    "init_resolution2",
    metavar="FLOAT",
    type=float,
    help="Initial resolution of the automatic search resolution in second partition"
    ", only used when `--resolution2=-1`.",
    default=0.8,
    show_default=True
)
@click.option(
    "-mw",
    "--min-weight",
    "min_weight",
    help="Minimum weight of graph",
    type=float,
    default=0.1,
    show_default=True
)
@click.option(
    "-mcw",
    "--min-cis-weight",
    "min_cis_weight",
    help="Minimum conitg's cis weight of graph",
    type=float,
    default=5.0,
    show_default=True
)
@click.option(
    "-q1",
    "--min-quality1",
    "min_quality1",
    help="Minimum quality of hyperedge in first cluster.",
    type=click.IntRange(0, 60),
    default=1,
    show_default=True, 
)
@click.option(
    "-q2",
    "--min-quality2",
    "min_quality2",
    help="Minimum quality of hyperedge in second cluster.",
    type=click.IntRange(0, 60),
    default=2,
    show_default=True,
)
@click.option(
    "-ms",
    "--min-scaffold-length",
    "min_scaffold_length",
    help="The minimum length of the output scaffolding.",
    type=float,
    default=5e5,
    show_default=True
)
@click.option(
    "--enable-misassembly-remove",
    "enable_misassembly_remove",
    help="ensable misassembly after second round partition.",
    is_flag=True,
    default=True,
    show_default=True
)
@click.option(
    "--is-recluster-contigs",
    "--recluster-contigs",
    is_flag=True,
    default=False,
    show_default=True,
    hidden=True
)
@click.option(
    "--threshold",
    metavar="FLOAT",
    help="Threshold of reweight.",
    type=float,
    default=0.01,
    show_default=True,
    hidden=True,
)
@click.option(
    "--max-round",
    "max_round",
    help="Maximize round of IRMM algorithms",
    metavar="INT",
    type=int,
    default=1,
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=1,
    metavar='INT',
    show_default=True,
)
# @click.option(
#     '-cs',
#     '--chunksize',
#     help='Chunk size of edges',
#     type=float,
#     default=None,
#     metavar='Float',
#     show_default=True
# )
def hyperpartition(hypergraph, 
                    contigsizes, 
                    output,
                    pairs,
                    porec,
                    contacts,
                    ultra_long,
                    ul_weight,
                    mode,
                    n,
                    fasta,
                    alleles_kmer_size,
                    alleles_window_size,
                    alleles_minimum_similarity,
                    alleles_diff_thres,
                    alleletable,
                    prunetable,
                    normalize,
                    allelic_factor,
                    cross_allelic_factor,
                    allelic_similarity,
                    min_allelic_overlap,
                    incremental,
                    kprune_norm_method,
                    # ultra_complex,
                    merge_cluster,
                    first_cluster,
                    disable_merge_in_first,
                    exclude_group_to_second,
                    exclude_group_from_first,
                    whitelist,
                    blacklist,
                    hcr_bed,
                    hcr_invert,
                    min_contacts,
                    min_length, 
                    resolution1,
                    resolution2,
                    init_resolution1,
                    init_resolution2,
                    min_weight,
                    min_cis_weight,
                    min_quality1,
                    min_quality2,
                    min_scaffold_length,
                    enable_misassembly_remove,
                    is_recluster_contigs,
                    threshold, 
                    max_round, 
                    threads,
                    # chunksize
                    ):
    """
    Separate contigs into groups based on hypergraph (a: hp).

        Hypergraph : Path of the hypergraph.

        Contig_sizes : Path of contig sizes.

        Output : Path of output clusters.

    """
    import msgspec 
    import re
    from .hypergraph import HyperExtractor, Extractor
    from .hyperpartition import HyperPartition
    from .algorithms.hypergraph import HyperEdges, merge_hyperedges
    from .utilities import read_chrom_sizes, is_file_changed
    
    assert not all([porec, pairs]), "confilct parameters, only support one type data"

    ultra_complex = None
    # mode = "basal" if mode == "haploid" else mode 
    n = re.split(":|x|\|", n) if n else None
   
    mode = None if mode == "None" else mode
    
    if mode == "phasing" and n is not None:
        if len(n) == 1 and incremental is False:
            if alleletable is None and prunetable is None:
                mode = "haploid"
            else:
                mode = "basal_withprune"
                logger.info(f"Mode was set to `phasing` mode, but the group number of second cluster was not be specified, set the mode to {mode}")
            
    if mode == "basal":
        incremental = False
        logger.info("Running hyperpartition with `basal(haploid)` mode.")
        if alleletable or prunetable:
            logger.warning("allelic information will not be used in basal mode")
            alleletable = None
            prunetable = None  
    elif mode == "phasing":
        incremental = True 
        assert alleletable or prunetable or fasta,\
            "phasing mode must add `--fasta` or `-pt` or `-at` param"
        logger.info("Running hyperpartition with `phasing` mode.")

    elif mode == "basal_withprune":
        incremental = False
        assert alleletable or prunetable, \
            "basal_withprune modemust add `-pt` or `-at` param"
        logger.info("Running hyperpartition with `basal_withprune` mode.")
    else:
        incremental = incremental 
        if incremental is False: 
            if (alleletable is not None or prunetable is not None):
                if n:
                    if len(n) >= 2:
                        logger.info(f"You have specified two group numbers of `{':'.join(n)}`, "
                                    "the mode was set to `phasing` to run two-layer cluster.")
                        incremental = True
                    else:
                        logger.warning("Mode not be specified, running `basal_withprune` mode, "
                                "if you want to phase the diploid or polyploid, please "
                                "set the mode to `phasing`, or add `-inc` parameter.")
                else:
                    logger.warning("Mode not be specified, running basal_withprune mode, "
                                "if you want to phase the diploid or polyploid, please "
                                "set the mode to `phasing`, or add `-inc` parameter.")
            else:
                logger.warning("Mode not be specified, running `basal(haploid)` mode")
        else:
            logger.info("`-inc` be specified, will run in `phasing` mode")

    if kprune_norm_method == "auto":
        if pairs:
            kprune_norm_method = "cis"
        elif porec:
            kprune_norm_method = "cis"
        else:
            kprune_norm_method = "cis"
    
    if n is not None:
        # n = re.split(":|x|\|", n) 
        if len(n) <= 1 and incremental is True:
            logger.warning("Second round partition will not be run, if you want to run second round partition the `-inc` parameters must be added")
    else:
        n = [None, None]

    if (n[0] is None or n[0] == '0') and resolution1 < 0:
        logger.info(f"The group number of first cluster was not be specified, set the resolution1 from {resolution1} to 1.0.")
        resolution1 = 1.0 

    if (len(n) == 2):
        if (n[1] is None or n[1] == '0') and resolution2 < 0:
            logger.info(f"The group number of second cluster was not be specified, set the resolution1 from {resolution2} to 1.0.")
            resolution2 = 1.0

    for i, v in enumerate(n):
        if v:
            if i == 1:
                try:
                    n[i] = int(v)
                except ValueError:
                    if Path(n[i]).exists():
                        logger.info(f"Load 2nd group number list f`{n[i]}`")
                        tmp_df = pd.read_csv(n[i], sep='\s+', header=None, index_col=None)

                        n[i] = tmp_df.to_dict()[1]

                    else:
                        assert ValueError, "n2 must be a int or a exists file"
            else:
                n[i] = int(v)
        else:
            n[i] = 0

    contigsizes = read_chrom_sizes(contigsizes)
    if contigsizes.empty:
        logger.error("The contigsizes file is empty. please input correct file")
    prefix = Path(hypergraph).stem
    if porec:
        contigs = contigsizes.index.values.tolist()
        contigs = list(map(str, contigs))

        # contigs = natsorted(contigs)
        contig_idx = defaultdict(None, dict(zip(contigs, range(len(contigs)))))
        if is_file_changed(hcr_bed) or is_file_changed(hypergraph) or not Path(f"{prefix}.q{min_quality1}.hg").exists():
            logger.info(f"Load raw hypergraph from porec table `{hypergraph}`")
            
            he = HyperExtractor(hypergraph, contig_idx, contigsizes.to_dict()['length'], 
                                min_quality=min_quality1, hcr_bed=hcr_bed, 
                                hcr_invert=hcr_invert, threads=threads)
            he.save(f"{prefix}.q{min_quality1}.hg")
            hypergraph = he.edges
        else:
            if hcr_bed:
                logger.warning(f"Load raw hypergraph from existed file of `{prefix}.q{min_quality1}.hg`, if the {hcr_bed} changed, you should remove this existing hg.")
            else:
                logger.warning(f"Load raw hypergraph from existed file of `{prefix}.q{min_quality1}.hg`.")
            hypergraph = msgspec.msgpack.decode(open(f"{prefix}.q{min_quality1}.hg", 'rb').read(), type=HyperEdges)

    elif pairs:
        contigs = contigsizes.index.values.tolist()
        contigs = list(map(str, contigs))
        # contigs = natsorted(contigs)
        contig_idx = defaultdict(None, dict(zip(contigs, range(len(contigs)))))
        if is_file_changed(hcr_bed) or is_file_changed(hypergraph) or not Path(f"{prefix}.q{min_quality1}.hg").exists():
            logger.info(f"Load raw hypergraph from pairs file `{hypergraph}`")
            he = Extractor(hypergraph, contig_idx, contigsizes.to_dict()['length'], 
                           min_quality=min_quality1, hcr_bed=hcr_bed, 
                           hcr_invert=hcr_invert, threads=threads)
            he.save(f"{prefix}.q{min_quality1}.hg")
            hypergraph = he.edges
        else:
            logger.warning(f"Load raw hypergraph from exists file of `{prefix}.q{min_quality1}.hg`, if the input pairs changed, you should remove this existing hg.")
            hypergraph = msgspec.msgpack.decode(open(f"{prefix}.q{min_quality1}.hg", 'rb').read(), type=HyperEdges)
        
        
    else:
        logger.info(f"Load raw hypergraph from `{hypergraph}")
        try:
            hypergraph = msgspec.msgpack.decode(open(hypergraph, 'rb').read(), type=HyperEdges)
        except msgspec.ValidationError:
            raise msgspec.ValidationError(f"`{hypergraph}` is not the hypergraph. Please input correct hypergraph format")
            
    if ultra_long and ul_weight:
        logger.info("Load raw hypergraph from ultra-long hypergraph: `{ultra_long}`")
        try:
            ultra_long_hypergraph = msgspec.msgpack.decode(open(ultra_long, 'rb').read(), type=HyperEdges)
        except msgspec.ValidationError:
            raise msgspec.ValidationError(f"`{ultra_long}` is not the hypergraph. Please input correct hypergraph format")
        
        HE_list = [hypergraph]
        for i in range(ul_weight):
            HE_list.append(ultra_long_hypergraph)
   
        hypergraph = merge_hyperedges(HE_list)
    
    hypergraph.to_numpy()
    


    if whitelist:
        whitelist_file = whitelist 
        whitelist = [i.strip().split()[0] for i in open(whitelist) if i.strip()]
        logger.info(f"Only `{len(whitelist)}` contigs will be load to cluster, which specified in `{whitelist_file}`.")
    if blacklist:
        blacklist_file = blacklist
        blacklist = [i.strip().split()[0] for i in open(blacklist) if i.strip()]

        logger.info(f"`{len(blacklist)}` contigs will be removed to cluster, which specified in `{blacklist_file}`..")

    if merge_cluster:
        assert n[0] != 0, "parameter `k` must add to run merge cluster"
        logger.info("Only run the algorithm of clusters merging.")
        _merge_cluster = ClusterTable(merge_cluster)
        if whitelist:
            whitelist = list(set(whitelist) & set(_merge_cluster.contigs))
        else:
            whitelist = _merge_cluster.contigs


    assert whitelist is None or blacklist is None, \
        "Only support one list of whitelist or blacklist"
    
    if exclude_group_to_second:
        exclude_group_to_second = exclude_group_to_second.split(",")
        try:
            exclude_group_to_second = list(map(int, exclude_group_to_second))
        except ValueError:
            logger.warning("The exclude groups must be several numbers with comma seperated")

    if exclude_group_from_first:
        exclude_group_from_first = exclude_group_from_first.split(",")
        try:
            exclude_group_from_first = list(map(int, exclude_group_from_first))
        except ValueError:
            logger.warning("The exclude groups must be several numbers with comma seperated")

    if alleletable:
        if prunetable:
            logger.warning("The allele table will not be used, becauese prunetable parameter added.")

    if incremental is False and first_cluster is not None:
        logger.warning("First cluster only support for incremental method, will be not used.")
    
    if not prunetable and not alleletable and not fasta:
        logger.info("Not inplement the allelic and cross-allelic reweight algorithm")
    
    is_remove_misassembly = False if enable_misassembly_remove else True
    mapq_filter_1_to_2 = True if first_cluster else False
    hp = HyperPartition(hypergraph, 
                            contigsizes,
                            None,
                            ul_weight,
                            n,
                            fasta,
                            alleles_kmer_size,
                            alleles_window_size,
                            alleles_minimum_similarity,
                            alleles_diff_thres,
                            alleletable,
                            prunetable,
                            normalize,
                            contacts,
                            kprune_norm_method,
                            allelic_factor,
                            cross_allelic_factor,
                            allelic_similarity,
                            min_allelic_overlap,
                            ultra_complex,
                            disable_merge_in_first,
                            exclude_group_to_second,
                            exclude_group_from_first,
                            whitelist,
                            blacklist,
                            min_contacts,
                            min_length,
                            resolution1, 
                            resolution2,
                            init_resolution1,
                            init_resolution2,
                            min_weight,
                            min_cis_weight,
                            min_quality1,
                            min_quality2,
                            mapq_filter_1_to_2,
                            min_scaffold_length,
                            is_remove_misassembly,
                            is_recluster_contigs,
                            threshold, 
                            max_round, 
                            threads, 
                            # chunksize
                            )
    
  
    if incremental:
        if first_cluster and op.exists(first_cluster):
            first_cluster = ClusterTable(first_cluster)
            if n[0] != len(first_cluster.groups) and n[0]:
                logger.warning("The group number of first cluster is conflicted with the specified group number specified in `-n`.")
        else:
            first_cluster = None 
        
        if merge_cluster:
            logger.error("Imcomplete function of merge cluster in phasing mode.")
            sys.exit(-1)
            hp.phased_merge_cluster(_merge_cluster)
        else:
            hp.incremental_partition(n, first_cluster)
    else:
        if merge_cluster:
            groups = list(_merge_cluster.data.values())
            hp.merge_cluster(groups)
        else:
            hp.single_partition(int(n[0]))
    
    hp.to_cluster(output)

@cli.command(cls=RichCommand, hidden=True, epilog=__epilog__)
@click.argument(
    "hypergraph",
    metavar="HyperGraph",
    type=click.Path(exists=True)
)
@click.argument(
    "contigsizes",
    metavar="ContigSizes",
    type=click.Path(exists=True),
)
@click.argument(
    "clustertable",
    metavar="ClusterTable",
    type=click.Path(exists=True),
)
@click.argument(
    "collapsed_contigs",
    metavar="Collapsed_Contigs_Table",
    type=click.Path(exists=True),
)
@click.option(
    "-n",
    "--ploidy",
    help="ploidy level of genome.",
    metavar="INT",
    default=2,
    type=int
)
@click.option(
    "-at",
    "--alleletable",
    metavar="PATH",
    help="""
    Path to allele table that is generated from `alleles`. 
        Skip when prunetable parameter added. [default: None]
    """,
    default=None,
    show_default=True,
    type=click.Path(exists=True),
)
@click.option(
    "-as",
    "--allelic-similarity",
    "allelic_similarity",
    metavar="FLOAT",
    help="The similarity of allelic, which used to merge clusters",
    default=0.85,
    type=click.FloatRange(0.0, 1.0),
    show_default=True
)
def collapsed_rescue(hypergraph, contigsizes, clustertable, 
                    collapsed_contigs, ploidy,
                    alleletable, allelic_similarity):
    import msgspec
    from .collapse import CollapsedRescue
    from .core import AlleleTable, ClusterTable 
    from .algorithms.hypergraph import HyperEdges, HyperGraph
    
    hyperedge = msgspec.msgpack.decode(open(hypergraph, 'rb').read(), type=HyperEdges)
    hyperedge.to_numpy()
    hypergraph = HyperGraph(hyperedge, min_quality=2)
    
    at = AlleleTable(alleletable, sort=False, fmt="allele2")
    ct = ClusterTable(clustertable)
    collapsed_contigs = pd.read_csv(collapsed_contigs, sep='\t', index_col=0, header=None,
                                        names=["contig", "depth", "CN"])

    collapsed_contigs = collapsed_contigs[collapsed_contigs["CN"] <= (ploidy)]

    cr = CollapsedRescue(
        HG=hypergraph,
        clustertable=ct,    
        alleletable=at,
        collapsed_contigs=collapsed_contigs,
        allelic_similarity=allelic_similarity,
    )

    cr.rescue()


@cli.command(cls=RichCommand, hidden=True, epilog=__epilog__)
@click.argument(
    "hypergraph",
    metavar="HyperGraph",
    type=click.Path(exists=True)
)
@click.argument(
    "fasta",
    metavar="Fasta",
    type=click.Path(exists=True),
)
@click.argument(
    "agp_file",
    metavar="AGP",
    type=click.Path(exists=True),
)
@click.argument(
    "collapsed_contigs",
    metavar="Collapsed_Contigs_Table",
    type=click.Path(exists=True),
)
@click.option(
    "-n",
    "--ploidy",
    help="ploidy level of genome.",
    metavar="INT",
    default=2,
    type=int
)
@click.option(
    "-at",
    "--alleletable",
    metavar="PATH",
    help="""
    Path to allele table that is generated from `alleles`. 
        Skip when prunetable parameter added. [default: None]
    """,
    default=None,
    show_default=True,
    type=click.Path(exists=True),
)
@click.option(
    "-c",
    "-s",
    "-sc",
    "--split-contacts",
    "split_contacts",
    metavar="SPLIT_CONTACTS",
    help="Split contacts from `cphasing-rs pairs2clm`, only used in `precision` or `fast` mode.",
    type=click.Path(exists=True),
    default=None,
    show_default=True
)
@click.option(
    "-as",
    "--allelic-similarity",
    "allelic_similarity",
    metavar="FLOAT",
    help="The similarity of allelic, which used to merge clusters",
    default=0.85,
    type=click.FloatRange(0.0, 1.0),
    show_default=True
)
def collapsed_rescue2(hypergraph, fasta, agp_file, 
                    collapsed_contigs, ploidy,
                    alleletable, split_contacts, allelic_similarity):
    import msgspec
    from .agp import import_agp
    from .collapse import CollapsedRescue2
    from .core import AlleleTable, ClusterTable 
    from .algorithms.hypergraph import HyperEdges, HyperGraph
    from .utilities import read_chrom_sizes
    hyperedge = msgspec.msgpack.decode(open(hypergraph, 'rb').read(), type=HyperEdges)
    hyperedge.to_numpy()
    hypergraph = HyperGraph(hyperedge, min_quality=2)
    
    at = AlleleTable(alleletable, sort=False, fmt="allele2")
    
    
    collapsed_contigs = pd.read_csv(collapsed_contigs, sep='\t', index_col=0, header=None,
                                        names=["contig", "depth", "CN"])

    collapsed_contigs = collapsed_contigs[collapsed_contigs["CN"] <= (ploidy)]

    cr = CollapsedRescue2(
        HG=hypergraph,
        agp=agp_file,    
        fasta=fasta,
        alleletable=at,
        split_contacts=split_contacts,
        collapsed_contigs=collapsed_contigs,
        allelic_similarity=allelic_similarity,
    )

    cr.run()

@cli.command(cls=RichCommand, epilog=__epilog__)
@click.argument(
    "clustertable",
    metavar="ClusterTable",
    type=click.Path(exists=True),
)
@click.argument(
    "count_re", 
    metavar="Count_RE",
    type=click.Path(exists=True),
)
@click.argument(
    "clm",
    metavar="CLM",
    type=click.Path(exists=True),
)
@click.option(
    "-c",
    "-s",
    "-sc",
    "--split-contacts",
    "split_contacts",
    metavar="SPLIT_CONTACTS",
    help="Split contacts from `cphasing-rs pairs2clm`, only used in `precision` or `fast` mode.",
    type=click.Path(exists=True),
    default=None,
    show_default=True
)
@click.option(
    "-a",
    "-at",
    "--allele-table",
    "allele_table",
    metavar="AlleleTable",
    help="Input allele table from `cphasing alleles` to adjust the haplotype to parallel",
    default=None,
    type=click.Path(exists=True)
)
@click.option(
    "-f",
    "--fasta",
    metavar="Fasta",
    help="Input contig-level fasta to build the assembly result."
    "If None, tour files will be generated.",
    default=None,
    type=click.Path(exists=True)
)
@click.option(
    '--corrected',
    is_flag=True,
    default=False,
    help="Contigs were corrected, used for building agp and fasta."
)
@click.option(
    "-m",
    "--method",
    "-scaf-method",
    "--scaf-method",
    "--scaffolding-method",
    "method",
    metavar="STR",
    help="The method of scaffolding, `['precision', 'allhic', 'fast']`."
    "precision: haphic_fastsort + allhic, which will quicker than allhic only. "
    "It is worth noting that `allhic` in `C-Phasing` is parameterized to "
    "achieve better results than the previous version",
    default="precision",
    type=click.Choice(["precision", "allhic", "fast"]),
    show_default=True
)
@click.option(
    "--keep-temp",
    "keep_temp",
    help="Dont not delete the tempdir.",
    default=False,
    is_flag=True,
    show_default=True
)
@click.option(
    "-o",
    "-oa",
    "--output",
    metavar="STR",
    help="Output agp file, only valid when fasta file exists.",
    default="groups.agp",
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
# @click.argument(
#     'coolfile',
#     metavar="INPUT_COOL_PATH",
#     type=click.Path(exists=True)
# )
# @click.argument(
#     'output',
#     metavar="OUTPUT_SCORE_PATH",
# )
def scaffolding(clustertable, count_re, clm, 
                split_contacts, 
                allele_table,
                fasta,
                corrected,
                method,
                keep_temp,
                output, threads):
    """
    Ordering and orientation the contigs (a: sf).

        ClusterTable: Path of cluster table from hyperpartition.

        Count_RE: Path of counts RE file.

        CLM: Path of clm file.

    """

    from .algorithms.scaffolding import AllhicOptimize, HapHiCSort
    
    if split_contacts is None and (method == "precision" or method == "fast"):
        logger.warning("The split contacts not be specified, change the method to `allhic`.")
        method = "allhic"

    if fasta is None:
        logger.warning("Will not generate the agp file, because the fasta not be specified")                               

    if method == "allhic":
        ao = AllhicOptimize(clustertable, count_re, clm, allele_table=allele_table, corrected=corrected,
                                fasta=fasta, keep_temp=keep_temp, output=output, threads=threads)
        ao.run()
    elif method == "precision":
        assert split_contacts is not None, "split_contacts file must specified by `-sc` parameters"
        hs = HapHiCSort(clustertable, count_re, clm, split_contacts, corrected=corrected,
                                allele_table=allele_table, keep_temp=keep_temp,
                                fasta=fasta, output=output, threads=threads)
        hs.run()

    else:
        assert split_contacts is not None, "split_contacts file must specified by `-sc` parameters"
        hs = HapHiCSort(clustertable, count_re, clm, split_contacts, skip_allhic=True, corrected=corrected,
                        allele_table=allele_table, fasta=fasta, output=output, threads=threads)
        hs.run()

@cli.command(hidden=True)
@click.argument(
    "hypergraph",
    metavar="HyperGraph",
    type=click.Path(exists=True)
)
def hyperoptimize(hypergraph):
    import msgspec
    from cphasing.algorithms.hypergraph import HyperEdges, HyperGraph
    from cphasing.algorithms.scaffolding import HyperOptimize
    he = msgspec.msgpack.decode(open(hypergraph,  'rb').read(), type=HyperEdges)
    HG = HyperGraph(he)
    ho = HyperOptimize(HG, 2)
    order = np.arange(len(HG.nodes))
    order2 = ho.shuffle()
    print(ho.fitness(order))
    
    ho.evaluate()


@cli.command(cls=RichCommand)
@click.argument(
    "fasta",
    metavar="Fasta",
    type=click.Path(exists=True)
)
@click.option(
    '--corrected',
    is_flag=True,
    help="Is corrected contigs in tours, which used when fasta and tour inconsistent",
    default=False, 
)
@click.option(
    "-o",
    "--output",
    metavar="STR",
    help="Output genome fasta, '.gz' supported.",
    default="groups.asm.fasta",
    show_default=True
)
@click.option(
    "-oa",
    "--output-agp",
    "output_agp",
    metavar="STR",
    help="Output agp file",
    default="groups.agp",
    show_default=True
)
@click.option(
    "--only-agp",
    "only_agp",
    help="Only output the agp file.",
    is_flag=True,
    default=False,
    show_default=True
)
def build(fasta, corrected, output, output_agp, only_agp):
    """
    Build genome release from tour files.

    Fasta : contig-level fasta file
    """
    from .build import Build
    Build(fasta, output,  corrected=corrected,
            output_agp=output_agp, only_agp=only_agp)


@cli.command(cls=RichCommand, short_help="Rename and orient the groups according to a reference.")
@click.option(
    '-r',
    '--ref',
    metavar="REFERENCE FASTA",
    help="Genome of monoploid or closely related species.",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '-f',
    '--fasta',
    metavar="DRAFT FASTA",
    help="Path to draft assembly",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '-a',
    '--agp',
    help="agp file.",
    metavar='AGP',
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '--unphased',
    help="Haplotypes are not parallel align.",
    is_flag=True,
    default=False,
)
@click.option(
    '-o',
    '--output',
    metavar="OUTPUT_AGP",
    default="groups.renamed.agp",
    help="Output path of renamed agp file"
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
def rename(ref, fasta, agp, 
                unphased, output, threads):
    """
    Rename and orientation the groups according to a refernce.

        To speed up this function, we only align the first haplotype to reference, 
            because the orientation among different haplotypes has been paralleled in `scaffolding`. 
        If you want to orient all of the haplotypes you can specify the `--unphased` parameters.
    """
    
    from .algorithms.scaffolding import Rename 

    r = Rename(ref=ref,
                    fasta=fasta,
                    agp=agp,
                    hap_aligned=False if unphased else True,
                    output=output,
                    threads=threads)
    r.run()

@cli.command(cls=RichCommand, hidden=HIDDEN, short_help="Calculate the size of all contigs. (hidden)")
@click.argument(
    "fasta",
    metavar="INPUT_FASTA_PATH",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    metavar="OUTPUT",
    help="output file path.",
    default="-",
    show_default=True
)
def chromsizes(fasta, output):
    cmd = ["cphasing-rs", "chromsizes", f"{fasta}", "-o", f"{output}"]
    flag = run_cmd(cmd)
    assert flag == 0, "Failed to execute command, please check log."
    

@cli.command(hidden=HIDDEN, short_help="Convert paf to pairs. (hidden)")
@click.argument(
    "paf",
    metavar="INPUT_PAF_PATH",
    type=click.Path(exists=True)
)
@click.argument(
    "contigsizes",
    metavar="CONTIGSIZES",
    type=click.Path(exists=True)
)
@click.option(
    "-b",
    "--bed",
    default=None,
    metavar="STR",
    help="bed file, which only retain the fragment start and end located in it.",
    type=click.Path(exists=True),
    show_default=True,
)
@click.option(
    "-q",
    "--min-mapq",
    "min_mapq",
    default=1,
    metavar="INT",
    help="Minimum mapping quality of alignments",
    type=click.IntRange(0, 60),
    show_default=True,
)
@click.option(
    "-p",
    "--min-identity",
    "min_identity",
    default=0.75,
    metavar="FLOAT",
    help="Minimum identity of alignments",
    type=click.FloatRange(0, 1.0),
    show_default=True,
)
@click.option(
    "-l",
    "--min-length",
    "min_length",
    default=30,
    metavar="FLOAT",
    help="Minimum length of alignments",
    type=int,
    show_default=True,
)
@click.option(
    '-e',
    '--max-edge',
    'max_edge',
    help='Maximum length of fragment located in the edge of contigs.',
    metavar='INT',
    default=0,
    show_default=True
)
@click.option(
    "--min-order",
    "min_order",
    metavar="INT",
    help="Minimum order of porec reads.",
    default=2,
    type=int,
    show_default=True
)
@click.option(
    "--max-order",
    "max_order",
    metavar="INT",
    help="Maximum order of porec reads.",
    default=50,
    type=int,
    show_default=True
)
@click.option(
    "-o",
    "--output",
    metavar="OUTPUT",
    default="-",
    show_default=True
)
def paf2pairs(paf, contigsizes, 
              bed, min_mapq, 
              min_identity, min_length, 
              max_edge,
              min_order, max_order, output):
    if not bed:
        cmd = ["cphasing-rs", "paf2pairs", paf, 
             "-q", str(min_mapq), "-p", str(min_identity), 
            "-l", str(min_length), "-m", str(min_order), 
            "-e", str(max_edge),
            "-M", str(max_order), contigsizes, "-o", output]
    else:
        cmd = ["cphasing-rs", "paf2pairs", paf, 
            "-b", bed, "-q", str(min_mapq), "-p", str(min_identity), 
            "-e", str(max_edge),
            "-l", str(min_length), "-m", str(min_order), 
            "-M", str(max_order), contigsizes, "-o", output]


    flag = run_cmd(cmd)
    assert flag == 0, "Failed to execute command, please check log."


@cli.command(cls=RichCommand, hidden=HIDDEN, short_help="Convert paf to porec table. (hidden)")
@click.argument(
    "paf",
    metavar="INPUT_PAF_PATH",
    type=click.Path(exists=True)
)
@click.option(
    "-b",
    "--bed",
    default=None,
    metavar="STR",
    help="bed file, which only retain the fragment start and end located in it.",
    type=click.Path(exists=True),
    show_default=True,
)
@click.option(
    "-q",
    "--min-mapq",
    "min_mapq",
    default=1,
    metavar="INT",
    help="Minimum mapping quality of alignments",
    type=click.IntRange(0, 60),
    show_default=True,
)
@click.option(
    "-p",
    "--min-identity",
    "min_identity",
    default=0.75,
    metavar="FLOAT",
    help="Minimum identity of alignments",
    type=click.FloatRange(0, 1.0),
    show_default=True,
)
@click.option(
    "-l",
    "--min-length",
    "min_length",
    default=30,
    metavar="FLOAT",
    help="Minimum length of alignments",
    type=int,
    show_default=True,
)
@click.option(
    '-e',
    '--max-edge',
    'max_edge',
    help='Maximum length of fragment located in the edge of contigs.',
    metavar='INT',
    default=2000,
    show_default=True
)
@click.option(
    "--max-order",
    "max_order",
    metavar="INT",
    help="Maximum order of porec reads.",
    default=50,
    type=int,
    show_default=True
)
@click.option(
    "-o",
    "--output",
    metavar="OUTPUT",
    help="output path",
    default="-",
    show_default=True
)
def paf2porec(paf, bed, min_quality, 
              min_identity, min_length, 
              max_edge,
              max_order, output):
    if not bed:
        cmd = ["cphasing-rs", "paf2porec",  "-p", str(min_identity), "-M", str(max_order),
                "-q", str(min_quality), "-l", str(min_length),  paf, 
                "-e", str(max_edge), "-o", output]
    else:
        cmd = ["cphasing-rs", "paf2porec",  "-p", str(min_identity), "-b",  bed, 
               "-M", str(max_order), "-q", str(min_quality), "-l", str(min_length), 
               "-e", str(max_edge), "-o", output]

    flag = run_cmd(cmd)
    assert flag == 0, "Failed to execute command, please check log."


@cli.command(cls=RichCommand, hidden=HIDDEN, short_help="Convert porec to pairs. (hidden)")
@click.argument(
    "porec",
    metavar="INPUT_POREC_PATH",
    type=click.Path(exists=True)
)
@click.argument(
    "contigsizes",
    metavar="CONTIGSIZES",
    type=click.Path(exists=True)
)
@click.option(
    "-q",
    "--min-mapq",
    "min_mapq",
    default=1,
    metavar="INT",
    help="Minimum mapping quality of alignments",
    type=click.IntRange(0, 60),
    show_default=True,
)
@click.option(
    "--min-order",
    "min_order",
    metavar="INT",
    help="Minimum order of porec reads.",
    default=2,
    type=int,
    show_default=True
)
@click.option(
    "--max-order",
    "max_order",
    metavar="INT",
    help="Maximum order of porec reads.",
    default=50,
    type=int,
    show_default=True
)
@click.option(
    "-o",
    "--output",
    metavar="OUTPUT",
    default="-",
    show_default=True
)
def porec2pairs(porec, min_mapq, 
                min_order, max_order,
                contigsizes, output):
    cmd = ["cphasing-rs", "porec2pairs", porec, contigsizes, "-q", str(min_mapq),
           "--min-order", str(min_order), "--max-order", str(max_order), "-o", output]

    flag = run_cmd(cmd)
    assert flag == 0, "Failed to execute command, please check log."


@cli.command(cls=RichCommand, hidden=HIDDEN, short_help="Convert pairs to bam file. (hidden)")
@click.argument(
    "bam",
    metavar="INPUT_BAM_PATH",
    type=click.Path(exists=True)
)
@click.option(
    '-q',
    '--min_quality',
    help='Minimum quality of mapping [0, 255].',
    metavar='INT',
    type=click.IntRange(0, 255, clamp=True),
    default=1,
    show_default=True
)
@click.option(
    "-o",
    "--output",
    metavar="OUTPUT",
    default="-",
    show_default=True
)
def bam2pairs(bam, min_quality, output):
    cmd = ["cphasing-rs", "bam2pairs", f"{bam}",  
           
            "-q", f"{min_quality}", "-o", f"{output}"]
    flag = run_cmd(cmd)
    assert flag == 0, "Failed to execute command, please check log."

@cli.command(cls=RichCommand, hidden=HIDDEN, short_help="Convert pairs to bam file. (hidden)")
@click.argument(
    "pairs",
    metavar="INPUT_PAIRS_PATH",
    type=click.Path(exists=True)
)
@click.option(
    '-q',
    '--min_quality',
    help='Minimum quality of mapping [0, 60].',
    metavar='INT',
    type=click.IntRange(0, 60, clamp=True),
    default=1,
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads. (unused)',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
@click.option(
    "-o",
    "--output",
    metavar="OUTPUT",
    default="-",
    show_default=True
)
def pairs_filter(pairs, min_quality, threads, output):
    cmd = ["cphasing-rs", "pairs-filter", f"{pairs}", "-t", f"{threads}",
            "-q", f"{min_quality}", "-o", f"{output}"]
    flag = run_cmd(cmd)
    assert flag == 0, "Failed to execute command, please check log."




@cli.command(cls=RichCommand, hidden=HIDDEN, short_help="Convert pairs to clm. (hidden)")
@click.argument(
    "pairs",
    metavar="INPUT_PAIRS_PATH",
    type=click.Path(exists=True)
)
@click.option(
    '--min-contacts',
    'min_contacts',
    help='Minimum contacts for contig pair',
    default=1,
    show_default=True,
    type=int,
)
@click.option(
    '-q',
    '--min_quality',
    help='Minimum quality of mapping [0, 255].',
    metavar='INT',
    type=click.IntRange(0, 255, clamp=True),
    default=1,
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads. (unused)',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
@click.option(
    "-o",
    "--output",
    metavar="OUTPUT",
    default="-",
    show_default=True
)
def pairs2clm(pairs, min_contacts, min_quality, threads, output):
    cmd = ["cphasing-rs", "pairs2clm", f"{pairs}", "-c", f"{min_contacts}", 
            "-q", f"{min_quality}", "-t", f"{threads}", "-o", f"{output}"]
    flag = run_cmd(cmd)
    assert flag == 0, "Failed to execute command, please check log."

@cli.command(cls=RichCommand, hidden=HIDDEN, short_help="Convert pairs to clm. (hidden)")
@click.argument(
    "pairs",
    metavar="INPUT_PAIRS_PATH",
    type=click.Path(exists=True)
)
@click.option(
    '-b',
    '--binsize',
    'binsize',
    help='Binsize of the depth calculation',
    default=10000,
    show_default=True,
    type=int,
)
@click.option(
    '-q',
    '--min_quality',
    help='Minimum quality of mapping [0, 255].',
    metavar='INT',
    type=click.IntRange(0, 255, clamp=True),
    default=0,
    show_default=True
)
@click.option(
    "-o",
    "--output",
    metavar="OUTPUT",
    default="-",
    show_default=True
)
def pairs2depth(pairs, binsize, min_quality, output):
    cmd = ["cphasing-rs", "pairs2depth", f"{pairs}", "-b", f"{binsize}", 
            "-q", f"{min_quality}", "-o", f"{output}"]
    flag = run_cmd(cmd)
    assert flag == 0, "Failed to execute command, please check log."

@cli.command(hidden=True, short_help="Convert pairs to contacts. (hidden)")
@click.argument(
    "pairs",
    metavar="INPUT_PAIRS_PATH",
    type=click.Path(exists=True)
)
@click.option(
    '--min-contacts',
    'min_contacts',
    help='Minimum contacts for contig pair',
    metavar="INT",
    default=1,
    show_default=True,
    type=int,
)
@click.option(
    '-n',
    '-s',
    '--split-num',
    'split_num',
    metavar="INT",
    help="Split contigs into several bins",
    type=int,
    default=1,
    show_default=True,
)
@click.option(
    '-q',
    '--min_quality',
    help='Minimum quality of mapping [0, 255].',
    metavar='INT',
    type=click.IntRange(0, 255, clamp=True),
    default=1,
    show_default=True
)
@click.option(
    "-o",
    "--output",
    metavar="OUTPUT",
    default="-",
    show_default=True
)
def pairs2contacts(pairs, min_contacts, split_num, 
                   min_quality, output):
    cmd = ["cphasing-rs", "pairs2contacts", f"{pairs}", 
           "-c", f"{min_contacts}", "-q", f"{min_quality}",
            "-n", f"{split_num}", "-o", f"{output}"]
    flag = run_cmd(cmd)
    assert flag == 0, "Failed to execute command, please check log."


@cli.command(cls=RichCommand, hidden=HIDDEN, short_help="Convert pairs to mnd file. (hidden)")
@click.argument(
    "pairs",
    metavar="INPUT_PAIRS_PATH",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    metavar="OUTPUT",
    default="-",
    show_default=True
)
@click.option(
    "-q",
    "--min-mapq",
    "min_mapq",
    default=1,
    metavar="INT",
    help="Minimum mapping quality of alignments",
    type=click.IntRange(0, 60),
    show_default=True,
)
def pairs2mnd(pairs, output, min_mapq):

    cmd = ["cphasing-rs", "pairs2mnd", pairs, "-o", output, "-q", str(min_mapq)]
    flag = run_cmd(cmd)
    assert flag == 0, "Failed to execute command, please check log."


@cli.command(cls=RichCommand, epilog=__epilog__)
@click.argument(
    "pairs",
    metavar="INPUT_PAIRS_PATH",
    type=click.Path(exists=True)
)
@click.argument(
    "chromsize",
    metavar="CHROM_SIZE",
    type=click.Path(exists=True)
)
@click.argument(
    "outcool",
    metavar="OUT_COOL_PATH"
)
@click.option(
    "-bs",
    "--binsize",
    help="Bin size in bp. Enabled with suffix of [k, m]",
    type=str,
    default="10000",
    show_default=True
)
@click.option(
    "-q",
    "--min-mapq",
    "min_mapq",
    default=1,
    metavar="INT",
    help="Minimum mapping quality of alignments",
    type=click.IntRange(0, 60),
    show_default=True,
)
@click.option(
    '--fofn',
    help="""
    If this flag is set then the pairs is a file of
    filenames corresponding to the 4DN pairs you want to merge.
    This is workaround for when the command line gets too long.
    """,
    is_flag=True,
    default=False,
    show_default=True,
    hidden=True
)
@click.option(
    '--low-memory',
    help="Reduce memory usage.",
    is_flag=True,
    default=False,
)
@click.option(
    '-cs',
    '--chunksize',
    help='Chunk size of loading pairs.',
    type=float,
    default=100000000,
    metavar='Float',
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=10,
    metavar='INT',
    show_default=True,
)
def pairs2cool(pairs, chromsize, outcool,
               binsize, min_mapq, fofn, 
               low_memory, chunksize, threads):
    """
    Convert pairs file into a specified resolution cool file.

        INPUT_PAIRS_PATH : Path of pairs file, can be compressed.

        CHROM_SIZE : Two columns of chromosomes or contigs size.

        OUT_COOL_PATH : Output path of cool file.
    """
    from subprocess import PIPE, Popen
    from cooler.cli.cload import pairs as cload_pairs
    from .core import Pairs2
    

    Path("logs").mkdir(exist_ok=True)

    binsize = humanized2numeric(binsize)
    logger.info(f"Load pairs: `{pairs}`.")
    logger.info(f"Matrix's bin size: {to_humanized2(binsize)}")
    
    if not low_memory:
        available_memory = psutil.virtual_memory().available / 1024 / 1024 / 1024
        pairs_size = op.getsize(pairs) / 1024 / 1024 / 1024

        if pairs_size * 20 > available_memory:
            logger.warning(f"Memory usage is too high, use low memory mode.")
            low_memory = True

    if not low_memory:
        p = Pairs2(pairs, min_mapq=min_mapq, threads=threads, chunksize=chunksize)
        p.to_cool(outcool, binsize=binsize)
        return 

    if fofn:
        pairs_files = [i.strip() for i in open(pairs) if i.strip()]
        pid = os.getpid()
        if pairs_files[0].endswith(".gz"):
            os.system(f'zgrep "^#" {pairs_files[0]} > temp.{pid}.header')
            pairs_files = ' '.join(pairs_files)
            os.system(f'zcat {pairs_files} | grep -v "^#" > temp1.{pid}.pairs')
            os.system(f'cat temp.{pid}.header  temp1.{pid}.pairs > temp.{pid}.pairs')
            os.remove(f'temp1.{pid}.pairs')
        else:
            os.system(f'grep "^#" {pairs_files[0]} > temp.{pid}.header')
            pairs_files = ' '.join(pairs_files)
            os.system(f'cat {pairs_files} | grep -v "^#" > temp1.{pid}.pairs')
            os.system(f'cat temp.{pid}.header  temp1.{pid}.pairs > temp.{pid}.pairs')
            os.remove(f'temp1.{pid}.pairs')
        pairs = f'temp.{pid}.pairs'

    if min_mapq == 0:
        try:
            cload_pairs.main(args=[
                            f"{chromsize}:{binsize}",
                            pairs, 
                            outcool, 
                            "-c1", 2,
                            "-p1", 3,
                            "-c2", 4,
                            "-p2", 5],
                            prog_name='cload')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
    else:
        logger.info(f"Only load pair that minimum quality >= {min_mapq}.")
        cmd1 = ["cphasing-rs", "pairs-filter", 
                pairs, "-q", f"{min_mapq}", "2>", "logs/pairs-filter.log"]
        
        cmd2 = ["cooler", "cload", "pairs", 
                f"{chromsize}:{binsize}", 
                "-", outcool, 
                "-c1", "2", "-c2", "4",
                "-p1", "3", "-p2", "5", "2>", "logs/pairs2cool.log"]
        
        os.system(f"{' '.join(cmd1)} | {' '.join(cmd2)}")

        
    logger.info(f'Output binning contact matrix into `{outcool}`')
    if fofn:
        if op.exists(f'temp.{pid}.pairs'):
            os.remove(f'temp.{pid}.pairs')
        if op.exists(f'temp.{pid}.header'):
            os.remove(f'temp.{pid}.header')


@cli.command(cls=RichCommand, epilog=__epilog__, hidden=True)
@click.argument(
    "pairs",
    metavar="INPUT_PAIRS_PATH",
    type=click.Path(exists=True)
)
@click.argument(
    "outcool",
    metavar="OUT_COOL_PATH"
)
@click.option(
    "-bs",
    "--binsize",
    help="Bin size in bp. Enabled with suffix of [k, m]",
    type=str,
    default="10000",
    show_default=True
)
@click.option(
    "-q",
    "--min-mapq",
    "min_mapq",
    default=1,
    metavar="INT",
    help="Minimum mapping quality of alignments",
    type=click.IntRange(0, 60),
    show_default=True,
)
@click.option(
    '-cs',
    '--chunksize',
    help='Chunk size of loading pairs.',
    type=float,
    default=100000000,
    metavar='Float',
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=10,
    metavar='INT',
    show_default=True,
)
def pairs2cool2(pairs, outcool,
               binsize, min_mapq,
               chunksize, threads):
    """
    Convert pairs file into a specified resolution cool file.

        INPUT_PAIRS_PATH : Path of pairs file, can be compressed.

        OUT_COOL_PATH : Output path of cool file.
    """
    from .core import Pairs2 

    logger.info(f"Load pairs: `{pairs}`.")
    logger.info(f"Output cool's bin size: {binsize}")
    binsize = humanized2numeric(binsize)
    p = Pairs2(pairs, min_mapq=min_mapq, 
               chunksize=chunksize, threads=threads)
    p.to_cool(outcool, binsize=binsize)


@cli.command(cls=RichCommand, epilog=__epilog__)
@click.option(
    '-m',
    '--matrix',
    metavar='COOL',
    help="Contacts matrix stored by Cool format.",
    required=True,
    type=click.Path(exists=True)
)
@click.option(
    '-a',
    '--agp',
    metavar='AGP',
    type=click.Path(exists=True)
)
# @click.option(
#     '--factor',
#     '-k',
#     help='Factor of plot matrix. '
#             'If you input 10k matrix and want to plot heatmap at 500k, '
#             'factor should be set with 50.',
#     type=int,
#     default=50,
#     show_default=True
# )
@click.option(
    '-bs',
    '--binsize',
    metavar='STR',
    help='Bin size of the heatmap you want to plot. '
            ' Enabled suffix with k or m. [defalt: input matrix binsize]',
    default=None,
)
@click.option(
    '--no-adjust',
    'no_adjust',
    help='Matrix is chromosome-level, no need adjust.',
    default=False,
    is_flag=True,
    show_default=True,
    hidden=True,
)
@click.option(
    '--only-adjust',
    'only_adjust',
    help='Only adjust the matrix by agp, do not need plot the heatmap.',
    default=False,
    is_flag=True,
    show_default=True,
)
# @click.option(
#     '--coarsen',
#     'coarsen',
#     help='Need coarsen matrix to base bin size * k resolution, when only plot mode (agp not specified).',
#     default=False,
#     is_flag=True,
#     show_default=True
# )
# @click.option(
#     '--no-coarsen',
#     'no_coarsen',
#     help='The resolution of matrix is already for plotting, no need coarsen, when after adjust mode (agp specified).',
#     default=False,
#     is_flag=True,
#     show_default=True
# )
@click.option(
    '--only-coarsen',
    'only_coarsen',
    help='Only coarsen the input matrix, do not need plot the heatmap.',
    default=False,
    is_flag=True,
    show_default=True
)
@click.option(
    '--only-plot',
    'only_plot',
    help='Only plot the matrix.',
    default=False,
    is_flag=True,
    show_default=True,
    hidden=True,
)
@click.option(
    '-o',
    '--output',
    help='Output path of file.',
    default="plot.heatmap.png",
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads. (unused)',
    type=int,
    default=8,
    metavar='INT',
    show_default=True,
)
@click.option(
    "-s",
    "--scale",
    metavar="STR",
    help="Method of contact normalization",
    default="log1p",
    type=click.Choice(["log1p", "log", "none"]),
    show_default=True,
)
@click.option(
    "--triangle",
    is_flag=True,
    help="Plot the heatmap in triangle",
    default=False,
)
@click.option(
    "-vmin",
    "--vmin",
    metavar="FLOAT",
    default=None,
    type=float,
    show_default=True,
)
@click.option(
    "-vmax",
    "--vmax",
    metavar="FLOAT",
    default=None,
    type=float,
    show_default=True,
)
@click.option(
    '-c',
    '--chromosomes',
    help='Chromosomes and order in which the chromosomes should be plotted. '
            'Comma seperated. or a one column file',
    default=''
)
@click.option(
    "-dns",
    "--disable-natural-sort",
    is_flag=True,
    help="Disable natural sort of chromosomes, only used for `--chromosomes` or `--only-chr`",
    default=False,
)
@click.option(
    '-pc',
    '--per-chromosomes',
    '--per-chromosome',
    'per_chromosomes',
    help='Instead of plotting the whole matrix, '
            'each chromosome is plotted next to the other.',
    is_flag=True,
    default=False
)
@click.option(
    '-oc',
    '--only-chr',
    help='Only plot the chromosomes that ignore unanchored contigs. '
         'When `--chromosomes` specifed, this parameter will be ignored. '
         'The default use prefix of `Chr` to find the chromosomes. '
        '`--chr-prefix` can be used to change this.',
    is_flag=True,
    default=False
)
@click.option(
    '-cp',
    '--chr-prefix',
    metavar="STR",
    help="Prefix of the chromosomes, only used for `--only-chr`",
    default='Chr',
    show_default=True
)
@click.option(
    '-cpr',
    '--chrom-per-row',
    'chrom_per_row',
    help='Number of chromosome plot in each row',
    type=int,
    default=4,
    show_default=True
)
@click.option(
    '--fontsize',
    metavar="INT",
    help="Fontsize of the ticks, default is auto",
    type=int,
    default=None,
    show_default=True
)
@click.option(
    '-dpi',
    '--dpi',
    help='Resolution for the image.',
    default=600,
    show_default=True
)
@click.option(
    '-cmap',
    "-cm",
    '--cmap',
    "--colormap",
    help="""
    Colormap of heatmap. 
    Available values can be seen : 
    https://pratiman-91.github.io/colormaps/ 
    and http://matplotlib.org/examples/color/colormaps_reference.html and `whitered` .
    """,
    default='redp1_r',
    show_default=True
)
@click.option(
    '-wr',
    '--whitered',
    help="""
    Preset of `--scale none --colormap whitered`
    """,
    is_flag=True,
    default=False,
    show_default=True
)
@click.option(
    '--balance',
    help="""
    balance the matrix.
    """,
    is_flag=True,
    default=False,
    show_default=True
)
@click.option(
    '--balanced',
    help="""
    Plot balanced values, which need cool have weights columns in bins.
    """,
    is_flag=True,
    default=False,
    show_default=True
)
@click.option(
    '-nl',
    '--no-lines',
    'no_lines',
    help="""
    Don't add dash line in chromosome boundaries.
    """,
    default=False,
    show_default=True,
    is_flag=True,
)
@click.option(
    '-nt',
    '--no-ticks',
    'no_ticks',
    help="""
    Don't add ticks both in x axis and y axis.
    """,
    default=False,
    show_default=True,
    is_flag=True,
)
@click.option(
    '-rx',
    '--rotate-xticks',
    'rotate_xticks',
    help="Rotate the x ticks",
    is_flag=True,
    default=True, 
    show_default=True
)
@click.option(
    '-ry',
    '--rotate-yticks',
    'rotate_yticks',
    help="Rotate the x ticks",
    is_flag=True,
    default=False, 
    show_default=True
)
def plot(matrix, 
            agp, 
            binsize,
            no_adjust, 
            only_adjust, 
            # coarsen,
            # no_coarsen,
            only_coarsen,
            only_plot, 
            output,
            threads, 
            scale,
            triangle,
            vmin,
            vmax,
            chromosomes, 
            disable_natural_sort,
            per_chromosomes,
            only_chr,
            chr_prefix,
            chrom_per_row, 
            fontsize,
            dpi, 
            cmap,
            whitered,
            balance,
            balanced,
            no_lines,
            no_ticks,
            rotate_xticks,
            rotate_yticks,):
    """
    Adjust or Plot the contacts matrix after assembling.

    > **Usage:**
    > - adjust the matrix by agp and plot a heatmap 
    ```bash
    cphasing plot -a groups.agp -m sample.10000.cool -o groups.500k.wg.png
    ```
    > - adjust the matrix by agp and plot a 100k resolution heatmap
    ```bash
    cphasing plot -a groups.agp \\
        -m sample.10000.cool \\
        -o groups.100k.wg.png \\
        -bs 100k
    ```
    > - only plot a heatmap
    ```bash
    cphasing plot -m sample.100k.cool -o sample.100k.png
    ```
    > - Plot some chromosomes 
    ```bash
    cphasing plot -m sample.100k.cool -c Chr01,Chr02 -o Chr01_Chr02.100k.png
    ```
    """
    import cooler
    from .plot import (
        adjust_matrix,
        coarsen_matrix, 
        balance_matrix,
        plot_heatmap
        
    )
    coarsen = False
    no_coarsen = False
    if agp is None:
        only_plot = True 
        if not only_coarsen:
            logger.warning( "Only plot the matrix. "
                    "If you want to adjust matrix to chromosome-level, please provide agp file. ")
        else:
            logger.warning("Only coarsen the input matrix.")

    if not cooler.fileops.is_cooler(matrix):
        logger.error(f"Input file `{matrix}` is not a cool file.")
        sys.exit(-1)

    cool = cooler.Cooler(matrix)
    cool_binsize = cool.info['bin-size']
    bins = cool.bins()[:]
    if cool_binsize is None:
        cool_binsize = np.argmax(np.bincount(bins['end'] - bins['start']))
   
        
    if binsize is None and agp is not None:
        binsize = "500k"

    if binsize is not None:
        binsize = humanized2numeric(binsize)
        if binsize < cool_binsize:
            logger.warning(f"The matrix's binsize is larger than the specified binsize, "
                        f"set the heatmap's `{binsize}` to `{to_humanized2(cool_binsize)}`.")
            binsize = cool_binsize
            factor = 1
        elif binsize > cool_binsize:
            factor = binsize // cool_binsize
            
            
            if cool_binsize * factor < binsize:
                logger.warning(f"The input matrix's binsize is smaller than the heatmap's {to_humanized2(binsize)}, "
                                f"the specified binsize ({to_humanized2(binsize)}) should be a `factor * {to_humanized2(cool_binsize)}`, "
                               f" automaticly set the binsize to {to_humanized2(cool_binsize * factor)}.")
            else:
                logger.info(f"The input matrix's binsize is smaller than the heatmap's {to_humanized2(binsize)}, "
                        f"coarsen the matrix's binsize to {to_humanized2(cool_binsize * factor)}.")
            coarsen = True
    else:
        factor = 1
        binsize = cool_binsize


    if not only_plot:
        if not no_adjust:
            if which("bedtools") is None:
                raise ValueError(f"bedtools: command not found.")
            matrix = adjust_matrix(matrix, agp, threads=threads)

        if only_adjust:
            sys.exit()
        
        if not no_coarsen and factor > 1:
            matrix = coarsen_matrix(matrix, factor, None, threads)   
    
    else:
        if only_coarsen or coarsen:
            matrix = coarsen_matrix(matrix, factor, None, threads)   

    if balance:
        cool = cooler.Cooler(matrix)
        try:
            cool.bins()[:]['weight']
            if balanced:
                logger.warning("Weight existed in matrix, "
                            "the `--balanced --balance` parameter added, force execute balance.")
                balance_matrix(matrix, force=True, threads=threads)
            else:
                balanced = True
                logger.warning("Weight existed in matrix, skip the balance step, "
                        "or add `--balanced --balance` to force rerun balance.")
        except KeyError:

            balance_matrix(matrix, threads)
            balanced = True

    if op.exists(chromosomes):
        logger.info(f"Load chromosomes list from the file `{chromosomes}`.")
        chromosomes = [i.strip().split("\t")[0] for i in open(chromosomes) if i.strip()]
        
    else:
        
        chromosomes = chromosomes.strip().strip(",").split(',') if chromosomes else None
    
    if not chromosomes and only_chr:
        regex = re.compile(f"^{chr_prefix}")
        chroms = cooler.Cooler(matrix).chromnames 
        chromosomes = list(filter(lambda x: regex.findall(x), chroms))

    if chromosomes and only_chr:
        if not disable_natural_sort:
            chromosomes = natsorted(chromosomes)

    
    if no_ticks:
        xticks = False 
        yticks = False
    else: 
        xticks = True 
        yticks = True

    if balanced:
        cool = cooler.Cooler(matrix)
        try:
            cool.bins()[:]['weight']
            logger.info("Plot the balanced matrix.")
        except KeyError:
            logger.warning("There is not column of `bins/weight` in cool file, set `balanced=False`")
            balanced = False
    
    if whitered:
        logger.info("You have specified `--whitered`, set `--scale none --cmap whitered`.")
        scale = "none"
        cmap = "whitered"

    if not only_coarsen:
        plot_heatmap(matrix,
                 output,
                 chromosomes=chromosomes,
                 per_chromosomes=per_chromosomes,
                 chrom_per_row=chrom_per_row,
                 fontsize=fontsize,
                 dpi=dpi,
                 cmap=cmap, 
                 scale=scale,
                 triangle=triangle,
                 vmin=vmin,
                 vmax=vmax,
                 balanced=balanced,
                 xticks=xticks, 
                 yticks=yticks,
                 rotate_xticks=rotate_xticks,
                 rotate_yticks=rotate_yticks,
                 add_lines=False if no_lines else True,
                 threads=threads)



## hic subcommand
from .hic.cli import hic

@cli.group(cls=CommandGroup, hidden=True, short_help='Evaluation tools.')
@click.pass_context
def evaluator(ctx):
    """
    **Lack of validation**, please dont use these function.
    """
    pass

@evaluator.command(short_help="Estimate the highly homologous ratio")
@click.option(
    "-f",
    "--fasta",
    metavar="STR",
    help="Fasta file",
    required=True,
)
@click.option(
    '-ms',
    '--min-similarity',
    metavar='FLOAT',
    default=0.999,
    type=click.FloatRange(0, 1.0),
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=10,
    metavar='INT',
    show_default=True,
)
@click.option(
    '--aligner',
    metavar="STR",
    help="alinger",
    type=click.Choice(['wfmash', 'minigraph']),
    default='minigraph',
    show_default=True,
)
def homo(fasta, min_similarity, threads, aligner):
    from .evaluator.seqeval import HomoAnalysis
    ha = HomoAnalysis(fasta=fasta,
                      min_similarity=min_similarity,
                      threads=threads,
                      aligner=aligner)
    ha.run()

@evaluator.command(short_help='Calculate the allelic error rate', hidden=True)
@click.argument(
    "cluster",
    metavar='Cluster',
    type=click.Path(exists=True)
)
@click.argument(
    'alleletable',
    metavar='AlleleTable',
    type=click.Path(exists=True)
)
@click.argument(
    "contigsizes",
    metavar="contigsizes",
    type=click.Path(exists=True),
)
def allelic_error(cluster, alleletable, contigsizes):
    from .evaluate import allelic_error
    allelic_error(cluster, alleletable, contigsizes)


@evaluator.command(short_help="Evaluate the switch errors by ultra-long reads")
@click.option(
    "-f",
    "--fasta",
    metavar="STR",
    help="Fasta file",
    required=True,
)
@click.option(
    "-i",
    "--reads",
    metavar="STR",
    help="Reads file path, multiple file use comma separated",
    required=True,
)
@click.option(
    '-w',
    '--window',
    help='window size',
    default=5000,
    show_default=True,
    type=int
)
@click.option(
    '-m',
    '--min-windows',
    'min_windows',
    help='minimum windows of read',
    default=3,
    show_default=True,
    type=int,
)
@click.option(
    '-g',
    '--max-gap',
    'max_gap',
    help='maximum gap of two adjency alignments',
    default=1000000,
    show_default=True,
    type=int,
)
@click.option(
    '--hifi',
    help="Input hifi data.",
    default=False,
    is_flag=True,
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=10,
    metavar='INT',
    show_default=True,
)
@click.option(
    '-o',
    '--output',
    help='output file prefix',
    default='output',
    show_default=True
)
def sw(fasta, reads, window, min_windows, 
       max_gap, hifi, threads, output):
    from .evaluator.onteval import SwitchError
    se = SwitchError(fasta, reads,
                    window=window,
                    min_windows=min_windows,
                    maximum_gap=max_gap,
                    is_hifi=hifi,
                    threads=threads,
                    outprefix=output)
    se.run()


@evaluator.command(short_help='Analysis sequences compoent by hifi or ont reads ')
@click.option(
    "-f",
    "--fasta",
    metavar="STR",
    help="Fasta file",
    required=True,
)
@click.option(
    "-i",
    "--reads",
    metavar="STR",
    help="Reads file path, multiple file use comma separated",
    required=True,
)
@click.option(
    '-w',
    '--window',
    help='window size',
    default=5000,
    show_default=True,
    type=int
)
@click.option(
    '-s',
    '--step',
    help='step size',
    default=1000,
    show_default=True,
    type=int
)
@click.option(
    '-m',
    '--min-windows',
    'min_windows',
    help='minimum windows of read',
    default=2,
    show_default=True,
    type=int,
)
@click.option(
    '--is-ont',
    help="Input ont data.",
    default=False,
    is_flag=True,
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=10,
    metavar='INT',
    show_default=True,
)
@click.option(
    '-o',
    '--output',
    help='output file prefix',
    default='output',
    show_default=True
)
def ca(fasta, reads, window,
        step, min_windows,
       is_ont, threads, output):
    from .evaluator.onteval import ComponentAnalysis
    ca = ComponentAnalysis(fasta, reads, 
                           window=window, step=step,
                           min_windows=min_windows,
                           is_ont=is_ont,
                           threads=threads,
                           output=output)
    ca.run()


@cli.group(cls=CommandGroup, short_help='Misc tools.', epilog=__epilog__)
@click.pass_context
def utils(ctx):
    pass

@utils.command(cls=RichCommand, short_help='Convert agp to assembly file.')
@click.argument(
    "agpfile",
    metavar="AGP",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout
)
@click.option(
    "--add_gap",
    metavar="Bool",
    is_flag=True,
    default=False,
    show_default=True,
    help="add gap into assembly."
)
def agp2assembly(agpfile, output, add_gap):
    """
    Convert agp to cluster file.

    AGP : Path to agp file.
    """
    from .agp import agp2assembly
    agp2assembly(agpfile, output, add_gap)


@utils.command(cls=RichCommand, short_help='Convert agp to cluster file.')
@click.argument(
    "agpfile",
    metavar="AGP",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout
)
def agp2cluster(agpfile, output):
    """
    Convert agp to cluster file.

    AGP : Path to agp file.
    """
    from .agp import agp2cluster
    agp2cluster(agpfile, output)


@utils.command(cls=RichCommand, short_help='Convert agp to fasta file.')
@click.argument(
    "agpfile",
    metavar="AGP",
    type=click.Path(exists=True)
    )
@click.argument(
    "fasta",
    metavar="Fasta",
    type=click.Path(exists=True)

)
@click.option(
    "-o",
    "--output",
    help="Output of results",
    type=click.File('w'),
    default=sys.stdout
)
@click.option(
    "-c",
    "--contig",
    "--contigs",
    help="Outptu contig-level fasta",
    is_flag=True,
    default=False,
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
    hidden=True

)
def agp2fasta(agpfile, fasta, output, contig, threads):
    """
    Convert agp to fasta

        AGP: Path of AGP file

        FASTA: Path of draft assembly
    """
    from .agp import agp2fasta
    agp2fasta(agpfile, fasta, output, contig)


@utils.command(cls=RichCommand)
@click.argument(
    "agp",
    metavar="AGP",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--outdir",
    metavar='STR',
    help="Output directory of tour.",
    default="tour",
    show_default=True
)
@click.option(
    "-f",
    "--force",
    default=False,
    is_flag=True,
    show_default=True,
    help="Force output."
)
def agp2tour(agp, outdir, force):
    """
    Convert agp to several tour file.

    AGP : Path to agp file
    """
    from .agp import agp2tour

    agp2tour(agp, outdir, force)

@utils.command(cls=RichCommand, short_help='Rename dupliacted contig in agp file')
@click.argument(
    "agp",
    metavar="AGP",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Output of results",
    type=click.File('w'),
    default=sys.stdout
)
def agp_dup(agp, output):
    """
    Rename duplicated contig to another name.
    Note: this function conflict with break contig with raw name, which mean breaked contig name should be renamed.
    """

    from .agp import agp_dup
    agp_dup(agp=agp, output=output)


@utils.command(cls=RichCommand, short_help="Convert assembly to agp file and contigs file")
@click.argument(
    "assembly",
    metavar="Assembly",
    type=click.Path(exists=True)
)
# @click.argument(
#     "fasta",
#     metavar="Fasta",
#     type=click.Path(exists=True)

# )
@click.option(
    "-n",
    help="""
    Number of groups. If set to 0 or None, the partition will be run automatically.
    If in phasing mode, you can set to n1:n2, 
    which meaning that generate n1 groups in the first round partition
    and generate n2 groups in the second round partition. 
    Also, you can set `-n 8:0` to set the second round partition to automatical mode.
     The n2 parameter also can be set a file, which containing two columns, 
      the firse column repersent the group index from `first.clusters.txt`, and the 
       second column repersent the n2 groups in each n1 group. [default: 0]
    """,
    default=None,
    show_default=True,
)
# @click.option(
#     "-m",
#     "--max-chrom",
#     metavar="INT",
#     help="maximum number of chromosome or scaffold, and remain scaffolds will output into contig-level",
#     type=int,
#     default=None,
#     show_default=True
# )
@click.option(
    "-cp",
    "--chr-prefix",
    help="The prefix of chromosome",
    metavar="STR",
    default="Chr",
    show_default=True
)
@click.option(
    "-p",
    "--phased",
    help="Output phased chromosome name, {prefix}g{num}.",
    is_flag=True,
    default=False,
)
@click.option(
    "-s",
    "--sort",
    "sort_by_length",
    help="Sort scaffolds by scaffold-length",
    default=False,
    is_flag=True
)
@click.option(
    '-o',
    '--outprefix',
    help='output prefix, if none use the prefix of assembly',
    default=None,
    show_default=True
)
def assembly2agp(assembly, n, 
                chr_prefix, phased, sort_by_length, outprefix):
    from .agp import assembly2agp 

    n = re.split(":|x|\|", n) if n else None
    if n is not None:
        if len(n) == 1:
            n = [int(n[0])]
        elif len(n) >= 2:
            n = list(map(int, n))
    
    assembly2agp(assembly, n, 
                    chrom_prefix=chr_prefix, phased=phased,
                    sort_by_length=sort_by_length, outprefix=outprefix)


@utils.command(cls=RichCommand, short_help='Statistics of contigsizes')
@click.argument(
    'contigsizes',
    metavar="ContigSizes",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Output of results",
    type=click.File('w'),
    default=sys.stdout
)
def stat_contigsizes(contigsizes, output):
    from .utilities import stat_contigsizes, read_chrom_sizes

    contigsizes = read_chrom_sizes(contigsizes)
    df = stat_contigsizes(contigsizes)
    df.to_csv(output, sep='\t', index=True, header=None)


@utils.command(short_help='Statistics of AGP.')
@click.argument(
    "agp",
    metavar='AGP',
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Output of results",
    type=click.File('w'),
    default=sys.stdout
)
def statagp(agp, output):
    """
    Statistics of AGP.

    AGP : Path to agp file.

    """
    from .agp import statagp
    
    statagp(agp, output)


@utils.command(cls=RichCommand, short_help='Statistics of ClusterTable.')
@click.argument(
    "cluster",
    metavar='Cluster',
    type=click.Path(exists=True)
)
@click.argument(
    "contigsizes",
    metavar="contigsizes",
    type=click.Path(exists=True),
)
@click.option(
    "-o",
    "--output",
    help="Output of results",
    type=click.File('w'),
    default=sys.stdout
)
def statcluster(cluster, contigsizes, output):
    from .core import ClusterTable
    from .utilities import read_chrom_sizes
    
    ct = ClusterTable(cluster)
    df = read_chrom_sizes(contigsizes)

    for group in ct.groups:
        _contigs = ct.data[group]
        print(group, df.loc[_contigs]['length'].sum(), file=output)

@utils.command(cls=RichCommand)
@click.argument(
    "cluster",
    metavar='Cluster',
    type=click.Path(exists=True)
)
@click.argument(
    "contigsizes",
    metavar="Contigsizes",
    type=click.Path(exists=True),
)
@click.option(
    "--init-agp",
    metavar="STR",
    help="input a init agp to set the contig orientation, default set all to +",
    default=None,
    show_default=True,
)
@click.option(
    "-o",
    "--output",
    help="Output of results",
    type=click.File('w'),
    default=sys.stdout
)
def cluster2agp(cluster, contigsizes, init_agp, output):
    """
    Convert cluster table to a pseudo agp file

        which the order and orientation of contig are pseudo

        Cluster: Path of the cluster table

        Contigsizes: Path of contig sizes 

    """
    from .core import ClusterTable
    from .agp import import_agp 
    orientation_db = {}
    ct = ClusterTable(cluster)
    if init_agp:
        agp_df, _ = import_agp(init_agp)
        tmp_df = agp_df[['id', 'orientation']]
        tmp_df.set_index('id', inplace=True)
        orientation_db = tmp_df.to_dict()['orientation']
    
    ct.to_agp(contigsizes=contigsizes, output=output, orientation_db=orientation_db)


@utils.command(cls=RichCommand, short_help='Convert cluster to several count RE files.')
@click.argument(
    "cluster",
    metavar='Cluster',
    type=click.Path(exists=True)
)
@click.argument(
    'count_re',
    metavar='CountRE',
    type=click.Path(exists=True)
)
def cluster2count(cluster, count_re):
    """
    Convert cluster to several count RE files.

    ClusterTable : Path to cluster table.

    CountRE : Path to countRE file.
    
    """
    from .core import ClusterTable
    ct = ClusterTable(cluster)
 
    ct.to_countRE(count_re)


@utils.command(cls=RichCommand, short_help='Extract fasta by cluster table')
@click.argument(
    "cluster",
    metavar='Cluster',
    type=click.Path(exists=True)
)
@click.argument(
    'fasta',
    metavar='Fasta',
    type=click.Path(exists=True)
)
def cluster2fasta(cluster, fasta):
    """
    Extract fasta by cluster table
    """
    from .core import ClusterTable
    ct = ClusterTable(cluster)
    ct.to_fasta(fasta=fasta)

@utils.command(cls=RichCommand, short_help='Convert cluster to pseudo assembly file.')
@click.argument(
    "cluster",
    metavar='Cluster',
    type=click.Path(exists=True)
)
@click.argument(
    'fasta',
    metavar='Fasta',
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout
)
def cluster2assembly(cluster, fasta, output):
    """
    Convert cluster to a pseudo assembly file.

    ClusterTable : Path to cluster table.

    Fasta : Path to fasta file.
    
    """
    from .core import ClusterTable
    ct = ClusterTable(cluster)
 
    ct.to_assembly(fasta, output)

@utils.command(short_help='Convert cluster to pseudo tour files.')
@click.argument(
    "cluster",
    metavar='Cluster',
    type=click.Path(exists=True)
)
def cluster2tour(cluster):
    """
    Convert cluster to several pseudo tour files.

    ClusterTable : Path to cluster table.
    
    """
    from .core import ClusterTable
    ct = ClusterTable(cluster)
 
    ct.to_tour()

@utils.command(cls=RichCommand)
@click.argument(
    "count_re",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout
)
@click.option(
    '--fofn',
    help="""
    If this flag is set then the countREs is a file of
    filenames corresponding to the contact tables you want to merge.
    This is workaround for when the command line gets too long.
    """,
    is_flag=True,
    default=False,
    show_default=True
)
def countRE2cluster(count_re, output, fofn):
    """
    Convert serveral count RE table to cluster file.
    """
    if fofn:
        countREs = [i.strip() for i in open(count_re) if i.strip()]
    else:
        countREs = [count_re]
    countREs = list(map(CountRE, countREs))

    group_name = list(map(lambda x: x.filename.rsplit(".", 1)[0], countREs))
    for group_name, cr in zip(group_name, countREs):
        print(f"{group_name}\t{cr.ncontigs}\t{' '.join(cr.contigs)}", 
                file=output) 

@utils.command(cls=RichCommand)
@click.argument(
    "coolfile",
    metavar="INPUT_COOL_PATH",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout
)
def cool2depth(coolfile, output):
    """
    Calculate the depth from the cool file.
    """
    from .utilities import cool2depth
    cool2depth(coolfile, output)     


@utils.command(cls=RichCommand)
@click.argument(
    "coolfile",
    metavar="INPUT_COOL_PATH",
    type=click.Path(exists=True)
)
@click.argument(
    "outcool",
    metavar="OUTPUT_COOL_PATH"
)
@click.option(
    '-b',
    '--balanced',
    help='Load the balanced matrix. If not provided, load unbalanced counts',
    default=False,
    is_flag=True,
    show_default=True,
)
@click.option(
    '--cis',
    help='output cis contig contacts',
    default=False,
    is_flag=True,
    show_default=True
)
@click.option(
    '--min-contacts',
    'min_contacts',
    help='Minimum contacts for contig pair',
    default=1,
    show_default=True,
    type=int,
)
@click.option(
    '--no-mask-nan',
    'no_mask_nan',
    help='Do not mask nan bins.',
    default=False,
    show_default=True
)
@click.option(
    '--symmetric-upper',
    'symmetric_upper',
    default=True,
    show_default=True,
)
def merge_cool(coolfile, 
            outcool, 
            balanced,
            cis,
            min_contacts, 
            no_mask_nan,
            symmetric_upper):
    """
    merge slidewindows matrix into whole contig matrix.
    
    INPUT_COOL_PATH : Path to COOL file.

    OUTPUT_COOL_PATH : Path to output COOL file.

    """
    from .utilities import merge_matrix
    no_dia = False if cis else True
    merge_matrix(coolfile, outcool, 
                    balance=balanced,
                    min_contacts=min_contacts, 
                    no_mask_nan=no_mask_nan, 
                    no_dia=no_dia,
                    symmetric_upper=symmetric_upper)
    
@utils.command(cls=RichCommand)
@click.argument(
    "coolfile",
    metavar="INPUT_COOL_PATH"
)
@click.argument(
    "chromlist",
    metavar="CHROM_LIST_PATH"
)
@click.argument(
    "outcool",
    metavar="OUTPUT_COOL_PATH"
)
def extract_matrix(coolfile, chromlist, outcool):
    """
    extract matrix from a cool by chrom/contig list

        INPUT_COOL_PATH: Path of cool file.

        CHROM_LIAT_PATH: Path of chrom/contig list, one column.

        OUTPUT_COOL_PATH: Path of output cool file.

    """
    from .utilities import extract_matrix
    chromlist = [i.strip() for i in open(chromlist) if i.strip()]
    extract_matrix(coolfile, chromlist, outcool)

@utils.command(hidden=True, cls=RichCommand)
@click.argument(
    "coolfile",
    metavar="INPUT_COOL_PATH"
)
@click.argument(
    "prunetable",
    metavar="PRUNE_TABLE_PATH"
)
@click.argument(
    "outcool",
    metavar="OUTPUT_COOL_PATH"
)
def prune_matrix(coolfile, prunetable, outcool):
    """
    prune matrix by a prunetable

        INPUT_COOL_PATH: Path of cool file.

        PRUNE_TABLE_PATH: Path of prune table.

        OUTPUT_COOL_PATH: Path of output cool file.

    """
    from .utilities import prune_matrix
    from .core import PruneTable
    pt = PruneTable(prunetable)
    prunepairs = pt.data[['contig1', 'contig2']].values.tolist()
    prunepairs = list(map(tuple, prunepairs))

    prune_matrix(coolfile, prunepairs, outcool)

@utils.command(cls=RichCommand)
@click.argument(
    "real_list",
    metavar="Real_List",
    type=click.Path(exists=True)
)
@click.argument(
    "contigsizes",
    metavar="ContigSizes",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout
)
def pseudo_agp(real_list, contigsizes, output):
    """
    create a pseudo agp from simulation data

        Real_list: Path of two columns file, which first is chrom and second is contig

        ContigSizes: Path of contig sizes file.
    """
    from .agp import pseudo_agp
    pseudo_agp(real_list, contigsizes, output)


@utils.command(cls=RichCommand)
@click.argument(
    'alleletable',
    metavar='AlleleTable',
    type=click.Path(exists=True)
)
@click.argument(
    'contacts',
    metavar='INPUT_CONTACTS_PATH',
    type=click.Path(exists=True)
)
@click.option(
    "-m",
    "--min-values",
    "min_values",
    metavar="FLOAT",
    default=1.0,
    type=float,
    show_default=True 
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout
)
def filter_high_similarity_contigs(alleletable, 
                                   contacts,
                                   min_values,
                                   output):
    
    """
    filter high similarity by h-trans contacts 
    """
    from .alleles import filter_high_similarity_contigs

    filter_high_similarity_contigs(alleletable,
                                   contacts,
                                   min_values,
                                   output)


cli.add_command(statagp)
cli.add_command(agp2fasta)

ALIASES = {
    "allele": alleles,
    "pipe": pipeline,
    "align": alignments,
    "alignment": alignments,
    "hg": hypergraph,
    "hp": hyperpartition,
    "partition": hyperpartition,
    "ho": hyperoptimize,
    "optimize": scaffolding,
    "sf": scaffolding,
    "scaf": scaffolding, 
    "scaffold": scaffolding,
    "scaffolds": scaffolding,
    "contigsizes": chromsizes,
    "pairs2contact": pairs2contacts,
    "pair2cool": pairs2cool,
    "p2c": pairs2cool,
    "agp2fa": agp2fasta,
    "eval": evaluator,
}