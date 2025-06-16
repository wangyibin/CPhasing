#!/usr/bin/env python

"""
cli for methylation align pipelines
"""
import argparse
import logging
import os
import os.path as op
import sys

import click
import logging
import sys
from rich_click import RichCommand
from pathlib import Path
from click_didyoumean import DYMGroup

from .. import __epilog__
from ..cli import CommandGroup
from ..cli import cli 
from ..utilities import run_cmd

logger = logging.getLogger(__name__)
class CommandGroup(DYMGroup, RichCommand):
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


@cli.group(cls=CommandGroup,
           context_settings={"help_option_names": ["-h", "--help", "-help"]},
           epilog=__epilog__,
           short_help='Using methylation information to refine the alignments. (nightly)')
@click.pass_context
def methalign(ctx):
    pass


@methalign.command(cls=RichCommand, short_help='Refine the alignments by methylation signal.')
@click.argument(
    'fasta',
    type=click.Path(exists=True),
)
@click.argument(
    'bedgraph',
    type=click.Path(exists=True),
)
@click.argument(
    'bam',
    nargs=-1,
    type=click.Path(exists=True),
)
@click.option(
    '-p',
    '--penalty',
    metavar="INT",
    default=2,
    type=click.IntRange(0, 20),
    help='Penalty for inconsistent 5mC ',
    show_default=True
)
@click.option(
    '-rc',
    '--ref-prob-cutoff',
    type=click.IntRange(0, 100),
    default=50,
    metavar="INT",
    help="probability cutoff of methylation of reference in bedgraph",
    show_default=True
)
@click.option(
    '-c',
    '--prob-cutoff',
    type=click.IntRange(0, 255),
    default=128,
    metavar="INT",
    help="probability cutoff of ML in porec bam",
    show_default=True
)
@click.option(
    '-dq',
    '--designate-mapq',
    type=click.IntRange(-1, 60),
    default=2,
    metavar="INT",
    help="designate MAPQ for the best alignments. " \
    "Set the value to -1 to disable the function, smaller than designate_mapq will be modified to designate_mapq, when refine steps", 
    show_default=True
)
@click.option(
    '-r',
    '--recalculate-all',
    is_flag=True,
    help="recalculate scores for all alignments although there may not be any secondary alignments",
    default=False,
    show_default=True,
)
@click.option(
    '-o',
    '--output',
    type=click.Path(),
    help='Output paf file name',
    default='methalign.refined.paf.gz',
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    type=int,
    default=8,
    help='Number of threads for single bam processing',
    show_default=True
)
@click.option(
    "--process",
    type=int,
    default=10,
    help="Number of processes to process bam",
    show_default=True
)
def refine(fasta, bedgraph, bam, penalty, ref_prob_cutoff,
            prob_cutoff, designate_mapq, recalculate_all,
            output, threads, process):
    """
    Refine the alignments by methylation information.


        FASTA: reference genome in fasta format

        BEDGRAPH: bedgraph file containing the methylation information, four columns are required: chromosome, start, end, and methylation level(0-100)

        BAM: bam file containing the alignments and MM/ML tags on the primary alignments
    """
    from joblib import Parallel, delayed
    from .filter_bam_methyl import (
        parse_bam, 
        parse_bedgraph
        )
    from ..utilities import read_fasta

    fa_dict = read_fasta(fasta)

    hifi_methylation_dict = parse_bedgraph(bedgraph, ref_prob_cutoff, threads=threads)

    if len(bam) == 0:
        logging.error("No bam files are provided.")
        sys.exit(1)
    
    if len(bam) > 1:
        def generate_args(bam):
            for bam in bam:
                yield (
                    bam,
                    fa_dict,
                    hifi_methylation_dict,
                    penalty,
                    prob_cutoff,
                    designate_mapq,
                    recalculate_all,
                    threads
                )
        
        out_bams = Parallel(n_jobs=process)(
            delayed(parse_bam)(*arg) for arg in generate_args(bam)
        )

    else:
        out_bam = parse_bam(
            bam[0],
            fa_dict,
            hifi_methylation_dict,
            penalty,
            prob_cutoff,
            designate_mapq,
            recalculate_all,
            threads
        )

        out_bams = [out_bam]

    def bam2paf(bam):
        cmd = ['cphasing-rs', 'bam2paf', bam, '-o', bam.replace('.bam', '.paf.gz')]

        run_cmd(cmd)

        return bam.replace('.bam', '.paf.gz')
    

    pafs = Parallel(n_jobs=process)(
            delayed(bam2paf)(bam) for bam in out_bams
        )
    
    if len(pafs) > 1:
        cmd = ['cat'] + pafs + ['>', output]
        flag = os.system(' '.join(cmd))
        assert flag == 0, "Failed to concatenate paf files."

        logger.info(f"Output merged paf file is saved to {output}")

        ## remove pafs file 
        for paf in pafs:
            if Path(paf).exists():
                os.remove(paf)  

    else:
        os.rename(pafs[0], output)
        logger.info(f"Output paf file is saved to {output}")

    paf2porec_cmd = ['cphasing-rs', 'paf2porec', output, '-q', '0',
                     '-o', output.replace('.paf.gz', '.porec.gz')]
    run_cmd(paf2porec_cmd)

    logger.info("All done.")


ALIASES = {
   
}