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
# from ..cli import CommandGroup
from ..cli import cli, LOG_LEVEL
from ..utilities import run_cmd

logger = logging.getLogger(__name__)


# class CommandGroup(DYMGroup, RichCommand):
#     """
#     List subcommand in the order there were added.
#     """
#     def list_commands(self, ctx):
#         return list(self.commands)

#     def get_command(self, ctx, cmd_name):
#         """
#         Alias command, https://stackoverflow.com/questions/46641928/python-click-multiple-command-names
#         """
#         try:
#             cmd_name = ALIASES[cmd_name].name 
#         except KeyError:
#             pass 
#         return super().get_command(ctx, cmd_name)


# @cli.group(cls=CommandGroup,
#            context_settings={"help_option_names": ["-h", "--help", "-help"]},
#            epilog=__epilog__,
#            short_help='Using methylation information to refine the alignments. (nightly)')
# @click.pass_context
# def methalign(ctx):
#     pass


@cli.command(cls=RichCommand, epilog=__epilog__, short_help='Refine the alignments by methylation signal.')
@click.help_option('-h', '--help', '-help')
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
    hidden=True,
)
@click.option(
    '-o',
    '--output',
    type=str,
    help='Output suffix of paf file name',
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
def methalign(fasta, bedgraph, bam, penalty, ref_prob_cutoff,
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
    from .refine import refine
    from ..utilities import read_fasta
    
    method = "fast"
    

    if len(bam) == 0:
        logging.error("No bam files are provided.")
        sys.exit(1)

    logger.info("Input files:")
    logger.info(f"Reference fasta: {fasta}")
    logger.info(f"Reference bedgraph: {bedgraph}")
    logger.info(f"Bam files: {', '.join(bam)}")

    if method == "v1":
        fa_dict = read_fasta(fasta)

        hifi_methylation_dict = parse_bedgraph(bedgraph, ref_prob_cutoff, threads=threads)
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
    
    else:
        if len(bam) > 1:
            args = []
            for b in bam:
                _output = Path(b).stem + '.refined.bam'
                args.append((
                    b,
                    fasta,
                    bedgraph,
                    penalty,
                    ref_prob_cutoff,
                    prob_cutoff,
                    designate_mapq,
                    _output,
                    threads
                ))
            out_bams = Parallel(n_jobs=process)(
                delayed(refine)(*arg) for arg in args
            )

        else:
            out_bam = refine(
                bam[0],
                fasta,
                bedgraph,
                penalty,
                ref_prob_cutoff,
                prob_cutoff,
                designate_mapq,
                output=Path(bam[0]).stem + '.refined.bam',
                threads=threads
            )

            out_bams = [out_bam]
    
    def bam2paf(bam):
        cmd = ['cphasing-rs', 'bam2paf', bam, '-o', bam.replace('.bam', '.paf.gz')]
        log_dir = Path("logs").mkdir(exist_ok=True)
        log = op.join("logs", f"{Path(bam).stem}.bam2paf.log")
        flag = run_cmd(cmd, log=log)
        assert flag == 0, f"Failed to convert {bam} to paf format."

        return bam.replace('.bam', '.paf.gz')
    
    output_prefix = Path.Path(out_bams[0]).stem.split(".")[0]


    if len(out_bams) == 1:
        pafs = [bam2paf(out_bams[0])]
    else:
        pafs = Parallel(n_jobs=process)(
                delayed(bam2paf)(bam) for bam in out_bams
            )

    if len(pafs) > 1:
     
        cmd = ['cat'] + pafs + ['>', f"{output_prefix}.{output}"]
        flag = os.system(' '.join(cmd))
        assert flag == 0, "Failed to concatenate paf files."

        logger.info(f"Output merged paf file is saved to {output_prefix}.{output}")

        ## remove pafs file 
        for paf in pafs:
            if Path(paf).exists():
                os.remove(paf)  
        

    else:
        os.rename(pafs[0], f"{output_prefix}.{output}")
        logger.info(f"Output paf file is saved to {output_prefix}.{output}")

    paf2porec_cmd = ['cphasing-rs', 'paf2porec', output, '-q', '0',
                     '-o', f"{output_prefix}.{output}".replace('.paf.gz', '.porec.gz')]
    logs = op.join("logs", f"{output_prefix}.paf2porec.log")
    flag = run_cmd(paf2porec_cmd, log=logs)
    assert flag == 0, "Failed to convert paf to porec format."
    output_porec = f"{output_prefix}.{output}".replace('.paf.gz', '.porec.gz')
    logger.info(f"Output porec file is saved to {output_porec}")

    logger.info("All done.")
    

ALIASES = {
   
}