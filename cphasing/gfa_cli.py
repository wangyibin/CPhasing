#!/usr/bin/env python

import logging
import sys 
import os
import os.path as op
import psutil 
import re

import rich_click as click
from rich_click import RichCommand
from click_didyoumean import DYMGroup, DYMCommandCollection

from . import __version__, __epilog__
from ._config import *

from .cli import cli

logger = logging.getLogger(__name__)

class CommandGroup(DYMGroup, click.RichGroup, RichCommand):
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


@cli.group(cls=CommandGroup, epilog=__epilog__, no_args_is_help=True, short_help="Gfa tools. (nightly)")
@click.pass_context
def gfa(ctx):
    pass 

@gfa.command(cls=RichCommand, epilog=__epilog__, no_args_is_help=True, short_help="Remove overlapped regions on utg")
@click.argument(
    "gfa",
    metavar="GFA",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout,
)
def remove_overlap(gfa, output):
    from .gfa import Gfa
    g = Gfa(gfa)
    overlap_df = g.overlap_to_bed(invert=True)
    overlap_df.set_index("Chromosome", inplace=True)

    for seqid, seq in g.seqment_seqs.items():
        if seqid in overlap_df.index:
            s, e = overlap_df.loc[seqid, ['Start', 'End']]
            new_seq = seq[s:e]
            print(f">{seqid}\n{new_seq}", file=output)
        else:
            print(f">{seqid}\n{seq}", file=output)


@gfa.command(cls=RichCommand, epilog=__epilog__, no_args_is_help=True, short_help="Get overlapped regions on utg")
@click.argument(
    "gfa",
    metavar="GFA",
    type=click.Path(exists=True)
)
@click.option(
    "-v",
    "--invert",
    help="Invert the overlapped regions, which mean that output non-overlapped regions",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout,
)
def overlap_bed(gfa, invert, output):
    from .gfa import Gfa
    g = Gfa(gfa)
    overlap_df = g.overlap_to_bed(invert=invert)
    overlap_df.set_index("Chromosome", inplace=True)

    overlap_df.to_csv(output, sep='\t', index=True, header=False)

@gfa.command(cls=RichCommand,
              epilog=__epilog__, no_args_is_help=True, 
              short_help="Get clm from gfa file")

@click.argument(
    "gfa",
    metavar="GFA",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout,
)
def gfa2clm(gfa, output):
    from .gfa import Gfa
    g = Gfa(gfa)
    g.generate_mock_hic_files()

@gfa.command(cls=RichCommand, epilog=__epilog__, no_args_is_help=True, short_help="Get contig sizes")
@click.argument(
    "gfa",
    metavar="GFA",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout,
)
def contigsizes(gfa, output):
    from .gfa import Gfa
    g = Gfa(gfa)
    for contig, length in g.segments.items():
        print(f"{contig}\t{length}", file=output)


@gfa.command(cls=RichCommand, epilog=__epilog__, no_args_is_help=True, short_help="Get overlapping region coverge of utg")
@click.argument(
    "gfa",
    metavar="GFA",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout,
)
def overlap_coverage(gfa, output):
    from .gfa import Gfa
    g = Gfa(gfa)
    res_df = g.calculate_overlap_percentage()
    res_df.to_csv(output, sep='\t', index=False, header=False)


def _process_one_gfa_summary(gfa_file):
    import pandas as pd
    from .gfa import Gfa
    g = Gfa(gfa_file)
    summary_dict = g.summary()
    return pd.Series(summary_dict, name=gfa_file)

@gfa.command(cls=RichCommand, epilog=__epilog__, no_args_is_help=True, short_help="Get summary of gfa")
@click.argument(
    "gfa",
    metavar="GFA",
    nargs=-1,
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout,
)
@click.option(
    "-t",
    "--threads",
    help="Number of threads for parallel processing.",
    type=int,
    default=4,
    show_default=True
)
def summary(gfa, output, threads):
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from .gfa import Gfa

    def _process_one(gfa_file):
        
        g = Gfa(gfa_file)
        summary_dict = g.summary()
        return pd.Series(summary_dict, name=gfa_file)
    
    summaries = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(_process_one_gfa_summary, gfa_file): gfa_file for gfa_file in gfa}
        for fut in as_completed(futures):
            try:
                summaries.append(fut.result())
            except Exception as e:
                logger.error(f"Failed to summarize {futures[fut]}: {e}")

    if summaries:
        summary_df = pd.DataFrame(summaries)
        summary_df = summary_df.T
        summary_df.to_csv(output, sep='\t', index=True, header=True)
    
@gfa.command(cls=RichCommand, epilog=__epilog__, no_args_is_help=True, short_help="Get redundant contigs")
@click.argument(
    "gfa",
    metavar="GFA",
    type=click.Path(exists=True)
)
@click.option(
    "-p",
    "--percent",
    help="The percentage of overlap to be considered as redundant contig",
    type=click.FloatRange(0.0, 100.0, clamp=True),
    default=100.0,
    show_default=True,
    
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout,
)
def redundant_contigs(gfa, percent, output):
    from .gfa import Gfa
    g = Gfa(gfa)
    res_df = g.calculate_overlap_percentage()

    res_df = res_df[res_df['Overlap_Percentage'] >= percent]
    logger.info("Number of redundant contigs: {}".format(len(res_df)))
    logger.info("Length of redundant contigs: {}".format(res_df['Length'].sum()))

    res_df['Chromosome'].to_csv(output, sep='\t', index=False, header=False)

@gfa.command(cls=RichCommand, epilog=__epilog__, no_args_is_help=True, short_help="Reassembe by rescue unanchored contigs and fill gaps in agp file with gfa")
@click.argument(
    "agp",
    metavar="AGP",
    type=click.Path(exists=True)
)
@click.argument(
    "gfa",
    metavar="GFA",
    type=click.Path(exists=True)
)
@click.option(
    "--disable-reassemble",
    help="Disable reassemble, which mean that program skip fill gap and rescue unanchor contigs",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    "--disable-remove-redundant",
    help="Disable remove redundant contigs, which mean that all contigs in agp will be rescued.",
    is_flag=True,
    default=False,
    show_default=True
)
@click.option(
    "-p",
    "--percent",
    help="The percentage of overlap to be considered as redundant contig",
    type=click.FloatRange(0.0, 100.0, clamp=True),
    default=100.0,
    show_default=True,
    
)
@click.option(
    "--chimeric-bed",
    help="The bed file of chimeric contigs.",
    type=click.Path(exists=True),
    show_default=True,
    default=None,
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout,
)
def reassemble(agp, gfa, disable_reassemble, disable_remove_redundant, percent, 
                chimeric_bed, output):
    from .gfa import Gfa
    g = Gfa(gfa)
    if chimeric_bed:
        g.apply_chimeric_bed(chimeric_bed)
    g.rescue_agp_with_gfa(agp, output,
                          redundant_perc=percent,
                          remove_redundant=not disable_remove_redundant,
                          rescue=not disable_reassemble)

@gfa.command(cls=RichCommand, epilog=__epilog__, no_args_is_help=True, short_help="Remove redundant contigs in agp file with gfa")
@click.argument(
    "agp",
    metavar="AGP",
    type=click.Path(exists=True)
)
@click.argument(
    "gfa",
    metavar="GFA",
    type=click.Path(exists=True)
)
@click.option(
    "-p",
    "--percent",
    help="The percentage of overlap to be considered as redundant contig",
    type=click.FloatRange(0.0, 100.0, clamp=True),
    default=100.0,
    show_default=True,
    
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout,
)
def remove_redundant(agp, gfa, percent, output):
    from .gfa import Gfa
    g = Gfa(gfa)
    g.rescue_agp_with_gfa(agp, output,
                          redundant_perc=percent,
                          remove_redundant=True,
                          rescue=False)


ALIASES = {
}