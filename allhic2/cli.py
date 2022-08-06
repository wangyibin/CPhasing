#!/usr/bin/env python

import click 
import sys 

from . import __version__

class CommandGroup(click.Group):
    """
    List subcommand in the order there were added.
    """
    def list_commands(self, ctx):
        return list(self.commands)

@click.version_option(__version__, "-V", "--version")
@click.group(context_settings={"help_option_names": ["-h", "--help"]},
            cls=CommandGroup)
def cli():
    pass

def auto():
    pass


@cli.command()
@click.option(
    '-r',
    'reference',
    metavar='STR',
    required=True
)
@click.option(
    '-1',
    'fastq1',
    metavar='STR',
    required=True
)
@click.option(
    '-2',
    'fastq2',
    metavar='STR',
    required=True
    )
@click.option(
    '-t',
    '--threads',
    help='Number of threads',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,

)
def mapper(
    reference,
    fastq1, 
    fastq2,
    threads

):
    """
    mapper of Hi-C reads
    """
    from .mapper import HisatMapper
    hm1 = HisatMapper(reference, fastq1)
    hm1.run()
    hm2 = HisatMapper(reference, fastq2)
    hm2.run()
    pass

def correct():
    pass

def pregroup():
    pass

def extract():
    pass

def normalize():
    pass

def partition():
    pass

@cli.command()
@click.argument(
    'clustertable',
    metavar='ClusterTable'
)
@click.argument(
    'alleletable',
    metavar='AlleleTable',
)
@click.argument(
    'count_re',
    metavar='CountRE'
)
@click.argument(
    'pairs',
    metavar='Pairs'
)
@click.argument(
    'ploidy',
    metavar='Ploidy',
    type=int
)
@click.option(
    '--method',
    help=('method of recluster, must in {greedy, strict}. '
        '`greedy`, greedy to rescue uncluster contigs, '
        'which contig rescued in first gap. '
        '`strict`, strict to rescue uncluster contigs, '
        'which contig must only rescued to one gap.'),
    type=click.Choice(['greedy', 'strict']),
    default='greedy',
    show_default=True,
    metavar='STR',
)
def recluster(
    clustertable, 
    alleletable, 
    count_re,
    pairs,
    ploidy,
    method
):
    """
    recluster partition results by allele table.

        ClusterTable : Path to cluster table.

        AlleleTable : Path to allele table.

        CountRE : Path to countRE file.

        Pairs : Path to pairs.

        Ploidy : Ploidy.

    """
    from .recluster import reCluster

    rc = reCluster(clustertable, alleletable, ploidy, count_re, pairs)
    rc.run(method=method)

def rescue():
    pass

def plot():
    pass

@cli.group(cls=CommandGroup, short_help='Misc tools.')
@click.pass_context
def utils(ctx):
    pass


@utils.command()
@click.argument(
    "agp",
)
@click.option(
    "-o",
    "--output",
    type=click.File('w'),
    default=sys.stdout
)
def agp2cluster(agp, output):
    from .agp import agp2cluster
    agp2cluster(agp, output)

@utils.command()
@click.argument(
    "agp",
    )
def agp2fasta(agp):
    pass

@utils.command()
@click.argument(
    "agp",
)
def agp2tour(agp):
    pass