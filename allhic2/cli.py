#!/usr/bin/env python

from allhic2.core import ClusterTable
import click
import logging
import sys 

from . import __version__

logger = logging.getLogger("allhic2")
class CommandGroup(click.Group):
    """
    List subcommand in the order there were added.
    """
    def list_commands(self, ctx):
        return list(self.commands)

@click.version_option(__version__, "-V", "--version")
@click.group(context_settings={"help_option_names": ["-h", "--help"]},
            cls=CommandGroup)
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
    if quiet:
        logger.setLevel(logging.CRITICAL)
    elif verbose > 0:
        LOG_LEVELS = [logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG]
        offset = 2 
        idx = min(len(LOG_LEVELS) - 1, offset + verbose)
        logger.setLevel(LOG_LEVELS[idx])
    else:
        logger.setLevel(logging.INFO)

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
    '-q',
    '--mapq',
    help='minimum quality of mapping',
    type=int,
    default=10,
    show_default=True
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
    Mapper for Hi-C reads mapping.
    """
    from .mapper import HisatMapper
    hm1 = HisatMapper(reference, fastq1, threads=threads)
    hm1.run()
    hm2 = HisatMapper(reference, fastq2, threads=threads)
    hm2.run()
    pass

def correct():
    pass

@cli.command()
@click.option(
    "-r",
    "--fasta",
    help="polyploid contig-level fasta",
    required=True
)
@click.option(
    "--cds"
)
@click.option(
    "--method",
    default='gmap',
    show_default=True,
    type=click.Choice(['gmap', 'synteny', 'seq_similarity']),
    
)
def alleles(fasta, cds, method):
    """
    Build allele table for prune.
    """
    logger.warning('Incomplete developing function')
    pass

def pregroup():
    pass

def prune():
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
    Recluster partition results by allele table.

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

@cli.command()
@click.argument(
    "fasta",
    metavar="Fasta"
)
@click.option(
    "-o",
    "--output",
    metavar="STR",
    help="Output genome fasta, '.gz' supported.",
    default="groups.asm.fasta",
    show_default=True
)
def build(fasta, output):
    """
    Build genome release.

    Fasta : contig-level fasta file
    """
    from .build import Build
    Build(fasta, output)

def plot():
    pass

def assess():
    pass

@cli.group(cls=CommandGroup, short_help='Misc tools.')
@click.pass_context
def utils(ctx):
    pass


@utils.command(short_help='Convert agp to cluster file.')
@click.argument(
    "agp",
    metavar="AGP"
)
@click.option(
    "-o",
    "--output",
    help="Output of results",
    type=click.File('w'),
    default=sys.stdout
)
def agp2cluster(agp, output):
    """
    Convert agp to cluster file.

    AGP : Path to agp file.
    """
    from .agp import agp2cluster
    agp2cluster(agp, output)

@utils.command(short_help='Convert agp to fasta file.')
@click.argument(
    "agp",
    )
@click.argument(
    "fasta",
    metavar="Fasta",

)
@click.option(
    "-o",
    "--output",
    help="Output of results",
    type=click.File('w'),
    default=sys.stdout
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
def agp2fasta(agp, fasta, output, threads):
    logger.warning('Incomplete developing function')
    pass

@utils.command()
@click.argument(
    "agp",
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

@utils.command(short_help='Statistics of AGP.')
@click.argument(
    "agp",
    metavar='AGP'
)
@click.option(
    "-o",
    "--output",
    help="Output of results",
    type=click.File('w'),
    default=sys.stdout
)
def agpstat(agp, output):
    """
    Statistics of AGP.

    AGP : Path to agp file.

    """
    from .agp import agpstat
    
    agpstat(agp, output)

@utils.command(short_help='Convert cluster to several count RE files.')
@click.argument(
    "cluster",
    metavar='Cluster'
)
@click.argument(
    'count_re',
    metavar='CountRE'
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
