#!/usr/bin/env python

import click
import logging
import sys 
import os
import os.path as op

import pandas as pd

from pathlib import Path

from . import __version__
from .core import (
    AlleleTable, 
    ClusterTable, 
    CountRE, 
    PairTable
)
from .utilities import run_cmd

logger = logging.getLogger("cphasing")

banner = """
   ____      ____  _               _             
  / ___|    |  _ \| |__   __ _ ___(_)_ __   __ _ 
 | |   _____| |_) | '_ \ / _` / __| | '_ \ / _` |
 | |__|_____|  __/| | | | (_| \__ \ | | | | (_| |
  \____|    |_|   |_| |_|\__,_|___/_|_| |_|\__, |
                                           |___/ 
"""
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
        LOG_LEVELS = [logging.CRITICAL, logging.ERROR, 
                        logging.WARNING, logging.INFO, 
                        logging.DEBUG]
        offset = 2 
        idx = min(len(LOG_LEVELS) - 1, offset + verbose)
        logger.setLevel(LOG_LEVELS[idx])
    else:
        logger.setLevel(logging.INFO)

@cli.command()
@click.argument(
    'fasta'
)
@click.argument(
    'data'
)
@click.option(
    '--method',
    type=click.Choice(['diploid', 'polyploid'])
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
def pipeline(fasta, data, method):
    """
    Developing function.
    """
    pass


@cli.group(cls=CommandGroup, short_help='Process Pore-C alignments.')
@click.pass_context
def alignments(ctx):
    pass

@alignments.command()
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
    show_default=True
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

@alignments.command()
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
def paf2table(paf, output, threads, 
                min_quality, min_identity, min_length):
    """
    Convert paf to pore_c_table.
    """
    from .core import PAFTable

    paf = PAFTable(paf, threads=threads, 
                    min_quality=min_quality, 
                    min_identity=min_identity, 
                    min_length=min_length)
    paf.filter()
    pore_c_table = paf.to_pore_c_table()

    pore_c_table.save(output, paf.tmpdir)
    paf.clean_tempoary()

@alignments.command()
@click.argument(
    "pore_c_table",
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
    
@alignments.command()
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
def pore_c_chrom2contig(
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
    df = pd.read_parquet(pore_c_table)
    #df = df.query('pass_filter == "True"')

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


@cli.command()
@click.argument(
    'fasta',
    type=click.Path(exists=True)
)
@click.argument(
    'pairs',
    type=click.Path(exists=True)
)
@click.option(
    '-p', 
    '--percent',
    help='Percent of the map to saturate.',
    type=float,
    default=.95,
    show_default=True
)
@click.option(
    '-s', 
    '--sensitive',
    help='Sensitivity to depletion score.',
    type=click.FloatRange(0, 1, clamp=True),
    default=.5,
    show_default=True
)
@click.option(
    '-r',
    '--resolutions',
    help='Comma-separated list of target resolutions. ',
    default='20000,10000,5000,2500',
    show_default=True,
)
@click.option(
    '-d',
    '--depletion',
    help='The number of bins to calculate the sat score.',
    type=int,
    metavar='INT',
    default=4,
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
@click.option(
    '-f',
    '--force',
    help='Remove exisiting contact files.',
    default=False,
    show_default=True,
    is_flag=True
)
@click.option(
    '-o',
    '--output',
    help='Output path of corrected fatsa.',
    metavar='FILE',
    default='corrected.fasta'
)
@click.option(
    '-ob',
    '--outbed',
    help='Output path of corrected positions.',
    metavar='FILE',
    default=None
)
@click.option(
    '-op',
    '--outpairs',
    help='Output path of corrected pairs.',
    metavar='FILE',
    default=None
)
def correct(
    fasta, 
    pairs, 
    percent,
    sensitive,
    resolutions,
    depletion,
    threads,
    force,
    output,
    outbed,
    outpairs
):
    """
    Correct chimeric contigs by contacts.

        Fasta : Path of raw assembly.

        Pairs : Path of contact pairs file.

    """
    from .correct import Corrector
    c = Corrector(fasta, pairs, 
                    threads=threads,
                    force=force,
                    percent=percent,
                    sensitive=sensitive,
                    resolutions=resolutions,
                    depletion=depletion,
                    output=output,
                    outbed=outbed,
                    outpairs=outpairs)
    c.run()


@cli.command()
@click.option(
    "-f",
    "--fasta",
    metavar="FILE",
    help="polyploid contig-level fasta.",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    "--method",
    help="method for allele table constructions.",
    default='similarity',
    show_default=True,
    type=click.Choice(['gene', 'similarity'])#, 'synteny']),
)
@click.option(
    "-c",
    "--cds",
    help="the cds file of reference.",
    metavar="FILE",
    type=click.Path(exists=True)
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
    "-m",
    "minimum_similarity",
    help="minimum k-mer similarity for similarity calculation.",
    metavar="FLOAT",
    type=float,
    default=.8,
    show_default=True
)
@click.option(
    "-b",
    "--bed",
    help="the four columns bed file of reference, "
    "(chrom, start, end, gene)",
    metavar="FILE",
    type=click.Path(exists=True)
)
@click.option(
    "-p",
    "--ploidy",
    help="ploidy of genome.",
    metavar="INT",
    type=int,
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads. Only use for gmap method.',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
@click.option(
    "--skip_gmap_index",
    help="gmap index already existed and named `DB`, skip.",
    default=False,
    is_flag=True,
    show_default=True
)
def alleles(fasta, method, 
                kmer_size, window_size, minimum_similarity,
                cds, bed, ploidy, skip_gmap_index, threads):
    """
    Build allele table for prune.
    """
    
    from .alleles import GmapAllele, PartigAllele
    
    if method == 'gene':
        ga = GmapAllele(fasta, cds, bed, ploidy, 
                        skip_index=skip_gmap_index,
                        threads=threads)
        ga.run()
    elif method == 'similarity':
        prefix = op.basename(fasta).rsplit(".", 1)[0] 
        output = f"{prefix}.allele.table"
        pa = PartigAllele(fasta, kmer_size, window_size, 
                                minimum_similarity, output)
        pa.run()
    
    else:
        logger.warning('Incompletely developing function')


@cli.command()
@click.argument(
    'infile',
    metavar='InFile',
    type=click.Path(exists=True)

)
@click.argument(
    'fastafile',
    metavar='FastaFile',
    type=click.Path(exists=True)
)
@click.option(
    '-e',
    '--enzyme',
    metavar='ENZYME',
    help='Restriction enzyme name, e.g. MboI, HindIII, Arima.',
    default='MboI',
    show_default=True,
)
@click.option(
    '--minLinks',
    'minLinks',
    help='Minimum number of links for contig pair.',
    metavar='INT',
    default=3,
    type=int,
    show_default=True
)
def extract(infile, fastafile, enzyme, minLinks):
    """
    Extract countRE and pair table from 4DN pairs file.

    InFile : Path of 4DN pairs  file (Accepts compressed files).

    FastaFile : Path of fasta file.
    """
    from .utilities import restriction_site
    cmd = ['allhic', 'extract', infile, fastafile, 
            '--RE', ",".join(restriction_site(enzyme)), 
            '--minLinks', str(minLinks)]
    print(restriction_site(enzyme))
    run_cmd(cmd)

def normalize():
    pass


@cli.command()
@click.argument(
    'alleletable',
    metavar='AlleleTable',
    type=click.Path(exists=True)
)
@click.argument(
    'count_re',
    metavar='CountRE',
    type=click.Path(exists=True)
)
@click.argument(
    'pairtable',
    metavar='PairTable',
    type=click.Path(exists=True)
)
@click.option(
    '-n',
    '--normalize',
    help='If normalize the contacts.',
    is_flag=True,
    default=False,
    show_default=True
)
def prune(
    alleletable, 
    count_re,
    pairtable,
    normalize
):
    """
    Prune allelic signal by allele table.

        AlleleTable : Path to allele table.

        CountRE : Path to countRE file.

        PairTable : Path to allhic pairs table.

    """
    from .prune import Prune

    at = AlleleTable(alleletable, sort=False)
    cr = CountRE(count_re)
    pt = PairTable(pairtable, symmetric=False)

    Prune(at, cr, pt, normalize)
    pt.save(pt.filename.replace(".txt", ".prune.txt"))

@cli.command()
@click.argument(
    'alleletable',
    type=click.Path(exists=True),
    metavar='AlleleTable'
)
@click.argument(
    'count_re',
    metavar='CountRE',
    type=click.Path(exists=True)
)
@click.argument(
    'pairtable',
    metavar='PairTable',
    type=click.Path(exists=True)
)
@click.option(
    '-f',
    '--fasta',
    metavar='Fasta',
    help='Input fasta file to split by allele table.',
    default=None,
    show_default=True,
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--outdir",
    metavar='STR',
    help="Output directory of pregroup results.",
    default="outdir",
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
def pregroup(alleletable, count_re, pairtable, fasta, outdir, threads):
    """
    Pregroup countRE and pairs by homologous groups.
        

        AlleleTable : Path to allele table.

        CountRE : Path to countRE file.

        PairTable : Path to pair table.
    """
    from .pregroup import pregroup
    at = AlleleTable(alleletable, sort=False)
    cr = CountRE(count_re)
    pt = PairTable(pairtable)
    pregroup(at, cr, pt, fasta, outdir, threads)

@cli.command()
@click.argument(
    'count_re',
    metavar='CountRE',
    type=click.Path(exists=True)
)
@click.argument(
    'pairtable',
    metavar='PairTable',
    type=click.Path(exists=True)
)
@click.argument(
    'k',
    metavar='K',
    type=int
)
@click.option(
    '--maxLinkDensity',
    'maxLinkDensity',
    help='Density threshold before marking contig as repetive.'
            'If use adaptive mode must input a tuple of range.',
    default=2,
    type=str,
    show_default=True
)
@click.option(
    '--minREs',
    'minREs',
    help='Minimum number of RE sites in a contig to be clustered. '
            'If use adaptive mode must input a tuple of range.',
    default=10,
    type=str,
    show_default=True
)
@click.option(
    '--nonInformativeRatio',
    'nonInformativeRatio',
    help='Cutoff for recovering skipped contigs back into the clusters.',
    default=3,
    type=str,
    show_default=True
)
@click.option(
    '--adaptive',
    help='Adaptively looking for the best results',
    default=False,
    is_flag=True,
    show_default=True
    
)
@click.option(
    '--maxLinkDensity_step',
    'maxLinkDensity_step',
    metavar=True,
    help='The step of maxLinkDensity, only for adaptive mode.',
    type=int,
    default=2,
    show_default=True
)
@click.option(
    '--minREs_step',
    'minREs_step',
    metavar=True,
    help='The step of minREs, only for adaptive mode.',
    type=int,
    default=2,
    show_default=True
)
@click.option(
    '-t',
    '--threads',
    help='Number of threads, only for adaptive mode.',
    type=int,
    default=5,
    metavar='INT',
    show_default=True,
)
def partition(
    count_re, 
    pairtable, 
    k, 
    maxLinkDensity,
    minREs,
    nonInformativeRatio,
    adaptive,
    maxLinkDensity_step,
    minREs_step,
    threads):
    """
    Separate contigs into k groups.
        
        CountRE : Path to countRE file.

        PairTable : Path to pair table.

        K : Number or partitions. 

    """
    from .partition import Partitioner, AdaptivePartitioner
    
    if adaptive is False:
        
        Partitioner.partition(count_re, 
                                pairtable, 
                                k, 
                                maxLinkDensity, 
                                minREs, 
                                nonInformativeRatio)
    else:
        try: 
            minREs_tuple = eval(minREs)
        except:
            minREs_tuple = (50, 300)

        assert isinstance(minREs_tuple, tuple), \
                    "minREs must a tuple of range."
        
        try:
            maxLinkDensity_tuple = eval(maxLinkDensity)
        except:
            maxLinkDensity_tuple = (2, 10)

        assert isinstance(maxLinkDensity_tuple, tuple), \
                    "minREs must a tuple of range."

        ap = AdaptivePartitioner(count_re, 
                                pairtable, 
                                k,
                                maxLinkDensity_tuple,
                                maxLinkDensity_step,
                                minREs_tuple,
                                minREs_step,
                                threads=threads)
        ap.run()


@cli.command()
@click.argument(
    'clustertable',
    metavar='ClusterTable',
    type=click.Path(exists=True)
)
@click.argument(
    'alleletable',
    metavar='AlleleTable',
    type=click.Path(exists=True)
)
@click.argument(
    'count_re',
    metavar='CountRE',
    type=click.Path(exists=True)
)
@click.argument(
    'pairtable',
    metavar='PairTable',
    type=click.Path(exists=True)
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
)
def recluster(
    clustertable, 
    alleletable, 
    count_re,
    pairtable,
    ploidy,
    method
):
    """
    Recluster partition results by allele table.

        ClusterTable : Path to cluster table.

        AlleleTable : Path to allele table.

        CountRE : Path to countRE file.

        PairTable : Path to pair table.

        Ploidy : Ploidy.

    """
    from .recluster import reCluster

    rc = reCluster(clustertable, alleletable, ploidy, count_re, pairtable)
    rc.run(method=method)

@cli.command()
@click.argument(
    'clustertable',
    metavar='ClusterTable',
    type=click.Path(exists=True)
)
@click.argument(
    'count_re',
    metavar='CountRE',
    type=click.Path(exists=True)
)
@click.argument(
    'pairtable',
    metavar='PairTable',
    type=click.Path(exists=True)
)
@click.option(
    '-e',
    '--exclude',
    help='Groups list to exclude rescuing.',
    multiple=True,
    default=None,
    show_default=True
)
@click.option(
    '-m',
    '--min-score',
    'min_score',
    help='Minimum score for rescuing contigs.',
    type=float,
    metavar='FLOAT',
    default=0.01,
    show_default=True
)
@click.option(
    "-o",
    "--output",
    help="Output of results [default: stdout]",
    type=click.File('w'),
    default=sys.stdout
)
def rescue(
    clustertable, 
    count_re, 
    pairtable,
    exclude,
    min_score,
    output
    ):
    """
    Rescue uncluster contigs into already groups.

        ClusterTable : Path to cluster table.

        CountRE : Path to countRE file.

        PairTable : Path to pair table.
    """
    from .rescue import Rescuer

    r = Rescuer(clustertable, 
                count_re, 
                pairtable, 
                exclude, 
                min_score
                )
    r.rescue()
    r.save(output)

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

@cli.command()
@click.option(
    '-m',
    '--matrix',
    metavar='COOL',
    required=True
)
@click.option(
    '-a',
    '--agp',
    metavar='AGP'
)
@click.option(
    '--factor',
    '-k',
    help='Factor of plot matrix. '
            'If you input 10k matrix and want to plot heatmap at 500k,'
            'factor should be set 50.',
    type=int,
    default=50
)
@click.option(
    '--chromosomes',
    help='Chromosomes and order in which the chromosomes should be plotted. '
            'Comma seperated.',
    default=''
)
@click.option(
    '--per-chromosomes',
    'per_chromosomes',
    help='Instead of plotting the whole matrix, '
            'each chromosome is plotted next to the other.',
    is_flag=True,
    default=False
)
@click.option(
    '--no-adjust',
    'no_adjust',
    help='Matrix is chromosome-level, no need adjust.',
    default=False,
    is_flag=True,
    show_default=True,
)
@click.option(
    '--only-adjust',
    'only_adjust',
    help='Only adjust the matrix by agp.',
    default=False,
    is_flag=True,
    show_default=True,
)
@click.option(
    '--no-coarsen',
    'no_coarsen',
    help='The resolution of matrix is already for plotting, no need coarsen.',
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
)
@click.option(
    '-o',
    '--output',
    help='Output path of file.',
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
@click.option(
    '--dpi',
    help='Resolution for the image.',
    default=600,
    show_default=True
)
@click.option(
    '--cmap',
    help='Colormap of heatmap.',
    default='YlOrRd',
    show_default=True
)
def plot(matrix, agp, factor,
        chromosomes, per_chromosomes, 
        no_adjust, only_adjust, no_coarsen,
        only_plot, output,
        threads, dpi, cmap):
    """
    Adjust or Plot the contacts matrix after assembling.
    """
    from .plot import (
        adjust_matrix, 
        plot_matrix,
        coarsen_matrix,
    )
    if not only_plot:
        if not no_adjust:
            assert agp is not None, \
                "Must provide agp file for matrix adjustment."
            matrix = adjust_matrix(matrix, agp)

        if only_adjust:
            sys.exit()

        if not no_coarsen:
            matrix = coarsen_matrix(matrix, factor, None, threads)   
    
    chromosomes = chromosomes.strip().strip(",").split(',')

    plot_matrix(matrix, 
                output,
                chromosomes,
                per_chromosomes,
                dpi=dpi,
                cmap=cmap)


def assess():
    pass

@cli.group(cls=CommandGroup, short_help='Sub-command for Hi-C pipeline.')
@click.pass_context
def hic(ctx):
    pass

@hic.command()
@click.option(
    '-r',
    '--reference',
    help='Path of reference fasta file.',
    metavar='FILE',
    #type=click.Path(exists=True),
    required=True
)
@click.option(
    '-1',
    '--read1',
    help='Path of read 1.',
    metavar='FILE',
    type=click.Path(exists=True),
    multiple=True,
    required=True
)
@click.option(
    '-2',
    '--read2',
    help='Path of read 2.',
    metavar='FILE',
    type=click.Path(exists=True),
    multiple=True,
    required=True
    )
@click.option(
    '-e',
    '--enzyme',
    metavar='ENZYME',
    help='Restiction enzyme. i.e. MboI, HindIII.',
    default=None,
    show_default=True
)
@click.option(
    '-q',
    '--mapq',
    help='Minimum quality of mapping [0, 60].',
    metavar='INT',
    type=click.IntRange(0, 60, clamp=True),
    default=1,
    show_default=True
)
@click.option(
    '--aligner',
    help='Aligner executable.',
    type=click.Choice(['chromap']),#, 'hisat2']),
    default='chromap',
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
def mapper(
    reference,
    read1, 
    read2,
    enzyme,
    mapq,
    aligner,
    threads

):
    """
    Mapper for reads mapping.
    """
    assert len(read1) == len(read2), "reads must paired."

    if aligner == 'hisat2':
        if enzyme is None:
            logger.error('Missing option "-e"')
            sys.exit()
        from .mapper import HisatMapper
        hm1 = HisatMapper(reference, read1, enzyme, 
                            min_quality=mapq, threads=threads)
        hm1.run()
        hm2 = HisatMapper(reference, read2, enzyme, 
                            min_quality=mapq, threads=threads)
        hm2.run()
        
    elif aligner == 'chromap':
        from .mapper import ChromapMapper
        from .core import PairHeader
        res = []
        for r1, r2 in zip(read1, read2):
            cm = ChromapMapper(reference, r1, r2, 
                                min_quality=mapq, 
                                threads=threads)
            cm.run()
            res.append(cm.output_pairs)
        else:
            if len(res) > 1:
                header = PairHeader([]).from_file(res[0])
                header.save("temp.pairs.header")

                cmd = ['cat', 'temp.pairs.header'] + res
                cmd2 = ['LC_ALL=C', 'grep', '-v', '#']


@cli.group(cls=CommandGroup, short_help='Misc tools.')
@click.pass_context
def utils(ctx):
    pass

@utils.command(short_help='Convert agp to assembly file.')
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


@utils.command(short_help='Convert agp to cluster file.')
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

@utils.command(short_help='Convert agp to fasta file.')
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
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,

)
def agp2fasta(agpfile, fasta, output, threads):
    from .agp import agp2fasta
    agp2fasta(agpfile, fasta, output)

@utils.command()
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

@utils.command(short_help='Convert cluster to pesudo assembly file.')
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
    Convert cluster to a pesudo assembly file.

    ClusterTable : Path to cluster table.

    Fasta : Path to fasta file.
    
    """
    from .core import ClusterTable
    ct = ClusterTable(cluster)
 
    ct.to_assembly(fasta, output)

@utils.command()
@click.argument(
    "count_re"
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

@utils.command()
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

@utils.command()
@click.argument(
    "pairs",
    metavar="Pairs",
    type=click.Path(exists=True)
)
@click.argument(
    "output",
    metavar="Output"
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
def pairs2mnd(pairs, output, threads):
    """
    convert 4DN pairs to mnd file.

        Pairs : Path of pairs file.

        Output : Path of output mnd file.
    """
    from .core import Pairs 
    p = Pairs(pairs)
    p.to_mnd(output, threads) 

@utils.command()
@click.argument(
    "pore_c_tables",
    metavar="Pore_C_Tables",
    type=click.Path(exists=True)
)
@click.argument(
    "output",
    metavar="OUTPUT",
)
@click.option(
    "-k",
    help="Number of groups.",
    default=None,
    show_default=True,
)
@click.option(
    "--fasta",
    help="Path of fasta file.",
    default=None,
    show_default=True,
    type=click.Path(exists=True),
    required=True
)
@click.option(
    "--prune",
    metavar="STR",
    help="prune list from prune.",
    default=None,
    show_default=True,
    type=click.Path(exists=True)
)
@click.option(
    "--min-order",
    "min_order",
    help="Minimum contig order of pore-c reads",
    metavar="INT",
    type=int,
    default=2,
    show_default=True
)
@click.option(
    "--max-order",
    "max_order",
    help="Maximum contig order of pore-c reads",
    metavar="INT",
    type=int,
    default=15,
    show_default=True
)
@click.option(
    "--min-alignments",
    "min_alignments",
    help="Minimum length of alignments",
    metavar="INT",
    type=int,
    default=500,
    show_default=True
)
@click.option(
    "--min-length",
    "min_length",
    help="Minimum length of contigs",
    metavar="INT",
    type=int,
    default=10000,
    show_default=True
)
@click.option(
    "--threshold",
    metavar="FLOAT",
    help="Threshold of reweight",
    type=float,
    default=0.01,
    show_default=True
)
@click.option(
    "--max-round",
    "max_round",
    help="Maximize round of reweight",
    metavar="INT",
    type=int,
    default=50,
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
    '-t',
    '--threads',
    help='Number of threads.',
    type=int,
    default=4,
    metavar='INT',
    show_default=True,
)
def hyperpartition(pore_c_tables, output,
                    k, fasta, prune,
                    min_order, max_order,
                    min_alignments,
                    min_length, threshold, 
                    max_round, fofn, 
                    threads):
    """
    Separate contigs into several groups by hypergraph cluster.
        
        Pore_C_TABLE : Path of Pore-C table.

        OUTPUT: Path of output clusters.
    """
    from .partition import HyperPartition
    if fofn:
        pore_c_tables = [i.strip() for i in open(pore_c_tables)]

    print((pore_c_tables, fasta,
                            k,
                            prune, min_order, 
                            max_order, min_alignments, 
                            min_length,
                            threshold, max_round, 
                            threads))
    hp = HyperPartition(pore_c_tables, 
                            k, fasta,
                            prune, min_order, 
                            max_order, min_alignments, 
                            min_length,
                            threshold, max_round, 
                            threads)

    if not k:
        hp.single_partition()
    if k:
        if len(k.split(":")) > 1:
            hp.multi_partition()
        else:
            hp.single_partition()
            
    hp.to_cluster(output)