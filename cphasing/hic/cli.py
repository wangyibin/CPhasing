#!/usr/bin/env python

"""
cli for hic pipelines
"""

import click
import logging
import sys


from ..cli import CommandGroup
from ..cli import cli 
from ..core import (
    AlleleTable, 
    ClusterTable, 
    CountRE, 
    PairTable
)

from ..utilities import run_cmd

logger = logging.getLogger(__name__)

@cli.group(cls=CommandGroup, short_help='Sub-command for the legacy Hi-C pipeline.')
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
    help='Restiction enzyme. i.e. MboI, HindIII., (Deprecated)',
    default=None,
    show_default=True,
    hidden=True
)
@click.option(
    "-k",
    "kmer_size",
    help="kmer size for mapping.",
    metavar="INT",
    type=int,
    default=17,
    show_default=True
)
@click.option(
    "-w",
    "window_size",
    help="minimizer window size for mapping.",
    metavar="INT",
    type=int,
    default=7,
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
    '-a',
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
    kmer_size,
    window_size,
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
        from ..mapper import HisatMapper
        hm1 = HisatMapper(reference, read1, enzyme, 
                            min_quality=mapq, threads=threads)
        hm1.run()
        hm2 = HisatMapper(reference, read2, enzyme, 
                            min_quality=mapq, threads=threads)
        hm2.run()
        
    elif aligner == 'chromap':
        from ..mapper import ChromapMapper
        from ..core import PairHeader
        res = []
        for r1, r2 in zip(read1, read2):
            cm = ChromapMapper(reference, r1, r2, 
                                kmer_size=kmer_size,
                                window_size=window_size,
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


@hic.command()
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


@hic.command()
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
    help="path of output allele table [default: fasta_prefix]",
    default=None
)
@click.option(
    "-c",
    "--cds",
    help="the cds file of reference.",
    metavar="FILE",
    type=click.Path(exists=True)
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
    default=None,
    show_default=True,
)
@click.option(
    "--skip_gmap_index",
    help="gmap index already existed and named `DB`, skip.",
    default=False,
    is_flag=True,
    show_default=True
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
def alleles(fasta, output, cds, bed, ploidy, skip_gmap_index, threads):
    """
    Build allele table by gmap
    """
    from ..alleles import GmapAllele 

    assert ploidy is not None, "-p parameter should be provide"
    ga = GmapAllele(fasta, cds, bed, ploidy, 
                        skip_index=skip_gmap_index,
                        threads=threads)

    
@hic.command()
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
    from ..utilities import restriction_site
    cmd = ['allhic', 'extract', infile, fastafile, 
            '--RE', ",".join(restriction_site(enzyme)), 
            '--minLinks', str(minLinks)]
    print(restriction_site(enzyme))
    run_cmd(cmd)

@hic.command()
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

@hic.command()
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


@hic.command()
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

@hic.command()
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


@hic.command()
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