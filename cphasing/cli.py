#!/usr/bin/env python

import click
import logging
import sys 
import os
import os.path as op

import numpy as np
import pandas as pd

from collections import defaultdict
from pathlib import Path
from pytools import natsorted
from shutil import which

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
@click.option(
    '-f',
    '--fasta',
    metavar="FASTA",
    help="Path to draft assembly",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '-pcd',
    '--porec-data',
    'porec_data',
    metavar="PoreC_Data",
    help="",
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
    '-p',
    '--pattern',
    help='Pattern of restrict enzyme. Commam separate for multiple pattern',
    metavar="STR",
    default="AAGCTT",
    show_default=True,
)
@click.option(
    "-hcr",
    "--hcr",
    is_flag=True,
    default=False,
    help="Only retain high confidence regions to subsequence analysis.",
    show_default=True
)
@click.option(
    '--mode',
    metavar="STR",
    help="mode of hyperpartition",
    type=click.Choice(['basal', 'phasing', 'basal_withprune']),
    default='phasing',
    show_default=True,
)
@click.option(
    '-s',
    '--steps',
    metavar='STR',
    help="steps",
    default="1,2,3,4",
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
    default=1.0,
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
    '-fc',
    '--first-cluster',
    'first_cluster',
    metavar="PATH",
    help='Use the existing first cluster results to second round cluster',
    default=None,
    show_default=True
)
@click.option(
    "-as",
    "--allelic-similarity",
    "allelic_similarity",
    metavar="FLOAT",
    help="The similarity of allelic",
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
    '--factor',
    '-k',
    help='Factor of plot matrix. '
            'If you input 10k matrix and want to plot heatmap at 500k, '
            'factor should be set with 50.',
    type=int,
    default=50,
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
            porec_data, 
            porectable, 
            pairs, 
            pattern,
            hcr,
            mode, 
            steps,
            skip_steps, 
            n,
            resolution1,
            resolution2, 
            first_cluster,
            allelic_similarity,
            min_allelic_overlap,
            whitelist,
            factor,
            threads):
    """
    A pipeline of polyploid phaseing and scaffolding\n
    Steps:\n
        0. mapper;\t\t\t\t\t\t\t\t\t
        1. alleles;\t\t\t\t\t\t\t\t\t 
        2. prepare;\t\t\t\t\t\t\t\t\t
        3. hyperpartiton;\t\t\t\t\t\t\t\t\t
        4. scaffolding;\t\t\t\t\t\t\t\t\t
        5. plot;\t\t\t\t\t\t\t\t\t

    Usages:\n
        Input pore-c data\t\t\t\t\t\t\t\t\t
        $ cphasing pipeline -f contigs.fasta -pcd sample.fastq.gz -t 10 -s "all"
        Input pore-c table\t\t\t\t\t\t\t\t\t
        $ cphasing pipeline -f contigs.fasta -pct sample.porec.gz -t 10
        Input pairs\t\t\t\t\t\t\t\t\t
        $ cphasing pipeline -f contigs.fasta -prs sample.pairs.gz -t 10 

        Skip step:\t\t\t\t\t\t\t\t\t
        $ cphasing pipeline -f contigs.fasta -pct sample.porec.gz -t 10 -ss 1 
        Only run a step:\t\t\t\t\t\t\t\t\t
        $ cphasing pipeline -f contigs.fasta -pct sample.porec.gz -t 10 -s 1

        Basal mode:\t\t\t\t\t\t\t\t\t
        $ cphasing pipeline -f contigs.fasta -pct sample.porec.gz -t 10 --mode basal


        
    """
    from .pipeline.pipeline import run 
    
    assert any([(porec_data is not None), (porectable is not None), (pairs is not None)]), \
        "PoreC data or PoreC table or Pairs must specified"


    if steps:
        if steps == "all":
            steps = set(["0", "1", "2", "3", "4", "5"])

        else:
            steps = steps.strip().split(",")
    else:
        steps = set(map(str, [0, 1, 2, 3, 4]))



    if skip_steps:
        skip_steps = set(skip_steps.strip().split(","))
    else:
        skip_steps = set()
    
    run(fasta, porec_data,
         porectable, pairs, 
         pattern=pattern,
         hcr_flag=hcr,
        mode=mode, 
        steps=steps,
        skip_steps=skip_steps,
        n=n,
        resolution1=resolution1,
        resolution2=resolution2,
        first_cluster=first_cluster,
        whitelist=whitelist,
        allelic_similarity=allelic_similarity,
        min_allelic_overlap=min_allelic_overlap,
        factor=factor,
        threads=threads)
    

## Subcommand of UL ONT pipeline
from .hitig.cli import hitig 

@cli.command()
@click.argument(
    "reference",
    type=click.Path(exists=True)
)
@click.argument(
    "fastq",
    type=click.Path(exists=True)
)
@click.option(
    "-p",
    "--pattern",
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
    default=1,
    show_default=True
)
@click.option(
    '--realign',
    help="realign to rescue multiple alignments, only support for contig-level",
    is_flag=True,
    default=False,
    show_default=True
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
def mapper(reference, fastq, pattern, kmer_size, 
            window_size, mm2_params, mapq, force, 
            realign, outprefix, threads):
    """
    mapper for pore-c reads.

        REFERENCE: Path of reference

        FASTQ: Path of pore-c reads

    """
    from .mapper import PoreCMapper

    additional_arguments = mm2_params.strip().split()
 
    pcm = PoreCMapper(reference, fastq,
                        pattern=pattern,
                        k = kmer_size,
                        w = window_size,
                        force=force,
                        realign=realign,
                        min_quality=mapq,
                        additional_arguments=additional_arguments,
                        outprefix=outprefix,
                        threads=threads)
    pcm.run()


@cli.group(cls=CommandGroup, short_help='Process Pore-C alignments.')
@click.pass_context
def alignments(ctx):
    pass

@alignments.command(hidden=True, deprecated=True)
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


@alignments.command(hidden=True, deprecated=True)
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
    print(paf, output, 
                min_quality, min_identity, min_length, 
                chunksize, threads)
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

@alignments.command(hidden=True, deprecated=True)
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
    
@alignments.command(hidden=True)
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

@alignments.command(hidden=True)
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


@alignments.command()
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
    hidden=True
)
def pairs2mnd(pairs, output, threads):
    """
    convert 4DN pairs to mnd file.

        Pairs : Path of pairs file.

        Output : Path of output mnd file.
    """
    # from .core import Pairs 
    # p = Pairs(pairs)
    # p.to_mnd(output, threads) 
    cmd = ["cphasing-rs", "pairs2mnd", pairs, "-o", output]
    flag = run_cmd(cmd, log=os.devnull)
    assert flag == 0, "Failed to execute command, please check log."


@alignments.command()
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

@alignments.command()
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

@alignments.command()
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

@cli.command(short_help='Only retain the HCRs from Pore-C/Hi-C data')
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
    "-bs",
    "--binsize",
    help="Bin size in bp.",
    type=int,
    default=10000,
    show_default=True
)
@click.option(
    "-p",
    "--percent",
    help="Percentile of high confidence regions.",
    type=click.IntRange(0, 100),
    default=95,
    show_default=True,
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
def hcr(porectable, pairs, contigsize, binsize, percent, 
        fofn, output):
    """
    Only retain high confidence regions to subsequence analysis.
        High confidence regions are identified from contacts.

    """
    from .hitig.hcr.hcr_by_contacts import hcr_by_contacts

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
                prefix = Path(pore_c_table).with_suffix("")
                while prefix.suffix in {'.gz', 'gz', 'porec', ".porec"}:
                    prefix = prefix.with_suffix('')
                
                pairs = f"{prefix}.pairs.gz"
                pairs_files.append(pairs)
                if not Path(pairs).exists():
                    cmd = ["cphasing-rs", "porec2pairs", porectable, contigsize,
                    "-o", pairs]
                
                    flag = run_cmd(cmd, log=f"logs/porec2pairs.log")
                    assert flag == 0, "Failed to execute command, please check log."
                else:
                    logger.warn(f"Use exists 4DN pairs file of `{pairs}`.")
        
    else:
        if fofn:
            pairs_files = [i.strip() for i in open(pairs) if i.strip()]
        else:
            pairs_files = [pairs]
        
        for pairs in pairs_files:
            prefix = Path(pairs).with_suffix("")
            while prefix.suffix in {'.gz', 'gz', '.pairs'}:
                prefix = prefix.with_suffix('')


    out_small_cool = f"{prefix}.{binsize}.cool"
    if not Path(out_small_cool).exists():
        try: 
            pairs2cool.main(args=[
                                pairs,
                                contigsize,
                                out_small_cool
                            ],
                            prog_name='pairs2cool')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
        
    else:
        logger.warn(f"Use exists cool file of `{out_small_cool}`.")
    
    hcr_by_contacts(out_small_cool, f"{prefix}.{binsize}.hcr.bed" , percent)

    if porectable:
        cmd = ["cphasing-rs", "porec-intersect", porectable, f"{prefix}.{binsize}.hcr.bed",
               "-o", f"{prefix}_hcr.porec.gz"]
        flag = run_cmd(cmd, log=f"logs/hcr_intersect.log")
        assert flag == 0, "Failed to execute command, please check log."
        logger.info(f'Successful out high confidence high-orrder '
                    f'contacts into `{f"{prefix}_hcr.porec.gz"}`')
        cmd = ["cphasing-rs", "porec2pairs", f"{prefix}_hcr.porec.gz", contigsize,
               "-o", f"{prefix}_hcr.pairs.gz", "-q", "1"]
        
        flag = run_cmd(cmd, log=f"logs/hcr_porec2pairs.log")
        assert flag == 0, "Failed to execute command, please check log."
        logger.info(f'Successful out high confidence contacts into `{f"{prefix}_hcr.pairs.gz"}`')

    else:
        cmd = ["cphasing-rs", "pairs-intersect", pairs, f"{prefix}.{binsize}.hcr.bed",
               "-o", f"{prefix}_hcr.pairs.gz"]
        flag = run_cmd(cmd, log=f"logs/hcr_intersect.log")
        assert flag == 0, "Failed to execute command, please check log."
        logger.info(f'Successful out high confidence contacts into `{f"{prefix}_hcr.pairs.gz"}`')


@cli.command(short_help='Prepare data for subsequence analysis.')
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
    help='Pattern of restrict enzyme. Commam separate for multiple pattern',
    metavar="STR",
    default="AAGCTT",
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
def prepare(fasta, pairs, min_contacts, pattern, threads, outprefix):
    """

    """
    from .prepare import pipe 

    pipe(fasta, pairs, pattern, min_contacts, threads, outprefix)
    pass 


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
    "-o",
    "--output",
    metavar="PATH",
    help="path of output allele table [default: fasta_prefix]",
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
    "-m",
    "minimum_similarity",
    help="minimum k-mer similarity for similarity calculation.",
    metavar="FLOAT",
    type=float,
    default=.2,
    show_default=True
)
def alleles(fasta, output, 
                kmer_size, window_size, minimum_similarity):
    """
    Build allele table.
    """
    
    from .alleles import PartigAllele

    if not output:
        fasta_prefix = Path(fasta).with_suffix("")
        while fasta_prefix.suffix in {".fasta", "gz", "fa", '.fa', '.gz'}:
            fasta_prefix = fasta_prefix.with_suffix("")
    
        output = f"{fasta_prefix}.allele.table"
    
    pa = PartigAllele(fasta, kmer_size, window_size, 
                            minimum_similarity, output)
    pa.run()


@cli.command(short_help="Generate the allelic contig and cross-allelic contig table.")
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
    default='greedy',
    type=click.Choice(["fast", "greedy"]),
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
    help="Normalization method of contacts for kprune",
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


@cli.command()
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
            threads):
    """
    Construct hypergraph from contacts.

    The hypergraph or graph of pore-c or hic enabled 
        extract from pore-c table or 4DN pairs. 

        Contacts : Path of Pore-C table or 4DN pairs.
        
        Contig_sizes : Path of contig sizes.

        OUTPUT : Path of output hypergraph.

    """
    from .hypergraph import (
        HyperExtractor, 
        HyperExtractorSplit,
        Extractor
        )
    from .utilities import read_chrom_sizes

    contigsizes = read_chrom_sizes(contigsize)
    contigs = natsorted(contigsizes.index.values.tolist())

    if whitelist:
        whitelist = set([i.strip() for i in open(whitelist) if i.strip()])
        contigs = list(filter(lambda x: x in whitelist, contigs))

    
    contig_idx = defaultdict(None, dict(zip(contigs, range(len(contigs)))))

    if not pairs:
        if fofn:
            pore_c_tables = [i.strip() for i in open(contacts) if i.strip()]
        else:
            pore_c_tables = contacts

        if not split:
            he = HyperExtractor(pore_c_tables, contig_idx, contigsizes.to_dict()['length'], 
                                min_order, max_order, min_alignments, min_quality, threads)
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

        e = Extractor(pairs_files, contig_idx, contigsizes.to_dict()['length'], threads)
        e.save(output)


@cli.command()
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
    type=click.Path(exists=True)
)
@click.option(
    "-ul",
    "--ultra-long",
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
    default=1.0,
    type=float,
    show_default=True,
    hidden=True,
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
    "-at",
    "--alleletable",
    metavar="PATH",
    help="""
    Path to allele table that is generated from `alleles`. 
        Skip when prunetable parameter added.
    """,
    default=None,
    show_default=True,
    type=click.Path(exists=True),
    hidden=True
)
@click.option(
    "-pt",
    "--prunetable",
    metavar="PATH",
    help="Path to prune table that is generated from `kprune`.",
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
    help="The similarity of allelic",
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
    '--mode',
    metavar="STR",
    help="mode of hyperpartition, conflict of `-inc`",
    type=click.Choice(["basal", "phasing", "basal_withprune", None]),
    default=None,
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
    help="Don't run hyperpartition, only merge this cluster into k group",
    default=None,
    show_default=True,
)
@click.option(
    '-fc',
    '--first-cluster',
    'first_cluster',
    metavar="PATH",
    help='Use the existing first cluster results to second round cluster',
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
    """,
    default=None,
    show_default=True,
    type=click.Path(exists=True)
)
@click.option(
    "-mc",
    "--min-contacts",
    "min_contacts",
    help="Minimum contacts of contigs",
    metavar="INT",
    type=int,
    default=1,
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
    default=1.0,
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
    "-mw",
    "--min-weight",
    "min_weight",
    help="Minimum weight of graph",
    type=float,
    default=0.1,
    show_default=True
)
@click.option(
    "-ms",
    "--min_scaffold_length",
    help="The minimum length of the output scaffolding.",
    type=float,
    default=5e5,
    show_default=True
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
    help="Maximize round of reweight",
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
                    n,
                    alleletable,
                    prunetable,
                    normalize,
                    allelic_factor,
                    cross_allelic_factor,
                    allelic_similarity,
                    min_allelic_overlap,
                    mode,
                    incremental,
                    kprune_norm_method,
                    # ultra_complex,
                    merge_cluster,
                    first_cluster,
                    whitelist,
                    blacklist,
                    min_contacts,
                    min_length, 
                    resolution1,
                    resolution2,
                    min_weight,
                    min_scaffold_length,
                    threshold, 
                    max_round, 
                    threads,
                    # chunksize
                    ):
    """
    Separate contigs into groups based on hypergraph.

        Hypergraph : Path of the hypergraph.

        Contig_sizes : Path of contig sizes.

        Output : Path of output clusters.

    """
    import msgspec 
    from .hypergraph import HyperExtractor, Extractor
    from .hyperpartition import HyperPartition
    from .algorithms.hypergraph import HyperEdges
    from .utilities import read_chrom_sizes 
    
    assert not all([porec, pairs]), "confilct parameters, only support one type data"

    ultra_complex = None

    if mode == "basal":
        incremental = False
        if alleletable or prunetable:
            logger.warn("allelic information will not be used in basal mode")
            alleletable = None
            prunetable = None
    elif mode == "phasing":
        incremental = True 
        assert alleletable or prunetable,\
            "phasing mode must add `-pt` or `-at` param"
    elif mode == "basal_withprune":
        incremental = False
        assert alleletable or prunetable, \
            "basal_withprune modemust add `-pt` or `-at` param"
    else:
        incremental = incremental
        logger.warn("Mode not be specified, running basal_withprune mode")

    if kprune_norm_method == "auto":
        if pairs:
            kprune_norm_method = "cis"
        elif porec:
            kprune_norm_method = "none"
        else:
            kprune_norm_method = "none"
    
    if n is not None:
        if ":" in n and incremental is False:
            logger.warn("Second round partition will not be run, or `-inc` parameters must be added")

    n = n.split(":") if n else [None, None]

    for i, v in enumerate(n):
        if v:
            if i == 1:
                try:
                    n[i] = int(v)
                except ValueError:
                    if Path(n[i]).exists():
                        logger.info(f"Load 2nd group number list f`{n[i]}`")
                        tmp_df = pd.read_csv(n[i], sep='\t', header=None, index_col=None)

                        n[i] = tmp_df.to_dict()[1]

                    else:
                        assert ValueError, "n2 must be a int or a exists file"
            else:
                n[i] = int(v)
        else:
            n[i] = 0

    contigsizes = read_chrom_sizes(contigsizes)
    prefix = Path(hypergraph).stem
    if porec:
        contigs = contigsizes.index.values.tolist()
        contig_idx = defaultdict(None, dict(zip(contigs, range(len(contigs)))))
        if not Path(f"{prefix}.hg").exists():
            logger.info(f"Load raw hypergraph from porec table `{hypergraph}`")
            
            he = HyperExtractor(hypergraph, contig_idx, contigsizes.to_dict()['length'])
            he.save(f"{prefix}.hg")
            hypergraph = he.edges
        else:
            logger.warn(f"Load raw hypergraph from exists file of `{prefix}.hg`")
            hypergraph = msgspec.msgpack.decode(open(f"{prefix}.hg", 'rb').read(), type=HyperEdges)

        

    elif pairs:
        contigs = contigsizes.index.values.tolist()
        contig_idx = defaultdict(None, dict(zip(contigs, range(len(contigs)))))
        if not Path(f"{prefix}.hg").exists():
            logger.info(f"Load raw hypergraph from pairs file `{hypergraph}`")
            he = Extractor(hypergraph, contig_idx, contigsizes.to_dict()['length'])
            he.save(f"{prefix}.hg")
            hypergraph = he.edges
        else:
            logger.warn(f"Load raw hypergraph from exists file of `{prefix}.hg`")
            hypergraph = msgspec.msgpack.decode(open(f"{prefix}.hg", 'rb').read(), type=HyperEdges)
        
        
    else:
        logger.info(f"Load raw hypergraph from `{hypergraph}")
        hypergraph = msgspec.msgpack.decode(open(hypergraph, 'rb').read(), type=HyperEdges)
        

    if whitelist:
        whitelist = [i.strip() for i in open(whitelist) if i.strip()]
    if blacklist:
        blacklist = [i.strip() for i in open(blacklist) if i.strip()]

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
    

    if alleletable:
        if prunetable:
            logger.warn("The allele table will not be used, becauese prunetable parameter added.")

    if incremental is False and first_cluster is not None:
        logger.warn("First cluster only support for incremental method, will be not used.")
    
    if not prunetable and not alleletable:
        logger.info("Not inplement the allelic and cross-allelic reweight algorithm")
    
    hp = HyperPartition(hypergraph, 
                            contigsizes,
                            ultra_long,
                            ul_weight,
                            n,
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
                            whitelist,
                            blacklist,
                            min_contacts,
                            min_length,
                            resolution1, 
                            resolution2,
                            min_weight,
                            min_scaffold_length,
                            threshold, 
                            max_round, 
                            threads, 
                            # chunksize
                            )
    
  
    if incremental:
        if first_cluster and op.exists(first_cluster):
            first_cluster = ClusterTable(first_cluster)
        else:
            first_cluster = None 
        
        hp.incremental_partition(n, first_cluster)
    else:
        if merge_cluster:
            hp.merge_cluster(_merge_cluster)
        else:
            hp.single_partition(int(n[0]))
    
    hp.to_cluster(output)



@cli.command()
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
    "-sc",
    "--split-contacts",
    "split_contacts",
    metavar="SPLIT_CONTACTS",
    type=click.Path(exists=True),
    default=None,
    show_default=True
)
@click.option(
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
    "-m",
    "--method",
    metavar="STR",
    help="Method of scaffolding",
    default="haphic",
    type=click.Choice(["haphic", "allhic", "haphic_fastsort"]),
    show_default=True
)
@click.option(
    "-o",
    "--output",
    metavar="STR",
    help="Output agp file, only valid when fasta file exists",
    default="groups.agp",
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
                fasta, method,
                output, threads):
    """
    Ordering and orientation the contigs.

        ClusterTable: Path of cluster table from hyperpartition.

        Count_RE: Path of counts RE file.

        CLM: Path of clm file.

    """

    from .algorithms.scaffolding import AllhicOptimize, HapHiCSort
    
    if method == "allhic":
        ao = AllhicOptimize(clustertable, count_re, clm, allele_table=allele_table,
                                fasta=fasta, output=output, threads=threads)
        ao.run()
    elif method == "haphic":
        assert split_contacts is not None, "split_contacts file must specified by `-sc` parameters"
        assert fasta is not None, "fasta must specified by `-f` parameters"
        hs = HapHiCSort(clustertable, count_re, clm, split_contacts, 
                                allele_table=allele_table,
                                fasta=fasta, output=output, threads=threads)
        hs.run()

    else:
        assert split_contacts is not None, "split_contacts file must specified by `-sc` parameters"
        assert fasta is not None, "fasta must specified by `-f` parameters"
        hs = HapHiCSort(clustertable, count_re, clm, split_contacts, skip_allhic=True,
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


@cli.command()
@click.argument(
    "fasta",
    metavar="Fasta",
    type=click.Path(exists=True)
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
def build(fasta, output, output_agp, only_agp):
    """
    Build genome release.

    Fasta : contig-level fasta file
    """
    from .build import Build
    Build(fasta, output, 
            output_agp=output_agp, only_agp=only_agp)

@cli.command()
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
    help="Bin size in bp.",
    type=int,
    default=10000,
    show_default=True
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
    show_default=True
)
def pairs2cool(pairs, chromsize, outcool,
               binsize, fofn):
    """
    Convert pairs file into a specified resolution cool file.

        INPUT_PAIRS_PATH : Path of pairs file, can be compressed.

        CHROM_SIZE : Two columns of chromosomes or contigs size.

        OUT_COOL_PATH : Output path of cool file.
    """

    from cooler.cli.cload import pairs as cload_pairs
    from .utilities import merge_matrix

    logger.info(f"Load pairs: `{pairs}`.")
    logger.info(f"Bin size: {binsize}")
    
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
    
    logger.info(f'Output binning contact matrix into `{outcool}`')
    if fofn:
        if op.exists(f'temp.{pid}.pairs'):
            os.remove(f'temp.{pid}.pairs')
        if op.exists(f'temp.{pid}.header'):
            os.remove(f'temp.{pid}.header')

    # merge_matrix(outcool, outcool=f"{outcool.rsplit('.', 2)[0]}.whole.cool")


@cli.command()
@click.option(
    '-m',
    '--matrix',
    metavar='COOL',
    required=True,
    type=click.Path(exists=True)
)
@click.option(
    '-a',
    '--agp',
    metavar='AGP',
    type=click.Path(exists=True)
)
@click.option(
    '--factor',
    '-k',
    help='Factor of plot matrix. '
            'If you input 10k matrix and want to plot heatmap at 500k, '
            'factor should be set with 50.',
    type=int,
    default=50,
    show_default=True
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
# @click.option(
#     '-t',
#     '--threads',
#     help='Number of threads. (unused)',
#     type=int,
#     default=4,
#     metavar='INT',
#     show_default=True,
# )
@click.option(
    '-c',
    '--chromosomes',
    help='Chromosomes and order in which the chromosomes should be plotted. '
            'Comma seperated. or a one column file',
    default=''
)
@click.option(
    '-pc',
    '--per-chromosomes',
    'per_chromosomes',
    help='Instead of plotting the whole matrix, '
            'each chromosome is plotted next to the other.',
    is_flag=True,
    default=False
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
    '-dpi',
    '--dpi',
    help='Resolution for the image.',
    default=600,
    show_default=True
)
@click.option(
    '-cmap',
    '--cmap',
    help="""
    Colormap of heatmap. 
    Available values can be seen : 
    https://pratiman-91.github.io/colormaps/ 
    and http://matplotlib.org/examples/color/colormaps_reference.html
    """,
    default='redp1_r',
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
def plot(matrix, 
            agp, 
            factor,
            no_adjust, 
            only_adjust, 
            no_coarsen,
            only_plot, 
            output,
            # threads, 
            chromosomes, 
            per_chromosomes,
            chrom_per_row, 
            dpi, 
            cmap,
            no_lines,
            no_ticks,):
    """
    Adjust or Plot the contacts matrix after assembling.
    """
    from .plot import (
        adjust_matrix,
        coarsen_matrix, 
        plot_heatmap
        
    )
    threads = 1

    if agp is None:
        only_plot = True 
        logger.warning( "Only plot the matrix. "
                "If you want to adjust matrix to chromosome-level, please provide agp file. ")

    if not only_plot:
        if not no_adjust:
            if which("bedtools") is None:
                raise ValueError(f"bedtools: command not found.")
            matrix = adjust_matrix(matrix, agp)

        if only_adjust:
            sys.exit()

        if not no_coarsen and factor > 1:
            matrix = coarsen_matrix(matrix, factor, None, threads)   
    
    if op.exists(chromosomes):
        logger.info(f"Load chromosomes list from the file `{chromosomes}`.")
        chromosomes = [i.strip() for i in open(chromosomes) if i.strip()]
        
    else:
        
        chromosomes = chromosomes.strip().strip(",").split(',') if chromosomes else None
        

    if no_ticks:
        xticks = False 
        yticks = False
    else: 
        xticks = True 
        yticks = True

    plot_heatmap(matrix,
                 output,
                 chromosomes=chromosomes,
                 per_chromosomes=per_chromosomes,
                 chrom_per_row=chrom_per_row,
                 dpi=dpi,
                 cmap=cmap, 
                 xticks=xticks, 
                 yticks=yticks,
                 add_lines=False if no_lines else True,
                 threads=threads)





ALIASES = {
    "pipe": pipeline,
    "align": alignments,
    "hg": hypergraph,
    "hp": hyperpartition,
    "ho": hyperoptimize,
    "sf": scaffolding,
}

## hic subcommand
from .hic.cli import hic

@cli.group(cls=CommandGroup, short_help='Evaluation tools.', hidden=True)
@click.pass_context
def evaluate(ctx):
    pass

@evaluate.command(short_help='Calculate the allelic error rate')
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
    pass 



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
    """
    Convert agp to fasta

        AGP: Path of AGP file

        FASTA: Path of draft assembly
    """
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
def statagp(agp, output):
    """
    Statistics of AGP.

    AGP : Path to agp file.

    """
    from .agp import statagp
    
    statagp(agp, output)


@utils.command(short_help='Statistics of ClusterTable.')
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

@utils.command(short_help='Convert cluster to pseudo assembly file.')
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

@utils.command()
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
        


@utils.command()
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
    
@utils.command()
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

@utils.command(hidden=True)
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

@utils.command()
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