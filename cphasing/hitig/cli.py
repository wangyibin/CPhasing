#!/usr/bin/env python

"""
cli for ultra-long ont pipelines
"""

import rich_click as click
from rich_click import RichCommand
from click_didyoumean import DYMGroup

import logging
import sys

from rich.logging import Console, RichHandler

# from ..cli import CommandGroup
from ..cli import cli 
from ..utilities import read_chrom_sizes
from .. import __epilog__

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=Console(stderr=True))]
)


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
           short_help='Get high-quality contigs by ONT/HiFi data',
           epilog=__epilog__)
@click.pass_context
def hitig(ctx):
    """
    Higher quality of contig can help us to phasing and scaffolding assembly well.\n
    So, we can use Ultra-long or HiFi data to correct chimeric or identify high confidence regions (HCRs).

    """
    logger.setLevel(logging.INFO)
    pass



@hitig.command(cls=RichCommand, epilog=__epilog__)
@click.option(
    '-f',
    '--fasta',
    required=True,
    type=click.Path(exists=True)
)
@click.option(
    '-i', 
    '--fastq',
    required=True,
)
@click.option(
    '-a',
    '--min-as',
    'min_as',
    metavar='INT',
    help='Minimum alignment score (AS)',
    default=2000,
    show_default=True,
    type=int,
)
@click.option(
    '-q',
    '--min-mapq',
    'min_mapq',
    metavar='INT',
    help='minimum mapping quality',
    default=2,
    show_default=True,
    type=int
)
@click.option(
    '-n',
    '--nhap',
    help='maximum number of supplement alignment records',
    default=4,
    show_default=True,
    type=int,
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
    default=10,
    show_default=True,
    type=int,
)
@click.option(
    '--step-size',
    'step_size',
    help='step size',
    default=1000,
    show_default=True,
    type=int
)
@click.option(
    '-ms',
    '--min-sa',
    'min_sa',
    metavar='INT',
    help='Number of minmum split alignments in a window.',
    default=5,
    type=int,
    show_default=True,
)
@click.option(
    '-e',
    '--edge',
    metavar='INT',
    help='Minimum length of contigs, which mean too short contigs will be skipped.',
    default=20000,
    type=int,
    show_default=True,
)
@click.option(
    '-d',
    '--min-depth',
    'min_depth',
    metavar='INT',
    help='minimum depth of windows',
    type=int,
    default=1,
    show_default=True
)
@click.option(
    '-c', 
    '--cutoff',
    metavar='FLOAT',
    help='cutoff of identification chimeric contigs, '
    'which equal (count of splited alignment reads)/(avarage depth in chimeric region).',
    type=click.FloatRange(0, 1.0),
    default=0.75,
    show_default=True
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
@click.option(
    '-s',
    '--steps',
    metavar='STR',
    help="steps",
    default="1,2,3",
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
def pipeline(fasta, fastq, min_as, min_mapq, 
    nhap, window, min_windows, step_size, 
    min_sa, edge, min_depth,
    cutoff, hifi, 
    threads, output, 
    steps, skip_steps):
    """
    Pipeline of chimeric correct by hitig.

    > **Steps:**\n
        1. correct-alignments: mapping and correct alignments\n
        2. find-chimeric: find chimeric contigs \n
    
    > **Usages:**\n
    > - Input Ultra-long data\n
    ```bash
    $ hitig pipeline -f contigs.fasta -i ultra-long.fastq.gz -t 20 
    ```
    > - Input HiFi data\n
    ```bash
    $ hitig pipeline -f contigs.fasta -i hifi.fasta.gz -t 20 --hifi -m 2
    ```
    
    """
    # """
    # > - Input HiFi data\n
    # ```bash
    # $ hitig pipeline -f contigs.fasta -i hifi.fasta.gz -t 20 --hifi -w 2000 --step-size 500 -ms 2 -m 5
    # ```
    # """
    
    from .pipeline import run 

    if steps:
        if steps == "all":
            steps = set(["1", "2"])
        else:
            steps = steps.strip().split(",")
    else:
        steps = set(map(str, [1, 2]))
    
    if skip_steps:
        skip_steps = set(skip_steps.strip().split(","))
    else:
        skip_steps = set()
    
    run(
        fasta, fastq, min_as, min_mapq, 
    nhap, window, min_windows, step_size,
    min_sa, edge, min_depth,
    cutoff, hifi,
    threads, output, 
    steps, skip_steps
    )


@hitig.command(cls=RichCommand, hidden=True)
@click.option(
    '-i',
    '--fastq',
    metavar='FilePath',
    default=None,
    show_default=True,
    help='fastq file',
    required=True
)
@click.option(
    '-w', 
    '--window',
    metavar='INT',
    default=5000,
    type=int,
    help='window size',
    show_default=True
)
@click.option(
    '-p',
    '--npart',
    metavar='INT',
    help='number of part',
    default=10,
    type=int,
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
def split_reads(fastq, window, npart, threads):
    """
    Split reads by window size.

    This is the script for split fastq file into window reads. For accelerating this process, we split fastq into N parts at first, then slice reads with window.
    """
    from .split_reads import pipe 
    pipe(fastq, window, npart, threads)


@hitig.command(cls=RichCommand)
@click.option(
    '-f',
    '--fasta',
    required=True,
    type=click.Path(exists=True)
)
@click.option(
    '-i', 
    '--fastq',
    required=True,
    type=click.Path(exists=True)
)
@click.option(
    '-a',
    '--min-as',
    'min_as',
    metavar='INT',
    help='Minimum alignment score (AS)',
    default=2000,
    show_default=True,
    type=int,
)
@click.option(
    '-q',
    '--min-mapq',
    'min_mapq',
    metavar='INT',
    help='minimum mapping quality',
    default=2,
    show_default=True,
    type=int
)
@click.option(
    '-n',
    '--nhap',
    help='maximum number of supplement alignment records',
    default=4,
    show_default=True,
    type=int,
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
    default=10,
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
def correct_alignments(
    fasta, fastq, min_as, min_mapq, 
    nhap, window, min_windows, hifi,
    threads, output
):
    """
    Mapping and correct alignments.

        The pipeline including two steps. First, 
        the 5000bp window reads are compared to the draft assembly
        using minimap2. Because the short window reads are prone 
        to occur alignment errors, we further correct errors based on 
        the longest increasing subsequence (LIS) method.
    """
    from .correct_alignments import workflow
    workflow(fasta, fastq, threads, output, window, min_windows, nhap, min_as, min_mapq, hifi)


@hitig.command(cls=RichCommand)
@click.option(
    '-p',
    '--corrected-paf',
    'corrected_paf',
    metavar='PATH',
    help='the corrected UL-ONT reads/Hifi reads mapping result, must be .paf format.',
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    '-l', 
    '--lis',
    metavar='PATH',
    help='Corrected alignment result with mapping quality greater than thredshold (eg. outputmapq.LIS.gtf).',
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    '-f',
    '--fasta',
    metavar='PATH',
    help='the raw contigs file.',
    type=click.Path(exists=True),
    required=True
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
    '--step-size',
    'step_size',
    help='step size',
    default=1000,
    show_default=True,
    type=int
)
@click.option(
    '-m',
    '--min-sa',
    'min_sa',
    metavar='INT',
    help='Number of minmum split alignments in a window.',
    default=5,
    type=int,
    show_default=True,
)
@click.option(
    '-e',
    '--edge',
    metavar='INT',
    help='Minimum length of contigs, which mean too short contigs will be skipped.',
    default=20000,
    type=int,
    show_default=True,
)
@click.option(
    '-d',
    '--min-depth',
    'min_depth',
    metavar='INT',
    help='minimum depth of windows',
    type=int,
    default=1,
    show_default=True
)
@click.option(
    '-c', 
    '--cutoff',
    metavar='FLOAT',
    help='cutoff of identification chimeric contigs, '
    'which equal (count of splited alignment reads)/(avarage depth in chimeric region).',
    type=click.FloatRange(0, 1.0),
    default=0.75,
    show_default=True
)
@click.option(
    '-o',
    '--output',
    help='output file prefix',
    default='output',
    show_default=True
)
def find_chimeric(
    corrected_paf, 
    lis, 
    fasta, 
    window, 
    step_size,
    min_sa, 
    edge,
    min_depth,
    cutoff, 
    output):
    """
    Finding split alignments and correct the chimeric contigs.

    """
    from .find_chimeric import (
        find_split_alignment,
        paf2depth,
        norm_merge_bp,
        correct_fasta
    )
 
    break_point_file = find_split_alignment.workflow(lis, window, min_sa, edge, output)

    contigsizes, depth_file = paf2depth.workflow(corrected_paf, fasta, window, step_size, outPre=output)

    break_pos_file = norm_merge_bp.workflow(break_point_file, depth_file, contigsizes,
                                                window, min_depth, cutoff, edge, output)

    correct_fasta.workflow(fasta, break_pos_file, output)


@hitig.command(cls=RichCommand)
@click.option(
    '-f',
    '--fasta',
    metavar='PATH',
    help='the raw contigs file.',
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '-l',
    '--lis',
    metavar='PATH',
    help='Corrected alignment result with mapping quality greater than thredshold (eg. outputmapq.LIS.gtf).',
    required=True,
    type=click.Path(exists=True)
)
@click.option(
    '-sa',
    '--split-align',
    'split_align',
    metavar='PATH',
    help='split alignments file, '
        '5 cols: <chr> <start> <end> <count of split-align> <split-aligned reads>',
    required=True,
    type=click.Path(exists=True)
)
@click.option(
    '-p',
    '--paf',
    required=True,
    metavar='PATH',
    help='Ul-ONT/HiFi raw mapping results in paf format.',
    type=click.Path(exists=True),
)
@click.option(
    '-d',
    '--depth',
    metavar='PATH',
    help='Ul-ONT/HiFi reads depth file, 4 cols: <chr> <start> <end> <depth>.',
    type=click.Path(exists=True),
    hidden=True,
)
@click.option(
    '-b',
    '--break-pos',
    'break_pos',
    metavar='PATH',
    help="use break positions to correct HCRs.",
    type=click.Path(exists=True),
    default=None,
    show_default=True,
    hidden=True,
)
# @click.option(
#     '-c',
#     '--contig-sizes',
#     'contig_sizes',
#     metavar='PATH',
#     help="contig sizes file for HCRs correct by break positions.",
#     type=click.Path(exists=True),
#     default=None,
#     show_default=True
# )
@click.option(
    '-q',
    '--min-mapq',
    'min_mapq',
    metavar='INT',
    help='minimum mapping quality',
    default=0,
    show_default=True,
    type=int
)
@click.option(
    '-mc',
    '--min-count',
    'min_count',
    metavar='INT',
    help='minimum count of windows in LIS',
    default=3,
    type=int,
    show_default=True,
)
@click.option(
    '-ms',
    '--min-score',
    'min_score',
    metavar='FLOAT',
    help='minimum score of mapping'
        ' s = (count of windows with mapq10 / count of all windows)',
    default=0.6,
    type=float,
    show_default=True,
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
    '--step-size',
    'step_size',
    help='step size',
    default=1000,
    show_default=True,
    type=int
)
@click.option(
    '-M',
    '--max-depth',
    'max_depth',
    metavar='INT',
    help='maximum depth',
    type=int,
    default=None,
    show_default=True
)
@click.option(
    '-ucr',
    '--use-continous-regions',
    '--use-continous-region',
    "use_continous_region",
    help="only retain contious mapping regions of min count window to generate hcr",
    default=False,
    is_flag=True,
)
@click.option(
    '-o',
    '--output',
    help='output file prefix',
    default='output',
    show_default=True
)
def hcr(fasta, lis, split_align, paf, 
            depth, break_pos, min_mapq,
            min_count, min_score, window, 
            step_size, max_depth, 
            use_continous_region, output):
    """
    Identify high confidence regions (HCRs).

        Identify the high confidence regions of genome to help the subsequence analysis.

    """
    from .hcr import hcr, bed2depth
    from .hcr.clean import clean
    from .find_chimeric import paf2depth 

    contigsizes, depth = paf2depth.workflow(paf, fasta, window, step_size, output, min_mapq=min_mapq)
    if contigsizes and not break_pos:
        contig_sizes_df = read_chrom_sizes(contigsizes)
        contig_sizes_db = contig_sizes_df.to_dict()['length']

        all_contig_depth_df, high_coverage_df, low_coverage_df, hcr_from_depth_file =\
              bed2depth.workflow(depth, window, output, max_depth)
        
        # clean(fasta, low_coverage_df, f"{output}.cleaned.fasta")

    elif contigsizes and break_pos:
        contig_sizes_df = read_chrom_sizes(contigsizes)
        contig_sizes_db = contig_sizes_df.to_dict()['length']

        all_contig_depth_df, high_coverage_df, low_coverage_df, hcr_from_depth_file =\
              bed2depth.workflow(depth, window, output, max_depth)
        
        # clean(fasta, low_coverage_df, f"{output}.cleaned.fasta", break_pos)


    if break_pos and contigsizes:
        hcr.workflow(lis, split_align, hcr_from_depth_file, 
                        min_count, min_score, output, break_pos, contigsizes, 
                        use_continous_regions=use_continous_region)
    else:
        if break_pos:
            logger.warning("The file of break position and contig sizes should be input at the same time.")
        hcr.workflow(lis, split_align, hcr_from_depth_file, 
                        min_count, min_score, output, use_continous_regions=use_continous_region)
        

@hitig.command(cls=RichCommand, hidden=True)
@click.option(
    '-f',
    '--fasta',
    metavar='Fasta',
    help='the raw contigs file.',
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '-d',
    '--depth',
    metavar="Depth",
    help="Depth file of mapping quality 0",
    type=click.Path(exists=True),
    required=True
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
    "-oc",
    "--output-collapsed",
    is_flag=True,
    default=False,
    show_default=True
)
@click.option(
    '--force',
    help='Force run all the command, ignore existing results.'
    ' The index file also will be removed.',
    is_flag=True,
    default=False,
    show_default=True,
)
def clean(fasta, depth, threads, 
          output_collapsed, force):
    """
    Remove false duplicated and dup collapsed contigs.
    """
    from .hcr.clean import Clean 

    c = Clean(fasta, 
              depth, 
              output_collapsed=output_collapsed,
              threads=threads, 
              force=force)
    c.run()





ALIASES = {
    "sr": split_reads,
    "ca": correct_alignments,
    "fc": find_chimeric,
}