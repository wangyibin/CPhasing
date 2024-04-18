#!/usr/bin/env python

"""
cli for ultra-long ont pipelines
"""

import click
import logging
import sys


# from ..cli import CommandGroup
from ..cli import cli 

logger = logging.getLogger(__name__)


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

@cli.group(cls=CommandGroup, 
           context_settings={"help_option_names": ["-h", "--help", "-help"]},
           short_help='Get high-quality contigs by ONT/HiFi data')
@click.pass_context
def hitig(ctx):
    """
    Higher quality of contig can help us to phasing and scaffolding assembly well.
    So, we can use Ultra-long or HiFi data to correct chimeric or identify high confidence regions (HCRs).

    """
    logger.setLevel(logging.INFO)
    pass



@hitig.command()
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
    default=3,
    show_default=True
)
@click.option(
    '-c', 
    '--cutoff',
    metavar='FLOAT',
    help='cutoff of identification chimeric contigs, '
    'which equal (count of splited alignment reads)/(avarage depth in chimeric region).',
    type=click.FloatRange(0, 1.0),
    default=0.8,
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
    default="1,2",
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
    nhap, window, min_windows,
    min_sa, edge, min_depth,
    cutoff, threads, output,
    steps, skip_steps):
    """
    Pipeline of chimeric correct by hitig.
    """
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
    nhap, window, min_windows,
    min_sa, edge, min_depth,
    cutoff, threads, output,
    steps, skip_steps
    )


@hitig.command(hidden=True)
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


@hitig.command()
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
    nhap, window, min_windows,
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
    workflow(fasta, fastq, threads, output, window, min_windows, nhap, min_as, min_mapq)


@hitig.command()
@click.option(
    '-p',
    '--corrected-paf',
    'corrected_paf',
    metavar='PATH',
    help='the corrected UL-ONT reads/Hifi reads mapping result, must be .paf format.',
    required=True,
)
@click.option(
    '-l', 
    '--lis',
    metavar='PATH',
    help='Corrected alignment result with mapping quality greater than thredshold (eg. outputmapq.LIS.gtf).',
    required=True,
)
@click.option(
    '-f',
    '--fasta',
    metavar='PATH',
    help='the raw contigs file.',
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
    default=3,
    show_default=True
)
@click.option(
    '-c', 
    '--cutoff',
    metavar='FLOAT',
    help='cutoff of identification chimeric contigs, '
    'which equal (count of splited alignment reads)/(avarage depth in chimeric region).',
    type=click.FloatRange(0, 1.0),
    default=0.8,
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

    depth_file = paf2depth.workflow(corrected_paf, fasta, window, output)

    break_pos_file = norm_merge_bp.workflow(break_point_file, depth_file, 
                                                window, min_depth, cutoff, output)

    correct_fasta.workflow(fasta, break_pos_file, output)


@hitig.command()
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
    '-d',
    '--depth',
    metavar='PATH',
    help='Ul-ONT/HiFi reads depth file, 4 cols: <chr> <start> <end> <depth>.',
    required=True,
    type=click.Path(exists=True)
)
@click.option(
    '-b',
    '--break-pos',
    'break_pos',
    metavar='PATH',
    help="use break positions to correct HCRs.",
    type=click.Path(exists=True),
    default=None,
    show_default=True
)
@click.option(
    '-c',
    '--contig-sizes',
    'contig_sizes',
    metavar='PATH',
    help="contig sizes file for HCRs correct by break positions.",
    type=click.Path(exists=True),
    default=None,
    show_default=True
)
@click.option(
    '-mc',
    '--min-count',
    'min_count',
    metavar='INT',
    help='minimum count of windows in LIS',
    default=5,
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
    '-M',
    '--max-depth',
    'max_depth',
    metavar='INT',
    help='maximum depth',
    type=int,
    default=100,
    show_default=True
)
@click.option(
    '-o',
    '--output',
    help='output file prefix',
    default='output',
    show_default=True
)
def hcr(lis, split_align, depth, break_pos, 
            contig_sizes, min_count, min_score, window, max_depth, output):
    """
    Identify high confidence regions (HCRs).

        Identify the high confidence regions of genome to help the subsequence analysis.

    """
    from .hcr import hcr, bed2depth

    hcr_from_depth_file = bed2depth.workflow(depth, window, max_depth, output)
    
    if break_pos and contig_sizes:
        hcr.workflow(lis, split_align, hcr_from_depth_file, 
                        min_count, min_score, output, break_pos, contig_sizes)
    else:
        logger.warning("The file of break position and contig sizes should be input at the same time.")
        hcr.workflow(lis, split_align, hcr_from_depth_file, 
                        min_count, min_score, output)



ALIASES = {
    "sr": split_reads,
    "ca": correct_alignments,
    "fc": find_chimeric,
}