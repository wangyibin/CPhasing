#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
mapper of ALLHiC2
"""

import logging
import os
import os.path as op
import sys

from pathlib import Path
from shutil import which
from subprocess import Popen, PIPE

from .utilities import run_cmd

logger = logging.getLogger(__name__)
class HisatMapper(object):
    """
    single ends mapping by hisat2

    Params:
    --------

    Returns:
    --------

    Examples:
    --------
    >>> mapper = HisatMapper('reference.fasta', 'sample_R1.fastq.gz')
    """
    def __init__(self, index, fastq, min_quality=10, 
                    threads=4, 
                    additional_arguments=(),
                    hisat2_path='hisat2'):

        self.index = index
        self.fastq = Path(fastq)
        self.threads = threads
        self.min_quality = min_quality
        self.additional_arguments = additional_arguments
        self._path = hisat2_path
        if which(self._path) is None:
            raise ValueError(f"{self._path}: command not found")

        self.prefix = self.fastq.with_suffix('')
        while self.prefix.suffix in {'.fastq', 'gz', 'fq'}:
            self.prefix = self.prefix.with_suffix('')

        self.global_bam = Path(f'{self.prefix}.global.bam')
        self.unmap_fastq = Path(f'{self.prefix}.unmap.fastq')
        self.trimmed_fastq = Path(f'{self.prefix}.trimed.fastq')
        self.local_bam = Path(f'{self.prefix}.local.bam')
        self.merge_bam = Path(f'{self.prefix}.merge.bam')
        self.sorted_bam = Path(f'{self.prefix}.sorted.bam')

    def global_mapping(self):
       

        map_command = [f'{self._path}', '-x', self.index, '-U', str(self.fastq), 
                    '-k', '1', '--no-spliced-alignment', 
                    '--un', str(self.unmap_fastq),
                    '--no-softclip', '--threads', str(self.threads)]

        bam_command = ['samtools', 'view', '-bS', '-@', '4', '-F', '4', '-']
        
        logger.info('Running command:')
        logger.info('\t' + ' '.join(map_command) + ' | ' + ' '.join(bam_command)
                    + ' > ' + str(self.global_bam))

        pipelines = []
        try:
            pipelines.append(
                Popen(map_command, stdout=PIPE, 
                        stderr=open(os.devnull, 'w'), bufsize=-1)
            )

            pipelines.append(
                Popen(bam_command, stdin=pipelines[-1].stdout,
                        stdout=open(self.global_bam, 'wb'), 
                        bufsize=-1)
            )

            pipelines[-1].wait()
        
        finally:
            for p in pipelines:
                if p.poll() is None:
                    p.terminate()

    def trim_fastq(self):
        
        # command = ['../bin/cutsite_trimming', '--fastq', fastq, '--cutsite', 
        #             cutsite, '--out', self.trimed_fastq]
        
        # run_cmd(command)

        from .cutsite import cutsite_trimming
        cutsite_trimming(self.unmap_fastq, 'AAGCTAGCTT', self.trimmed_fastq)

    def trimmed_mapping(self):
        
        map_command = [self._path, '-x', self.index, '-U', str(self.trimmed_fastq), 
                    '-k', '1', '--no-spliced-alignment', '--no-softclip',
                    '--threads', str(self.threads)]
        bam_command = ['samtools', 'view', '-@', '4', '-bS', '-']

        logger.info('Running command:')
        logger.info('\t' + ' '.join(map_command) + ' | ' + ' '.join(bam_command)
                    + ' > ' + str(self.local_bam))

        pipelines = []
        try:
            pipelines.append(
                Popen(map_command, stdout=PIPE, 
                        stderr=open(os.devnull, 'w'), bufsize=-1)
            )

            pipelines.append(
                Popen(bam_command, stdin=pipelines[-1].stdout,
                        stdout=open(self.local_bam, 'wb'), 
                        bufsize=-1)
            )

            pipelines[-1].wait()
        
        finally:
            for p in pipelines:
                if p.poll() is None:
                    p.terminate()


    def combine(self):
        command = ['samtools', 'merge', '-f', '-@', str(self.threads), 
                    str(self.merge_bam), str(self.global_bam), str(self.local_bam)]
        

        run_cmd(command)

    def sort(self):
        command = ['samtools', 'sort', '-@', str(self.threads), 
                    '-n', str(self.merge_bam), '-o', str(self.sorted_bam)]
        
        
        run_cmd(command)

    def clean(self):
        self.global_bam.unlink()
        self.trimmed_fastq.unlink()
        self.unmap_fastq.unlink()
        self.local_bam.unlink()
        self.merge_bam.unlink()

    def run(self):
        self.global_mapping()
        self.trim_fastq()
        self.trimmed_mapping()
        self.combine()
        self.sort()

    @classmethod
    def pair():
        command = []
        

