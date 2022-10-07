#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
mapper of Hi-C data.
"""

import logging
import os
import os.path as op
import sys

from pathlib import Path
from shutil import which
from subprocess import Popen, PIPE

from .utilities import run_cmd, ligation_site

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
    def __init__(self, index, fastq, enzyme, 
                    min_quality=10, 
                    threads=4, 
                    additional_arguments=(),
                    hisat2_path='hisat2', 
                    log_dir='logs'):

        self.index = index
        self.fastq = Path(fastq)
        self.enzyme = enzyme
        self.ligation_site = ligation_site(self.enzyme)[0]
        self.threads = threads
        self.min_quality = min_quality
        self.additional_arguments = additional_arguments
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

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
                        stderr=open(f'{self.log_dir}/{self.prefix}'
                        '.global_mapping.log','w'), 
                        bufsize=-1)
            )

            pipelines.append(
                Popen(bam_command, stdin=pipelines[-1].stdout,
                        stdout=open(self.global_bam, 'wb'), 
                        stderr=open(f'{self.log_dir}/{self.prefix}'
                        '.global_mapping.log','w'),
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
        cutsite_trimming(self.unmap_fastq, self.ligation_site, self.trimmed_fastq)

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
                        stderr=open(f'{self.log_dir}/{self.prefix}'
                        '.trimmed_mapping.log','w'), 
                        bufsize=-1)
            )

            pipelines.append(
                Popen(bam_command, stdin=pipelines[-1].stdout,
                        stdout=open(self.local_bam, 'wb'), 
                        stderr=open(f'{self.log_dir}/{self.prefix}'
                        '.global_mapping.log','w'),
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
        

        run_cmd(command, log=f'{self.log_dir}/{self.prefix}.bam_merge.log')

    def sort(self):
        command = ['samtools', 'sort', '-@', str(self.threads), 
                    '-n', str(self.merge_bam), '-o', str(self.sorted_bam)]
        
        
        run_cmd(command, log=f'{self.log_dir}/{self.prefix}.bam_sort.log')

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

    @classmethod
    def create_index(self, reference):
        cmd = [f'{self._path}-build', '-p', str(self.threads),
                str(self.index.parent), str(self.index)]
        
        run_cmd(cmd)

class ChromapMapper:
    """
    Mapper for Hi-C reads using chromap.

    Params:
    --------
    reference: str
        contig-level assembly.
    read1: str
        Hi-C paired-end read1.
    read2: str
        Hi-C paired-end read2.
    min_quality: int, default 1
        minimum number of mapping quality
    threads: int, default 4
        number of threads.
    additional_arguments: tuple, default None
        additional arguments of chromap.
    path: str, default "chromap"
        Path of chromap.
    log_dir: str, default "logs"
        directory of logs.
    
    Returns:
    --------
    object:
        object of ChromapMapper
    
    Examples:
    --------
    >>> cm = ChromapMapper(reference, read1, read2, 
                min_quality=mapq, threads=threads)
    >>> cm.run()
    """
    def __init__(self, reference, read1, read2, min_quality=30, 
                    threads=4, additional_arguments=(), 
                    path='chromap', log_dir='logs'):
        self.reference = Path(reference)
        self.index_path = Path(f'{self.reference.stem}.index')
        self.read1 = Path(read1)
        self.read2 = Path(read2)
        self.threads = threads
        self.min_quality = min_quality
        self.additional_artuments=()
        
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self._path = path
        if which(self._path) is None:
            raise ValueError(f"{self._path}: command not found")

        self.prefix = Path(self.read1.stem).with_suffix('')
        while self.prefix.suffix in {'.fastq', 'gz', 'fq'}:
            self.prefix = self.prefix.with_suffix('')
        self.prefix = Path(str(self.prefix).replace('_R1', ''))

        self.output_pairs = Path(f'{self.prefix}.pairs')

    def index(self):
        """
        Create chromap index.
        """
        cmd = [self._path, '-t', str(self.threads), 
                '-i', '-r', str(self.reference), '-o', 
                str(self.index_path)]
        
        run_cmd(cmd, log=f'{str(self.log_dir)}/{self.index_path}.log')
    
    def mapping(self):
        cmd = [self._path, '-t', str(self.threads), 
                '--preset', 'hic', 
                '-q', str(self.min_quality),
                '-x', str(self.index_path), 
                '-r', str(self.reference), '-1', str(self.read1),
                '-2', str(self.read2), '-o', str(self.output_pairs)]
        
        run_cmd(cmd, log=f'{str(self.log_dir)}/{self.prefix}_mapping.log')

    def run(self):
        if not self.index_path.exists():
            self.index()
        else:
            logger.warning(f'The index of `{self.index_path}` was exisiting, skipped ...')
        
        self.mapping()


class ReadPair:
    def __init__(self):
        pass
    