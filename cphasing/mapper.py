#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
mapper of Pore-C/CiFi and Hi-C data.
"""

import logging
import os
import os.path as op
import sys
import subprocess
import pandas as pd 

from pathlib import Path
from shutil import which
from subprocess import Popen, PIPE

from .utilities import (
    is_compressed_table_empty,
    is_empty,
    run_cmd,
    ligation_site,
    to_humanized
)

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
                    # additional_arguments=("-B" "5" "-O" "2" "-E" "1"),
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
            raise ValueError(f"{self._path}: command not found.")

        self.prefix = self.fastq.with_suffix('')
        while self.prefix.suffix in {'.fastq', 'gz', 'fq', '.fq', '.gz'}:
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

        from .hic.cutsite import cutsite_trimming
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
    read1: list
        Hi-C paired-end read1.
    read2: list
        Hi-C paired-end read2.
    min_quality: int, default 0
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
    def __init__(self, reference, read1, read2, output_format='pairs.pqs',
                    kmer_size=17, window_size=7, min_quality=0, 
                    threads=4, additional_arguments=None, 
                    path='_chromap', log_dir='logs'):
        

        self.reference = Path(reference)
        self.index_path = Path(f'{self.reference.stem}.index')
        self.contigsizes = Path(f'{self.reference.stem}.contigsizes')
        self.read1 = read1
        self.read2 = read2
        self.output_format = output_format
        self.threads = threads
        self.kmer_size = kmer_size
        self.window_size = window_size
        self.min_quality = min_quality
        self.additional_artuments = additional_arguments
        
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self._path = path
        if which(self._path) is None:
            if which(self._path.rstrip("_")) is None:
                raise ValueError(f"{self._path}: command not found")
            else:
                self._path = self._path.rstrip("_")

        if which('crabz') is None and which('pigz') is None:
            raise ValueError(f"pigz: command not found")
    
        if which('crabz') is None:
            self.compressor = 'pigz'
        else:
            self.compressor = 'crabz'

        self.prefix = Path(Path(self.read1[0]).stem).with_suffix('')
        while self.prefix.suffix in {'.fastq', 'gz', 'fq', '.fq', '.gz', '_R1', '_1', '_2', '.fasta', 'fasta', '.fa'}:
            self.prefix = self.prefix.with_suffix('')
        
        if str(self.prefix).endswith("_R1"):
            self.prefix = Path(str(self.prefix)[:-3])
        elif str(self.prefix).endswith("_1"):
            self.prefix = Path(str(self.prefix)[:-2])

        self.prefix = Path(str(self.prefix).replace('_R1', '').replace('_1', ''))

        self.output_pairs = Path(f'{self.prefix}.pairs')

    def get_contig_sizes(self):
        cmd = ["cphasing-rs", "chromsizes", str(self.reference), 
                "-o", str(self.contigsizes)]
        flag = run_cmd(cmd, log=f"{self.log_dir}/{self.prefix}.contigsizes.log")
        assert flag == 0, "Failed to execute command, please check log."

    def get_genome_size(self):
        from .utilities import get_genome_size

        return get_genome_size(self.reference)

    def index(self):
        """
        Create chromap index.
        """
        cmd = [self._path, '-t', str(self.threads), 
                '-k', str(self.kmer_size), '-w', str(self.window_size),
                '-i', '-r', str(self.reference), '-o', 
                str(self.index_path)]
        with open("hicmapper.index.sh", "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")

        flag = run_cmd(cmd, log=f'{str(self.log_dir)}/{self.index_path}.log')
        assert flag == 0, "Failed to execute command, please check log."

        if not Path(self.index_path).exists():
            logger.error(f"Failed to create index of `{self.index_path}`.")
            sys.exit(1)

    def mapping(self):
        reads_1 = " ".join(self.read1)
        reads_2 = " ".join(self.read2)
        cmd = f"{self._path} -t {self.threads} " \
                f"-q {self.min_quality} "\
                f"--preset hic " \
                f"-x {str(self.index_path)} " \
                f"--remove-pcr-duplicates " \
                f"-r {str(self.reference)} " \
                f"-1 <(pigz -dc -p 8 {reads_1}) " \
                f"-2 <(pigz -dc -p 8 {reads_2}) " \
                f"-o {str(self.output_pairs)} " \
                f"2>{self.log_dir}/{self.prefix}.mapping.log"

        with open("hicmapper.cmd.sh", "w") as f:
            f.write("#!/bin/bash\n")
            f.write(cmd + "\n")
        logger.info('Running command:')
        logger.info('\t' + cmd)
        subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)


        if not Path(self.output_pairs).exists():
            logger.error(f"Failed to create pairs of `{self.output_pairs}`.")
            sys.exit(1)

        if is_empty(self.output_pairs):
            logger.error(f"Empty pairs of `{self.output_pairs}`.")
            sys.exit(1)    

        logger.info("Done.")

    def compress(self):
        if self.compressor == 'crabz':
            cmd = ['crabz', '-I', '--format', 'mgzip', '-p', str(self.threads), f'{str(self.output_pairs)}']
        else:
            cmd = ['pigz', '-f', '-p', str(self.threads), f'{str(self.output_pairs)}']

        with open("hicmapper.cmd.sh", "+a") as f:
            f.write(" ".join(cmd) + "\n")

        flag = run_cmd(cmd, log=f'{str(self.log_dir)}/{self.prefix}.compress.log')
        assert flag == 0, "Failed to execute command, please check log."
        logger.info("Done.")

    def pairs2pqs(self):
        from .cli import pairs2pqs
        
        logger.info('Converting pairs to pairs.pqs ...')

        args = [str(self.output_pairs) + ".gz", "-t", str(self.threads)]
        try:
            pairs2pqs.main(args=args, prog_name='pairs2pqs')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
        
        logger.info(f"Successfully converted pairs to `{str(self.output_pairs)}.pqs`")
            

    def run(self):
        self.get_contig_sizes()
        if not self.index_path.exists():
            self.index()
        else:
            logger.warning(f'The index of `{self.index_path}` was exisiting, skipped ...')
        
        self.mapping()
        self.compress()
        if self.output_format == 'pairs.pqs':
            self.pairs2pqs()

class MinimapMapper:
    """
    Mapper for Hi-C reads using minimap2.

    Params:
    --------
    reference: str
        contig-level assembly.
    read1: str
        Hi-C paired-end read1.
    read2: str
        Hi-C paired-end read2.
    min_quality: int, default 0
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
    def __init__(self, reference, read1, read2, 
                    kmer_size=21, window_size=11, min_quality=0, 
                    threads=4, additional_arguments=(), 
                    output_format='pairs.pqs',
                    path='minimap2', log_dir='logs'):
        

        self.reference = Path(reference)
        self.contigsizes = Path(f'{self.reference.stem}.contigsizes')
        self.read1 = read1 
        self.read2 = read2
        self.threads = threads
        self.kmer_size = kmer_size
        self.window_size = window_size
        self.min_quality = min_quality
        self.output_format = output_format
        self.additional_artuments=()
        
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self._path = path

        if self.get_genome_size() > 4e9:
            self.batchsize = to_humanized(self.get_genome_size())
            self.batchsize = str(int(self.batchsize[:-1]) + 1) + self.batchsize[-1]
        else:
            self.batchsize = "4g"


        if which('crabz') is None and which('pigz') is None:
            raise ValueError(f"pigz: command not found")
    
        if which('crabz') is None:
            self.compressor = 'pigz'
        else:
            self.compressor = 'crabz'

        self.prefix = Path(Path(self.read1[0]).stem).with_suffix('')
        while self.prefix.suffix in {'.fastq', 'gz', 'fq', '.fq', '.gz', '_R1', '_1', '_2', '.fasta', 'fasta', '.fa'}:
            self.prefix = self.prefix.with_suffix('')
        
        if str(self.prefix).endswith("_R1"):
            self.prefix = Path(str(self.prefix)[:-3])
        elif str(self.prefix).endswith("_1"):
            self.prefix = Path(str(self.prefix)[:-2])

        self.prefix = Path(str(self.prefix).replace('_R1', '').replace('_1', ''))

        self.output_pairs = Path(f'{self.prefix}.{self.output_format}')

    def get_contig_sizes(self):
        cmd = ["cphasing-rs", "chromsizes", str(self.reference), 
                "-o", str(self.contigsizes)]
        flag = run_cmd(cmd, log=f"{self.log_dir}/{self.prefix}.contigsizes.log")
        assert flag == 0, "Failed to execute command, please check log."

    def get_genome_size(self):
        from .utilities import get_genome_size

        return get_genome_size(self.reference)

   
    def mapping(self):
        if self.compressor == 'crabz':
            decompress_cmd = f'crabz -d -p 8 2>{self.log_dir}/{self.prefix}.hic.mapping.decompress.log'
            compress_cmd = f'crabz -p 8 --format mgzip 2>{self.log_dir}/{self.prefix}.hic.mapping.compress.log'
        else:
            decompress_cmd = f'pigz -dc -p 8 2>{self.log_dir}/{self.prefix}.hic.mapping.decompress.log'
            compress_cmd = f'pigz -c -p 8 2>{self.log_dir}/{self.prefix}.hic.mapping.compress.log'
        self.secondary = "no"
        reads_1 = " ".join(self.read1)
        reads_2 = " ".join(self.read2)
        cmd = f"{self._path} -t {self.threads} -c -x sr " \
                f"-I {self.batchsize} " \
                f"-k {self.kmer_size} -w {self.window_size} {str(self.reference)} " \
                f"<({decompress_cmd} {reads_1}) <({decompress_cmd} {reads_2}) " \
                f"--secondary={self.secondary} " \
                f"2>{self.log_dir}/{self.prefix}.hic.mapping.log | " \
                f"{compress_cmd} > {str(self.prefix)}.paf.gz"
        
        self.output_paf = f"{str(self.prefix)}.paf.gz"
        with open("hicmapper.cmd.sh", "w") as f:
            f.write("#!/bin/bash\n")
            f.write(cmd + "\n")
        logger.info('Running command:')
        logger.info('\t' + cmd)

        subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)


    def paf2pairs(self):
        cmd = ['cphasing-rs', 'paf2pairs',
               "-l", str(30),
               "-p", str(0.99),
               '-q', str(self.min_quality),
               str(self.output_paf), str(self.contigsizes),
               '-o', str(self.output_pairs)]
        
        with open("hicmapper.cmd.sh", "a") as f:
            f.write(" ".join(cmd) + "\n")

        flag = run_cmd(cmd, log=f'{str(self.log_dir)}/{self.prefix}.paf2pairs.log')
        assert flag == 0, "Failed to execute command, please check log."
        self.output_porec = f'{self.prefix}.porec.gz'

        if not Path(self.output_pairs).exists():
            logger.error(f"Failed to create pairs of `{self.output_pairs}`.")
            sys.exit(1)

        if is_empty(self.output_pairs):
            logger.error(f"Empty pairs of `{self.output_pairs}`.")
            sys.exit(1)    
        # if Path(self.output_paf).exists():
        #     Path(self.output_paf).unlink()
            
        #     logger.info(f"Removed intermediate file: `{self.output_paf}`.")

        if Path(self.output_porec).exists():
            Path(self.output_porec).unlink()
            if Path(f"{self.prefix}.porec.read.summary").exists():
                Path(f"{self.prefix}.porec.read.summary").unlink()
            if Path(f"{self.prefix}.pairs.concatemer.summary").exists():
                Path(f"{self.prefix}.pairs.concatemer.summary").unlink()
            logger.info(f"Removed intermediate file: `{self.output_porec}`.")

        logger.info("Done.")         

    def run(self):
        self.get_contig_sizes()
        self.mapping()
        self.paf2pairs()

class BwaMapper:
    """
    Mapper for Hi-C reads using bwa-mem2.

    Params:
    --------
    reference: str
        contig-level assembly.
    read1: list
        Hi-C paired-end read1.
    read2: list
        Hi-C paired-end read2.
    min_quality: int, default 0
        minimum number of mapping quality
    threads: int, default 4
        number of threads.
    path: str, default "bwa-mem2"
        Path of bwa-mem2.
    log_dir: str, default "logs"
        directory of logs.
    """
    def __init__(self, reference, read1, read2, 
                    min_quality=0, threads=4, 
                    additional_arguments=(), 
                    output_format='pairs.pqs',
                    path='bwa-mem2', log_dir='logs'):
        
        self.reference = Path(reference)
        self.contigsizes = Path(f'{self.reference.stem}.contigsizes')
        self.read1 = read1 
        self.read2 = read2
        self.threads = threads
        self.min_quality = min_quality
        self.output_format = output_format
        self.additional_arguments = additional_arguments
        
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self._path = path
        if which(self._path) is None:
             raise ValueError(f"{self._path}: command not found")
        
        if which('samtools') is None:
            raise ValueError("samtools: command not found")

        if which('crabz') is None and which('pigz') is None:
            raise ValueError(f"pigz: command not found")
    
        if which('crabz') is None:
            self.compressor = 'pigz'
        else:
            self.compressor = 'crabz'

        self.prefix = Path(Path(self.read1[0]).stem).with_suffix('')
        while self.prefix.suffix in {'.fastq', 'gz', 'fq', '.fq', '.gz', '_R1', '_1', '_2', '.fasta', 'fasta', '.fa'}:
            self.prefix = self.prefix.with_suffix('')
        
        if str(self.prefix).endswith("_R1"):
            self.prefix = Path(str(self.prefix)[:-3])
        elif str(self.prefix).endswith("_1"):
            self.prefix = Path(str(self.prefix)[:-2])

        self.prefix = Path(str(self.prefix).replace('_R1', '').replace('_1', ''))

        self.output_bam = Path(f'{self.prefix}.bam')
        self.output_pairs = Path(f'{self.prefix}.{self.output_format}')

    def get_contig_sizes(self):
        cmd = ["cphasing-rs", "chromsizes", str(self.reference), 
                "-o", str(self.contigsizes)]
        flag = run_cmd(cmd, log=f"{self.log_dir}/{self.prefix}.contigsizes.log")
        assert flag == 0, "Failed to execute command, please check log."

    def index(self):
        """
        Create bwa-mem2 index.
        """
        if Path(f"{self.reference}.0123").exists():
            logger.warning(f'The index of `{self.reference}` was existing, skipped ...')
            return

        cmd = [self._path, 'index', str(self.reference)]
        
        with open("hicmapper.index.sh", "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")

        flag = run_cmd(cmd, log=f'{str(self.log_dir)}/{self.reference.name}.index.log', out2err=True)
        assert flag == 0, "Failed to execute command, please check log."

    def mapping(self):
        if self.compressor == 'crabz':
            decompress_cmd = f'crabz -d -p 8 2>{self.log_dir}/{self.prefix}.hic.mapping.decompress.log'
        else:
            decompress_cmd = f'pigz -dc -p 8 2>{self.log_dir}/{self.prefix}.hic.mapping.decompress.log'
        
        reads_1 = " ".join(self.read1)
        reads_2 = " ".join(self.read2)
        
        cmd = f"{self._path} mem -t {self.threads} -5SPM " \
              f"{self.reference} " \
              f"<({decompress_cmd} {reads_1}) <({decompress_cmd} {reads_2}) " \
              f"2>{self.log_dir}/{self.prefix}.hic.mapping.log | " \
              f"samblaster 2>{self.log_dir}/{self.prefix}.hic.samblaster.log | " \
              f"samtools view -@ 12 -bS -F 3340 - > {self.output_bam}"

        with open("hicmapper.cmd.sh", "w") as f:
            f.write("#!/bin/bash\n")
            f.write(cmd + "\n")
        
        logger.info('Running command:')
        logger.info('\t' + cmd)

        subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)

        if not self.output_bam.exists():
            logger.error(f"Failed to create bam of `{self.output_bam}`.")
            sys.exit(1)

    def bam2pairs(self):
        cmd = ['cphasing-rs', 'bam2pairs',
               '-q', str(self.min_quality),
               str(self.output_bam),
               ]

        cmd2 = ['cphasing-rs', 'pairs2pqs', 
                '-', 
                '-o', f'{str(self.output_pairs)}',
               ]
        
        cmd = cmd + [f"2>{str(self.log_dir)}/{str(self.prefix)}.bam2pairs.log"] + ['|'] + cmd2 + [f"2>>{str(self.log_dir)}/{str(self.prefix)}.bam2pairs.log"]
        
        with open("hicmapper.cmd.sh", "a") as f:
            f.write(" ".join(cmd) + "\n")

        logger.info('Running command:')
        logger.info('\t' + " ".join(cmd))
    
        flag = subprocess.run(" ".join(cmd), shell=True, executable='/bin/bash')
        assert flag.returncode == 0, "Failed to execute command, please check log."
        if not Path(self.output_pairs).exists():
            logger.error(f"Failed to create pairs of `{self.output_pairs}`.")
            sys.exit(1)

        if is_empty(self.output_pairs):
            logger.error(f"Empty pairs of `{self.output_pairs}`.")
            sys.exit(1)

        logger.info("Done.")

    def run(self):
        self.get_contig_sizes()
        self.index()
        self.mapping()
        self.bam2pairs()
class PoreCMapper:
    def __init__(self, reference, reads, pattern="GATC", 
                    k=15, w=10, min_quality=1, 
                    min_identity=0.8, min_length=150, 
                    max_edge=0, secondary=False, realign=False,
                    additional_arguments="-x map-ont", outprefix=None,
                    threads=4, path='minimap2', log_dir='logs',
                    force=False):
        self.reference = Path(reference)
        self.index_path = Path(f'{self.reference.stem}.index')
        self.contigsizes = Path(f'{self.reference.stem}.contigsizes')
        self.reads = reads
        self.pattern = pattern
        self.additional_arguments = additional_arguments
        self.secondary = secondary
        self.realign = realign
        if self.realign:
            self.secondary = True


        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self._path = path
        if which(self._path) is None:
            raise ValueError(f"{self._path}: command not found")
        
        if which('crabz') is None and which('pigz') is None:
            raise ValueError(f"pigz: command not found")
    
        if which('crabz') is None:
            self.compressor = 'pigz'
        else:
            self.compressor = 'crabz'

        if which('cphasing-rs') is None:
            raise ValueError(f"cphasing-rs: command not found")


        self.prefix = Path(Path(self.reads[0]).stem).with_suffix('')
        while self.prefix.suffix in {'.fastq', 'gz', 'fq', "fasta", 'fa'}:
            self.prefix = self.prefix.with_suffix('')
       
        self.prefix = Path(str(self.prefix)) if not outprefix else outprefix


        if self.get_genome_size() > 4e9:
            self.batchsize = to_humanized(self.get_genome_size())
            self.batchsize = str(int(self.batchsize[:-1]) + 1) + self.batchsize[-1]
        else:
            self.batchsize = "4g"

        self.threads = threads
        self.min_quality = min_quality
        self.min_identity = min_identity
        self.min_length = min_length
        self.max_edge = max_edge
        self.k = k 
        self.w = w
        self.force = force
  
        self.outpaf = f'{self.prefix}.paf.gz'
        self.realign_outpaf = f'{self.prefix}.realign.paf.gz'
        if realign:
            self.outporec = f'{self.prefix}.realign.porec.gz'
            self.outpairs = f'{self.prefix}.realign.pairs.pqs'
        else:
            self.outporec = f'{self.prefix}.porec.gz'
            self.outpairs = f'{self.prefix}.pairs.pqs'
        

    def get_genome_size(self):
        from .utilities import get_genome_size

        return get_genome_size(self.reference)

    def get_contig_sizes(self):
        cmd = ["cphasing-rs", "chromsizes", str(self.reference), 
                "-o", str(self.contigsizes)]
        flag = run_cmd(cmd, log=f"{self.log_dir}/{self.prefix}.contigsizes.log")
        assert flag == 0, "Failed to execute command, please check log."

    def get_digest_bed(self):
        cmd = ["cphasing-rs", "digest", str(self.reference), "-p", str(self.pattern), 
                    "-o", f"{self.reference.stem}.slope.{self.pattern}.bed" ]
        flag = run_cmd(cmd, log=f"{self.log_dir}/{self.prefix}.digest.log")
        assert flag == 0, "Failed to execute command, please check log."

    def index(self):
        """
        Create minimap2 index
        """
        cmd = [self._path, "-t", str(self.threads), 
                "-k", str(self.k),
                "-w", str(self.w),
                "-I", self.batchsize,
                "-d", str(self.index_path),
                str(self.reference)]

        flag = run_cmd(cmd, log=f'{str(self.log_dir)}/{self.index_path}.log')
        assert flag == 0, "Failed to execute command, please check log."


    def mapping(self):

        secondary = self.secondary 
        
        secondary = "--secondary=yes" if secondary else "--secondary=no"
        if "secondary" in self.additional_arguments:
            secondary = ""
            
        if self.compressor == 'crabz':
            decompress_cmd = f'crabz -d -p 8 2>{self.log_dir}/{self.prefix}.mapping.decompress.log'
            compress_cmd = f'crabz -p 8 --format mgzip 2>{self.log_dir}/{self.prefix}.mapping.compress.log'
        else:
            decompress_cmd = f'pigz -dc -p 8 2>{self.log_dir}/{self.prefix}.mapping.decompress.log'
            compress_cmd = f'pigz -c -p 8 2>{self.log_dir}/{self.prefix}.mapping.compress.log'

        reads = ' '.join(self.reads)
        cmd = f"{self._path} -t {self.threads} " \
                f"-I {self.batchsize} " \
                f"-c {secondary} " \
                f"{self.additional_arguments} " \
                f"{str(self.reference)} " \
                f"<({decompress_cmd} {reads}) " \
                f" 2> {self.log_dir}/{self.prefix}.mapping.log " \
                f"| {compress_cmd} > {str(self.outpaf)}"

        with open("mapper.cmd.sh", "w") as f:
            f.write("#!/bin/bash\n")
            f.write(cmd + "\n")
        logger.info('Running command:')
        logger.info('\t' + cmd)
        subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)
    
    def run_realign(self):
        cmd = ["cphasing-rs", "realign", f"{self.outpaf}", "-o", f"{self.realign_outpaf}"]
        flag = run_cmd(cmd, log=f'{self.log_dir}/{self.prefix}.realign.log')
        assert flag == 0, "Failed to execute command, please check log."

    def paf2porec(self):
        paf = self.realign_outpaf if self.realign else self.outpaf
        if self.pattern:
            cmd = ["cphasing-rs", "paf2porec", f"{paf}", "-b",
                   f"{self.reference.stem}.slope.{self.pattern}.bed", "-q", 
                    f"{self.min_quality}", "-l", f"{self.min_length}", 
                    "-e", f"{self.max_edge}", "-p", f"{self.min_identity}",
                    "-o", f"{self.outporec}",
                ]
        else:
            cmd = ["cphasing-rs", "paf2porec", f"{paf}", "-q", 
                    f"{self.min_quality}", "-l", f"{self.min_length}",
                    "-e", f"{self.max_edge}", "-p", f"{self.min_identity}",
                      "-o", f"{self.outporec}"]

        with open("mapper.cmd.sh", "a") as f:
            f.write(" ".join(cmd) + "\n")
        
        flag = run_cmd(cmd, log=f'{self.log_dir}/{self.prefix}.paf2porec.log')
        assert flag == 0, "Failed to execute command, please check log."

    def porec2pairs(self):
        cmd = ["cphasing-rs", "porec2pairs", f"{self.outporec}", 
               str(self.contigsizes), "-q", f"{self.min_quality}",
                "-o", f"{self.outpairs}"]
        with open("mapper.cmd.sh", "a") as f:
            f.write(" ".join(cmd) + "\n")
        flag = run_cmd(cmd, log=f'{self.log_dir}/{self.prefix}.porec2pairs.log')
        assert flag == 0, "Failed to execute command, please check log."
    
    def run(self):
        
        if not Path(f"{self.prefix}.contigsizes").exists() or self.force:
            self.get_contig_sizes()
        else:
            logger.warning(f"The contigsizes of `{self.prefix}.contigsizes "
                           "existing, skipped ...")
        
        # if not self.index_path.exists() or self.force:
        #     self.index()
        # else:
        #     logger.warning(f'The index of `{self.index_path}` was exisiting, skipped ...')
        
        if self.pattern:
            if not op.exists(f"{self.prefix}.slope.{self.pattern}.bed") or self.force:
                self.get_digest_bed()
            else:
                logger.warning(f"The digest bed of `{self.prefix}.slope.{self.pattern}.bed "
                               "existing, skipped ...")

        if not op.exists(self.outpaf) or self.force:
            self.mapping()
            if is_compressed_table_empty(self.outpaf):
                logger.error(f"Empty mapping result: `{self.outpaf}`, please check log.")
                sys.exit(1)
        else:
            if is_compressed_table_empty(self.outpaf):
                logger.warning(f"Empty existing mapping result: `{self.outpaf}`, force rerun mapper.")
                self.mapping()
            else:
                logger.warning(f"The paf of `{self.outpaf} existing, skipped `reads mapping` ...")

        if self.realign:
            if not op.exists(self.realign_outpaf) or self.force:
                self.run_realign()
            else:
                logger.warning(f"The realign paf of `{self.realign_outpaf} existing, skipped `realign` ...")

        if not op.exists(self.outporec) or self.force:
            self.paf2porec()
        else:
            logger.warning(f"The porec table of {self.outporec} existing, skipped `paf2porec`...")
        
        if not op.exists(self.outpairs) or self.force:
            self.porec2pairs()
        else:
            if is_compressed_table_empty(self.outpairs):
                logger.warning(f"Empty existing pairs result: `{self.outpairs}`, force rerun porec2pairs.")
                self.porec2pairs()
            else:
                logger.warning(f"The pairs of {self.outpairs} existing, skipped `porec2pairs` ...")

        logger.info("Mapping done.")