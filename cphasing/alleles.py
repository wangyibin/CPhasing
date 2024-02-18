#!/usr/bin/env python
# -*- coding:utf-8 -*-


import logging
import os
import os.path as op
import sys
import re

import numpy as np

from Bio import SeqIO
from collections import defaultdict, OrderedDict
from io import StringIO
from subprocess import Popen
from pathlib import Path
from pyfaidx import Fasta

from .utilities import xopen


from .utilities import (
    cmd_exists, 
    choose_software_by_platform,
    run_cmd, 
    get_genome_size
)

logger = logging.getLogger(__name__)


class PartigLine:
    def __init__(self, line):
        line_list = line.strip().split()
        if line_list[0] == "C":
            self.type = line_list[0]
            self.seqName = line_list[1]
            self.seqLen = line_list[2]
            self.minimizerConsidered = line_list[3]
            self.mzUnique = line_list[4]
        else:
            self.type = line_list[0]
            self.seqName1 = line_list[1]
            self.seqName2 = line_list[2]
            self.strand = line_list[3]
            self.mzConsidered1 = line_list[4]
            self.mzConsidered2 = line_list[5]
            self.mzShared = line_list[6]
            self.kmerSimilarity = line_list[7]

class PartigRecords:
    def __init__(self, infile, symmetric=True):
        self._file = infile
        logger.info(f'Load file `{self._file}`.')
        self.symmetric = symmetric
        self.parse()

    def parse(self):
        self.C = []
        self.S = []
        with open(self._file) as fp:
            for line in fp:
                line = PartigLine(line)
                if line.type == 'C':
                    self.C.append(line)
                else:
                    self.S.append(line)
    
    @property
    def nSeq(self):
        return len(self.C)

    @property
    def seqNames(self):
        names = [name.seqName for name in self.C]

        return names

    @property
    def pairs(self):
        if self.symmetric:
            return [(i.seqName1, i.seqName2) for i in self.S]
        else:
            return [(i.seqName1, i.seqName2) for i in self.S 
                                if i.seqName1 < i.seqName2]
    
    @property
    def orientation(self):
        if self.symmetric:
            return {(i.seqName1, i.seqName2): i.strand for i in self.S}
        else:
            return {(i.seqName1, i.seqName2): i.strand for i in self.S 
                                            if i.seqName1 < i.seqName2}

    def convert(self, fasta):
        """
        convert ID to contig name.

        """
        if fasta[-3:] == ".gz":
            handle = xopen(fasta)
            fasta = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        else:
            fasta = Fasta(fasta)

        fastadb = dict(zip(self.seqNames, fasta.keys()))
    
        for i in range(len(self.C)):
            self.C[i].seqName = fastadb[self.C[i].seqName]
        
        for i in range(len(self.S)):
            self.S[i].seqName1 = fastadb[self.S[i].seqName1]
            self.S[i].seqName2 = fastadb[self.S[i].seqName2]

    def convertToIndex(self):
        """
        convert ID to index.
        """
        db = {}
        for i in range(len(self.C)):
            db[self.C[i].seqName] = i  
            self.C[i].seqName = i 

        for i in range(len(self.S)):
            self.S[i].seqName1 = db[self.S[i].seqName1]
            self.S[i].seqName2 = db[self.S[i].seqName2]
    
    def get_minimizerConsidered_dict(self):
        """
        get a dictionary that contains the minimizerConsidered informations.

        Returns:
        --------
        db: dict

        Examples:
        --------
        >>> pr.convert('sample.fasta')
        >>> pr.get_minimizerConsidered_dict()
        Ordereddict(('utg0000001l', 33333))
        
        """
        db = OrderedDict()
        for i in self.C:
            db[i.seqName] = int(i.minimizerConsidered)
        
        return db

    def to_alleletable(self, fasta, output, fmt="allele2"):
        """
        convert partig table to allele table

        Params:
        --------
        fasta: str
            path of fasta file
        output: _io.TextIOWrapper
            writable object of output file
        fmt: str, [default: 'allele2']
            the format of alleletable, must in {'allele1', 'allele2'}

        Examples:
        --------
        >>> out = open('output.allele.table', 'w')
        >>> pr.to_alleletable('sample.fasta', out)
        """
        self.convert(fasta)
        i = 0

        if fmt == "allele2":
            for record in self.C:
                print(f"#{record.seqName} "
                      f"{record.seqLen} "
                      f"{record.minimizerConsidered} "
                      f"{record.mzUnique}", file=output) 
            
        
        for record in self.S:
            i += 1
            
            if fmt == "allele1":
                print(i, i, record.seqName1, 
                            record.seqName2, sep='\t', file=output)
            elif fmt == "allele2":
                strand = 1 if record.strand == "+" else -1
                print(i, i, record.seqName1, 
                            record.seqName2, 
                            record.mzConsidered1,
                            record.mzConsidered2,
                            record.mzShared,
                            record.kmerSimilarity, 
                            strand, sep='\t', file=output)


class PartigAllele:
    def __init__(self,
                    fasta,
                    k=19,
                    w=19,
                    c=100,
                    n=5,
                    m=0.8,
                    d=0.2,
                    output='Allele.ctg.table',
                    log_dir='logs'):
        self.fasta = fasta
        fasta_prefix = Path(fasta).with_suffix("")
        while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
            fasta_prefix = fasta_prefix.with_suffix("")
        self.prefix = fasta_prefix
        self.partig_res = f"{self.prefix}.similarity.res"

        ## parameters for partig
        self.k = k
        self.w = w
        self.c = c 
        self.n = n
        self.m = m 
        self.d = d

        self.output = output
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self.path = choose_software_by_platform('partig')
        if not cmd_exists(self.path):
            logger.error(f'No such command of `{self.path}`.')
            sys.exit()

    def partig(self):
        """
        get partig record
        """
        cmd = [self.path, f'-k{self.k}', f'-w{self.w}',
                f'-m{self.m}', f'-c{self.c}', f'-n{self.n}',
                f'-d{self.d}', self.fasta]
        
        logger.info('Calculating the similarity of sequences ...')
        pipelines = []
        try:
            pipelines.append(
                Popen(cmd, stdout=open(self.partig_res, 'w'),
                stderr=open(f"{self.log_dir}/alleles.core.log", 'w'),
                bufsize=-1)
            )
            pipelines[-1].wait()
        except:
            raise Exception('Failed to execute command:' 
                                f'\t{cmd}.')
        finally:
            for p in pipelines:
                if p.poll() is None:
                    p.terminate()
                else:
                    assert pipelines != [], \
                        "Failed to execute command, please check log."
                if p.returncode != 0:
                    raise Exception('Failed to execute command:' 
                                        f'\t{" ".join(cmd)}.')

    def get_pr(self):
        self.pr = PartigRecords(self.partig_res)
        self.pr.convert(self.fasta)

    def to_alleletable(self, fmt="allele2"):
        with open(self.output, 'w') as output:
            self.pr.to_alleletable(self.fasta, output, fmt)
        logger.info(f'Successful output allele table in `{self.output}`.')

    def run(self):
        self.partig()
        self.get_pr()
        self.to_alleletable()

    
class GmapAllele:
    def __init__(self, 
                    fasta,
                    cds, 
                    bed,
                    ploidy,
                    skip_index=False,
                    output='Allele.ctg.table',
                    output_gff3='gmap.gff3',
                    threads=4, 
                    log_dir='logs'):
    
        self.cds = cds
        self.bed = bed
        self.fasta = fasta
        self.ploidy = ploidy
        self.skip_index = skip_index
        self.output = output
        self.output_gff3 = output_gff3
        self.threads = threads
        self.dbname = 'DB'

        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self.genome_size = get_genome_size(self.fasta)
        if self.genome_size < 4e9:
            self.gmap_cmd = 'gmap'
        else:
            self.gmap_cmd = 'gmapl'
        
        if not cmd_exists('gmap'):
            logger.error('No such command of `gmap`.')
            sys.exit()
        
    def build_index(self):
        cmd = ['gmap_build', '-D', '.', '-db', 
                self.dbname,
                self.fasta]
        
        run_cmd(cmd, log=f'{str(self.log_dir)}/gmap_build.log', out2err=True)
    
    def gmap(self):
        cmd = [self.gmap_cmd, 
                '-d', self.dbname,
                '-D', '.', 
                '-f', '2', 
                '-n', str(self.ploidy),
                '-t', str(self.threads),
                self.cds
                ]
        
        errout = open(f'{self.log_dir}/gmap.log', 'w')
        out = open(f'{self.output_gff3}', 'w')

        logger.info("Running command:")
        logger.info("\t" + " ".join(cmd))
        pipelines = []
        try:
            pipelines.append(
                Popen(cmd, stderr=errout,
                stdout=out, bufsize=-1)
            )
            pipelines[-1].wait()
        except:
            raise Exception('Failed to execute command:' 
                                f'\t{" ".join(cmd)}.')
        finally:
            for p in pipelines:
                if p.poll() is None:
                    p.terminate()
                
                if p.returncode != 0:
                    raise Exception('Failed to execute command:' 
                                        f'\t{" ".join(cmd)}.')
    
    def to_alleletable(self):
        def get_gene_id(attributes, key='Name'):
            db = dict(map(lambda x: x.split('='), 
                        [i for i in attributes.split(';') if i]))
            return db[key]
        
        
        db = defaultdict(list)
        with open(self.output_gff3) as fp:
            for line in fp:
                if line[0] == '#':
                    continue
                line_list = line.strip().split()
                if line_list[2] != 'gene':
                    continue
                
                gene = get_gene_id(line_list[8])

                db[gene].append(line_list[0])
        
        with open(self.bed) as fp, open(self.output, 'w') as out:
            for line in fp:
                if not line.strip():
                    continue
                line_list = line.strip().split()

                chrom, start, end, gene = line_list[:4]
                if gene not in db:
                    continue
            
                print('\t'.join([chrom, gene] + db[gene]), 
                        file=out)
        
        logger.info(f"Successful output allele table in `{self.output}`")

    def run(self):
        if self.skip_index:
            if not Path(self.dbname).exists():
                logger.error(f'No such file of `{self.dbname}` for gmap.')
                sys.exit()
        else:
            self.build_index()      
        self.gmap()
        self.to_alleletable()
