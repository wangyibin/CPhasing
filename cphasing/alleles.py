#!/usr/bin/env python
# -*- coding:utf-8 -*-


import logging
import os
import os.path as op
import sys
import re

import numpy as np

from collections import defaultdict, OrderedDict
from subprocess import Popen
from pathlib import Path
from pyfaidx import Fasta


from .utilities import (
    cmd_exists, 
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
        logger.info(f'Load file {self._file}.')
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
    
    def convert(self, fasta):
        """
        convert ID to contig name.

        """
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
    

class PartigAllele:
    def __init__(self,
                    fasta,
                    k=19,
                    w=19,
                    m=0.8,
                    output='Allele.ctg.table',
                    log_dir='logs'):
        self.fasta = fasta
        self.prefix = op.basename(fasta).rsplit(".", 1)[0] 
        self.partig_res = f"{self.prefix}.partig.res"

        ## parameters for partig
        self.k = k
        self.w = w
        self.m = m 

        self.output = output
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        if not cmd_exists('partig'):
            logger.error('No such command of `partig`.')
            sys.exit()

    def partig(self):
        """
        get partig record
        """
        cmd = ['partig', f'-k{self.k}', f'-w{self.w}',
                f'-m{self.m}', '-c 100', self.fasta]
        
        logger.info('Calculating the similarity of sequences ...')
        pipelines = []
        try:
            pipelines.append(
                Popen(cmd, stdout=open(self.partig_res, 'w'),
                stderr=open(f"{self.log_dir}/partig.log", 'w'),
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
                                        
        self.pr = PartigRecords(self.partig_res)
        self.pr.convert(self.fasta)

    def to_alleletable(self):
        with open(self.output, 'w') as output:
            i = 0
            for record in self.pr.S:
                i += 1
                print(i, i, record.seqName1, record.seqName2, 
                        sep='\t', file=output)

        logger.info(f'Successful output allele table in `{self.output}`.')

    def run(self):
        self.partig()
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
