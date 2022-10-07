#!/usr/bin/env python
# -*- coding:utf-8 -*-


import logging
import os
import os.path as op
import sys
import re

from collections import defaultdict
from subprocess import Popen
from pathlib import Path

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
    def __init__(self, infile):
        self._file = infile
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
    def seqNames(self):
        names = []
        for name in self.C:
            names.append(name.seqName)
        return names

    @property
    def pairs(self):
        return [(i.seqName1, i.seqName2) for i in self.S]
    
    

class PartigAllele:
    def __init__(self,
                    fasta,
                    output='Allele.ctg.table',
                    log_dir='logs'):
        self.fasta = fasta 
        self.output = output
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        if not cmd_exists('partig'):
            logger.error('No such command of `partig`.')
            sys.exit()
        
    

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
        finally:
            for p in pipelines:
                if p.poll() is None:
                    p.terminate()
    
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
