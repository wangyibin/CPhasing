#!/usr/bin/env python
# -*- coding:utf-8 -*-


import logging
import os
import os.path as op
import sys
import re

import pandas as pd
import numpy as np

from Bio import SeqIO
from collections import defaultdict, OrderedDict
from joblib import Parallel, delayed
from itertools import combinations, permutations
from io import StringIO
from subprocess import Popen
from shutil import which
from pathlib import Path
from pyfaidx import Fasta

from .core import AlleleTable
from .utilities import xopen, list_flatten, read_fasta


from .utilities import (
    cmd_exists, 
    choose_software_by_platform,
    run_cmd, 
    get_genome_size,
    humanized2numeric
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

    def to_alleletable(self, fasta, output, fmt="allele2", k=0, 
                       is_filter=True,
                       filter_value=5000):
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
                # if k:
                #     if int(record.mzShared) * k * float(record.kmerSimilarity) < filter_value:
                #         continue 
                if is_filter:
                    record.kmerSimilarity = float(record.kmerSimilarity)
                    similarity1 = int(record.mzShared) / int(record.mzConsidered1)
                    similarity2 = int(record.mzShared) / int(record.mzConsidered2)
                    if similarity1 > record.kmerSimilarity or similarity2 > record.kmerSimilarity:
                        record.kmerSimilarity = max([similarity1, similarity2])
                        if record.kmerSimilarity  > 1.0:
                            record.kmerSimilarity = 1.0

                    if record.kmerSimilarity < 0.80:
                        continue

                    record.kmerSimilarity = f"{record.kmerSimilarity:.4}"

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
                    filter_value=1000,
                    output='Allele.ctg.table',
                    threads=4,
                    log_dir='logs'):
        self.fasta = fasta
        fasta_prefix = Path(Path(fasta).name).with_suffix("")
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
        self.filter_value = filter_value
        self.slide_window = 100_000_000

        self.output = output
        self.threads = threads
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self.path = choose_software_by_platform('partig')
        if not cmd_exists(self.path):
            logger.error(f'No such command of `{self.path}`.')
            sys.exit()

    def split_contig(self):
        """
        Split contig when contig > 130 Mb to fix error of partig
        """
        from .utilities import get_contig_length

        contigsize = get_contig_length(self.fasta)
        if any(list(map(lambda x: x > 130_000_000, contigsize.values()))):
            
            cmd = ["cphasing-rs", "slidefastq", "-w", str(self.slide_window), str(self.fasta),
                    "-o", f"tmp.{self.prefix}.slide.fasta"]
            flag = run_cmd(cmd, log=f"logs/{self.prefix}.fasta.slide.log")
            assert flag == 0, "Failed to execute command, please check log."
            
            self.fasta = f"tmp.{self.prefix}.slide.fasta"

            return True 
    
        else:
            return False
        
    def partig(self):
        """
        get partig record
        """
        cmd = [self.path, f'-k{self.k}', f'-w{self.w}',
                f'-m{self.m}', f'-c{self.c}', f'-n{self.n}',
                f'-d{self.d}', f'-t{self.threads}', self.fasta]
        
        logger.info('Calculating the similarity of sequences ...')
        pipelines = []
        try:
            pipelines.append(
                Popen(cmd, stdout=open(self.partig_res, 'w'),
                stderr=open(f"{self.log_dir}/{self.prefix}.alleles.core.log", 'w'),
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
        
        if self.is_split:
            with open(f"tmp.{self.prefix}.allele.table", 'w') as output:
                self.pr.to_alleletable(self.fasta, output, fmt, self.k, 
                                       is_filter=True)
            
            header_lines = []
            lines = []
            with open(f"tmp.{self.prefix}.allele.table", 'r') as fp:
                for line in fp:
                    if line.strip().startswith("#"):
                        header_lines.append(line.strip().split())
                    else:
                        lines.append(line.strip().split())
            header = pd.DataFrame(header_lines).astype({1: np.int64, 2: np.int64, 3: np.int64})
            if len(lines) != 0:
                data = pd.DataFrame(lines).astype({4: np.int64, 5: np.int64,
                                                6: np.int64, 7: np.float64,
                                                8: np.int8})
                
                data[2] = data[2].str.rsplit('_', n=1).map(lambda x: x[0])
                data[3] = data[3].str.rsplit('_', n=1).map(lambda x: x[0])
                data = data[(data[4] > 1000) & (data[5] > 1000)]
                max_idx = data.groupby([2, 3, 8], as_index=False).agg({6: "sum"})
                
                max_idx = max_idx.loc[max_idx.groupby([2, 3])[6].idxmax()].drop(6, axis=1)
                max_idx.set_index([2, 3], inplace=True)
                data = data.groupby([2, 3]).agg(
                                                {4: "sum", 5: "sum",
                                                6: "sum", 7: "mean", 
                                                })

                data = pd.concat([data, max_idx], axis=1).reset_index()
                data = data[data[6] > 500]
                data = data[data[7] > 0.8]
                data[7] = data[7].map(lambda x: f'{x: .4f}')
                header[0] = header[0].str.rsplit('_', n=1).map(lambda x: x[0])
                header = header.groupby([0]).agg("sum").reset_index()

                header.to_csv(self.output, sep=' ', header=None, index=None)
                data.reset_index().reset_index().to_csv(self.output, mode='a', sep='\t', header=None, index=None)    
            else:
                with open(self.output, 'w') as out:
                    pass 

            if Path(f"tmp.{self.prefix}.allele.table").exists():
                    os.remove(f"tmp.{self.prefix}.allele.table")
        else:
            with open(self.output, 'w') as output:
                self.pr.to_alleletable(self.fasta, output, fmt, self.k, 
                                       is_filter=True, filter_value=self.filter_value)
        logger.info(f'Successful output allele table in `{self.output}`.')

    def run(self):
        self.is_split = self.split_contig()
        self.partig()
        self.get_pr()
        self.to_alleletable()

        if self.is_split:
            if Path(f"tmp.{self.prefix}.slide.fasta").exists():
                os.remove(f"tmp.{self.prefix}.slide.fasta")
            if Path(f"tmp.{self.prefix}.slide.fasta.fai").exists():
                os.remove(f"tmp.{self.prefix}.slide.fasta.fai")
    
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
        cmd = ['gmap_build', '-D', '.', '-d', 
                self.dbname,
                self.fasta]
        
        flag = run_cmd(cmd, log=f'{str(self.log_dir)}/gmap_build.log', out2err=True)
        assert flag == 0, "Failed to execute command, please check log."

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
        logger.info("\t" + " ".join(cmd) + f"> {self.output_gff3}" )
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

class AlignmentAlleles:
    PAF_HADER = ["contig1", "length1", "start1", "end1", "strand",
                 "contig2", "length2", "start2", "end2", "mismatch",
                 "matches", "mapq", "identity"]
    PAF_HADER2 = ["contig2", "length2", "start2", "end2", "strand",
                 "contig1", "length1", "start1", "end1", "mismatch",
                 "matches", "mapq", "identity"]
    def __init__(self, fasta, ploidy, k=19, w=19,
                  s="5k", l="10k", p=98, 
                 H=10.0, no_raw_table=True,
                 output=None, log_dir="logs", threads=4):
        self.file = fasta 
        self.fasta = fasta 
        
        self.ploidy = ploidy 
        self.k = k
        self.w = w 
        self.s = s 
        self.l = l 
        self.p = p
        self.H = H 
        
        self.no_raw_table = no_raw_table
        self.threads = threads
        
        self.n = self.ploidy - 1 

        if not cmd_exists("wfmash"):
            logger.error(f'No such command of `wfmash`.')
            sys.exit()

        fasta_prefix = Path(Path(fasta).name).with_suffix("")
        while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
            fasta_prefix = fasta_prefix.with_suffix("")

        self.prefix = fasta_prefix

        self.paf = f"{self.prefix}.selfalign.paf"

        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        if not output:
            self.output = f"{self.prefix}.allele.table"
        else:
            self.output = output

    def align(self):
        # cmd = ["wfmash", self.fasta, "-m", "-t", str(self.threads),
        #        "-H", "0.001", "-n", f"{self.n}", "-l", f"{self.l}",
        #         "-L",
                #  "-Y", "#", "-s", f"{self.s}", "-p", f"{self.p}" ]
        cmd = ["wfmash", self.fasta, "-m", "-t", str(self.threads),
               "-H", f"{self.H}", "-n", f"{self.n}", "-l", f"{self.l}",
                "-k", f"{self.k}", 
                  "-Y", "'#'", "-s", f"{self.s}", "-p", f"{self.p}",
                   ">", self.paf, f"2>{str(self.log_dir)}/alleles.align.log" ]
        
        logger.info("Self mapping ...")
        os.system(" ".join(cmd))

        
        # pipelines = []
        # try:
        #     pipelines.append(
        #         Popen(cmd, stdout=open(self.paf, "w"),
        #               stderr=open(f"{self.log_dir}/alleles.align.log", "w"),
        #               bufsize=-1)
        #     )
        #     pipelines[-1].wait()
        # except:
        #     raise Exception('Failed to execute command:' 
        #                         f'\t{cmd}.')
        # finally:
        #     for p in pipelines:
        #         if p.poll() is None:
        #             p.terminate()
        #         else:
        #             assert pipelines != [], \
        #                 "Failed to execute command, please check log."
        #         if p.returncode != 0:
        #             raise Exception('Failed to execute command:' 
        #                                 f'\t{" ".join(cmd)}.')
                
    def read_paf(self):
        logger.info(f"Load alignments results `{self.paf}`")
        df = pd.read_csv(self.paf, sep='\t', header=None, usecols=range(13),
                         names=self.PAF_HADER, index_col=None)
        df['identity'] = df['identity'].map(lambda x: x.replace("id:f:", ""))
        df1 = df[df['contig1'] < df['contig2']]
        df2 = df[df['contig1'] > df['contig2']]
        df2.columns = self.PAF_HADER2
        df = pd.concat([df1, df2], axis=0)

        # df = df[df['matches'] >= 50000]
        df = df.sort_values(['contig1', 'start1', 'contig2', 'start2'])
        self.paf_df = df 

        return df 
        
    def remove_subset(self, allele_lines):
        allele_lines = [frozenset(allele_line) for allele_line in allele_lines]
        new_allele_lines = set(allele_lines)
        
        for allele_line in allele_lines:
            if len(allele_line) < self.ploidy:
                subsets = {other for other in allele_lines if other != allele_line and allele_line.issubset(other)}
                if subsets:
                    new_allele_lines.difference_update({allele_line})

        return [list(allele_line) for allele_line in new_allele_lines]

    def pairwise(self, allele_lines):
        
        pairs = set()
        for allele_line in allele_lines:
            pairs.update(combinations(allele_line, 2))
        
        return pairs

    def paf2allele(self):
        logger.info("Convert alignments to allele table ...")
        allele_lines = []
        ss = humanized2numeric(self.s)

        

        def func(contig, tmp_df):

            def process_row(row, ss, i):
                start1, end1 = row['start1'], row['end1']
                contig2 = row['contig2']
                a.loc[(start1 + ss//5) / ss: (end1 - ss/5)/ss, i] = contig2


            tmp_allele_lines = []
            if len(tmp_df) <= self.n:
                tmp_allele_lines.append([contig, *tmp_df['contig2']])
            else:
                length = tmp_df.head(1)['length1'].values[0]
                
                a = np.zeros((length // ss , 1))
                a = pd.DataFrame(a)
                a[0] = contig 
                i = 0

                for idx, row in tmp_df.iterrows():
                    i += 1
                    start1, end1 = row['start1'], row['end1']
                    contig2 = row['contig2']
                    a.loc[(start1 + ss//5) / ss: (end1 - ss/5)/ss, i] = contig2
                # tmp_df.apply(process_row, axis=1, args=(ss, i))
                a = a.drop_duplicates()

                res = a.apply(lambda x: x.dropna().values.tolist(), axis=1)
                res = list(filter(lambda x: len(x) > 1, res))
                tmp_allele_lines.extend(res)

            return tmp_allele_lines
        

        args = [(contig, tmp_df) for contig, tmp_df in self.paf_df.groupby('contig1')]
        allele_lines = Parallel(n_jobs=min(self.threads, len(args)))(
            delayed(func)(i, j) for i, j in args
        )
        allele_lines = list_flatten(allele_lines)
        allele_lines = self.remove_subset(allele_lines)
        # self.contig_pairs = self.pairwise(allele_lines)
         

        res_df = pd.DataFrame(allele_lines)
        res_df = res_df.drop_duplicates().sort_values(by=0).reset_index(drop=True)
        
        return res_df

    def export_allele2(self):
        df = self.paf_df[['contig1', 'contig2', 'length1', 'length2', 'matches', 'identity', 'strand']]
        df = df.drop_duplicates(['contig1', 'contig2'])
        df2 = df.copy()
        df2.columns = ['contig2', 'contig1', 'length1', 'length2', 'matches', 'identity',  'strand']
        df = pd.concat([df, df2], axis=0)
        df = df.dropna()
        df['matches'] = df['matches'].astype(int)

        # df = df.set_index(['contig1', 'contig2'])
        
        # unmap_contig_pairs = set(self.contig_pairs).difference(set(map(tuple, df.index.tolist())))
        # self.supplementary_align(unmap_contig_pairs)

        # df = df.reindex(self.contig_pairs)

        df = df.reset_index(drop=True).reset_index().reset_index()
        
        df['strand'] = df['strand'].map(lambda x: 1 if x == "+" else -1)
        if not self.output:
            output = f"{self.prefix}.allele.table"
        else:
            output = self.output
        df.to_csv(output, sep='\t', 
                  index=None, header=None)
        logger.info(f"Export the allele table `{output}`")
        
        return df 

    def supplementary_align(self, unmap_contig_pairs):
        fasta_db = read_fasta(self.fasta)
        with open(f"{self.prefix}.unmap.contig1.fasta", "w") as out1, \
                open(f"{self.prefix}.unmap.contig2.fasta", "w") as out2:
        
            for contig1, contig2 in unmap_contig_pairs:
                if contig1 > contig2:
                    continue
                seq1, seq2 = fasta_db[contig1], fasta_db[contig2]
                print(f">{contig1}\n{seq1}", file=out1)
                print(f">{contig2}\n{seq2}", file=out2)


    def save(self, output):
        self.allele_table.reset_index().reset_index().to_csv(
            output, header=None, index=None, sep='\t'
        )
        logger.info(f"Successful output raw allele table in `{self.output}`")

    def run(self):
        # if not Path(self.paf).exists():
        #     logger.info(f"No such file of `{self.paf}`")
        self.align()

        self.paf_df = self.read_paf()

        if not self.no_raw_table:
            self.allele_table = self.paf2allele()
            output = self.output.replace(".allele.table", ".raw.allele.table")
            self.save(output)
       
        self.export_allele2()


class AlignmentAlleles2:
    PAF_HADER = ["contig1", "length1", "start1", "end1", "strand",
                 "contig2", "length2", "start2", "end2", "matches",
                 "alignment", "mapq", "dv"]
    META = {
        "contig1": "object",
        "length1": int,
        "start1": int,
        "end1": int,
        "strand": "category",
        "contig2": "object",
        "length2": int,
        "start2": int,
        "end2": int,
        "matches": int,
        "aligns": int,
        "mapq": int,
        "dv": "object"
    }
    def __init__(self, fasta, k=19, w=19, c=100,
                 output=None, log_dir="logs", threads=4):
                 
        self.file = fasta
        self.fasta = fasta

        self.k = k
        self.w = w

        self.c = c 

        self.threads = threads
        if not cmd_exists("minigraph"):
            logger.error(f'No such command of `minigraph`.')
            sys.exit()

        fasta_prefix = Path(Path(fasta).name).with_suffix("")
        while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
            fasta_prefix = fasta_prefix.with_suffix("")

        self.prefix = fasta_prefix

        self.paf = f"{self.prefix}.selfalign.paf"
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        if not output:
            self.output = f"{self.prefix}.allele.table"
        else:
            self.output = output

    def align(self):
        cmd = ["minigraph", "-x", "asm", "-k", str(self.k), "-w", str(self.w),
               "-N", str(self.c), "-t", str(self.threads), "--secondary=yes",
                 self.fasta, self.fasta,
              f"2>{str(self.log_dir)}/alleles.align.log", ">", self.paf ]
        
        logger.info("Self mapping ...")
        os.system(" ".join(cmd))
    
    def read_paf(self):
        logger.info(f"Load alignments results `{self.paf}`")
        df = pd.read_csv(self.paf, sep='\t', header=None, usecols=list(range(12)) + [16],
                         names=self.PAF_HADER, index_col=None, dtype=self.META)

        df = df[df['contig1'] < df['contig2']]
        df = df[df['matches'] > 5000]
        df['dv'] = df['dv'].str.split(":").map(lambda x: x[-1]).astype(np.float64)
        df['identity'] = df['matches'] / df['alignment']
        df.reset_index(drop=True, inplace=True)


        def func(df):
            strand = df.groupby(['strand']).apply(lambda x: x['matches'].sum()).idxmax()
            tmp_df = df.sort_values('matches', ascending=False).head(1)
            tmp_df['dv'] = df['dv'].mean()
            tmp_df['matches'] = df['matches'].sum()
            tmp_df['alignment'] = df['alignment'].sum()
            tmp_df['start1'] = df['start1'].min()
            tmp_df['end1'] = df['end1'].max()
            tmp_df['mapq'] = df['mapq'].max()
            tmp_df['start2'] = df['start2'].min()
            tmp_df['end2'] = df['end2'].max()
            tmp_df['identity'] = tmp_df['matches'] / tmp_df['alignment']
            tmp_df['strand'] = strand
            
            return tmp_df

        df = df.groupby(['contig1', 'contig2']).apply(func).reset_index(drop=True)
        df['identity'] = 1- df['dv']

        return df 
    

    def export_allele(self):
        
        df = self.paf_df[['contig1', 'contig2', 'length1', 'length2', 'matches', 'strand', 'identity']]
        df2 = df.copy()
        df2.columns = ['contig2', 'contig1', 'length2', 'length1', 'matches', 'strand', 'identity']
        df = pd.concat([df, df2], axis=0)
        df = df.dropna()


        df = df.reset_index(drop=True).reset_index().reset_index()
        df['strand'] = df['strand'].map(lambda x: 1 if x == "+" else -1)
        if not self.output:
            output = f"{self.prefix}.allele.table"
        else:
            output = self.output

        df.to_csv(output, sep='\t',
                    index=None, header=None)
        logger.info(f"Export the allele table `{output}`")


    def run(self):
        self.paf_df = self.read_paf()
        self.export_allele() 


class ANIAlleles:
    def __init__(self, fasta, ploidy, percent,
                 output=None, log_dir="logs", threads=4):
        self.file = fasta 
        self.fasta = fasta 
        
        self.ploidy = ploidy 
        self.percent = percent  
        if not cmd_exists("skani"):
            logger.error(f'No such command of `skani`.')
            sys.exit()

        self.threads = threads

        fasta_prefix = Path(Path(fasta).name).with_suffix("")
        while fasta_prefix.suffix in {".fasta", "gz", "fa", ".fa", ".gz"}:
            fasta_prefix = fasta_prefix.with_suffix("")

        self.prefix = fasta_prefix

        self.ani_table = f"{self.prefix}.ani"

        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        if not output:
            self.output = f"{self.prefix}.allele.table"
        else:
            self.output = output


    def align(self):
        cmd = ["skani", "dist", "-t", str(self.threads), 
               "--ri", "-r", self.fasta, "--qi", "-q", self.fasta,
               "-s", str(self.percent), "-n", str(self.ploidy * 10),
               "-o", self.ani_table]
        
        logger.info("Self comparison ...")
        flag = run_cmd(cmd, log=f'{str(self.log_dir)}/alleles.ani.log', out2err=True)
        assert flag == 0, "Failed to execute command, please check log."

    def read_ani_table(self):
        logger.info(f"Load alignments results `{self.ani_table}`")
        df = pd.read_csv(self.ani_table, sep='\t', header=None, usecols=range(2, 7),
                         skiprows=1,
                         names=['identity', 'frac1', 'frac2', 'contig1', 'contig2',], index_col=None)

        df = df[df['contig1'] != df['contig2']]
        df2 = df.copy()
        df2.columns = ['identity', 'frac2', 'frac1', 'contig2', 'contig1']
        df = pd.concat([df, df2], axis=0)
        df.drop_duplicates(inplace=True)
        df.query('frac1 > 1.0 and frac2 > 1.0', inplace=True)

        df = df.sort_values(['contig1', 'contig2'])
        self.ani_df = df 


        return df
    
    def to_alleletable(self):
        
        df = self.ani_df[['contig1', 'contig2', 'frac1', 'frac2', 'identity']]
        df = df.dropna()
        df['identity'] = df['identity'].map(lambda x: f"{x/100:.4}")
       
        df2 = df.copy()
        df2.columns = ['contig2', 'contig1', 'frac2', 'frac1', 'identity']
        df = pd.concat([df, df2], axis=0)
        df = df[['contig1', 'contig2', 'frac1', 'frac2', 'identity', 'identity']]
        df['strand'] = 1
        df = df.reset_index(drop=True).reset_index().reset_index()
        if not self.output:
            output = f"{self.prefix}.allele.table"
        else:
            output = self.output
        df.to_csv(output, sep='\t', 
                  index=None, header=None)
        logger.info(f"Export the allele table `{output}`")
        
        return df

    def run(self):
        self.align()
        self.read_ani_table()
        self.to_alleletable()


def filter_high_similarity_contigs(alleletable, contacts, min_values, output):
    at = AlleleTable(alleletable, fmt='allele2', sort=False)
    tmp_data = at.data[at.data[1] < at.data[2]]

    db = pd.read_csv(contacts, sep='\t', index_col=(0, 1), header=None,
                     names=['contig1', 'contig2', 'count'], 
                     dtype={'contig1': 'U', 'contig2': 'U'}).query('contig1 <= contig2').to_dict()['count']
    
    def func(row):
        contig1, contig2 = row[1], row[2]
        try:
            cis1 = db[(contig1, contig1)]
            cis2 = db[(contig2, contig2)]
            trans = db[(contig1, contig2)]
        except KeyError:
            return None

        value = trans / ((cis1 * cis2) ** (1/2))

        return value 

    tmp_res = tmp_data.apply(func, axis=1).dropna() >= min_values
    res = tmp_data.loc[tmp_res.index[tmp_res]][[1, 2]]
  
    res = set(res.values.flatten().tolist())

    print("\n".join(res), file=output)
