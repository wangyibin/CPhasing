#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the algorithms of travelling salesman problem
"""

import argparse
import logging
import glob
import os
import os.path as op
import sys
import tempfile
import warnings
warnings.simplefilter('ignore')

import cooler
import random
import igraph as ig
import gc
import numpy as np
import networkx as nx 
import networkx.algorithms.approximation as nx_app 
import multiprocessing
import igraph as ig 
import pandas as pd
import polars as pl
import re
import scipy
import shutil

from collections import defaultdict, OrderedDict, Counter
from joblib import Parallel, delayed, cpu_count
from itertools import combinations, permutations, product
from pathlib import Path
from pytools import natsorted
from string import ascii_uppercase, ascii_lowercase
from scipy.sparse import tril, coo_matrix
from subprocess import Popen

from ..agp import agp2fasta, agp2tour
from ..algorithms.hypergraph import IRMM
from .._config import *
from ..core import (
                    AlleleTable, 
                    CountRE, 
                    ClusterTable, 
                    Tour
                    )
from ..utilities import (
    decompress_cmd,
    cmd_exists,
    choose_software_by_platform, 
    run_cmd,
    read_fasta
    )


logger = logging.getLogger(__name__)

class HaplotypeAlign:
    """

    Align haplotype into parallel line

    Params:
    --------
    at: AlleleTable
        AlleleTable with `allele2` format
    hap_tour_db: dict
        haplotype tour dict 
        {"Chr01": [Tour, Tour, Tour ...],
            "Chr02": [Tour, Tour, Tour ...]}
    """
    def __init__(self, at: str, tour_list: list,
                 workdir=".", threads: int = 4):
        at = AlleleTable(at, sort=False, fmt='allele2')
        self.tour_list = tour_list 
        self.allele_data = at.data.set_index([1, 2])
        self.hap_tour_db = self.get_hap_tour_db(self.tour_list) 
        self.workdir = workdir
        self.threads = threads 

    def get_hap_tour_db(self, tour_list):
        db = defaultdict(list)
        for tour in tour_list:
            try:
                hap, idx = tour.rsplit("g", 1)
            except:
                if "single" not in db:
                    db['single'] = []    
                db['single'].append(Tour(tour))
                continue

            db[hap].append(Tour(tour))

        return db 


    @staticmethod
    def align(tour1: Tour, tour2: Tour, data: pd.DataFrame, workdir: str):
        """
        align two haplotype contigs by allele table data

        Params:
        ----------
        tour1: Tour

        tour2: Tour

        data: pd.DataFrame
            1       2  mzShared  strand
            ctg1    ctg3  200   1
            ctg2    ctg4  20    -1

        workdir: str
        
            
        """
        hap1, hap2 = tour1.to_dict(1), tour2.to_dict(1)

        contig_pairs = {(contig1, contig2): hap1[contig1] * hap2[contig2]
                            for contig1, contig2 in product(hap1.keys(), hap2.keys())}

        # contig_pairs_df = pd.DataFrame(list(contig_pairs.items()))
        # contig_pairs_df.columns = [0, "value"]
        # contig_pairs_df[[1, 2]] = pd.DataFrame(contig_pairs_df[0].tolist())
        # contig_pairs_df = contig_pairs_df.drop(0, axis=1).set_index([1, 2])
        contig_pairs_df = pd.DataFrame.from_dict(contig_pairs, orient='index', columns=['value'])
        contig_pairs_df.index = pd.MultiIndex.from_tuples(contig_pairs_df.index)
        tmp_df = data.reindex(list(contig_pairs.keys())).dropna()

        tmp_df = tmp_df[tmp_df['similarity'] > 0.85]
        tmp_df['strand'] = tmp_df['strand'].astype(int)
        tmp_df['mzShared'] = tmp_df['mzShared'].astype(int)

        tmp_df['strand2'] = contig_pairs_df.loc[tmp_df.index]
        # tmp_df = tmp_df[tmp_df['mzShared'] > 2000]
        # tmp_df = tmp_df.assign(value=lambda x: x['mzShared'] * x['strand'] * x['strand2'])
        tmp_df['value'] = tmp_df['mzShared'] * tmp_df['strand'] * tmp_df['strand2']
        # tmp_df = tmp_df.drop(['mzShared', 'strand', 'strand2'], axis=1)

        if tmp_df['value'].sum() < 0:
            tour2.reverse()
            # tour2.rm()
            # tour2.backup("before_reverse")
            tour2.save(f"{workdir}/{tour2.filename}")    
    
    def run(self):
        logger.info("Adjust the tours to parallel among different haplotypes.")

        args = [(tours[0], tours[i], self.allele_data, self.workdir) 
                    for hap, tours in self.hap_tour_db.items() 
                    for i in range(1, len(tours))]
        
        # for item in args:
        #     self.align(*item)
        total_machine_cpu = cpu_count()
        
        try:
            Parallel(n_jobs=min(len(args), min(total_machine_cpu, self.threads)), backend="multiprocessing")(
                        delayed(self.align)(*a) for a in args
            )
        except ValueError:
            logger.warning("Failed to run HaplotypeAlign. skipped")

class Rename:
    PAF_HADER = ["contig1", "length1", "start1", "end1", "strand",
                 "contig2", "length2", "start2", "end2", "mismatch",
                 "matches", "mapq", "identity"]

    def __init__(self, ref, fasta, agp, 
                    hap_pattern=r'(Chr\d+)g(\d+)',
                    hap_aligned=True,
                    output="groups.renamed.agp",
                    suffix_style="number",
                    threads=4, log_dir="logs",
                    tmp_dir="rename_tmp"):
        

        self.ref = Path(ref).absolute()
        self.fasta = Path(fasta).absolute()
        self.agp = Path(agp).absolute()
        self.paf = "ref.align.paf"
        if hap_pattern:
            self.hap_pattern = r"{}".format(hap_pattern)
        else:
            self.hap_pattern = r'(Chr\d+)g(\d+)'
        self.hap_aligned = hap_aligned
        self.suffix_style = suffix_style
        self.output = Path(output).absolute()
        self.output_prefix = Path(self.output.stem).absolute()
        self.threads = threads 
        self.tmp_dir = tmp_dir
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        if not cmd_exists("wfmash"):
            logger.error(f'No such command of `wfmash`.')
            sys.exit()
    
    def get_contigs(self):
        agp_df = agp2tour(self.agp, "raw_tour")
        fasta_db = read_fasta(self.fasta)
        agp_df = agp_df.loc[agp_df['chrom'] != agp_df['id']]
        agp_df = agp_df.reset_index()
        agp_df['chrom'] = agp_df['chrom'].astype('object')

        hap_db = OrderedDict()

        contigs_db = OrderedDict()
        for chrom, cluster in agp_df.groupby('chrom', sort=False):
            if cluster.empty:
                continue
            if chrom == cluster['id'].values[0]:
                continue
            
            try:
                if self.hap_aligned:
                 
                    hap = re.match(self.hap_pattern, chrom)
                    hap = hap.groups() if hap is not None else [chrom]
                else:
                    raise AssertionError
            except:
                hap = [chrom]
            
            if hap[0] not in hap_db:
                hap_db[hap[0]] = []
            
            hap_db[hap[0]].append(chrom)
            contigs_db[chrom] = cluster['id'].values.tolist()


        first_contigs = OrderedDict()
        first_chrom_to_other_chroms = OrderedDict()
        with open("tmp.first.contigs.fasta", "w") as out:
            for hap in hap_db:
                chroms = hap_db[hap]
                chrom = chroms[0]
                first_chrom_to_other_chroms[chrom] = chroms
                first_contigs[chrom] = contigs_db[chrom]
                
                for contig in first_contigs[chrom]:
                    print(f">{contig}\n{fasta_db[contig]}", file=out)
        
        self.first_contig = "tmp.first.contigs.fasta"
        
        return "tmp.first.contigs.fasta", first_contigs, first_chrom_to_other_chroms

    def read_paf(self, first_contigs_db):
        logger.info(f"Load alignments results `{self.paf}`")
        df = pd.read_csv(self.paf, sep='\t', header=None, usecols=range(13),
                         names=self.PAF_HADER, index_col=None)
        df['identity'] = df['identity'].map(lambda x: x.replace("id:f:", ""))
    
        df = df.sort_values(['contig2', 'start2'])

        contig_to_chrom_db = dict()
        for chrom, contigs in first_contigs_db.items():
            for contig in contigs:
                contig_to_chrom_db[contig] = chrom
        
        df['group'] = df['contig1'].map(contig_to_chrom_db.get)
        df['strand'] = df['strand'].map(lambda x: 1 if x == "+" else -1)
        self.paf_df = df 

        return df 

    def global_align(self):
        logger.info("Mapping ...")
        cmd = ["wfmash", str(self.ref), self.first_contig, 
                    "-m", "-s", "50k", "-l", "250k", 
                    "-p", "90",
                    "-t", str(self.threads)]
        
        pipelines = []
        try:
            pipelines.append(
                Popen(cmd, stdout=open(self.paf, "w"),
                      stderr=open(f"../{self.log_dir}/ref.align.log", "w"),
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

    def align(self, paf, first_chrom_to_other_chroms, workdir):
        res = OrderedDict()
        for group, tmp_df in paf.groupby('group'):
            assign_ref_group = tmp_df.groupby('contig2')['matches'].sum().idxmax()
            tmp_df2 = tmp_df[tmp_df['contig2'] == assign_ref_group]

            res[group] = (assign_ref_group, tmp_df2[['contig1', 'strand', 'matches']])

        chrom_to_first_chrom = {}
        for chrom in first_chrom_to_other_chroms:
            for chrom2 in first_chrom_to_other_chroms[chrom]:
                chrom_to_first_chrom[chrom2] = chrom 

        tours = glob.glob("raw_tour/*.tour")
        tours = list(map(Tour, tours))

        renamed_chrom = []
        hap_idx_db = {}
        tour_results = {}
        for chrom1 in first_chrom_to_other_chroms:
            try:
                tmp_res = res[chrom1]
            except KeyError:
                for chrom2 in first_chrom_to_other_chroms[chrom1]:
                    tour_results[chrom2] = (f"unrenamed_{chrom2}", f"unrenamed_{chrom2}", Tour(f"raw_tour/{chrom2}.tour"))
                    renamed_chrom.append(f"unrenamed_{chrom2}")
                continue
            
            tour = Tour(f"raw_tour/{chrom1}.tour")
            
            tour_dict = tour.to_dict()
            tmp_df = tmp_res[1]
            tmp_df['tour_strand'] = tmp_df['contig1'].map(tour_dict.get)
            tmp_df['tour_strand'] = tmp_df['tour_strand'].map(lambda x: 1 if x == "+" else -1)

            v = tmp_df.eval('strand * matches * tour_strand').sum()

            for chrom2 in first_chrom_to_other_chroms[chrom1]:
                tour = Tour(f"raw_tour/{chrom2}.tour")
                if v < 0:
                    logger.debug(f"Reverse {tour.group}")
                    tour.reverse()
                # hap = tour.group.rsplit("g", 1)
                hap = re.match(self.hap_pattern, tour.group)
                hap = hap.groups() if hap is not None else [] 

                if len(hap) > 1 and self.hap_aligned:
                    
                    if self.suffix_style == "lowerletter":
                        new_chrom = f"{tmp_res[0]}{ascii_lowercase[int(hap[1]) - 1]}"
                    elif self.suffix_style == "upperletter":
                        new_chrom = f"{tmp_res[0]}{ascii_uppercase[int(hap[1]) - 1]}"
                    else:
                        new_chrom = f"{tmp_res[0]}g{hap[1]}"
                    tour_results[chrom2] = (tmp_res[0], new_chrom, tour)
                    
                    renamed_chrom.append(f"{tmp_res[0]}")
                else:
                    if f"{tmp_res[0]}" not in hap_idx_db:
                        hap_idx_db[tmp_res[0]] = 1
                    
                        renamed_chrom.append(f"{tmp_res[0]}")
                    else:
                        hap_idx_db[tmp_res[0]] += 1
                        renamed_chrom.append(f"{tmp_res[0]}")
    
                    if self.suffix_style == "lowerletter":
                        new_chrom = f"{tmp_res[0]}{ascii_lowercase[int(hap_idx_db[tmp_res[0]]) - 1]}"
                    elif self.suffix_style == "upperletter":
                        new_chrom = f"{tmp_res[0]}{ascii_uppercase[int(hap_idx_db[tmp_res[0]]) - 1]}"
                    else:
                        new_chrom = f"{tmp_res[0]}g{int(hap_idx_db[tmp_res[0]])}"

                    tour_results[chrom2] = (tmp_res[0], new_chrom, tour)

        
        renamed_counts = Counter(renamed_chrom)
        with open(f"{self.output_prefix}.rename.list", "w") as out:
            for chrom, (ref_chrom, new_chrom, tour) in tour_results.items():
                if renamed_counts[ref_chrom] == 1:
                    logger.debug(f"Rename {chrom} to {ref_chrom}")
                    tour.save(f"{workdir}/{ref_chrom}.tour")
                    print(f"{chrom}\t{ref_chrom}", file=out)
                else:
                    logger.debug(f"Rename {chrom} to {new_chrom}")
                    tour.save(f"{workdir}/{new_chrom}.tour")
                    print(f"{chrom}\t{new_chrom}", file=out)


        logger.info("Rename done, output the rename list into "
                    f"{self.output_prefix}.rename.list")

    def run(self):
        from ..cli import build
        tmpDir =  tempfile.mkdtemp(prefix=self.tmp_dir, dir='./')
        logger.info('Working on temporary directory: {}'.format(tmpDir))
    
        os.chdir(tmpDir)
        workdir = os.getcwd()
        
        _, first_contigs_db, first_chrom_to_other_chroms = self.get_contigs()
        self.global_align()

        ref_align_df = self.read_paf(first_contigs_db)

        self.align(ref_align_df, first_chrom_to_other_chroms, workdir)

        try:
            build.main(args=[str(self.fasta), "-oa", self.output, 
                             "-o", f"{self.output_prefix}.fasta"], 
                        prog_name="build")
            
        except SystemExit as e:
            exc_info = sys.exc_info() 
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e

        os.chdir("../")


        shutil.rmtree(tmpDir)


class AllhicOptimize:

    def __init__(self, clustertable, count_re, clm, 
                    allele_table=None,
                    fasta=None, 
                    corrected=False,
                    output="groups.agp", 
                    tmp_dir='scaffolding_tmp', 
                    keep_temp=False,
                    log_dir="logs",
                    threads=4):
        self.clusterfile = str(Path(clustertable).absolute())
        self.clustertable = ClusterTable(clustertable)
        
        self.count_re = CountRE(count_re, minRE=1)
        self.clm_file = str(Path(clm).absolute())
        self.allele_table = str(Path(allele_table).absolute()) if allele_table else None
        self.fasta = Path(fasta).absolute() if fasta else None
        self.corrected = corrected
        self.output = output
        self.tmp_dir = tmp_dir 
        self.delete_temp = False if keep_temp else True 
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.threads = threads 

        self.allhic_path = choose_software_by_platform("allhic")
        os.environ["POLARS_MAX_THREADS"] = str(self.threads)
        
    @property
    def clm(self):
        # df = pd.read_csv(self.clm_file, sep='\t', header=None, index_col=0)
        
        df = pl.scan_csv(self.clm_file, separator='\t', 
                         has_header=False, low_memory=True,
                         dtypes={"column_2": pl.datatypes.UInt32})       

        return df

    @staticmethod
    def extract_count_re(group, contigs, count_re):
        tmp_df = count_re.data.reindex(contigs)
        tmp_df.dropna(inplace=True, axis=0)
        tmp_df = tmp_df.astype({'RECounts': int, 'Length': int})
        tmp_df.to_csv(f"{group}.txt", sep='\t', header=None)

        return f"{group}.txt"

    @staticmethod
    def extract_clm_rust(clm_file, cluster_file, log_dir):
        if clm_file.endswith(".gz"):
            cmd0 = decompress_cmd(clm_file, threads=10)
            cmd = ['cphasing-rs', 'splitclm', '-', cluster_file]
            logger.info("Splitting clm file ...")
            flag = os.system(
                " ".join(cmd0)
                + f" 2>../{log_dir}/clm_decompress.log"
                + " | "
                + " ".join(cmd)
                + f" 2>../{log_dir}/splitclm.log"
            )
        else:
            cmd = ['cphasing-rs', 'splitclm', clm_file, cluster_file]

            flag = run_cmd(cmd, log=f"../{log_dir}/splitclm.log", out2err=True)

        assert flag == 0, "Failed to execute command, please check log."


    @staticmethod
    def extract_clm(group, contigs, clm):
        if len(contigs) > 1:
            contig_pairs = list(permutations(contigs, 2))
            contig_with_orientation_pairs = set(
                f"{pair[0]}{strand1} {pair[1]}{strand2}"
                for pair in contig_pairs
                for strand1, strand2 in [('+', '+'), ('+', '-'), ('-', '+'), ('-', '-')]
            )
            # clm = (pl.scan_csv(clm_file, separator='\t', 
            #                             has_header=False, low_memory=True,
            #                             dtypes={"column_2": pl.datatypes.UInt32})
            #                     .filter(pl.col("column_1")
            #                     .is_in(contig_with_orientation_pairs)).collect()
            #                     .write_csv(f"{group}.clm", separator='\t', include_header=False)   
            #                 )

            # tmp_df = clm.reindex(contig_with_orientation_pairs).dropna().astype({1: int})
            # tmp_df.to_csv(f"{group}.clm", sep='\t', header=None)
            # left_df = pl.DataFrame(pl.Series("column_1", contig_with_orientation_pairs))
            # tmp_df = left_df.join(clm, on=["column_1"], how='left').drop_nulls()

            tmp_df = clm.filter(pl.col("column_1").is_in(contig_with_orientation_pairs)).collect()
            tmp_df.write_csv(f"{group}.clm", separator='\t', include_header=False)
        else:
            with open(f"{group}.clm", 'w') as out:
                pass 

        return f"{group}.clm"

    @staticmethod
    def run_allhic_optimize(allhic_path, count_re, clm):
        cmd = [allhic_path, "optimize", count_re, clm]
        run_cmd(cmd, log=os.devnull, out2err=True)
        return count_re.replace(".txt", ".tour")

    @staticmethod
    def _run(allhic_path, count_re, clm, workdir):
        os.chdir(workdir)
        tmp_res = AllhicOptimize.run_allhic_optimize(allhic_path, count_re, clm)

        return tmp_res
    

    @staticmethod
    def _process_group(args):
        group, contigs, count_re, allhic_path, workdir = args
        tmp_count_re = AllhicOptimize.extract_count_re(group, contigs, count_re)
        tmp_clm = f"{group}.clm"
        return (allhic_path, tmp_count_re, tmp_clm, workdir)

    def run(self):
        from ..cli import build
    
        tmpDir =  tempfile.mkdtemp(prefix=self.tmp_dir, dir='./')
        logger.info('Working on temporary directory: {}'.format(tmpDir))
        os.chdir(tmpDir)
        workdir = os.getcwd()
        args = []

        groups = list(self.clustertable.data.keys())
        args = []
        # for group in groups:
        #     contigs = self.clustertable.data[group]
        #     args.append((group, contigs, clm))
        
        # clms = Parallel(n_jobs=min(len(args), self.threads))(delayed(
        #     AllhicOptimize.extract_clm)(i, j, k) for i, j, k in args
        # )
        args = []
        for group in groups:
            contigs = self.clustertable.data[group]
            args.append((group, contigs, self.count_re, self.allhic_path, workdir))
        
        with multiprocessing.Pool(processes=min(10, self.threads)) as pool:
            args = pool.map(AllhicOptimize._process_group, args)
 
        AllhicOptimize.extract_clm_rust(self.clm_file, self.clusterfile, self.log_dir)
        
        total_cpu_of_machine = cpu_count()
        if self.threads * 2 > total_cpu_of_machine:
            threads = self.threads // 2
            logger.info(f"Use {threads} threads to run optimize.")
        else:
            threads = self.threads

        logger.info("Running scaffolding in each group ...")
        tour_res = Parallel(n_jobs=min(len(args), threads))(delayed(
                        self._run)(i, j, k, l) for i, j, k, l in args)
        
        os.chdir(workdir)
        if self.allele_table and len(tour_res) > 1:
            hap_align = HaplotypeAlign(self.allele_table, tour_res, workdir, self.threads)
            hap_align.run()


        if not self.fasta:
            for file in glob.glob("*.tour"):
                shutil.copy(file, "../")
        else:
            try:
                if self.corrected:
                    build.main(args=[str(self.fasta), "--only-agp", "-oa", self.output, "--corrected"], 
                                prog_name='build')
                else:
                    build.main(args=[str(self.fasta), "--only-agp", "-oa", self.output], prog_name='build')
            except SystemExit as e:
                exc_info = sys.exc_info()
                exit_code = e.code
                if exit_code is None:
                    exit_code = 0
                
                if exit_code != 0:
                    raise e
                
            shutil.copy(f"{self.output}", "../")
            if self.corrected:
                if Path(self.output.replace("agp", "corrected.agp")).exists():
                    shutil.copy(f"{self.output.replace('agp', 'corrected.agp')}", "../")
             
        os.chdir("../")
        if self.delete_temp is True:
            shutil.rmtree(tmpDir)
            logger.info("Removed temporary directory.")

        logger.info("Scaffolding done.")


class AllhicReorient:
    def __init__(self, clm, tour):
        self.clm = clm
        self.tour = tour
    


    def reorient(self):
        pass

    

class HapHiCSort:
    from cphasing.algorithms.HapHiC_sort import parse_arguments, run

    def __init__(self, clustertable,
                    count_re, clm, 
                    split_contacts, 
                    skip_allhic=False,
                    allele_table=None,
                    fasta=None, 
                    corrected=False,
                    output="groups.agp", 
                    tmp_dir='scaffolding_tmp', 
                    keep_temp=False,
                    log_dir="logs",
                    threads=4):
        
        self.clusterfile = str(Path(clustertable).absolute())
        self.clustertable = ClusterTable(clustertable)
        self.count_re_path = Path(count_re).absolute()
        self.count_re = CountRE(count_re, minRE=1)
        
        self.clm_file = str(Path(clm).absolute())
        self.split_contacts = Path(split_contacts).absolute()
        self.skip_allhic = skip_allhic
        self.allele_table = str(Path(allele_table).absolute()) if allele_table else None
        self.fasta = Path(fasta).absolute() if fasta else None
        self.corrected = corrected
        self.output = output
        self.delete_temp = False if keep_temp else True
        self.tmp_dir = tmp_dir 
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.threads = threads 

        self.log_dir = "logs"
        Path(self.log_dir).mkdir(exist_ok=True)

        self.allhic_path = choose_software_by_platform("allhic")

        os.environ["POLARS_MAX_THREADS"] = str(self.threads)
        
    @property
    def clm(self):
        # df = pd.read_csv(self.clm_file, sep='\t', header=None, index_col=0)
        df = pl.scan_csv(self.clm_file, separator='\t', 
                         has_header=False, low_memory=True,
                         dtypes={"column_2": pl.datatypes.UInt32})   
               
        return df

    @staticmethod
    def extract_count_re(group, contigs, count_re):
        tmp_df = count_re.data.reindex(contigs)
        tmp_df.dropna(inplace=True, axis=0)
        tmp_df = tmp_df.astype({'RECounts': int, 'Length': int})
        tmp_df.to_csv(f"{group}.txt", sep='\t', header=None)

        return f"{group}.txt"

    @staticmethod
    def extract_clm_rust(clm_file, cluster_file, log_dir):
        if clm_file.endswith(".gz"):
            cmd0 = decompress_cmd(clm_file, threads=10)
            cmd = ["cphasing-rs", "splitclm", "-", cluster_file]
            logger.info("Splitting clm file ...")
            flag = os.system(
                " ".join(cmd0)
                + f" 2>../{log_dir}/clm_decompress.log"
                + " | "
                + " ".join(cmd)
                + f" 2>../{log_dir}/splitclm.log"
            )
        else:
            cmd = ['cphasing-rs', 'splitclm', clm_file, cluster_file]

            flag = run_cmd(cmd, log=f"../{log_dir}/splitclm.log", out2err=True)
        assert flag == 0, "Failed to execute command, please check log."

    @staticmethod
    def extract_clm(group, contigs, clm):
        if len(contigs) > 1:
            contig_pairs = list(permutations(contigs, 2))
            contig_with_orientation_pairs = set(
                    f"{pair[0]}{strand1} {pair[1]}{strand2}"
                    for pair in contig_pairs
                    for strand1, strand2 in [('+', '+'), ('+', '-'), ('-', '+'), ('-', '-')]
                )
            tmp_df = clm.filter(pl.col("column_1").is_in(contig_with_orientation_pairs)).collect()


            tmp_df.write_csv(f"{group}.clm", separator='\t', include_header=False)
        else:
            with open(f"{group}.clm", 'w') as out:
                pass 

        return f"{group}.clm"

    @staticmethod
    def run_haphic_optimize(total_count_re, split_contacts, tmp_dir, skip_allhic, threads=4, log_dir="logs", ):
        script_realpath = os.path.dirname(os.path.realpath(__file__))
        haphic_sort = f"{script_realpath}/HapHiC_sort.py"
        
        txt = glob.glob(f"{tmp_dir}/*.txt")
        
        if skip_allhic:
            cmd = ["python", 
                    haphic_sort, 
                    "--skip_allhic",
                    "--processes",
                    str(threads),
                    str(total_count_re),
                    str(split_contacts), 
                    tmp_dir]
        else:
            cmd = ["python", 
                haphic_sort, 
                "--processes",
                str(threads),
                str(total_count_re),
                str(split_contacts), 
                tmp_dir]

        
        cmd.extend(txt)
        
        run_cmd(cmd, log=f"../{log_dir}/HapHiC_sort.log", out2err=True)
      
    @staticmethod
    def _run(total_count_re, split_contacts, workdir, skip_allhic=True, threads=4, log_dir="logs"):
        HapHiCSort.run_haphic_optimize(total_count_re, split_contacts, workdir, skip_allhic, threads, log_dir)
    
    @staticmethod
    def _process_group(args):
        group, contigs, count_re, allhic_path, workdir = args
        tmp_count_re = HapHiCSort.extract_count_re(group, contigs, count_re)
        tmp_clm = f"{group}.clm"

        return (allhic_path, tmp_count_re, tmp_clm, workdir)
    
    def run(self):
        from ..cli import build
        tmpDir = tempfile.mkdtemp(prefix=self.tmp_dir, dir='./')
        logger.info('Working on temporary directory: {}'.format(tmpDir))
        os.chdir(tmpDir)
        workdir = os.getcwd()
      
        args = []
        for group in self.clustertable.data.keys():
            contigs = self.clustertable.data[group]
            args.append((group, contigs, self.count_re, self.allhic_path, workdir))
        
        with multiprocessing.Pool(processes=min(4, self.threads)) as pool:
            pool.map(HapHiCSort._process_group, args)
        
        HapHiCSort.extract_clm_rust(self.clm_file, self.clusterfile, self.log_dir)

        total_cpu_of_machine = cpu_count()
        if self.threads * 2 > total_cpu_of_machine:
            threads = total_cpu_of_machine // 2
            logger.info(f"Use {threads} threads to run optimize.")
        else:
            threads = self.threads
        HapHiCSort._run(self.count_re_path, self.split_contacts, "./", 
                            skip_allhic=self.skip_allhic, threads=threads,
                            log_dir=self.log_dir)
        
        os.chdir(workdir)

        tour_res = natsorted(glob.glob("./*.tour"))
        
        if len(tour_res) == 0:
            logger.warning("Failed to run HapHiC sort, please check log file in logs.")
            sys.exit(-1)

  
        if self.allele_table and len(tour_res) >= 2:
            hap_align = HaplotypeAlign(self.allele_table, tour_res, workdir, self.threads)
            hap_align.run()
        
        if not self.fasta:
            for file in glob.glob("*.tour"):
                shutil.copy(file, "../")
        else:
            try:
                if self.corrected:
                    build.main(args=[str(self.fasta), "--only-agp", "-oa", self.output, "--corrected"], 
                                prog_name='build')
                else:
                    build.main(args=[str(self.fasta), "--only-agp", "-oa", self.output], 
                                prog_name='build')
            except SystemExit as e:
                exc_info = sys.exc_info()
                exit_code = e.code
                if exit_code is None:
                    exit_code = 0
                
                if exit_code != 0:
                    raise e
                
            # os.system(f"cp {self.output} ../")
            shutil.copy(self.output, "../")
            if self.corrected:
                if Path(self.output.replace("agp", "corrected.agp")).exists():
                    shutil.copy(f"{self.output.replace('agp', 'corrected.agp')}", "../")
    
        os.chdir("../")

        if self.delete_temp:
            shutil.rmtree(tmpDir)
            logger.info("Removed temporary directory.")

        logger.info("Scaffolding done.")

##Deprecated
class OldOptimize0:
    
    def __init__(self, contigs, clm, threads=10):
        self.contigs = natsorted(contigs)
        self.contig_idx = dict(zip(self.contigs, range(len(self.contigs))))

        self.clm = clm

        self.data = self.parse()

    def parse(self):
        df = self.clm.dk_df 

        res = dict(zip(df.to_dict('split')['index'], 
                        df.to_dict('split')['data']))
        
        return res
    
    def filter(self, as_idx=False):
        """
        only retain the max score in different orientation pairs

        """
        score_res = {}
        orientation_res = {}
        for pair in self.data:
            scores = self.data[pair]
            max_index = np.argmax(scores)
            max_score = scores[max_index]
            if max_score == 0:
                continue

            score_res[pair] = max_score
            orientation_res[pair] = max_index
        
        graph_df = pd.DataFrame(score_res, index=['weight']).T
        graph_df = graph_df.reset_index()
        graph_df.columns = ['source', 'target', 'weight']
        
        if as_idx:
            graph_df['source'] = graph_df['source'].map(self.contig_idx.get)
            graph_df['target'] = graph_df['target'].map(self.contig_idx.get)
            graph_df = graph_df.dropna()
            graph_df['source'] = graph_df['source'].astype('int')
            graph_df['target'] = graph_df['target'].astype('int')

        return graph_df, orientation_res
    
    def save(self, output):
        """
        save the score data.

        Params:
        --------
        output: str 

        Examples:
        --------
        >>> so2.save("out.score.txt")

        """
        with open(output, "w") as out:
            for pair in self.data:
                for i, s in enumerate(self.data[pair]):
                    _pair = list(zip(pair, self.orientations[i]))
                    _pair = list(map(lambda x: "".join(x), _pair))
                    print(f"{_pair[0]} {_pair[1]}\t{s:f}", file=out)
        logger.info(f"Successful output the score of contig pairs into {output}")


    def graph(self):
        """
        construct a graph
        """
        graph_df, orientation_res = self.filter(as_idx=True)
        # graph_df['weight'] = 1 / graph_df['weight']
        G = ig.Graph.DataFrame(graph_df, directed=False)
        
        return G

    def _dfs(self, start):
        
        pass 

    def dfs(self):
        pass

    def minimum_spanning_tree(self):
        pass

class OldOptimize:
    orientations = ["++", "+-", "-+", "--"]
    
    def __init__(self, contigs, clm, threads=10):
        self.contigs = natsorted(contigs)
        self.contig_idx = dict(zip(self.contigs, range(len(self.contigs))))

        self.clm = clm

        self.data = self.parse()

    def parse(self):
        df = self.clm.dk_df 

        res = dict(zip(df.to_dict('split')['index'], 
                        df.to_dict('split')['data']))
        
        return res
    
    def filter(self, as_idx=False):
        """
        only retain the max score in different orientation pairs

        """
        score_res = {}
        orientation_res = {}
        for pair in self.data:
            scores = self.data[pair]
            max_index = np.argmax(scores)
            max_score = scores[max_index]
            if max_score == 0:
                continue

            score_res[pair] = max_score
            orientation_res[pair] = max_index
        
        graph_df = pd.DataFrame(score_res, index=['weight']).T
        graph_df = graph_df.reset_index()
        graph_df.columns = ['source', 'target', 'weight']
        
        if as_idx:
            graph_df['source'] = graph_df['source'].map(self.contig_idx.get)
            graph_df['target'] = graph_df['target'].map(self.contig_idx.get)
            graph_df = graph_df.dropna()
            graph_df['source'] = graph_df['source'].astype('int')
            graph_df['target'] = graph_df['target'].astype('int')

        return graph_df, orientation_res
    
    def save(self, output):
        """
        save the score data.

        Params:
        --------
        output: str 

        Examples:
        --------
        >>> so2.save("out.score.txt")

        """
        with open(output, "w") as out:
            for pair in self.data:
                for i, s in enumerate(self.data[pair]):
                    _pair = list(zip(pair, self.orientations[i]))
                    _pair = list(map(lambda x: "".join(x), _pair))
                    print(f"{_pair[0]} {_pair[1]}\t{s:f}", file=out)
        logger.info(f"Successful output the score of contig pairs into {output}")


    def graph(self):
        """
        construct a graph
        """
        graph_df, orientation_res = self.filter(as_idx=True)
        # graph_df['weight'] = 1 / graph_df['weight']
        G = ig.Graph.DataFrame(graph_df, directed=False)
        
        return G

    def _dfs(self, start):
        
        pass 

    def dfs(self):
        pass

    def minimum_spanning_tree(self):
        mst_tree = self.G.spanning_tree()


class SimpleOptimize:
    """

    """
    orientations = ["++", "+-", "-+", "--"]
    
    def __init__(self, contigs, cool, threads=10):
        self.cool = cool
        self.contigs = sorted(contigs, key=lambda x: self.cool.chromnames.index(x))
        self.contig_idx = dict(zip(self.contigs, range(len(self.contigs))))

        self.matrix = cool.matrix(balance=False, sparse=True)
        
        self.threads = threads 
        
        self.data = self.parse()

    @staticmethod
    def _parse(matrix, pair):
        """
        single pair parser
        """
        res = []
        contig1, contig2 = pair 
        sub_matrix = matrix.fetch(contig1, contig2)
        if sub_matrix.getnnz() == 0:
            return

        sub_matrix = sub_matrix.tocsr()
        l1, l2 = sub_matrix.shape
        d1 = float(int(np.ceil(l1 / 2)))
        d2 = float(int(np.ceil(l2 / 2)))

        k = l2 - l1 if l1 < l2 else 0 

        ## contig+ contig+
        c = tril(sub_matrix[d1:, :d2], k).sum()
        s = c / (d1 * d2)
        res.append(s)
        # print(f"{contig1}+ {contig2}+\t{s}")
        
        ## contig+ contig-
        c = tril(sub_matrix[:, ::-1][d1:, :d2], k).sum()
        s = c / (d1 * d2)
        res.append(s)
        # print(f"{contig1}+ {contig2}-\t{s}")

        ## contig- contig+
        c = tril(sub_matrix[::-1, :][d1:, :d2], k).sum()
        s = c / (d1 * d2)
        res.append(s)
        # print(f"{contig1}- {contig2}+\t{s}")

        ##contig- contig-
        c = tril(sub_matrix[::-1, ::-1][d1:, :d2], k).sum()
        s = c / (d1 * d2)
        res.append(s)
        # print(f"{contig1}- {contig2}-\t{s}")

        return pair, res 

    def parse(self):
        """
        calculate the score between different contig pairs

        |------||--------|        
        """
        
        res = defaultdict(list)

        args = []
        for pair in combinations(self.contigs, 2):
            
            args.append((self.matrix, pair))
        
        res = Parallel(n_jobs=self.threads)(
                delayed(SimpleOptimize._parse)(i, j) for i, j in args )
        
        ## remove None value
        res = dict(filter(lambda x: x is not None, res))
           
        return res

    def filter(self, as_idx=False):
        """
        only retain the max score in different orientation pairs

        """
        score_res = {}
        orientation_res = {}
        for pair in self.data:
            scores = self.data[pair]
            max_index = np.argmax(scores)
            max_score = scores[max_index]
            if max_score == 0:
                continue

            score_res[pair] = max_score
            orientation_res[pair] = max_index
        
        graph_df = pd.DataFrame(score_res, index=['weight']).T
        graph_df = graph_df.reset_index()
        graph_df.columns = ['source', 'target', 'weight']
        
        if as_idx:
            graph_df['source'] = graph_df['source'].map(self.contig_idx.get)
            graph_df['target'] = graph_df['target'].map(self.contig_idx.get)
            

        return graph_df, orientation_res
    
    def save(self, output):
        """
        save the score data.

        Params:
        --------
        output: str 

        Examples:
        --------
        >>> so2.save("out.score.txt")

        """
        with open(output, "w") as out:
            for pair in self.data:
                for i, s in enumerate(self.data[pair]):
                    _pair = list(zip(pair, self.orientations[i]))
                    _pair = list(map(lambda x: "".join(x), _pair))
                    print(f"{_pair[0]} {_pair[1]}\t{s}", file=out)
        logger.info(f"Successful output the score of contig pairs into {output}")


    def graph(self):
        """
        construct a graph
        """
        graph_df, orientation_res = self.filter(as_idx=True)
        # graph_df['weight'] = 1 / graph_df['weight']
        G = ig.Graph.DataFrame(graph_df, directed=False)
        
        return G

    def _dfs(self, start):
        
        pass 

    def dfs(self):
        pass

    def minimum_spanning_tree(self):
        pass

##Deprecated
class SimpleOptimize2:
    """

    """
    orientations = ["++", "+-", "-+", "--"]
    
    def __init__(self, contigs, cool, method="so", threads=10):
        self.cool = cool
        self.contigs = contigs #sorted(contigs, key=lambda x: self.cool.chromnames.index(x))
        self.contig_idx = dict(zip(self.contigs, range(len(self.contigs))))
        self.idx_to_contig = dict(zip(range(len(self.contigs)), self.contigs))

        self.matrix = cool.matrix(balance=False, sparse=True)
        
        self.method = method
        self.threads = threads 
        
        self.data = self.parse()
        self.score_df, self.orientation_res = self.filter(mode='score')


    @staticmethod
    def _parse_by_so(matrix, pair):
        """
        single pair parser
        """
        res = []
        contig1, contig2 = pair 
        sub_matrix = matrix.fetch(contig1, contig2)
        if sub_matrix.getnnz() == 0:
            return
        sub_matrix = sub_matrix.tocsr()
        l1, l2 = sub_matrix.shape
        d1 = int(np.ceil(l1 / 2))
        d2 = int(np.ceil(l2 / 2))

        k = l2 - l1 if l1 < l2 else 0 

        ## contig+ contig+
        c = tril(sub_matrix[d1:, :d2], k).sum()
        s = c / (d1 * d2)
        res.append(s)
        # print(f"{contig1}+ {contig2}+\t{s}")
        
        ## contig+ contig-
        c = tril(sub_matrix[:, ::-1][d1:, :d2], k).sum()
        s = c / (d1 * d2)
        res.append(s)
        # print(f"{contig1}+ {contig2}-\t{s}")

        ## contig- contig+
        c = tril(sub_matrix[::-1, :][d1:, :d2], k).sum()
        s = c / (d1 * d2)
        res.append(s)
        # print(f"{contig1}- {contig2}+\t{s}")

        ##contig- contig-
        c = tril(sub_matrix[::-1, ::-1][d1:, :d2], k).sum()
        s = c / (d1 * d2)
        res.append(s)
        # print(f"{contig1}- {contig2}-\t{s}")

        return pair, res 

    @staticmethod
    def _parse_by_so2(matrix, pair):
        """
        single pair parser
        """
        res = []
        contig1, contig2 = pair 
        sub_matrix = matrix.fetch(contig1, contig2)
        if sub_matrix.getnnz() == 0:
            return
        sub_matrix = sub_matrix.tocsr()
        l1, l2 = sub_matrix.shape
        d1 = int(np.ceil(l1 / 2))
        d2 = int(np.ceil(l2 / 2))

        d_matrix = np.zeros((l1, l2))
        for i in range(l1):
            for j in range(l2):
                d_matrix[i, j] = i + j + 1
        d_matrix = d_matrix[::-1]
        
    
        ## contig+ contig+
        s = (sub_matrix[d1:, :d2] / d_matrix[d1:, :d2]).sum()
        
        res.append(s)
        # print(f"{contig1}+ {contig2}+\t{s}")
        
        ## contig+ contig-
        s = (sub_matrix[:, ::-1][d1:, :d2] / d_matrix[d1:, :d2]).sum()
        res.append(s)
        # print(f"{contig1}+ {contig2}-\t{s}")

        ## contig- contig+
        s = (sub_matrix[::-1, :][d1:, :d2] / d_matrix[d1:, :d2]).sum()
        res.append(s)
        # print(f"{contig1}- {contig2}+\t{s}")

        ##contig- contig-
        s = (sub_matrix[::-1, ::-1][d1:, :d2] / d_matrix[d1:, :d2]).sum()
        res.append(s)
        # print(f"{contig1}- {contig2}-\t{s}")

        return pair, res 

    def parse(self):
        """
        calculate the score between different contig pairs

        |------||--------|        
        """
        
        res = defaultdict(list)

        args = []
        for pair in combinations(self.contigs, 2):
            
            args.append((self.matrix, pair))
        
        _parse = SimpleOptimize2._parse_by_so if self.method == "so" else SimpleOptimize2._parse_by_so2

        res = Parallel(n_jobs=self.threads)(
                delayed(_parse)(i, j) for i, j in args )
        
        ## remove None value
        res = dict(filter(lambda x: x is not None, res))
           
        return res

    def filter(self, mode='graph', as_idx=False):
        """
        only retain the max score in different orientation pairs

        """
        score_res = {}
        orientation_res = {}
        for pair in self.data:
            scores = self.data[pair]
            max_index = np.argmax(scores)
            max_score = scores[max_index]
            if max_score == 0:
                continue

            score_res[pair] = max_score
            orientation_res[pair] = max_index
        
        if mode == 'graph':
            graph_df = pd.DataFrame(score_res, index=['weight']).T
            graph_df = graph_df.reset_index()
            graph_df.columns = ['source', 'target', 'weight']
            
            if as_idx:
                graph_df['source'] = graph_df['source'].map(self.contig_idx.get)
                graph_df['target'] = graph_df['target'].map(self.contig_idx.get)
       
            return graph_df, orientation_res
        
        elif mode == 'score':
            score_df = pd.DataFrame(score_res, index=['score']).T
            score_df = score_df.reset_index()
            score_df.columns = ['contig1', 'contig2', 'score']
            score_df2 = score_df.rename(columns={"contig1": "contig2",
                                                    "contig2": "contig1"})
            score_df = pd.concat([score_df, score_df2], axis=0)

            score_df.set_index(["contig1", "contig2"], inplace=True)

            return score_df, orientation_res
        else:
            raise ValueError("mode must in {'graph', 'score'}")
        
        
    
    def save_score(self, output):
        """
        save the score data.

        Params:
        --------
        output: str 

        Examples:
        --------
        >>> so2.save("out.score.txt")

        """
        with open(output, "w") as out:
            for pair in self.data:
                for i, s in enumerate(self.data[pair]):
                    _pair = list(zip(pair, self.orientations[i]))
                    _pair = list(map(lambda x: "".join(x), _pair))
                    print(f"{_pair[0]} {_pair[1]}\t{s}", file=out)
        logger.info(f"Successful output the score of contig pairs into {output}")

    def graph(self):
        """
        construct a graph
        """

        graph_df, self.orientation_info = self.filter(as_idx=True)
    
        graph_df2 = graph_df[['target', 'source', 'weight']]
        graph_df2.columns = ['source', 'target', 'weight']

        graph_df = pd.concat([graph_df, graph_df2], axis=0)
        graph_df['weight'] = 1 / graph_df['weight']
        # graph_df['weight'] = max(graph_df['weight']) - graph_df['weight']
        # graph_df['weight'] = graph_df['weight'].astype(int)
        G = ig.Graph.DataFrame(graph_df, directed=False)
        

        return G, graph_df

    def get_global_score(self, order_list):
        order_pairs = [(i, j) for i, j in zip(order_list[:-1], order_list[1:])]
    
        scores = self.graph_df.reindex(order_pairs, fill_value=0)
        
        score = sum(scores['weight'])

        return score
    
    def bfs(self, start):
        res = self.G.bfs(start)[0]
        
        if len(res) == 1:
            return None
        # res = list(map(self.idx_to_contig.get, res[0]))
        return  self.get_global_score(res), res

    def dfs(self, start):
        
        res = self.G.dfs(start)[0]
        
        if len(res) == 1:
            return None
        # res = list(map(self.idx_to_contig.get, res[0]))
        return  self.get_global_score(res), res

    def order(self, algorithm='dfs'):
        
        args = []
        for idx in range(len(self.contigs)):
            args.append(idx)
        
        res = Parallel(n_jobs=self.threads)(
                delayed(getattr(self, algorithm))(i) for i in args)

        res = list(filter(lambda x: x is not None, res))
        res = sorted(res, key=lambda x: x[0])

        print(list(map(self.idx_to_contig.get, res[0][1])))
        return list(map(self.idx_to_contig.get, res[0][1]))

    def tsp_order(self):
        import networkx as nx
        tsp = nx.approximation.traveling_salesman_problem

        G = self.G.to_networkx()
        # SA_tsp = nx.approximation.simulated_annealing_tsp
        # method = lambda G, wt: SA_tsp(G, "greedy", weight=wt, move='1-0', temp=500)
        # path = tsp(G, cycle=False, method=method)
        
        path = tsp(G, cycle=False)

        return list(map(self.idx_to_contig.get, path))
    
    def lkh_order(self):
        import elkai
        import networkx as nx 
        G = self.G.to_networkx()
        M = nx.adjacency_matrix(G).todense()
        print(M)
        path = elkai.solve_int_matrix(M)

        return list(map(self.idx_to_contig.get, path))
    

    def nn_tsp(self, contigs, score_df, start=0):

        start_contig = "Chr5.ctg1" #contigs[start]
        candicate_contigs = set(set(contigs) - {start_contig})
        tour = [start_contig]
        print(start_contig)
        added_contigs = set() 
        added_contigs.add(start_contig)

        while candicate_contigs:
            next_contig = score_df.loc[tour[-1]].idxmax().values[0]
            tmp_df = score_df.loc[tour[-1]].sort_values(by=['score'], ascending=False)
            if next_contig in added_contigs:
                # print(score_df.loc[tour[-1]]['score'].nlargest(2, keep='all'))
                next_contig = score_df.loc[tour[-1]]['score'].nlargest(2, keep='all').index[1]
                i = 1
                while next_contig in added_contigs:
                    i += 1
                    next_contig = tmp_df.index[i]

            if next_contig not in contigs:
                continue
            
            print(next_contig)
            added_contigs.add(next_contig)
            tour.append(next_contig)

            candicate_contigs.remove(next_contig)

        return tour

    def orientation(self):
        
        for i, j in zip(self.ordering[:-1], self.ordering[1:]):
            try:
                print(i, j, self.orientations[self.orientation_info[(i, j)]])
            except:
                print(i, j)

    def save(self, output):
        with open(output, "w") as out:
            print("\n".join(self.ordering), file=out)

    def minimum_spanning_tree(self):
        pass


##Deprecated
class SAOptimizer:
    """
    Simulated Annealing
    """
    def __init__(self) -> None:
        pass 

    def cooling_schedule(self, t):
        return 0.99 * t
    
    def acceptance_probability(self, energy, new_energy, temperature):
        if new_energy < energy:
            return 1.0
        return np.exp((energy - new_energy) / temperature)
    
    def simulated_annealing(self, graph, start, iterations=1000):
        """
        Simulated Annealing

        Params:
        -------
        graph: np.array
            adjacency matrix
        start: list 
            start node
        iterations: int
            number of iterations
        
        Returns:
        --------
        current_path: list
        current_cost: int

        Example:
        --------
        >>> graph = np.array([[0, 1, 2, 3],
                                [1, 0, 4, 5],
                                [2, 4, 0, 6],  
                                [3, 5, 6, 0]])
        >>> start = 0
        >>> iterations = 1000
        >>> sa = SA()
        >>> current_path, current_cost = sa.simulated_annealing(graph, start, iterations)
        >>> print(current_path, current_cost)
        [0, 1, 2, 3] 12
        """
        # initial state
        current_path = [start]
        current_cost = 0
        unvisited_nodes = set(range(len(graph))) - {start}
    
        # iteration
        for i in range(iterations):
            # cooling
            T = self.cooling_schedule(iterations / i)
    
            # random neighbour
            current_node = current_path[-1]
            next_node = random.sample(unvisited_nodes, 1)[0]
    
            # calculate current and neighbour path cost
            current_cost += graph[current_node][next_node]
            new_cost = current_cost - graph[current_node][next_node] + graph[next_node][current_node]
    
            # accept neighbour or not
            if self.acceptance_probability(current_cost, new_cost, T) > random.random():
                current_cost = new_cost
                current_path.append(next_node)
                unvisited_nodes.remove(next_node)
    
        return current_path, current_cost


class HyperOptimize:
    def __init__(self, HG, split_num=2, mutapb=0.2, npop=100, ngen=5000):
        self.H = HG.incidence_matrix().T
        # print(self.H.shape)
        # print(self.H.sum(axis=0)[self.H.sum(axis=0) == 3].T.shape)
        # l = list(set(list(map(tuple, self.H.T.tolil().rows))))
        # rows, cols, vals = [], [], []
        # for row_idx, row in enumerate(l):
        #     rows.extend([row_idx] * len(row))
        #     cols.extend(row)
        #     vals.extend([1] * len(row))
        # self.H = coo_matrix((vals, (rows, cols))).tocsr()
        

        self.split_num = split_num 
        self.nodes = HG.nodes 
        self.contigs = pd.DataFrame(self.nodes)[0].str.rsplit("_").map(
                                        lambda x: x[0]).astype("category")
        self.contigsizes = HG.edges.contigsizes

        self.contig_idx = self.contigs.cat.codes

        self.mutapb = mutapb
        self.npop = npop 
        self.ngen = ngen
        self.order = self.init_order()

    def init_order(self):
        order = np.arange(len(self.nodes))

        return order

    def shuffle(self):
        split_num = self.split_num
        length = len(self.order)
        reshaped_order = np.reshape(self.order, (length//split_num, split_num))
    
        np.random.shuffle(reshaped_order)

        return reshaped_order.reshape(length)

    def init_population(self):
        population = [self.shuffle() for i in range(self.npop)]
        
        return population
    
    def fitness(self, order):
        a = self.H[:, order].tolil().rows
        a = list(map(lambda x: tuple(x), a))
        a = set(a)
        a = list(map(lambda x: np.array(x, dtype=np.float32), a))
        value = sum(list(map(lambda x: 1/(
                            np.prod(x[1:] - x[:-1])), a)))

        return value 
    
    @staticmethod
    def splice(order, split_num):
        length = len(order)
        pos = random.randint(0, length) 
        pos = pos - pos % split_num
        
        a, b = list(range(pos)), list(range(pos, length))
        
        return order[b + a]

    @staticmethod
    def permute(order, split_num):
        length = len(order)
        reshaped_order = np.reshape(order, (length//split_num, split_num))
        p, q = random.randint(0, length//split_num - 1), random.randint(0, length//split_num - 1)
    
        if p == q:
            return order 

        reshaped_order[[p, q]] = reshaped_order[[q, p]]

        return reshaped_order.reshape(length)

    @staticmethod
    def insertion(order, split_num):
        length = len(order)

        reshaped_order = np.reshape(order, (length//split_num, split_num))
        p, q = random.randint(0, length//split_num - 1), random.randint(0, length//split_num - 1 )
    
        if p == q:
            return order 

        if random.random() < 0.5:
            temp = reshaped_order[q].copy()
            reshaped_order[q] = reshaped_order[p]
            
            reshaped_order = np.vstack((np.delete(reshaped_order, p, axis=0), temp))
        else:
            temp = reshaped_order[p].copy()
            reshaped_order[p] = reshaped_order[q]

            reshaped_order = np.vstack((np.delete(reshaped_order, q, axis=0), temp))

        return reshaped_order.reshape(length)

    @staticmethod
    def inversion(order, split_num):
        length = len(order)
        
        p, q = random.randint(0, length), random.randint(0, length)

        reshaped_order = np.reshape(order, (length//split_num, split_num))
        p, q = random.randint(0, length//split_num - 1), random.randint(0, length//5 - 1)
    
        if p == q:
            return order 
        
        if p > q:
            p, q = q, p

        reshaped_order = np.vstack((reshaped_order[:p], 
                                    reshaped_order[p: q][::-1], 
                                    reshaped_order[q:]))

        return reshaped_order.reshape(length)

    def mutate(self, order):
        split_num = self.split_num
        random_num = random.random()

        if random_num < 0.2:
            new_order = self.permute(order, split_num)
        elif random_num < 0.4:
            new_order = self.splice(order, split_num)
        elif random_num < 0.7:
            new_order = self.insertion(order, split_num)
        else:
            new_order = self.inversion(order, split_num)
        
        return new_order 

    def crossover(self, order1, order2):
        split_num = self.split_num
        length = len(order1)
        
        order1 = np.reshape(order1, (length // split_num, split_num))
        order2 = np.reshape(order2, (length // split_num, split_num))
        
        offspring = np.ones(shape=(length // split_num, split_num))
        offspring = -offspring 
        random_num_1 = random.randint(0, length - 1)
        random_num_2 = random.randint(random_num_1, length)

        offspring[random_num_1: random_num_2] = order2[random_num_1: random_num_2]

        return offspring.reshape(length).astype(int)
                

    def select():
        pass 

    def evaluate(self):
        population = self.init_population()

        total_population = 0
        generation = 0
        stable_round = 0

        while True:
            generation += 1
            total_population += len(population)

            if generation != 1:
                scores = [self.fitness(order) for order in population]
                index = list(range(len(scores)))
                scores = dict(zip(index, scores))
                best = max(scores.items(), key=lambda x: x[1])
                print(best)
                selected = []
                while len(selected) < self.npop:
                    tmp = max([(i, scores[i]) for i in random.sample(index, 3)], 
                                key=lambda x: x[1])
                    selected.append(tmp[0])

                selected = [population[i] for i in selected]
                
                population = selected

            new_population = []
            for i in range(0, len(population) - 1, 2):
                try:
                    order1 = population[i]
                    order2 = population[i+1]
                    if random.random() < self.mutapb:
                        new_population.append(self.mutate(order1))
                        new_population.append(self.mutate(order2))
                    else:
                        new_population.append(order1)
                        new_population.append(order2)
                        # new_population.append(self.crossover(order1, order2))
                        # new_population.append(self.crossover(order1, order2))
                except:
                    new_population.append(population[-1])
           
            population = new_population
            
            if generation > self.ngen:
                break
            
            if generation > 2:
                if best[1] - previous_best[1] < 0.01:
                    stable_round += 1
                else:
                    stable_round = 0
            if generation > 2:        
                print(best)
                print(population[best[0]])
            if stable_round >= 10:
                print(stable_round)
                print(best)
                print(population[best[0]])
                break

            if generation > 1:
                previous_best = best 
              

def test(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('count_re', 
            help='contig group in countRE table')
    pReq.add_argument('cool', 
            help='Path to cool file')
    
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    contigs = CountRE(args.count_re).contigs
    cool = cooler.Cooler(args.cool)
    so =  SimpleOptimize(contigs, cool)
    so.parse()


class HyperScaffolding:
    def __init__(self, HG, contigs, resolution=1.0, min_contacts=1):
        self.min_contacts = min_contacts
        self.HG = HG 
        HG.extract_rows(contigs)
        self.resolution = resolution

        self.H = self.HG.incidence_matrix(min_contacts=self.min_contacts)
        self.length_array = np.array([self.HG.contigsizes[x] for x in self.HG.nodes])
        self.NW = self.get_normalize_weight()


    def get_normalize_weight(self):
        
        NW = np.outer(self.length_array, self.length_array)
        
        return NW 
    
    def cluster(self):

        self.A, cluster_assignments, self.K = IRMM(self.H, resolution=self.resolution)
        self.NA = np.array(self.A / self.NW)
        print(len(self.K))
        logger.debug(f"Partition contigs into {len(self.K)} group.")
 

    def orientation(self):
        pathes = []
        for sub_group in self.K:
            if len(sub_group) <= 2:
                pathes.append(list(sub_group))
                continue
            sub_group_idx = dict(zip(range(len(sub_group)), list(sub_group)))
            sub_NA = self.NA[list(sub_group), :][:, list(sub_group)]
            # sub_NA[np.where(sub_NA < 200)] = 0.0
            
            G = ig.Graph.Weighted_Adjacency(1 / sub_NA, mode='undirected', loops=False).to_networkx()
            # tree = nx.minimum_spanning_tree(G)
            ring_path = nx_app.christofides(G, weight='weight')
           
            print(self.HG.nodes[ring_path])
            ## break ring path
            
            tmp_res = {}
            for i in range(0, len(ring_path) - 1):
                j = i + 1
                i_idx, j_idx = ring_path[i], ring_path[j]
                tmp_res[(i, j)] = sub_NA[i_idx, j_idx]

            print(tmp_res)
            break_pos = sorted(tmp_res, key=lambda x: tmp_res[x])
            print(break_pos)
            break_pos = break_pos[0]
            logger.debug(f"Break ring path at {break_pos}.")

            path = []
            path.extend(ring_path[break_pos[1]:])
            path.extend(ring_path[:break_pos[1]][1:])
            new_path = [sub_group_idx[x] for x in path]

            pathes.append(new_path)

        if len(pathes) <= 1:
            return pathes[0]
        
        split_length_array = np.zeros(len(pathes) * 2)
        split_pathes = []
        print(pathes)
        for i in range(len(pathes)):
            split_pos = 0 
            length = self.length_array[pathes[i]].sum()
            cum_length = 0
            for idx in pathes[i]:
                cum_length += self.length_array[idx]
                split_pos += 1 
                if cum_length > length // 2:
                    if (cum_length - length // 2) > (self.length_array[idx] - (cum_length - length // 2 )):
                        split_pos -= 1 
                        cum_length -= self.length_array[idx]
                    break 
                   
            split_length_array[i*2] = cum_length 
            split_length_array[i*2 + 1] = length - cum_length 
            split_pathes.append((pathes[i][: split_pos], pathes[i][split_pos:]))


        split_group_A = np.zeros(shape=(len(pathes) * 2, len(pathes) * 2))
        for i in range(0, len(pathes) - 1):
            for j in range(i, len(pathes)):
                if i == j:
                    continue

                i1, j1, i2, j2 = i * 2, j * 2, i * 2 + 1, j * 2 + 1
                pathes_2, pathes_1 = split_pathes[i], split_pathes[j]
                # print(list(map(lambda x: self.HG.nodes[x], pathes_1)), list(map(lambda x: self.HG.nodes[x], pathes_2)))
                v1 = self.A[pathes_1[0], :][:, pathes_2[0]].sum() 
                v3 = self.A[pathes_1[0], :][:, pathes_2[1]].sum()
                v4 = self.A[pathes_1[1], :][:, pathes_2[0]].sum()
                v2 = self.A[pathes_1[1], :][:, pathes_2[1]].sum()
                a = np.array([[v1, v3], [v4, v2]])
                
                split_group_A[i1, j1] = v1  / (split_length_array[i1] * split_length_array[j1])
                split_group_A[i2, j2] = v2 / (split_length_array[i2] * split_length_array[j2])
                split_group_A[i2, j1] = v3  / (split_length_array[i2] * split_length_array[j1])
                split_group_A[i1, j2] = v4  / (split_length_array[i1] * split_length_array[j2])

        split_group_A += split_group_A.T - np.diag(split_group_A.diagonal())

        G = ig.Graph.Weighted_Adjacency(1 / split_group_A, mode='undirected', loops=False).to_networkx()
        # tree = nx.minimum_spanning_tree(G)
        ring_path = nx_app.christofides(G, weight='weight')

    
        group_A = np.zeros(shape=(len(pathes), len(pathes)))
        for i in range(len(pathes) - 1):
            for j in range(i + 1, len(pathes)):
                if i == j:
                    continue
                idx1 = pathes[i]
                idx2 = pathes[j]
                length1 = self.length_array[idx1].sum()
                length2 = self.length_array[idx2].sum()
                group_A[i, j] = self.NA[idx1, :][:, idx2].sum() / ( length1 * length2)

        group_A += group_A.T - np.diag(group_A.diagonal())
        print(group_A)
        G = ig.Graph.Weighted_Adjacency(1 / group_A, mode='undirected', loops=False).to_networkx()
        ring_path = nx_app.christofides(G, weight='weight')
        tmp_res = {}
        for i in range(0, len(ring_path) - 1):
            j = i + 1
            i_idx, j_idx = ring_path[i], ring_path[j]
            tmp_res[(i, j)] = group_A[i_idx, j_idx]

        break_pos = min(tmp_res, key=lambda x: tmp_res[x])
        logger.debug(f"Break ring path at {break_pos}.")

        path = []
        path.extend(ring_path[break_pos[1]:])
        path.extend(ring_path[:break_pos[1]][1:])

        
        new_path = []
        for idx1 in range(len(path) - 1):
           
            idx2 = idx1 + 1
            a = np.array([[split_group_A[idx1*2, idx2*2], 
                           split_group_A[idx1*2 +  1, idx2*2]], 
                        [split_group_A[idx1*2, idx2*2 + 1], 
                        split_group_A[idx1*2 + 1, idx2*2 + 1]]])
            
            max_idx = a.argmax()
            
            if max_idx == 0:
                if idx1 == 0:
                    new_path.extend(pathes[idx1][::-1])
                
                new_path.extend(pathes[idx2])
            elif max_idx == 1:
                if idx1 == 0:
                    new_path.extend(pathes[idx1])
                    
                new_path.extend(pathes[idx2])
            elif max_idx == 3:
                if idx1 == 0:
                    new_path.extend(pathes[idx1])
                new_path.extend(pathes[idx2][::-1])
            else:
                if idx1 == 0:
                    new_path.extend(pathes[idx1])
                    
                new_path.extend(pathes[idx2])

    
        return new_path


    def orient(self):
        pass 


def raw_sort(K, A, length_array, threads=4):

    def func(sub_group, A, length_array):
        sub_group = list(sub_group)
        if len(sub_group) <= 2:
            return sub_group
        NW = np.outer(length_array[sub_group], length_array[sub_group])
        
        sub_group_idx = dict(zip(range(len(sub_group)), list(sub_group)))
        
        sub_A = A[sub_group, :][:, sub_group]
        sub_NA = sub_A / NW
        try:
            sub_NA2 = 1 / sub_NA
        except TypeError:
            sub_NA = sub_NA.toarray()
            sub_NA2 = 1 / sub_NA

        try:
            G = ig.Graph.Weighted_Adjacency(sub_NA2, mode='undirected', loops=False).to_networkx()

            ring_path = nx_app.christofides(G, weight='weight')
        except:
            return sub_group
        tmp_res = {}
        for i in range(0, len(ring_path) - 1):
            j = i + 1
            i_idx, j_idx = ring_path[i], ring_path[j]
            tmp_res[(i, j)] = sub_NA[i_idx, j_idx]
        
        break_pos = sorted(tmp_res, key=lambda x: tmp_res[x])
        break_pos = break_pos[0]
        path = []
        path.extend(ring_path[break_pos[1]:])
        path.extend(ring_path[:break_pos[1]][1:])
        new_path = [sub_group_idx[x] for x in path]

        return new_path
    
    args = [(k, A, length_array) for k in K]
    # res = []
    # for i, j, k in args:
    #     res.append(func(i, j, k))
    if len(args) < 1:
        return K
    res = Parallel(n_jobs=min(threads, len(args)))(delayed(
        func)(i, j, k) for i, j, k in args 
    )
    

    return res 


def test_hyperscaffold(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('hg')
    pReq.add_argument('contig_list', 
            help='contig list')
    
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    import msgspec 
    from cphasing.algorithms.hypergraph import HyperEdges, HyperGraph

    contigs = [i.strip() for i in open(args.contig_list)]
    
    he = msgspec.msgpack.decode(open(args.hg, 'rb').read(), type=HyperEdges)
    he.to_numpy()
    hg = HyperGraph(he)

    hs = HyperScaffolding(hg, contigs)
    hs.cluster()
    res = hs.orientation()
    print(hs.HG.nodes[res])

if __name__ == "__main__":
    test_hyperscaffold(sys.argv[1:])