#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
correct chimeric contig by Hi-C signal.
"""

import logging
import os
import os.path as op
import sys
import shutil

import cooler
import numpy as np
import pandas as pd 


from collections import OrderedDict
from cooler.cli.cload import pairs as pairs_cmd
from joblib import Parallel, delayed

from pathlib import Path
from pyfaidx import Fasta
from pyranges import PyRanges
from scipy.signal import find_peaks
from scipy.sparse import tril

from .core import Pairs
from .utilities import read_fasta, xopen

logger = logging.getLogger(__name__)

def _zero_diags(chunk, n_diags):
    if n_diags > 0:
        if n_diags > 0:
            mask = np.abs(chunk['pixels']['bin1_id'] 
                            - chunk['pixels']['bin2_id']) < n_diags
            chunk['pixels']['count'][mask] = 0

    return chunk

class Corrector:
    """
    Params:
    --------
    fasta: str
        Draft assembly
    pairs: str
        4DN pairs
    
    Returns:
    --------

    Examples:
    --------
    >>> 
    """
    def __init__(self, fasta, pairs, 
                    percent=.95,
                    sensitive=.5,
                    resolutions='20000,10000,5000,1000',
                    depletion=4,
                    temp_dir='tmp', 
                    delete_temp=True, 
                    threads=1,
                    force=False,
                    output=sys.stdout,
                    outbed=None,
                    outpairs=None):
        self.fasta = Path(fasta)
        self.pairs = Path(pairs)
        self.prefix = Path(self.pairs.stem)
        while self.prefix.suffix in {'.gz', '.pairs'}:
            self.prefix = self.prefix.with_suffix('')
        
        try:
            self.resolutions = list(map(int, resolutions.split(',')))
        except ValueError:
            logger.Error('Resolutions must be in numberic.')
            sys.exit()

        self.threads = threads
        self.force = force
        self.percent = percent
        self.sensitive = sensitive
        
        self.depletion = depletion
        self.chrom_sizes_path = self.get_chrom_sizes()
        self.output = output
        self.outbed = outbed
        self.outpairs = outpairs
        

    def get_chrom_sizes(self):
        fasta = Fasta(self.fasta)
        chrom_sizes_path = f'{fasta.filename}.chromsizes'
        chrom_sizes_df = pd.read_csv(fasta.faidx.indexname, sep='\t',
                                        header=None, usecols=[0, 1])
        chrom_sizes_df.to_csv(chrom_sizes_path, sep='\t', 
                                index=False, header=False)
        chrom_sizes_df = chrom_sizes_df.set_index(0, drop=True)

        self.chrom_sizes = chrom_sizes_df

        return chrom_sizes_path    

    def get_contacts(self, resolution):
        res_cool_path = f'{self.prefix}.{resolution}.cool'
        
        if Path(res_cool_path).exists() and not self.force:
            logger.warning(f'Using existing contact file of `{res_cool_path}`')
            return resolution, res_cool_path
        try:
            pairs_cmd.main(args=[f'{self.chrom_sizes_path}:{resolution}', 
                                    str(self.pairs), 
                                    res_cool_path, 
                                    '-c1', '2', '-p1', '3', 
                                    '-c2', '4', '-p2', '5'], 
                                    prog_name='cooler')
        except SystemExit as e:
            exc_info = sys.exc_info()
            exit_code = e.code
            if exit_code is None:
                exit_code = 0
            
            if exit_code != 0:
                raise e
        
        logger.info(f'Output cool of resolution `{resolution}`.')
 
        return resolution, res_cool_path

    @staticmethod
    def _get_sat_level(matrix, percent):
        contact_array = matrix.data

        if contact_array.size == 0:
            return np.nan

        return np.percentile(contact_array, percent)

    @staticmethod 
    def _get_sat_score(matrix, bins, depletion, sat_level):
        
        scores = np.zeros(len(bins))
        
        s = matrix.multiply(matrix > sat_level)
        matrix[s.tocoo().col, s.tocoo().row] = sat_level

        for i in range(len(bins.index)):
            h = (i - depletion, i)
            v = (i + 1, i + depletion + 1)

            if h[0] < 0 or v[1] > len(bins.index):
                scores[i] = 0
            else:
                scores[i] = tril(matrix[h[0]:h[1], v[0]:v[1]]).sum()

        return scores
    
    @staticmethod
    def _get_sat_exp(sat_level, depletion):
        return depletion * (depletion + 1) * sat_level / 2
    
    @staticmethod
    def _get_wide_mismatch(scores, sensitive):
        troughs, _ = find_peaks(-np.log10(scores), prominence=sensitive)

        return troughs 

    def get_wide_mismatch(self, matrix, bins, contig):

        _m = matrix.fetch(contig)
        _bins = bins.fetch(contig)
        
        if len(_bins) <= 2 * self.depletion:
            return contig, None

        ## set diagonal to zero 
        _m.setdiag(0)
        _m.eliminate_zeros()
        _m = _m.tocsr() 
        if _m.nnz == 0:
            return contig, None
        ## get sat level and sat exp value
        sat_level = self._get_sat_level(_m, self.percent * 100)
        sat_exp = self._get_sat_exp(sat_level, self.depletion)
        
        scores = self._get_sat_score(_m, _bins, self.depletion, sat_level)
        
        troughs = self._get_wide_mismatch(scores, self.sensitive)
        
        return contig, _bins.iloc[troughs]
    
    def _fine_location(self, contig, contig_df_list):
    
        if contig_df_list[0] is None:
            return contig, None

        contig_df_list = filter(lambda x: x is not None, contig_df_list)
        if not contig_df_list:
            return contig, None
        contig_gr_list = list(map(lambda x: PyRanges(
                    chromosomes=x.chrom, starts=x.start, ends=x.end), 
                    contig_df_list))
        try:
            gr_init = contig_gr_list[0]
            gr = contig_gr_list[0]
           
        except IndexError:
            return contig, None
        if len(contig_gr_list) > 1:
            for i in range(1, len(contig_gr_list)):
                
                gr_init = gr_init.slack(self.resolutions[i]).intersect(contig_gr_list[i])
                
                if gr_init.empty is not True:
                    # gr_init.slack(self.resolutions[i])
                    gr = gr.join(gr_init, how='left').new_position("swap")
                    starts = gr.Start
                    ends = gr.End
                    idx1 = np.where(starts == -1)
                    idx2 = np.where(ends == -1)
                    idx = np.unique(np.r_[idx1, idx2])
                    df = gr.df
                
                
                    df.loc[idx, 'Start'] = df.loc[idx, 'Start_b']
                    df.loc[idx, 'End'] = df.loc[idx, 'End_b']
                
                    gr = PyRanges(df[['Chromosome', 'Start', 'End']])
                    
        return (contig, None) if gr.empty else (contig, gr.df)
            
    def fine_location(self, matrix, bins, contig, coarsen_position):
        print(coarsen_position)
        _m = matrix.fetch(contig)
        _bins = bins.fetch(contig)
        _bins_length = len(_bins)
        if len(_bins) <= 2 * self.depletion:
            return
        
        ## set diagonal to zero 
        _m.setdiag(0)
        _m.eliminate_zeros()
        _m = _m.tocsr() 
        if _m.nnz == 0:
            return 
        ## get sat level and sat exp value
        coarsen_position = coarsen_position.assign(
                                            start=lambda x: x["start"] // 1000,
                                            end=lambda x: x["end"] // 1000
        )
        sat_level = self._get_sat_level(_m, self.percent * 100)
        sat_exp = self._get_sat_exp(sat_level, self.depletion)

        for i, item in coarsen_position.iterrows():
            s, e = item.start, item.end
            s, e = s - 2*self.depletion, e + 2*self.depletion
            if s < 0:
                s = 0
            if e > _bins_length:
                e = _bins_length

            tmp_bins = _bins.iloc[s: e]
            tmp_matrix = _m[s: e, s: e]

            scores = self._get_sat_score(tmp_matrix, tmp_bins, self.depletion, sat_level)
            
            pos = self._get_wide_mismatch(scores, self.sensitive)
            
            pos = np.array(pos)
            pos = pos[(pos >= 2*self.depletion) & (pos < (len(scores) - 2*self.depletion))]
            #scores = scores[2*self.depletion: -2*self.depletion]
            print(sat_level, sat_exp, contig, scores, pos)

        #return scores, sat_exp

    @property
    def corrected_positions(self):
        self.uncorrected_count = 0
        self.corrected_count = 0
        db = []
        for contig, pos in self.chimeric_positions:
            length = self.chrom_sizes.loc[contig, 1]
            if pos is None:
                db.append((contig, 0, length, contig, 0))
                self.uncorrected_count += 1
            else:
                pos = pos[['Start', 'End']].values.flatten()
                pos = sorted(pos)
                start = np.r_[[0], pos]
                end = np.r_[pos, [length]]
                pos_list = list(zip(start, end))

                for start, end in pos_list:
                    db.append((contig, start, end, f'{contig}|{start}|{end}', 1))
                self.corrected_count += 1
        df = pd.DataFrame(db, columns=['chrom', 'start', 'end', 'name', 'flag'])
        
        return df
    
    def write_bed(self):
        """
        write corrected posistion as a bed file.
        """
        (self.corrected_positions[['chrom', 'start', 'end', 'name']]
        .to_csv(self.outbed, sep='\t', header=None, index=None)
        )
        logger.info(f'Written corrected position information into `{self.outbed}`.')

    def write_fasta(self):
        fasta_db = read_fasta(self.fasta)
        output = xopen(self.output, 'w')
        corrected_positions = self.corrected_positions
        for _, item in corrected_positions.iterrows():
            contig, start, end,name, flag = item
            seq = fasta_db[contig]
            if flag == 0:
                print(f'>{contig}\n{str(seq)}', file=output)
            else:
                print(f'>{name}\n{str(seq[start: end])}', 
                        file=output)
        
        logger.info(f'Corrected: {self.corrected_count}')
        logger.info(f'Written corrected fastat into `{self.output}`.')
    
    def write_pairs(self):
        corrected_positions = self.corrected_positions[['chrom', 'start', 'end', 'name']]
        corrected_positions.columns = ['Chromosome', 'Start', 'End', 'Name']
       
        
        pairs = Pairs(self.pairs)
        pairs.chrom2contig(corrected_positions, output=self.outpairs, 
                            threads=self.threads)

    def run(self):
        
        logger.info(f'Generating contacts on multi resolutions...')
        cool_path_db = Parallel(n_jobs=min(self.threads, len(self.resolutions)))(
                            delayed(self.get_contacts)(r) 
                            for r in self.resolutions)
        cool_path_db = dict(cool_path_db)
        
        res = {}
        for resolution in self.resolutions:
            cool = cooler.Cooler(cool_path_db[resolution])
            bins = cool.bins()
            matrix = cool.matrix(balance=False, sparse=True)
    
            args = []
            for contig in cool.chromnames:
                args.append((matrix, bins, contig))

            scores = Parallel(n_jobs=self.threads)(
                        delayed(self.get_wide_mismatch)(i, j, k)
                            for i, j, k in args)
            
            db = dict(scores)
            del args, scores
            res[resolution] = db 

        args = []
        for contig in cool.chromnames:
            tmp = [res[r][contig] for r in self.resolutions]
            args.append((contig, tmp))

        self.chimeric_positions = Parallel(n_jobs=self.threads)(
                                        delayed(self._fine_location)
                                        (i, j) for i, j in args)
        
        self.write_fasta()
        if self.outbed:
            self.write_bed()
        if self.outpairs:
            if self.corrected_count != 0:
                self.write_pairs()
            else:
                shutil.copy(self.pairs, self.outpairs)