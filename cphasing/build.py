#!/usr/bin/env python


import logging
import os
import os.path as op
import sys

import pandas as pd
import portion 

from Bio import SeqIO
from pathlib import Path
from pyfaidx import Fasta
from pytools import natsorted

from .agp import agp2fasta
from .core import Tour
from .utilities import xopen, list_flatten

logger = logging.getLogger(__name__)

GAP = 'N'*100

def Build(fasta_file, output='groups.asm.fasta', corrected=False,
            output_agp='groups.agp', only_agp=False):
    if fasta_file[-3:] == ".gz":
        handle = xopen(fasta_file)
        fasta = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    else:
        fasta = Fasta(fasta_file)

    p = Path('./')

    tour_list = list(p.glob('*.tour'))
    assert tour_list, "No such of tour files in currerent path."
    
    tour_list = natsorted(tour_list, key=lambda x: x.name)
    tour_list = map(Tour, tour_list)
    

    agp = output_agp

    agp_res = []
    corrected_contigs = []
    for tour in tour_list:
        _tmp_corrected_contigs, _tmp_agp = tour.to_agp(fasta, corrected=corrected)
        agp_res.append(_tmp_agp)
        if _tmp_corrected_contigs:
            corrected_contigs.append(_tmp_corrected_contigs)
    else:
        agp_df = pd.concat(agp_res, axis=0)
        if corrected_contigs:
            corrected_contigs = list_flatten(corrected_contigs)
            corrected_contigs_id = set(map(lambda x: x[0], corrected_contigs))
            anchor_contigs_df = agp_df[agp_df[5].isin(corrected_contigs_id)]
            
            contig_ranges = {}
            for contig in corrected_contigs_id:
                try:
                    contig_length = fasta.faidx.index[contig].rlen
                except KeyError:
                    logger.warning(f"Skip contig `{contig}`")
                    continue
                except AttributeError:
                    try:
                        length = len(fasta[contig])
                    except KeyError:
                        continue 
                
                contig_range = portion.closed(1, contig_length)
                contig_ranges[contig] = contig_range 
          
            for contig, tmp_df in anchor_contigs_df.groupby(5):
                tmp_range = tmp_df[[6, 7]].values.tolist()
                for i in tmp_range:
                    r = portion.closed(i[0], i[1])
                    
                    contig_ranges[contig] -= r 

            unanchor_tigs = set(fasta.keys()) - set(agp_df[agp_df[4] == 'W'][5].values.tolist())

        else:
            unanchor_tigs = set(fasta.keys()) - set(agp_df[agp_df[4] == 'W'][5].values.tolist())

        unanchor_res = []
        for unanchor_tig in unanchor_tigs:
            try:
                length = fasta.faidx.index[unanchor_tig].rlen
            except KeyError:
                logger.warning(f"Skip contig `{unanchor_tig}`")
                continue 
            except AttributeError:
                try:
                    length = len(fasta[unanchor_tig])
                except KeyError:
                    logger.warning(f"Skip contig `{unanchor_tig}`")
                    continue 
            unanchor_res.append([unanchor_tig, 1, length, 1, 'W', unanchor_tig, 
                1, length, '+'])
        
        
        if corrected_contigs:
            for contig in contig_ranges:
                for r in contig_ranges[contig]:
                    if r.upper - r.lower > 1:
                        tig_start = r.lower if r.left == portion.CLOSED else r.lower + 1
                        tig_end = r.upper if r.right == portion.CLOSED else r.upper - 1
                        tig_length = tig_end - tig_start + 1
                        unanchor_res.append([f"{contig}:{tig_start}-{tig_end}", 1, 
                                         tig_length, 1, 'W', contig, 
                                         tig_start, tig_end, '+'])
                      
        unanchor_df = pd.DataFrame(unanchor_res)

    agp_df = pd.concat([agp_df, unanchor_df], axis=0)
    agp_df.reset_index(drop=True, inplace=True)
    agp_df.to_csv(agp, sep='\t', header=False, index=False)
    logger.info(f"Output agp file `{agp}`")
  
    if not only_agp:
        with xopen(output, 'w') as out:
            try:
                agp2fasta(agp, fasta.filename, out)
            except AttributeError:
                agp2fasta(agp, fasta_file, out)
    
    if corrected and corrected_contigs:
        tmp_agp_df = agp_df
        x = tmp_agp_df[(tmp_agp_df[4] == "W") & (tmp_agp_df[5].isin(corrected_contigs_id))].index

        tmp_agp_df.loc[x, 5] = tmp_agp_df.loc[x, 5] + ":" + tmp_agp_df.loc[x, 6].astype(str) + "-" + tmp_agp_df.loc[x, 7].astype(str)
        tmp_agp_df.loc[x, 7] = tmp_agp_df.loc[x, 7].astype(int) - tmp_agp_df.loc[x, 6].astype(int) + 1
        tmp_agp_df.loc[x, 6] = 1
        

        tmp_agp_df.to_csv(agp.replace(".agp", ".corrected.agp"), sep='\t', 
                            header=False, index=False)
        logger.info(f"Output chimeric corrected agp file `{agp.replace('.agp', '.corrected.agp')}`")
  

    
