#!/usr/bin/env python


import logging
import os
import os.path as op
import sys

import pandas as pd

from Bio import SeqIO
from pathlib import Path
from pyfaidx import Fasta
from pytools import natsorted

from .agp import agp2fasta
from .core import Tour
from .utilities import xopen 

logger = logging.getLogger(__name__)

GAP = 'N'*100

def Build(fasta, output='groups.asm.fasta', 
            output_agp='groups.agp', only_agp=False):
    if fasta[-3:] == ".gz":
        handle = xopen(fasta)
        fasta = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    else:
        fasta = Fasta(fasta)
    p = Path('./')

    tour_list = list(p.glob('*.tour'))
    assert tour_list, "No such of tour files in currerent path."
    
    tour_list = natsorted(tour_list, key=lambda x: x.name)
    tour_list = map(Tour, tour_list)
    

    agp = output_agp

    agp_res = []
    for tour in tour_list:
        _tmp_agp = tour.to_agp(fasta)
        agp_res.append(_tmp_agp)
    else:
        agp_df = pd.concat(agp_res, axis=0)
        unanchor_tigs = set(fasta.keys()) - set(agp_df[agp_df[4] == 'W'][5].values.tolist())
        unanchor_res = []
        for unanchor_tig in unanchor_tigs:
            length = len(fasta[unanchor_tig])
            unanchor_res.append([unanchor_tig, 1, length, 1, 'W', unanchor_tig, 
                1, length, '+'])
        unanchor_df = pd.DataFrame(unanchor_res)
    
    agp_df = pd.concat([agp_df, unanchor_df], axis=0)
    agp_df.to_csv(agp, sep='\t', header=False, index=False)
    
    if not only_agp:
        with xopen(output, 'w') as out:
            agp2fasta(agp, fasta.filename, out)
    
    
