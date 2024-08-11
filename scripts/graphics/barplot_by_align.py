#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot a barplot to show phasing accuracy
"""

import argparse
import logging
import os
import os.path as op
import sys


import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import colormaps as cmaps
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Patch

import pandas as pd
import numpy as np

from collections import defaultdict, OrderedDict
from natsort import natsort_keygen
from pathlib import Path

from cphasing.utilities import list_flatten
from cphasing.utilities import to_humanized, chrom_ticks_convert


PAF_HADER = ["contig1", "length1", "start1", "end1", "strand",
                 "contig2", "length2", "start2", "end2", "mismatch",
                 "matches", "mapq", "identity"]


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('paf', 
            help='')
    pOpt.add_argument('--paf2', default=False, help="paf by mapping to haplotype resolved assembly")
    pOpt.add_argument('--cen', help="Cen bed")
    pOpt.add_argument('-o', '--output', type=str,
            default=None, help='output file [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    df = pd.read_csv(args.paf, sep='\t', header=None, usecols=range(13),
                         names=PAF_HADER, index_col=None)
    df['identity'] = df['identity'].map(lambda x: x.replace("id:f:", "")).astype(float)
    df = df.loc[df['identity'] >= 0.99]
    df = df.sort_values(by=['contig2', 'start2'], key=natsort_keygen())

    df = df.loc[~df['contig1'].str.startswith('utg')]
    df = df.loc[~df['contig1'].str.startswith('utig')]
    df = df.loc[~df['contig1'].str.startswith('ctg')]
    df = df.loc[~df['contig1'].str.startswith('contig')]
    
    df['align_length2'] = df['end2'] - df['start2']
    df['align_length1'] = df['end1'] - df['start1']   


    if args.paf2:
        paf2 = pd.read_csv(args.paf2, sep='\t', header=None, usecols=range(13),
                         names=PAF_HADER, index_col=None)
        duplicates = paf2[paf2.duplicated(subset=['contig1'], keep=False)]

        paf2 = paf2.drop_duplicates(subset=['contig1'], keep=False)
        paf2['identity'] = paf2['identity'].map(lambda x: x.replace("id:f:", "")).astype(float)
        paf2 = paf2.loc[paf2['identity'] >= 0.99]
        paf2['trio'] = paf2['contig2'].str.split("_").map(lambda x: x[1])
        paf2['trio'] = paf2['trio'].map(lambda x: 1 if x== 'PATERNAL' else 0)
        trio_info = paf2[['contig1', 'trio']]
        trio_info = trio_info.set_index('contig1').to_dict()['trio']
    
    if args.cen:
        cen_bed = pd.read_csv(args.cen, sep='\t', header=None, 
                              index_col=None, usecols=range(0, 4),
                              names=['chrom', 'start', 'end', 'name'])
        cen_bed = cen_bed.loc[~cen_bed['name'].str.contains('arm')]
        
        cen_bed = cen_bed.set_index('chrom')
    

    ref_chromsizes = df.groupby('contig2', sort=False)['length2'].unique().map(lambda x: x[0])

    
    bluered_12 = list(map(lambda x: mpl.colors.rgb2hex(x.colors),  list(cmaps.bluered_12)))
    rdbu = list(map(lambda x: mpl.colors.rgb2hex(x.colors), list(cmaps.rdbu)))
    puor = list(map(lambda x: mpl.colors.rgb2hex(x.colors), list(cmaps.puor)))
    prgn = list(map(lambda x: mpl.colors.rgb2hex(x.colors), list(cmaps.prgn)))
    brbg = list(map(lambda x: mpl.colors.rgb2hex(x.colors), list(cmaps.brbg)))
    piyg = list(map(lambda x: mpl.colors.rgb2hex(x.colors), list(cmaps.piyg))) 
    gy = list(map(lambda x: mpl.colors.rgb2hex(x.colors), list(cmaps.rdgy)))[6:]
    amwg_blueyellowred = list(map(lambda x: mpl.colors.rgb2hex(x.colors), list(cmaps.amwg_blueyellowred)))
    cmp_flux = list(map(lambda x: mpl.colors.rgb2hex(x.colors), list(cmaps.cmp_flux)))[::-1]
    rdylbu = list(map(lambda x: mpl.colors.rgb2hex(x.colors), list(cmaps.rdylbu))) 
    rdbu.pop(4)
    puor.pop(4)
    prgn.pop(4)
    brbg.pop(4)
    piyg.pop(4)

    rdylbu.pop(4)

    color_db = list_flatten([rdbu, puor, prgn, piyg, brbg, gy, cmp_flux, amwg_blueyellowred, rdylbu, bluered_12])

    ref_chromsizes_db = dict(zip(ref_chromsizes.index, range(len(ref_chromsizes.index))) ) 

    patches = []
    cen_patches = []
    colors = []
    ii = 0
    idx_db = defaultdict(int)
    legends = []
    ploted_idx = OrderedDict()
    for i, tmp_df in df.groupby('contig1', sort=False):
        
        idx_length_count = defaultdict(int)
        legends.append(Patch(facecolor=color_db[ii], label=i, edgecolor='none', linewidth=0))
        for j, row in tmp_df.iterrows():
            
            if args.paf2:
                try:
                    trio_idx = trio_info[j]
                except:
                    continue
                idx = ref_chromsizes_db[row.contig2]
                colors.append(color_db[ii])
                patches.append(Rectangle((row.start2, idx + trio_idx * 0.3), row.end2 - row.start2, 0.2, edgecolor='k', linewidth=0.05))
                if row.contig2 not in ploted_idx:
                    ploted_idx[row.contig2] = set()
                
                ploted_idx[row.contig2].add(idx + trio_idx * 0.3)

            else:
                colors.append(color_db[ii])
                idx = ref_chromsizes_db[row.contig2]
                idx_length_count[idx] += abs(row.end2 - row.start2)
                patches.append(Rectangle((row.start2, idx + idx_db[idx]), row.end2 - row.start2, 0.2, edgecolor='k', linewidth=0.05))

        if not args.paf2:
            idx = max(idx_length_count, key=lambda x: idx_length_count[x])
            idx_db[idx] += 0.3

        ii += 1
    
    
    for chrom in ploted_idx:
        for idx in ploted_idx[chrom]:
            for cen_chrom, cen_row in cen_bed.loc[chrom].iterrows():
                cen_patches.append(Rectangle((cen_row.start, idx - 0.04), cen_row.end - cen_row.start, 0.28, 
                                                    linewidth=0, edgecolor="none", alpha=0.1))

    plt.rcParams['font.family'] = 'Arial'
    
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    
    

    patch_collection = ax.add_collection(PatchCollection(patches, edgecolor='none', linewidth=0))
    patch_collection.set_facecolor(colors)
    cen_patch_collection = ax.add_collection(PatchCollection(cen_patches,  edgecolor='none', linewidth=0))
    cen_patch_collection.set_facecolor('#bcbcbc')
    cen_patch_collection.set_alpha(0.5)
    ax.autoscale()
    ax.invert_yaxis()
    _xticks = np.r_[ax.get_xticks()[1:-1]]
    ax.set_xticks(_xticks)
    xticklabels = chrom_ticks_convert(_xticks)
    ax.set_xticklabels(xticklabels, fontsize=8)
    ax.tick_params(axis='y', which='major', width=0)
    plt.xlim(0)
    plt.yticks(list(map( lambda x: x+0.2, range(len(ref_chromsizes)))), 
               ref_chromsizes_db.keys(), fontsize=8)
    ax.legend(handles=legends, loc=(0.25, 0.3), fontsize=4,
                ncol=2,  frameon=False, bbox_to_anchor=(0.5, 0., 0.5, 0.5))
    sns.despine( bottom=True, left=True)
    output = Path(args.paf).stem if not args.output else args.output
    plt.savefig(f"{output}.barplot.png", dpi=600, bbox_inches='tight')
    plt.savefig(f"{output}.barplot.pdf", dpi=300, bbox_inches='tight')
        

if __name__ == "__main__":
    main(sys.argv[1:])