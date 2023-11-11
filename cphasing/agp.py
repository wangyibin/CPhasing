#!/usr/bin/env python


import logging
import os
import os.path as op
import shutil
import sys

import pandas as pd 

from collections import OrderedDict
from natsort import natsort_keygen
from joblib import Parallel, delayed
from pathlib import Path

from .core import Tour
from .utilities import get_contig_length

logger = logging.getLogger(__name__)

AGP_NAMES_tig = ['chrom', 'start', 'end', 'number',
                 'type', 'id', 'tig_start', 'tig_end', 'orientation']
AGP_NAMES_gap = ['chrom', 'start', 'end', 'number',
                 'type', 'length', 'type', 'linkage', 'evidence']
def import_agp(agpfile, split=True):
    """
    import agp file and return a dataframe
    """
    df = pd.read_csv(agpfile, sep='\t', comment='#',
                     header=None, index_col=None,)
    logger.info('load agp file: `{}`'.format(agpfile))
    
    if split:
        tig_df = df[df[4] == 'W']
        gap_df = df[df[4] == 'U']
        tig_df.columns = AGP_NAMES_tig
        gap_df.columns = AGP_NAMES_gap
        tig_df = tig_df.astype(
            {'start': 'int64', 
            'end': 'int64', 'tig_start': 'int64', 
            'tig_end': 'int64',
            'orientation': 'category'})
        gap_df = gap_df.astype({'start': 'int64',
                                'end': 'int64',
                                'number': 'int64',
                                'length': 'int64'})

        tig_cat_dtype = pd.CategoricalDtype(categories=tig_df['chrom'].unique(), ordered=True)
        tig_df['chrom'] = tig_df['chrom'].astype(tig_cat_dtype)

        gap_cat_dtype = pd.CategoricalDtype(categories=gap_df['chrom'].unique(), ordered=True)
        gap_df['chrom'] = gap_df['chrom'].astype(gap_cat_dtype)

        tig_df.set_index('chrom', inplace=True)
        gap_df.set_index('chrom', inplace=True)
        
        return tig_df, gap_df
    else:
        return df

def agp2assembly(agpfile, output, add_gap=False):
    """
    convert agp to assembly.
    """
   
    tig_df, gap_df = import_agp(agpfile)
    gap_length = gap_df.iloc[0]['length']
    chrom_matrix_db = OrderedDict()
    tig_list = []
    for i, item in enumerate(tig_df.iterrows(), 1):
        chrom = item[0]
        tig = item[1].id
        tig_length = item[1].tig_end
        orientation = item[1].orientation

        tig_list.append([tig, i, tig_length])
        if chrom not in chrom_matrix_db:
            chrom_matrix_db[chrom] = []
        orientation = orientation if orientation == '-' else ""
        chrom_matrix_db[chrom].append(orientation + str(i))
    else:
        hic_gap_number = i + 1

    for item in tig_list:
        print(">{}".format(" ".join(map(str, item))), 
                            file=output)
    else:
        if add_gap:
            print(">{}".format(" ".join(['hic_gap_{}'.format(
                hic_gap_number), str(hic_gap_number), str(gap_length)])), 
                file=output)
            _gap = " {} ".format(hic_gap_number)
        else: 
             _gap = " "
        for chrom in chrom_matrix_db:
            print(_gap.join(
                map(str, chrom_matrix_db[chrom])), file=output)
    
    logger.info(f"Successful output assembly file into `{output.name}`.")



def agp2cluster(agp, store=None):
    """
    Convert agp to cluster file.

    Params:
    --------
    agp: str
    store: _io.TextWrapper
        output handle
    
    Returns:
    --------
    pd.DataFrame:
        Cluster dataframe
    
    Examples:
    --------
    >>> agp2cluster('groups.agp')

    """
    agp_df, _ = import_agp(agp)
    
    # remove contig
    agp_df.reset_index(inplace=True)
    agp_df = agp_df[agp_df['chrom'] != agp_df['id']]
    agp_df['chrom'] = agp_df['chrom'].astype('category')
    cluster_df = agp_df.groupby('chrom')['id'].apply(lambda x: list(x))
    
    if store:
        for i, cluster in cluster_df.iteritems():
            if not cluster:
                continue
            print("{}\t{}\t{}".format(i, len(cluster), " ".join(cluster)), 
                        file=store)

    return cluster_df

def agp2fasta(agp, fasta, output=sys.stdout, threads=1):
    """
    Convert agp to chromosome-level fasta file.

    Params:
    --------
    agp: str
        Path to agp file.
    fasta: str
        Path to fasta file.
    output: _io.TextIOWrapper
        Output handle
    threads: int
        Number of threads
    
    Returns:
    --------
    None

    Examples:
    --------
    >>> agp2fasta('groups.agp', 'contigs.fasta')
    """
    # from pyfaidx import Fasta
    # def get_seqs(chrom, cluster):
    #     seqs = []
    #     id_orient = cluster.loc[:, ['id', 'orientation']].values.tolist()
    #     for contig, orient in id_orient:
    #         if orient == '+':
    #             seqs.append(str(fasta[contig]))
    #         else:
    #             seqs.append(fasta[contig][::-1].complement.seq)
        
    #     out_seq = GAP.join(seqs)

    #     return chrom, out_seq
    
    # if isinstance(fasta, str):
    #     fasta = Fasta(fasta)
    # elif isinstance(fasta, Fasta):
    #     pass
    from .utilities import read_fasta
    seq_db = read_fasta(fasta)
    agp_df, gap_df = import_agp(agp)
   
    GAP = 'N' * gap_df.length[0]

    cluster_df = agp_df.groupby('chrom')
    for chrom, cluster in cluster_df:
        seqs = []
        id_orient = cluster.loc[:, ['id', 'orientation']].values.tolist()
    #outs = Parallel(n_jobs=threads, verbose=10)(delayed(get_seqs)(chrom, cluster) for chrom, cluster in cluster_df)
    #print(outs)
        for contig, orient in id_orient:
            if orient == '+':
                seqs.append(str(seq_db[contig]))
            else:
                seqs.append(str(seq_db[contig].reverse_complement()))
        
        out_seq = GAP.join(seqs)
        print(f'>{chrom}', file=output)
        print(out_seq, file=output)

    logger.info(f"Output chromosome-level fasta into `{output.name}`.")

def agp2tour(agp, outdir, force=False):
    """
    Convert agp to several tour files.

    Params:
    --------
    agp: str
        Path to agp file.
    outdir: str
        output directory of tour results.
    force: bool
        whether remove exists directory.
    
    Returns:
    --------
    None

    Examples:
    --------
    >>> agp2tour('groups.agp', 'tour')
    """
    agp_df, _ = import_agp(agp)

    # remove contig
    agp_df.reset_index(inplace=True)
    agp_df = agp_df[agp_df['chrom'] != agp_df['id']]
    agp_df['chrom'] = agp_df['chrom'].astype('category')

    cluster_df = agp_df.groupby('chrom')
    tour = Tour
    
    outdir = Path(outdir)
    if outdir.exists():
        if force:
            logger.info(f'Force output results, removed `{outdir}`')
            shutil.rmtree(outdir)
        else:
            logger.warn(f'The output directory of `{outdir}` exists.')
    outdir.mkdir(parents=True, exist_ok=True)

    for group, cluster in cluster_df:
        if cluster.empty is True:
            continue
        
        tour.from_tuples(cluster[['id', 'orientation']].values.tolist())
        with open(f'{str(outdir)}/{group}.tour', 'w') as out:
            print(' '.join(map(str, tour.data)), file=out)
            logger.info(f'Output tour: `{out.name}`')

    logger.info('ALL done.')

def statagp(agp, output):
    """
    Statistics of AGP.

    Params:
    --------
    agp: str
        Path to agp file
    output: _io.TextIOWrapper
        Output handle
    
    Returns:
    --------
    None

    Examples:
    --------
    >>> agp('groups.agp', sys.stdout)
    ChrID   Anchored_ctg    Length
    test    5       54
    Total Number of contigs:        6
    Total Length of contigs (bp):   155
    Total Number of anchored contigs:       5
    Total Length of chromosome level assembly (bp): 454
    Number of unanchored contigs:   1
    Length of unanchored contigs (bp):      101
    Anchor rate (%):        34.84
    """
    agp_df, gap_df = import_agp(agp)
    
    ## get length of each contigs
    agp_df['length'] = agp_df['tig_end'] - agp_df['tig_start'] + 1
    
    ## remove contig
    rm_tig_agp_df = agp_df.reset_index()
    rm_tig_agp_df = rm_tig_agp_df[rm_tig_agp_df['chrom'] != rm_tig_agp_df['id']]
    rm_tig_agp_df['chrom'] = rm_tig_agp_df['chrom'].astype('object')
    
    ## get anchored contigs dataframe
    tig_agp_df = agp_df.reset_index()
    tig_agp_df = tig_agp_df[tig_agp_df['chrom'] == tig_agp_df['id']]
    tig_agp_df['chrom'] = tig_agp_df['chrom'].astype('object')

    chrom_lengths = rm_tig_agp_df.groupby('chrom')['length'].sum()
    chrom_contigs_counts = rm_tig_agp_df.groupby('chrom')['id'].count()
    chrom_infos = pd.concat([chrom_contigs_counts, chrom_lengths], axis=1)
    chrom_infos = chrom_infos.sort_values('chrom', key=natsort_keygen())

    total_number_of_contigs = agp_df['id'].count()
    total_length_of_contigs = agp_df['length'].sum()
    total_number_of_anchored_contigs = rm_tig_agp_df['id'].count()

    gap_length = gap_df['length'].sum()
    chrom_length = rm_tig_agp_df['length'].sum()
    total_length_of_chromosome_level_assembly = chrom_length + gap_length

    number_of_unanchored_contigs = tig_agp_df['id'].count()
    length_of_unanchored_contigs = tig_agp_df['length'].sum()

    anchor_rate = chrom_length / (chrom_length + length_of_unanchored_contigs)
    
    ## output
    chrom_infos.reset_index(inplace=True)
    chrom_infos.columns = ['ChrID', 'Anchored_ctg', 'Length']
    chrom_infos.to_csv(output, sep='\t', index=False)
    print(f'Total Number of contigs:\t{total_number_of_contigs}',
        f'Total Length of contigs (bp):\t{total_length_of_contigs}', 
        f'Total Number of anchored contigs:\t{total_number_of_anchored_contigs}',
        f'Total Length of chromosome level assembly (bp):\t{total_length_of_chromosome_level_assembly}',
        f'Number of unanchored contigs:\t{number_of_unanchored_contigs}',
        f'Length of unanchored contigs (bp):\t{length_of_unanchored_contigs}',
        f'Anchor rate (%):\t' + f'{anchor_rate:.2%}'.replace("%", ""), 
        sep='\n', file=output)
    
    
def pseudo_agp(real_list, contigsizes, output):
    """
    Create a pseudo agp file from a list and contigsizes.
    """
    real_list = pd.read_csv(real_list, sep='\t', header=None, index_col=None)
    contig_sizes = dict(i.strip().split() for i in open(contigsizes) if i.strip())

    idx = 0
    for chrom, contig_df in real_list.groupby(0):
        start = 1 
        end = 1
        idx += 1 
        for j, contig in enumerate(contig_df[1].values.tolist()):
            start = end
            contig_length = int(contig_sizes[contig])
            end = start + contig_length - 1
            print("\t".join(map(str, [chrom, start, end, idx, "W", contig, 1, contig_length, "+"])), file=output)
            if j < len(contig_df) - 1:
                start = end
                end = start + 100
                print("\t".join(map(str, [chrom, start, end, idx, "U", 100, "contig", "yes", "map",])), file=output)

        
