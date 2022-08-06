#!/usr/bin/env python


import logging
import os
import os.path as op
import sys

import pandas as pd 


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
    logging.info('load agp file: `{}`'.format(agpfile))
    
    if split:
        tig_df = df[df[4] == 'W']
        gap_df = df[df[4] == 'U']
        tig_df.columns = AGP_NAMES_tig
        gap_df.columns = AGP_NAMES_gap
        tig_df = tig_df.astype(
            {'chrom': 'category', 'orientation': 'category'})
        gap_df = gap_df.astype({'chrom': 'category'})

        tig_df.set_index('chrom', inplace=True)
        gap_df.set_index('chrom', inplace=True)
        
        return tig_df, gap_df
    else:
        return df


def agp2cluster(agp, store=True):
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
            print("{}\t{}\t{}".format(i, len(cluster), " ".join(cluster)), file=store)

