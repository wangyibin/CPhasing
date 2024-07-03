#!/usr/bin/env python


import logging
import os
import os.path as op
import shutil
import sys

import pandas as pd 

from collections import OrderedDict, defaultdict
from natsort import natsort_keygen
from joblib import Parallel, delayed
from pathlib import Path

from .core import Tour
from .utilities import get_contig_length

logger = logging.getLogger(__name__)

AGP_NAMES_tig = ['chrom', 'start', 'end', 'number',
                 'type', 'id', 'tig_start', 'tig_end', 'orientation']
AGP_NAMES_gap = ['chrom', 'start', 'end', 'number',
                 'type', 'length', 'gap_type', 'linkage', 'evidence']
def import_agp(agpfile, split=True):
    """
    import agp file and return a dataframe
    """
    df = pd.read_csv(agpfile, sep='\t', comment='#',
                     header=None, index_col=None,)
    logger.info('Load agp file: `{}`'.format(agpfile))

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
        for i, cluster in cluster_df.items():
            if not cluster:
                continue
            print("{}\t{}\t{}".format(i, len(cluster), " ".join(cluster)), 
                        file=store)

    return cluster_df

def agp2fasta(agp, fasta, output=sys.stdout, output_contig=False, threads=1):
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

    if not output_contig:
        cluster_df = agp_df.groupby('chrom')
        for chrom, cluster in cluster_df:
            seqs = []
            id_orient = cluster.loc[:, ['id', 'tig_start', 'tig_end', 'orientation']].values.tolist()
        
            for contig, contig_start, contig_end, orient in id_orient:
                if orient == '+':
                    seqs.append(str(seq_db[contig][contig_start - 1: contig_end]))
                else:
                    seqs.append(str(seq_db[contig][contig_start - 1: contig_end].reverse_complement()))
            
            out_seq = GAP.join(seqs)
            print(f'>{chrom}', file=output)
            print(out_seq, file=output)

        logger.info(f"Output chromosome-level fasta into `{output.name}`.")
    
    else:
        tmp_data = agp_df.reset_index()[['chrom', 'id', 'tig_start', 'tig_end']].values.tolist()
        for record in tmp_data:
            output_contig, raw_contig, start, end = record 
            seq = seq_db[raw_contig]
            seq_length = len(seq)
            if start == 1 and end == seq_length:
                seq = str(seq)
                output_contig = raw_contig
            else:
                seq = str(seq[start-1: end])
                output_contig = f"{raw_contig}:{start}-{end}"

            print(f'>{output_contig}', file=output)
            print(seq, file=output)

        logger.info(f"Output contig-level fasta into `{output.name}`.")

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

        
def assembly2agp(assembly, 
                 chrom_num=None, sort_by_length=False,
                 gap_len=100,
                 phased=False, chrom_prefix="Chr",
                 outprefix=None):
    
    logger.info(f"Converting assembly to agp ...")
    if not outprefix:
        outprefix = Path(assembly).absolute().stem 

    contig_db = OrderedDict()
    contig_idx_db = OrderedDict()
    contig_idx_length_db = OrderedDict()
    fragment_db = OrderedDict()
    scaffolds = []
    with open(assembly) as fp:
        for line in fp:
            if line.strip().startswith(">"):
                contig, contig_idx, contig_len = line[1:].split()
                contig_db[contig] = (int(contig_idx), int(contig_len))
                contig_idx_db[int(contig_idx)] = contig 
                contig_idx_length_db[int(contig_idx)] = int(contig_len)

                if ":::" in contig:
                    raw_contig_data = contig.split(":::")
                    raw_contig = raw_contig_data[0]
                    raw_fragment_idx = int(raw_contig_data[1].replace("fragment_", ""))
        

                    if raw_contig not in fragment_db:
                        fragment_db[raw_contig] = []
                    fragment_db[raw_contig].append((contig, raw_fragment_idx, int(contig_len)))
            else:
                scaffolds.append(line.strip().split())

    
    fragments_range = OrderedDict()
    for raw_contig in fragment_db:
        fragments = sorted(fragment_db[raw_contig], key=lambda x: x[1])
        start = 1 
        for fragment in fragments:
            contig, fragment_idx, fragment_length = fragment 
            fragments_range[contig] = (raw_contig, start, start + fragment_length)
            start += fragment_length
    
    # print(sorted(fragments_range.items(), key=lambda x: x[1][1]))


    scaffolds_length_db = OrderedDict()
    raw_agp_records = OrderedDict()
    agp_records = OrderedDict()
    for i, group in enumerate(scaffolds):
        raw_record = []
        record = []
        record_num = 0
        start = 1 
        for idx in group:
            if idx.startswith("-"):
                orient = "-"
                idx = int(idx[1:])
            else:
                orient = "+"
                idx = int(idx )
            
            contig_length = contig_idx_length_db[idx]
            if i not in scaffolds_length_db:
                scaffolds_length_db[i] = 0

            scaffolds_length_db[i] += contig_length

            contig = contig_idx_db[idx]
            
           
    
            if contig in fragments_range:
                raw_contig, raw_contig_start, raw_contig_end = fragments_range[contig]
                # contig = f'{raw_contig}:{raw_contig_start}-{raw_contig_end}'

            else:
                raw_contig, raw_contig_start, raw_contig_end = contig, 1, contig_length
        
            record.append((i, start, start + contig_length, record_num,
                            "W", contig, 1, contig_length, orient))
            
            raw_record.append((
                i, start, start + contig_length - 1, record_num,
                "W", raw_contig, raw_contig_start, raw_contig_end, orient,
            ))
    

           
            record_num += 1 
            start += contig_length + gap_len
        
        agp_records[i] = record 
        raw_agp_records[i] = raw_record
        
        
    if sort_by_length:
        scaffolds_idx = sorted(scaffolds_length_db.items(), key=lambda x: x[1], reverse=True)
    else:
        scaffolds_idx = list(scaffolds_length_db.keys())
    

    raw_agp = f"{outprefix}.raw.agp"
    agp = f"{outprefix}.agp"

    with open(raw_agp, "w") as raw_out, open(agp, "w") as out:
        record_num = 1
        for group_idx in scaffolds_idx:
            records = agp_records[group_idx]
            raw_records = raw_agp_records[group_idx]
            unanchor_records = []
            if chrom_num:
                total_chrom_num = chrom_num[0] if len(chrom_num) == 1 else chrom_num[0] * chrom_num[1]

                for i, (record, raw_record) in enumerate(zip(records, raw_records)):
                    num = group_idx + 1
                    
                    if num > total_chrom_num:
                        unanchor_records.append(record)
                    else:
                        _, start, end, _, _type, contig, contig_start, contig_end, orient = record

                        if len(chrom_num ) == 1:
                            if phased:
                                chrom = f"{chrom_prefix}g{num}"
                            else:
                                chrom = f"{chrom_prefix}{num:0>2}"
                        else:
                            n1 = (num - 1) // chrom_num[1] + 1
                            n2 = (num - 1) % chrom_num[1] + 1
                           
                            chrom = f"{chrom_prefix}{n1:0>2}g{n2}"
                        
                        if contig in fragments_range:
                            raw_contig, raw_contig_start, raw_contig_end = fragments_range[contig]
                            contig = f"{raw_contig}:{raw_contig_start}-{raw_contig_end-1}"
                            raw_contig_end = raw_contig_end - 1
                            
                        print("\t".join(map(str, (chrom, start, end - 1, record_num, _type, 
                                                    contig, contig_start, contig_end, orient))), file=out)

                        _, raw_start, raw_end, _, _type, raw_contig, contig_start, contig_end, orient = raw_record
                        if raw_contig in fragment_db:
                            contig_end = contig_end - 1
                        print("\t".join(map(str, (chrom, raw_start, raw_end, record_num, _type, 
                                                    raw_contig, contig_start, contig_end, orient))), file=raw_out)
                        
                        if i < len(records) - 1:
                            record_num += 1
                            print("\t".join(map(str, (chrom, end, end + gap_len - 1, record_num, 
                                                        "U", gap_len, "contig", "yes", "map" ))), file=out)
                            print("\t".join(map(str, (chrom, raw_end + 1, raw_end+ gap_len, record_num, 
                                                        "U", gap_len, "contig", "yes", "map" ))), file=raw_out)
                        record_num += 1
                        
                else:
                    for record in unanchor_records:
                        _, start, end, _, _type, contig, contig_start, contig_end, orient = record
                        if contig in fragments_range:
                            raw_contig, raw_contig_start, raw_contig_end = fragments_range[contig]
                            chrom = f"{raw_contig}:{raw_contig_start}-{raw_contig_end - 1}"
                            contig = chrom
                            raw_contig_end = raw_contig_end - 1
                        else:
                            chrom = contig
                            raw_contig = contig
                            raw_contig_start, raw_contig_end = contig_start, contig_end

                        print("\t".join(map(str, (chrom, 1, raw_contig_end - raw_contig_start + 1 , record_num, _type, 
                                                    contig, 1, raw_contig_end - raw_contig_start + 1, '+'))), file=out)
                        

                        print("\t".join(map(str, (chrom, 1, raw_contig_end - raw_contig_start + 1, record_num, _type, 
                                                    raw_contig, raw_contig_start, raw_contig_end, '+'))), file=raw_out)
                        record_num += 1

            else:
                for i, (record, raw_record) in enumerate(zip(records, raw_records)):
                    chrom_num = group_idx + 1
                    _, start, end, _, _type, contig, contig_start, contig_end, orient = record

                    if phased:
                        chrom = f"{chrom_prefix}g{chrom_num}"
                    else:
                        chrom = f"{chrom_prefix}{chrom_num:0>2}"
                    
                    if contig in fragment_db:
                        raw_contig, raw_contig_start, raw_contig_end = fragments_range[contig]
                        contig = f"{raw_contig}:{raw_contig_start}-{raw_contig_end}"
                        print(raw_contig, raw_contig_start, raw_contig_end)
                    print("\t".join(map(str, (chrom, 1, contig_end - contig_start + 1, record_num, _type, 
                                                contig, contig_start, contig_end, orient))), file=out)
                    
                    _, raw_start, raw_end, _, _type, raw_contig, contig_start, contig_end, orient = raw_record

                    print("\t".join(map(str, (chrom, raw_start, raw_end, record_num, _type, 
                                                raw_contig, contig_start, contig_end, orient))), file=raw_out)
                    
                    if i < len(records) - 1:
                        record_num += 1
                        print("\t".join(map(str, (chrom, end + 1, end + gap_len, record_num, 
                                                    "U", 100, "contig", "yes", "map" ))), file=out)
                        print("\t".join(map(str, (chrom, raw_end +  1, raw_end + gap_len, record_num, 
                                                    "U", 100, "contig", "yes", "map" ))), file=raw_out)
                    record_num += 1

        logger.info(f"Successful output agp file `{agp}`.")
        logger.info(f"Successful output agp file `{raw_agp}`.")

    
    # if len(fragments_range) > 1:
        
    #     output_corrected_fasta = f"{outprefix}.corrected.fasta"
    #     with open(output_corrected_fasta, 'w') as output_fasta:
    #         agp2fasta(raw_agp, fasta, output=output_fasta, output_contig=True)
    # else:
    #     logger.info(f"Fragmented contig not found, the `{raw_agp}` equale to `{agp}`, "
    #                 f"and fasta can be used to subsequence analysis.")