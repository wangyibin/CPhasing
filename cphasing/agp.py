#!/usr/bin/env python


import logging
import os
import os.path as op
import shutil
import sys
import re

import io 

import pandas as pd 

from collections import OrderedDict, defaultdict
from natsort import natsort_keygen
from joblib import Parallel, delayed
from pathlib import Path

from .core import Tour
from .utilities import get_contig_length, read_fasta

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
    try:
        gap_length = gap_df.iloc[0]['length']
    except:
        gap_length = 100

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

def agp2fasta_v0(agp, fasta, output=sys.stdout, 
              output_contig=False, 
              skip_gap=False,
              threads=1):
    """
    Convert agp to chromosome-level fasta file.

    Params:
    --------
    agp: str
        Path to agp file.
    fasta: OrderedDict
        Fasta sequence database.
    output: _io.TextIOWrapper
        Output handle
    threads: int
        Number of threads
    
    Returns:
    --------
    None

    Examples:
    --------
    >>> fasta_db = read_fasta('contigs.fasta')
    >>> agp2fasta('groups.agp', fasta_db)
    """
    agp_df, gap_df = import_agp(agp)
    seq_db = fasta
   
    if skip_gap:
        GAP = ''
    else:
        try:
            GAP = 'N' * gap_df.length[0]
        except IndexError:
            GAP = 'N' * 100

    if not output_contig:
        cluster_df = agp_df.groupby('chrom')
        for chrom, cluster in cluster_df:
            seqs = []
            id_orient = cluster.loc[:, ['id', 'tig_start', 'tig_end', 'orientation']].values.tolist()
        
            for contig, contig_start, contig_end, orient in id_orient:
                fetch_contig = contig 
                if fetch_contig not in seq_db:
                    fetch_contig = re.sub(r'_d\d+$', '', contig)

                if orient == '+':
                    seqs.append(str(seq_db[fetch_contig][contig_start - 1: contig_end]))
                else:
                    seqs.append(str(seq_db[fetch_contig][contig_start - 1: contig_end].reverse_complement()))
            
            out_seq = GAP.join(seqs)
            print(f'>{chrom}', file=output)
            print(out_seq, file=output)

        logger.info(f"Output chromosome-level fasta into `{output.name}`.")
    
    else:
        tmp_data = agp_df.reset_index()[['chrom', 'id', 'tig_start', 'tig_end']].values.tolist()
        contig_idx = defaultdict(int)
        collapsed_rescue_flag = False 
        collapsed_rescued_contigs = []
        for record in tmp_data:
            output_contig, raw_contig, start, end = record 
            fetch_contig = raw_contig
            if fetch_contig not in seq_db:
                fetch_contig = re.sub(r'_d\d+$', '', contig)

            try:
                seq = seq_db[fetch_contig]
            except KeyError:
                logger.warning(f"Counld not found `{raw_contig}`, skipped.")
                continue
            seq_length = len(seq)
            if start == 1 and end == seq_length:
                seq = str(seq)
                output_contig = raw_contig
                contig_idx[output_contig] += 1
            else:
                seq = str(seq[start-1: end])
                output_contig = f"{raw_contig}:{start}-{end}"
                contig_idx[output_contig]

            if contig_idx[output_contig] > 1:
                logger.info(f"Duplicated contig `{output_contig}` find in agp records, suggest it is a collapsed contig, rename it to `{output_contig}_d{contig_idx[output_contig]}`")
                raw_contig = output_contig
                output_contig = f"{output_contig}_d{contig_idx[output_contig]}"
                collapsed_rescue_flag = True
                collapsed_rescued_contigs.append((raw_contig, output_contig))

            print(f'>{output_contig}', file=output)
            print(seq, file=output)

        logger.info(f"Output contig-level fasta into `{output.name}`.")

        if collapsed_rescue_flag:
            output_file = output.name.replace(".fasta", ".collapsed.contig.list")
            with open(f"{output_file}", "w") as out:
                for raw_contig, output_contig in collapsed_rescued_contigs:
                    out.write(f"{raw_contig}\t{output_contig}\n")
            logger.info(f"Output collapsed contigs: `{output_file}`")

def agp2fasta(agp, fasta, output=sys.stdout, 
              output_contig=False, 
              skip_gap=False,
              threads=1):
    """
    Convert agp to chromosome-level fasta file.

    Params:
    --------
    agp: str
        Path to agp file.
    fasta: OrderedDict
        Fasta sequence database.
    output: _io.TextIOWrapper
        Output handle
    output_contig: bool or str
        Output mode. 
        - False (or "chr"): Output chromosome-level sequences containing Ns (default).
        - True (or "raw"): Output raw individual sliced contigs.
        - "block": Join gapless adjacent contigs into blocks (scaffolds) and break output when hitting a gap.
    skip_gap: bool
        If True, ignore gap lines (U/N) and concatenate sequences without Ns.
    threads: int
        Number of threads
    
    Returns:
    --------
    None
    """
    agp_df = import_agp(agp, split=False)
    if isinstance(fasta, str):
        seq_db = read_fasta(fasta)
    else:
        seq_db = fasta
    
    if isinstance(output, str):
        output = open(output, "w")
        
   
    if output_contig is True or output_contig == 'block':
        mode = 'block'
    elif output_contig == 'raw':
        mode = 'raw'
    else:
        mode = 'chr'

    if mode == 'raw':
        agp_df_w = agp_df[agp_df[4].isin(['W', 'D', 'F'])]
        tmp_data = agp_df_w[[0, 5, 6, 7]].values.tolist()
        
        contig_idx = defaultdict(int)
        collapsed_rescue_flag = False 
        collapsed_rescued_contigs = []
        
        for record in tmp_data:
            chrom, raw_contig, start, end = record 
            start, end = int(start), int(end)
            
            fetch_contig = raw_contig
            if fetch_contig not in seq_db:
                fetch_contig = re.sub(r'_d\d+$', '', raw_contig)

            try:
                seq = seq_db[fetch_contig]
            except KeyError:
                logger.warning(f"Counld not found `{raw_contig}`, skipped.")
                continue
                
            seq_length = len(seq)
            if start == 1 and end == seq_length:
                seq_str = str(seq)
                output_name = raw_contig
            else:
                seq_str = str(seq[start-1: end])
                output_name = f"{raw_contig}:{start}-{end}"
                
            contig_idx[output_name] += 1
            if contig_idx[output_name] > 1:
                logger.info(f"Duplicated contig `{output_name}` found, rename to `{output_name}_d{contig_idx[output_name]}`")
                old_name = output_name
                output_name = f"{output_name}_d{contig_idx[output_name]}"
                collapsed_rescue_flag = True
                collapsed_rescued_contigs.append((old_name, output_name))

            print(f'>{output_name}', file=output)
            print(seq_str, file=output)

        logger.info(f"Output raw contig-level fasta into `{output.name}`.")

        if collapsed_rescue_flag:
            try:
                output_file = str(output.name).replace(".fasta", ".collapsed.contig.list")
                with open(f"{output_file}", "w") as out:
                    for raw_ctg, out_ctg in collapsed_rescued_contigs:
                        out.write(f"{raw_ctg}\t{out_ctg}\n")
                logger.info(f"Output collapsed contigs: `{output_file}`")
            except Exception:
                pass
                
    else:
        block_idx = 1

        for chrom, cluster_df in agp_df.groupby(0):
            block_seqs = []
        
            block_names = [] 
            
            for _, row in cluster_df.iterrows():
                comp_type = row[4]
                
                if comp_type in ['W', 'D', 'F']:
                    contig = row[5]
                    tig_start = int(row[6])
                    tig_end = int(row[7])
                    orient = row[8]
                    
                    fetch_contig = contig 
                    if fetch_contig not in seq_db:
                        fetch_contig = re.sub(r'_d\d+$', '', contig)

                    try:
                        tmp_seq = seq_db[fetch_contig][tig_start - 1: tig_end]
                    except KeyError:
                        logger.warning(f"Could not found `{contig}`, skipped.")
                        continue
                        
                    if orient == '-':
                        tmp_seq = tmp_seq.reverse_complement()
                        
                    block_seqs.append(str(tmp_seq))
                    block_names.append(f"{contig}:{tig_start}-{tig_end}{orient}")
                    
                elif comp_type in ['U', 'N']:
                    gap_length = int(row[5])
                    
                    if mode == 'block':
                        if not skip_gap and len(block_seqs) > 0:
                            out_seq = "".join(block_seqs)
                            block_header = "|".join(block_names)
                            block_name = f"unitig{block_idx:>07}l"
                            print(f'>{block_name} {block_header}', file=output)
                            print(out_seq, file=output)
                            block_idx += 1
                            block_seqs = []
                            block_names = [] 
                    elif mode == 'chr':
                        if not skip_gap and gap_length > 0:
                            block_seqs.append("N" * gap_length)
                            
            if len(block_seqs) > 0:
                out_seq = "".join(block_seqs)
                if mode == 'block':
                    block_header = "|".join(block_names)
                    block_name = f"unitig{block_idx:>07}l"
                    print(f'>{block_name} {block_header}', file=output)
                    block_idx += 1
                else:
                    print(f'>{chrom}', file=output)
                print(out_seq, file=output)
                
        if mode == 'block':
            logger.info(f"Output gapless block-level fasta into `{output.name}`.")
        else:
            logger.info(f"Output chromosome-level fasta into `{output.name}`.")

def agp2tour(agp, outdir="tour", force=False, store=True):
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
    store: bool
        whether to store to files

    Returns:
    --------
    dict

    Examples:
    --------
    >>> agp2tour('groups.agp', 'tour')
    """
    agp_df, _ = import_agp(agp)

    # remove contig
    agp_df.reset_index(inplace=True)
    agp_df = agp_df[agp_df['chrom'] != agp_df['id']]
    agp_df['chrom'] = agp_df['chrom'].astype('category')
    duplicated_df = agp_df[agp_df.duplicated(subset=['id'], keep=False)]
    duplicated_contigs = set(duplicated_df['id'].values.tolist())

    cluster_df = agp_df.groupby('chrom')
    tour = Tour
    
    
    db = OrderedDict()
    for group, cluster in cluster_df:

        if cluster.empty is True:
            continue
        data = []
        for i, row in cluster.iterrows():
            if row['id'] in duplicated_contigs: 
                data.append((f"{row.id}:{row.tig_start}-{row.tig_end}", row.orientation))
            else:
                data.append((row.id, row.orientation))

        db[group] = data
            
    outdir = Path(outdir)
    if outdir.exists():
        if force:
            logger.info(f'Force output results, removed `{outdir}`')
            shutil.rmtree(outdir)
        else:
            logger.warning(f'The output directory of `{outdir}` exists.')
    outdir.mkdir(parents=True, exist_ok=True)
        
    for group, data in db.items():  
        
        tour.from_tuples(data)
        with open(f'{str(outdir)}/{group}.tour', 'w') as out:
            print(' '.join(map(str, tour.data)), file=out)
            logger.debug(f'Output tour: `{out.name}`')

    logger.info('ALL done.')

    return agp_df

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
    tig_agp_df['chrom'] = tig_agp_df['chrom'].astype('object')
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
    print(f'Total number of contigs:\t{total_number_of_contigs}',
        f'Total length of contigs (bp):\t{total_length_of_contigs}', 
        f'Total number of anchored contigs:\t{total_number_of_anchored_contigs}',
        f'Total length of anchored contigs (bp):\t{chrom_length}',
        f'Total length of chromosome level assembly (bp):\t{total_length_of_chromosome_level_assembly}',
        f'Number of unanchored contigs:\t{number_of_unanchored_contigs}',
        f'Length of unanchored contigs (bp):\t{length_of_unanchored_contigs}',
        f'Anchor rate (%):\t' + f'{anchor_rate:.2%}'.replace("%", ""), 
        sep='\n', file=output)
    
    
def pseudo_agp(real_list, contigsizes, gap_size, output):
    """
    Create a pseudo agp file from a list and contigsizes.
    """
    real_list = pd.read_csv(real_list, sep='\t', header=None, index_col=None)
    contig_sizes = dict(i.strip().split()[:2] for i in open(contigsizes) if i.strip())

    idx = 0
    for chrom, contig_df in real_list.groupby(0):
        start = 1 
        end = 1
        idx += 1 
        for j, contig in enumerate(contig_df[1].values.tolist()):
            contig = str(contig)
            start = end
            contig_length = int(contig_sizes[contig])
            end = start + contig_length - 1
            print("\t".join(map(str, [chrom, start, end, idx, "W", contig, 1, contig_length, "+"])), file=output)
            if gap_size > 0:
                if j < len(contig_df) - 1:
                    start = end
                    end = start + gap_size
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
        scaffolds_idx = sorted(scaffolds_length_db, key=lambda x: scaffolds_length_db[x], reverse=True)
    else:
        scaffolds_idx = list(scaffolds_length_db.keys())
    
    raw_agp = f"{outprefix}.agp"

    agp = f"{outprefix}.corrected.agp"

    with open(raw_agp, "w") as raw_out, open(agp, "w") as out:
        record_num = 1
        for idx, group_idx in enumerate(scaffolds_idx):
            records = agp_records[group_idx]
            raw_records = raw_agp_records[group_idx]
            unanchor_records = []
            
            if chrom_num and isinstance(chrom_num, list):
                total_chrom_num = chrom_num[0] if len(chrom_num) == 1 else chrom_num[0] * chrom_num[1]

                for i, (record, raw_record) in enumerate(zip(records, raw_records)):
                    num = idx + 1
                    
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
                        # print(raw_contig, raw_contig_start, raw_contig_end)
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

def assembly2agp(assembly, 
                 chrom_num=None, sort_by_length=False,
                 gap_len=100,
                 phased=False, chrom_prefix="Chr",
                 outprefix=None,
                 original_agp=None):
    
    logger.info(f"Converting assembly to agp ...")
    if not outprefix:
        outprefix = Path(assembly).absolute().stem 

    orig_coords = {}
    no_gap_transitions = set()
    
    clean_ctg = lambda x: re.sub(r'_d\d+$', '', str(x))

    if original_agp:
        orig_df = import_agp(original_agp, split=False)
        for chrom, sub_df in orig_df.groupby(0):
            sub_df = sub_df.reset_index(drop=True)
            for idx in range(len(sub_df) - 1):
                row1 = sub_df.iloc[idx]
                row2 = sub_df.iloc[idx + 1]
                if row1[4] in ['W', 'D', 'F'] and row2[4] in ['W', 'D', 'F']:
                    ctg1, ori1 = clean_ctg(row1[5]), str(row1[8])
                    ctg2, ori2 = clean_ctg(row2[5]), str(row2[8])
                    no_gap_transitions.add(((ctg1, ori1), (ctg2, ori2)))
                    no_gap_transitions.add(((ctg2, '-' if ori2 == '+' else '+'), (ctg1, '-' if ori1 == '+' else '+')))
        
        orig_w_df = orig_df[orig_df[4].isin(['W', 'D', 'F'])]
        for _, row in orig_w_df.iterrows():
            comp_id = str(row[5])
            orig_coords[comp_id] = (int(row[6]), int(row[7]))

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

    scaffolds_length_db = OrderedDict()
    raw_agp_records = OrderedDict()
    agp_records = OrderedDict()
    scaffold_gaps = OrderedDict()

    for i, group in enumerate(scaffolds):
        raw_record = []
        record = []
        record_num = 0
        start = 1 
        gaps = []

        parsed_group = []
        for idx in group:
            if idx.startswith("-"):
                orient = "-"
                idx_val = int(idx[1:])
            else:
                orient = "+"
                idx_val = int(idx)
            parsed_group.append((contig_idx_db[idx_val], orient, idx_val))

        for g_idx, (contig, orient, idx_val) in enumerate(parsed_group):
            contig_length = contig_idx_length_db[idx_val]
            orig_start, orig_end = 1, contig_length
            
            if original_agp:
                clean_contig = clean_ctg(contig)
                if contig in orig_coords:
                    orig_start, orig_end = orig_coords[contig]
                    contig_length = orig_end - orig_start + 1
                elif clean_contig in orig_coords:
                    orig_start, orig_end = orig_coords[clean_contig]
                    contig_length = orig_end - orig_start + 1

            if i not in scaffolds_length_db:
                scaffolds_length_db[i] = 0

            scaffolds_length_db[i] += contig_length
    
            if contig in fragments_range:
                raw_contig, raw_contig_start, raw_contig_end = fragments_range[contig]
            else:
                raw_contig, raw_contig_start, raw_contig_end = contig, orig_start, orig_end
        
            record.append((i, start, start + contig_length, record_num,
                            "W", contig, orig_start, orig_end, orient))
            
            raw_record.append((
                i, start, start + contig_length - 1, record_num,
                "W", raw_contig, raw_contig_start, raw_contig_end, orient,
            ))
            
            current_gap_len = gap_len
            if g_idx < len(parsed_group) - 1:
                next_contig, next_orient, _ = parsed_group[g_idx + 1]
                if ((clean_ctg(contig), orient), (clean_ctg(next_contig), next_orient)) in no_gap_transitions:
                    current_gap_len = 0 

                gaps.append(current_gap_len)

            record_num += 1 
            start += contig_length + current_gap_len
        
        agp_records[i] = record 
        raw_agp_records[i] = raw_record
        scaffold_gaps[i] = gaps
        
        
    if sort_by_length:
        scaffolds_idx = sorted(scaffolds_length_db, key=lambda x: scaffolds_length_db[x], reverse=True)
    else:
        scaffolds_idx = list(scaffolds_length_db.keys())
    
    raw_agp = f"{outprefix}.agp"
    agp = f"{outprefix}.corrected.agp"
    has_corrected = len(fragment_db) > 0
    with open(raw_agp, "w") as raw_out:
        out = open(agp, "w") if has_corrected else io.StringIO()
        try:
            record_num = 1
            for idx, group_idx in enumerate(scaffolds_idx):
                records = agp_records[group_idx]
                raw_records = raw_agp_records[group_idx]
                unanchor_records = []
                
                if chrom_num and isinstance(chrom_num, list):
                    total_chrom_num = chrom_num[0] if len(chrom_num) == 1 else chrom_num[0] * chrom_num[1]

                    for i, (record, raw_record) in enumerate(zip(records, raw_records)):
                        num = idx + 1
                        
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
                                current_gap = scaffold_gaps[group_idx][i]
                                if current_gap > 0:
                                    record_num += 1
                                    print("\t".join(map(str, (chrom, end, end + current_gap - 1, record_num, 
                                                                "U", current_gap, "contig", "yes", "map" ))), file=out)
                                    print("\t".join(map(str, (chrom, raw_end + 1, raw_end + current_gap, record_num, 
                                                                "U", current_gap, "contig", "yes", "map" ))), file=raw_out)
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
                            # print(raw_contig, raw_contig_start, raw_contig_end)
                        print("\t".join(map(str, (chrom, start, end - 1, record_num, _type, 
                                                    contig, contig_start, contig_end, orient))), file=out)
                        
                        _, raw_start, raw_end, _, _type, raw_contig, contig_start, contig_end, orient = raw_record

                        print("\t".join(map(str, (chrom, raw_start, raw_end, record_num, _type, 
                                                    raw_contig, contig_start, contig_end, orient))), file=raw_out)
                        
                        if i < len(records) - 1:
                            current_gap = scaffold_gaps[group_idx][i]
                            if current_gap > 0:
                                record_num += 1
                                print("\t".join(map(str, (chrom, end, end + current_gap - 1, record_num, 
                                                            "U", current_gap, "contig", "yes", "map" ))), file=out)
                                print("\t".join(map(str, (chrom, raw_end + 1, raw_end + current_gap, record_num, 
                                                            "U", current_gap, "contig", "yes", "map" ))), file=raw_out)
                        record_num += 1
        finally:
            out.close()

        if has_corrected:
            logger.info(f"Successful output agp file `{agp}`.")
        logger.info(f"Successful output agp file `{raw_agp}`.")


def agp_dup(agp, output):
    """
    rename duplicated contigs to other name in agp file
    """

    # collapsed_contigs = {}
    # with open(collapsed_list) as fp:
    #     for line in fp:
    #         if not line.strip():
    #             continue
    #         line_list = line.strip().split()
    #         if len(line_list) < 2:
    #             continue

    #         contig1, contig2 = line_list 

    #         if contig1 not in collapsed_contigs:
    #             collapsed_contigs[contig1] = [contig1]
            
    #         collapsed_contigs[contig1].append(contig2)

    

    agp_df = import_agp(agp, split=False)
    agp_df.reset_index(inplace=True)
    
    duplicated_df = agp_df[agp_df[4] == 'W'][agp_df.duplicated(5)]

    duplicated_contigs = defaultdict(lambda :1)
    for idx, row in duplicated_df.iterrows():
        if row[4] == 'U':
            continue
        duplicated_contigs[row[5]] += 1
        suffix = duplicated_contigs[row[5]]
        if row[0] == row[5]:
            agp_df.loc[idx, 0] = f"{row[5]}_d{suffix}"

        agp_df.loc[idx, 5] = f"{row[5]}_d{suffix}"

    
    agp_df.drop('index', axis=1, inplace=True)

    agp_df.to_csv(output, sep='\t', header=None, index=None)


def split_agp(agp, hap_pattern=r"(Chr\d+)(?:g\d+)?",
              output="agp_split"):
    """
    Split agp file into different chro
    Params:
    --------
    agp: str
        Path to agp file.
    hap_pattern: str
        Regex pattern to identify haplotypes.

    outdir: str

        Output directory.

    Returns:
    --------
    None

    Examples:
    --------
    >>> split_agp('groups.agp', hap_pattern=r"(Chr\d+)g(\d+)")
    """
    agp_df = import_agp(agp, split=False)

    agp_df.reset_index(inplace=True)


    agp_df['hap'] = agp_df[0].apply(
        lambda x: re.match(hap_pattern, x).group(1) if re.match(hap_pattern, x) else "unhap"
    )
    agp_df = agp_df[agp_df['hap'] != "unhap"]
    # agp_df.drop(0, axis=1, inplace=True)
    hap_df = agp_df.groupby('hap')
    
 
    outdir = Path(output)
    outdir.mkdir(parents=True, exist_ok=True)
    for hap, df in hap_df:
        out_agp = outdir / f"{hap}.agp"
        # df.rename(columns={'hap': 0}, inplace=True)
        df = df[[0,1,2,3,4,5,6,7,8]]

        df.to_csv(out_agp, sep='\t', header=None, index=False)
        logger.info(f"Output split agp file: `{out_agp}`")

    