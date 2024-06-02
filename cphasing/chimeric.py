#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import logging
import re
import os
import os.path as op
import sys
import gc
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pandas as pd
import polars as pl 
import pyranges as pr  

from Bio.Seq import Seq
from collections import defaultdict
from joblib import Parallel, delayed
from pathlib import Path
from pyfaidx import Fasta
from scipy.signal import find_peaks, peak_widths
from sklearn.metrics import mean_squared_error

try: 
    from .core import PairHeader
    from .utilities import read_fasta, xopen, run_cmd
except ImportError:
    pass

logger = logging.getLogger(__name__)

def _calculate_depth(contig, contig_length, window_size, position_data):

    data = np.zeros(contig_length // window_size, dtype=np.int32)
    for pos1, pos2 in position_data:
        data[pos1: pos2 + 1] += 1

    return contig, data 

def import_pairs(pairs, window_size=500, min_mapq=0, threads=4):
    dtype = {'chrom1': pl.Categorical, 
            'pos1': pl.UInt32,
            'chrom2': pl.Categorical, 
            'pos2': pl.UInt32,
            'mapq': pl.Int8}
    
    ph = PairHeader([])
    ph.from_file(pairs)
    contigsizes = ph.chromsize

    new_columns = ['chrom1', 'pos1', 'chrom2', 'pos2', 'mapq'] if min_mapq > 0 else ['chrom1', 'pos1', 'chrom2', 'pos2']
    
    os.environ["POLARS_MAX_THREADS"] = str(threads)
    try:
        p = (pl.read_csv(pairs, separator='\t', has_header=False,
                            comment_prefix="#", columns=[1, 2, 3, 4, 7],
                            new_columns=new_columns,
                            dtypes=dtype, low_memory=True
                            ).filter(pl.col('chrom1') == pl.col('chrom2'))
                            .drop('chrom2')
                        )
        
        if min_mapq > 0:
            p = p.filter(pl.col('mapq') >= min_mapq).drop('map1')

    except pl.exceptions.OutOfBoundsError:
        p = (pl.read_csv(pairs, separator='\t', has_header=False,
                            comment_prefix="#", columns=[1, 2, 3, 4],
                            new_columns=['chrom1', 'pos1', 'chrom2', 'pos2'],
                            dtypes=dtype
                            ).filter(pl.col('chrom1') == pl.col('chrom2'))
                                .drop('chrom2')
                                )

    p = p.with_columns(
        (pl.when(pl.col('pos1') > pl.col('pos2'))
            .then(pl.col('pos2'))
            .otherwise(pl.col('pos1'))).alias('pos1'),
        (pl.when(pl.col('pos1') > pl.col('pos2'))
            .then(pl.col('pos1'))
            .otherwise(pl.col('pos2'))).alias('pos2')
    )

    p = p.with_columns(
        (pl.col('pos1') - 1) // window_size,
        (pl.col('pos2') - 1) // window_size,
    )
    
    p = p.to_pandas()

    return p, contigsizes

def calculate_depth(p, contigsizes, window_size=500, output_depth=False, threads=4):

    args = []
    for contig, tmp_df in p.groupby('chrom1', as_index=True):
        args.append((contig, contigsizes[contig], window_size, tmp_df[['pos1', 'pos2']].values))
    
    del p
    gc.collect()

    res = Parallel(n_jobs=threads)(
        delayed(_calculate_depth)(i, j, k, l)
            for i, j, k, l in args
    )

    res_db = dict(res)
    

    return res_db
    

def correct(depth_dict, window_size=500, min_windows=50, threads=4):
    args = []
    for contig in depth_dict:
        b = depth_dict[contig]
        if len(b) <= min_windows * 2:
            continue
        if b.mean() < 10:
            continue
        args.append((contig, b, min_windows))
    
    res = Parallel(n_jobs=threads)(
        delayed(find_break_point)(i, j, k)
            for i, j, k in args
    )
    
    res = list(filter(lambda x: x is not None, res))

    break_res = []
    for contig, break_points in res:
        for break_point in break_points:
            break_res.append((contig, break_point * (window_size) + window_size))
    
    break_points_df =  pd.DataFrame(break_res)
    try:
        contig_num = len(set(break_points_df[0].values.tolist()))
        break_points_df = break_points_df.sort_values(by=0)
    except KeyError:
        contig_num = 0

    logger.info(f"Identified `{len(break_points_df)}` "
                    f"break points in `{contig_num}` contigs.")

    return break_points_df

def edge_index(pos, length):
    """
    ------|-^----
    """
    mid_point = length // 2 
    if (length - pos) < mid_point:
        ## right 
        edge_index_value = (length - pos) / (length - mid_point)
    else:
        ## left 
        edge_index_value =  (mid_point - pos) / (mid_point - 0)
    
    return edge_index_value

def find_nearest(array, value):
    """
    https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()

    return idx, array[idx]

def find_break_point(contig, b, min_windows=50):
    b_len = len(b)
    a = np.arange(b_len)
    rmse = mean_squared_error(b, pd.DataFrame(b).shift(1).fillna(0))
    if rmse < 1.0 :
        return None
    
    try:
        z1 = np.polyfit(a, b, 50)
        p1 = np.poly1d(z1)
    except SystemError:
        return 
    yvalues = p1(a)

    peak_ind = find_peaks(yvalues, distance=10)[0]

    peak_ind = list(filter(lambda x: x >= min_windows and x <= (len(b)  - min_windows), list(peak_ind)))

    if len(peak_ind) <= 1:
        return None

    
    trough_ind = find_peaks(-yvalues, distance=10)[0]
    peak_widths_result = peak_widths(yvalues, peak_ind, rel_height=0.4)
    trough_widths_result = list(peak_widths(-yvalues, trough_ind, rel_height=.5))
    trough_widths_result[1] = -trough_widths_result[1] 
    trough_widths_x = list(zip(*trough_widths_result[2:]))
    trough_widths_y = list(trough_widths_result[1])
    peak_widths_x = list(zip(*peak_widths_result[2:]))   
    peak_widths_y = list(peak_widths_result[1])

    new_trough_ind = []
    suspicious_trough = []
    for idx, trough in enumerate(trough_ind):
        nearest_peak_idx, nearest_peak = find_nearest(peak_ind, trough)

        if (trough < min_windows) or ((len(b) - trough) < min_windows):
            continue

        if nearest_peak == 0 or nearest_peak == len(peak_ind) - 1:
            continue

        if trough > nearest_peak:
            try:
                nearest_peak2 = peak_ind[nearest_peak_idx + 1]
                nearest_peak2_idx = nearest_peak_idx + 1
            except IndexError:
                continue 
        else:
            try:
                nearest_peak2 = peak_ind[nearest_peak_idx - 1]
                nearest_peak2_idx = nearest_peak_idx - 1

                nearest_peak, nearest_peak2 = nearest_peak2, nearest_peak
                nearest_peak_idx, nearest_peak2_idx = nearest_peak2_idx, nearest_peak_idx
            except IndexError:
                continue 
        
        if nearest_peak_idx < 0 or nearest_peak2_idx < 0:
            continue
        
        # print(nearest_peak2_idx, nearest_peak_idx)
        if trough_widths_y[idx] < peak_widths_y[nearest_peak_idx] \
            and trough_widths_y[idx] < peak_widths_y[nearest_peak2_idx]:

            peak_width_x1 = list(map(int, peak_widths_x[nearest_peak_idx]))
            peak_width_x2 = list(map(int, peak_widths_x[nearest_peak2_idx]))
            if ((peak_width_x1[1] - peak_width_x1[0]) < 50) or ((peak_width_x2[1] - peak_width_x2[0]) < 50):
                continue
            peak1 = b[peak_width_x1[0]:peak_width_x1[1]].argmax() + peak_width_x1[0]
            peak2 = b[peak_width_x2[0]: peak_width_x2[1]].argmax() + peak_width_x2[0]
            peak1_edge = edge_index(peak1, b_len)
            peak2_edge = edge_index(peak2, b_len)
            peak1_y = b[peak_width_x1[0]:peak_width_x1[1]].mean()
            peak2_y = b[peak_width_x2[0]: peak_width_x2[1]].mean()


            if peak1_y < 10 or peak2_y < 10:
                continue
            # peak1_y = b[peak1]
            # peak2_y = b[peak2]
            trough_width_x = list(map(int, trough_widths_x[idx]))
            new_trough = b[trough_width_x[0]: trough_width_x[1]].argmin() + trough_width_x[0] 
            trough_y = b[new_trough]
            
            if trough_y < 10:
                continue

            # print(new_trough, trough_y, peak1_y, peak2_y, peak1, peak2, peak1_edge, peak2_edge)
            if (trough_y < (peak1_y / (1 + .5 * (peak1_edge )))) and (trough_y < (peak2_y / (1 + .5 * (peak2_edge)))):
                new_trough_ind.append(new_trough)
            # elif (trough_y > (peak1_y / 1.2)) and (trough_y > (peak2_y / 1.2)):
            #     suspicious_trough.append((new_trough, peak1, peak2))
    
    # if suspicious_trough:
    #     print(contig, suspicious_trough)
    
    new_new_trough_ind = []
    for new_trough in new_trough_ind:
        left_b = b[: new_trough]
        right_b = b[new_trough:]
        
        left_rmse = mean_squared_error(left_b, pd.DataFrame(left_b).shift(1).fillna(0))
        right_rmse = mean_squared_error(right_b, pd.DataFrame(right_b).shift(1).fillna(0))
        # print(left_rmse, right_rmse)
        if left_rmse < .95 or right_rmse < .95:
            continue 
        new_new_trough_ind.append(new_trough)


    if new_new_trough_ind:
        return contig, new_new_trough_ind
    else:
        return

def is_telomere(seq, motif="CCCTAA", window=2000):

    regex_right = re.compile(str(Seq(motif).complement()[::-1]))
    res_right = regex_right.finditer(seq[-window:])
    regex_left = re.compile(motif)
    res_left = regex_left.finditer(seq[0: window])

    res_left_pos = [x.start(0) for x in res_left]
    res_right_pos = [x.start(0) for x in res_right]
    
    if len(res_right_pos) > 30 or len(res_left_pos) > 30:
        return True
    
    else:
        return False 
    

def split_contigs(previous_break_points, contigsizes, pairs_df, split_num=2):
    split_contigsizes = {}
    split_contig_fragment_sizes = {}
    for contig in contigsizes:
        if len(previous_break_points) > 0:
            if contig in set(previous_break_points[0].values.tolist()): 
                continue 
        split_contig_fragment_size = contigsizes[contig] // split_num 
        split_contig_fragment_sizes[contig] = split_contig_fragment_size
        for i in range(split_num):
            new_contig = f"{contig}_{i}"
            split_contigsizes[new_contig] = split_contig_fragment_size
    
    pairs_df['split_length'] = pairs_df['chrom1'].map(split_contig_fragment_sizes.get)

    pairs_df['suffix_1'] = pairs_df['pos1'] // pairs_df['split_length']
    pairs_df['suffix_2'] = pairs_df['pos2'] // pairs_df['split_length']
    pairs_df = pairs_df.query('suffix_1 == suffix_2')
    pairs_df['pos1'] = (pairs_df['pos1'] - pairs_df['split_length'] * pairs_df['suffix_1']).astype(int)
    pairs_df['pos2'] = (pairs_df['pos2'] - pairs_df['split_length'] * pairs_df['suffix_2']).astype(int)
    pairs_df['chrom1'] = pairs_df['chrom1'].astype(str) + "_" + pairs_df['suffix_1'].astype(int).astype(str)
    pairs_df.drop(['split_length', 'suffix_1', 'suffix_2'], axis=1, inplace=True)
    
    
    return split_contigsizes, pairs_df


def break_points_to_regions(break_points_df, contigsizes):
    correct_contig_list = []
    for i, tmp_df in break_points_df.groupby(0):
        length = contigsizes[i]
        pos = tmp_df[1].values
        start = np.r_[[1], pos + 1]
        end = np.r_[pos, [length]]
        pos_list = list(zip(start, end))

        for start, end in pos_list:
            correct_contig_list.append((i, int(start), int(end), f"{i}:{start}-{end}"))

    correct_contig_list = pd.DataFrame(
        correct_contig_list, columns=['chrom', 'start', 'end', 'name'])

    return correct_contig_list

def write_fasta(fasta, corrected_positions, output):
    fasta_db = read_fasta(str(fasta))
    output = xopen(output, 'w')
    corrected_contigs = set(corrected_positions['chrom'].values.tolist())
    corrected_positions.set_index('chrom', inplace=True)
    for contig in fasta_db:
        seq = fasta_db[contig]
        if contig in corrected_contigs:
            for _, item in corrected_positions.loc[contig].iterrows():
                start, end, name = item
                print(f'>{name}\n{str(seq[start - 1: end])}',
                        file=output)

        else:
            print(f'>{contig}\n{str(seq)}', file=output)        

def corrected_pairs(pairs, break_bed, output):
    cmd = ["cphasing-rs", "pairs-break", pairs, break_bed, "-o", output]
    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    flag = run_cmd(cmd, log=f"logs/pairs-break.log")
    assert flag == 0, "Failed to execute command, please check log."


def run(fasta, pairs, window_size=500, 
        correct_round=1, telo_motif="CCCTAA",
        break_pairs=True, 
        output_depth=False, outprefix="output", threads=4):
    logger.info("Running chimeric correction ...")
    pairs_df, contigsizes = import_pairs(pairs, window_size=window_size, threads=threads)
    depth_dict = calculate_depth(pairs_df, contigsizes, window_size=window_size, 
                                        output_depth=output_depth, threads=threads)
    
    if output_depth:
        with open(f"{outprefix}.cis.w{window_size}.depth", 'w') as out:
            for contig in depth_dict:
                tmp_depths = depth_dict[contig]
                for i, tmp_depth in enumerate(tmp_depths):
                    print(f"{contig}\t{i*window_size}\t{(i+1)*window_size}\t{tmp_depth}", file=out)

    break_point_res = correct(depth_dict, window_size=window_size, threads=threads)
    
    del depth_dict  
    gc.collect()

    i = 1
    split_contigsizes = contigsizes.copy()
    new_break_point_res_df_list = []
    while correct_round - i:
        split_contigsizes, split_pairs_df = split_contigs(
                        break_point_res, split_contigsizes, pairs_df, split_num=2)
        split_depth_dict = calculate_depth(
                                split_pairs_df, split_contigsizes,
                                window_size=window_size, threads=threads)
        new_break_point_res = correct(
                                split_depth_dict, window_size=window_size, threads=threads)

        i += 1

        if len(new_break_point_res) > 0:
            # for i, tmp_df in new_break_point_res.iterrows():
            #     pass 
            new_break_point_res['length'] = new_break_point_res[0].map(split_contigsizes.get)
            new_break_point_res[['contig', 'suffix']] = new_break_point_res[0].str.split("_", expand=True)
            new_break_point_res[1] = new_break_point_res[1] + \
                    new_break_point_res['suffix'].astype(int) * new_break_point_res['length']
            new_break_point_res.drop(['length', 'suffix', 0], axis=1, inplace=True)
            new_break_point_res = new_break_point_res.rename(columns={'contig': 0})
            new_break_point_res_df_list.append(new_break_point_res)

    if len(new_break_point_res_df_list) > 0:
        break_point_res = pd.concat([break_point_res, 
                                     pd.concat(new_break_point_res_df_list, axis=0)], axis=0)

        del split_pairs_df, split_depth_dict
    del pairs_df
    gc.collect()

    if len(break_point_res) > 0:
        fasta_db = Fasta(fasta)
        corrected_positions = break_points_to_regions(break_point_res, contigsizes)
        filtered_pos = defaultdict(list)
        for i, tmp_df in corrected_positions.iterrows():
            contig = tmp_df['chrom']
            start = int(tmp_df['start'])
            end = int(tmp_df['end'])
            seq = str(fasta_db[contig][start: end])
            if is_telomere(seq, telo_motif):
                filtered_pos[contig].append((start, end))
        
        drop_idx = []
        for i, tmp_df in break_point_res.iterrows():
            contig = tmp_df[0]
            pos = tmp_df[1]
            filter_regions = filtered_pos[contig]
            for region in filter_regions:
                start, end = region
                if pos <= end and pos >= start - 1:
                    drop_idx.append(i)
        if len(drop_idx) > 0:
            logger.info(f"Filtered out {len(drop_idx)} break points, "
                            "because they contain telomere tandem repeats.")
        break_point_res.drop(drop_idx, axis=0, inplace=True) 
        corrected_positions = break_points_to_regions(break_point_res, contigsizes)
        break_point_res.to_csv(f"{outprefix}.breakPos.txt", 
                                        sep='\t', index=None, header=None)
        
        output_break_bed = f"{outprefix}.chimeric.contigs.bed"
        corrected_positions.to_csv(output_break_bed, 
                                    header=None, index=None, sep='\t')
        write_fasta(fasta, corrected_positions, f"{outprefix}.corrected.fasta")
        logger.info(f"Output corrected fasta in `{outprefix}.corrected.fasta`")

        if break_pairs:
            output_pairs_prefix = Path(Path(pairs).stem).with_suffix('')
            while output_pairs_prefix.suffix in {'.pairs', '.gz'}:
                output_pairs_prefix = output_pairs_prefix.with_suffix('')
            if_gz = ".gz" if "gz" in pairs else ""
            output_pairs = f"{output_pairs_prefix}.corrected.pairs{if_gz}"
            corrected_pairs(pairs, output_break_bed, output_pairs)
        else:
            output_pairs = pairs 
        return output_break_bed, f"{outprefix}.corrected.fasta", output_pairs
    
    else:
        logger.info(f"There is not chimeric contigs in `{fasta}`")
        return 

        
def main(args):
    from cphasing.core import PairHeader
    from cphasing.utilities import read_fasta, xopen
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('depth', 
            help='')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    df = pd.read_csv(args.depth, sep='\t', index_col=None, header=None)
    z1 = np.polyfit(df.index.values, df[3].values, 50)
    p1 = np.poly1d(z1)
    yvalues = p1(df.index.values)
    peak_ind = find_peaks(yvalues, distance=10)[0]
    trough_ind = find_peaks(-yvalues, distance=10)[0]

    peak_widths_result = peak_widths(yvalues, peak_ind, rel_height=0.2)
    trough_widths_result = peak_widths(-yvalues, trough_ind, rel_height=.5)
    print(trough_widths_result)
    print(trough_ind, peak_ind)
    plt.rcParams['font.family'] = 'Arial'
    plt.subplots(figsize=(5.5, 5))

    plt.plot(df.index.values, df[3].values, color="#209093")
    plt.plot(df.index.values, yvalues, color='#cb6e7f')
    for trough in trough_ind:
        plt.axvline(x=trough, linewidth=1, color='#bcbcbc', linestyle='--')

    plt.hlines(*peak_widths_result[1:], color="C2")
    plt.hlines(-trough_widths_result[1], *trough_widths_result[2:], colors="C3")

    plt.legend()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("Position (index)", fontsize=16)
    plt.ylabel("Depth", fontsize=16)

    plt.savefig('output.png', bbox_inches='tight')


if __name__ == "__main__":
    main(sys.argv[1:])