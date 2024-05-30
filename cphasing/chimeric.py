#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import logging
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

from joblib import Parallel, delayed
from pathlib import Path
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
    contigsizes_bed = pd.DataFrame(contigsizes, index=['length']).T.reset_index()
    contigsizes_bed['start'] = 0
    contigsizes_bed = contigsizes_bed[['index', 'start', 'length']]
    contigsizes_bed.columns = ['chrom', 'start', 'end']
    
    os.environ["POLARS_MAX_THREADS"] = str(threads)
    try:
        p = (pl.read_csv(pairs, separator='\t', has_header=False,
                            comment_prefix="#", columns=[1, 2, 3, 4, 7],
                            new_columns=['chrom1', 'pos1', 'chrom2', 'pos2', 'mapq'],
                            dtypes=dtype
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

def calculate_depth(p, contigsizes, window_size=500, threads=4):

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
            break_res.append((contig, break_point * (window_size) - window_size))
    
    break_points_df =  pd.DataFrame(break_res)
    try:
        contig_num = len(set(break_points_df[0].values.tolist()))
        break_points_df = break_points_df.sort_values(by=0)
    except KeyError:
        contig_num = 0

    logger.info(f"Identified `{len(break_points_df)}` "
                    f"break points in `{contig_num}` contigs.")

    return break_points_df



def find_nearest(array, value):
    """
    https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()

    return idx, array[idx]

def find_break_point(contig, b, min_windows=50):
    a = np.arange(len(b))

    rmse = mean_squared_error(b, pd.DataFrame(b).shift(1).fillna(0))
    if rmse < 1.0 :
        return None
    
    try:
        z1 = np.polyfit(a, b, 50)
        p1 = np.poly1d(z1)
    except:
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

            # print(new_trough, trough_y, peak1_y / 1.5, peak2_y / 1.5, peak1, peak2)
            if (trough_y < (peak1_y / 1.5)) and (trough_y < (peak2_y / 1.5)):
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

        if left_rmse < 1.0 or right_rmse < 1.0:
            continue 
        new_new_trough_ind.append(new_trough)

    if new_new_trough_ind:
        return contig, new_new_trough_ind
    else:
        return

def break_points_to_regions(break_points_df, contigsizes):
    correct_contig_list = []
    for i, tmp_df in break_points_df.groupby(0):
        length = contigsizes[i]
        pos = tmp_df[1].values
        start = np.r_[[1], pos + 1]
        end = np.r_[pos, [length]]
        pos_list = list(zip(start, end))

        for start, end in pos_list:
            correct_contig_list.append((i, start, end, f"{i}:{start}-{end}"))

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


def run(fasta, pairs, window_size=500, break_pairs=True, 
        outprefix="output", threads=4):
    logger.info("Running chimeric correction ...")
    pairs_df, contigsizes = import_pairs(pairs, window_size=window_size, threads=threads)
    depth_dict = calculate_depth(pairs_df, contigsizes, window_size=window_size, 
                                        threads=threads)
    break_point_res = correct(depth_dict, window_size=window_size, threads=threads)

    if len(break_point_res) > 0:
        break_point_res.to_csv(f"{outprefix}.breakPos.txt", 
                                        sep='\t', index=None, header=None)
        corrected_positions = break_points_to_regions(break_point_res, contigsizes)
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