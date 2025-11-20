#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Identify the high confidence regions of contigs by Pore-C/Hi-C contacts
"""

import argparse
import logging
import os
import os.path as op
import sys

import cooler 
import numpy as np 
import pandas as pd
import polars as pl
import pyranges as pr 

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.ticker import MaxNLocator
from pathlib import Path
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter1d

# from line_profiler import profile


logger = logging.getLogger(__name__)

# def plot(data, lower_value=0.1, upper_value=1.75, output="output"):
#     fig, ax = plt.subplots(figsize=(5.7, 5))
#     plt.rcParams['font.family'] = 'Arial'
#     data = data[data <= np.percentile(data, 98) * 1.5]
#     ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
#     ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
#     kdelines = sns.kdeplot(data, ax=ax, color='#253761', linewidth=2)
#     plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#     formatter = plt.gca().get_yaxis().get_major_formatter()
#     plt.gca().yaxis.set_major_formatter(formatter)
#     plt.gca().yaxis.get_offset_text().set_fontsize(14)
#     plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#     formatter = plt.gca().get_xaxis().get_major_formatter()
#     plt.gca().xaxis.set_major_formatter(formatter)
#     plt.gca().xaxis.get_offset_text().set_fontsize(14)
#     plt.xlim(0, np.percentile(data, 95) * 2)
#     plt.xticks(fontsize=18)
#     plt.yticks(fontsize=18)
#     plt.xlabel("Contacts", fontsize=24)
#     plt.ylabel("Density", fontsize=24)

#     x = kdelines.lines[0].get_xdata()
#     y = kdelines.lines[0].get_ydata()


#     peak_ind = find_peaks(y, distance=10)[0]

#     median_value = np.quantile(data, .3)
#     peak_ind = list(filter(lambda j: x[j] > median_value, peak_ind))
#     if len(peak_ind) == 0:
#         max_idx = np.argsort(x)[len(x)//2]
#     else:
#         max_idx = peak_ind[np.argmax(y[peak_ind])]
    
  
    
#     ax.fill_between((x[max_idx] * lower_value, x[max_idx] * upper_value), 
#                     0, ax.get_ylim()[1], alpha=0.5 , color='#bcbcbc')
#     ax.axvline(x[max_idx] * lower_value, linestyle='--', color='k')
#     ax.axvline(x[max_idx] * upper_value, linestyle='--', color='k')
#     ax.axvline(x[max_idx], linestyle='--', color='#cb6e7f')
#     ax.text(int(x[max_idx]), ax.get_ylim()[1] / 4, str(int(x[max_idx])), fontsize=10, color='#cb6e7f')    

#     # plt.plot(x[max_idx], y[max_idx], ms=10, color='r')
#     plt.savefig(f'{output}.kde.plot.png', dpi=600, bbox_inches='tight')
#     plt.savefig(f'{output}.kde.plot.pdf', dpi=600, bbox_inches='tight')
#     logger.info(f"Output kde plot of contacts distribution in `{output}.kde.plot.png`")
    
#     return int(x[max_idx]), x[max_idx] * lower_value, x[max_idx] * upper_value


def plot(data, lower_value=0.1, upper_value=1.75, output="output"):
    fig, ax = plt.subplots(figsize=(5.7, 5))
    plt.rcParams['font.family'] = 'Arial'
    data = data[data <= np.percentile(data, 98) * 1.5]
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    kdelines = sns.kdeplot(data, ax=ax, color='#253761', linewidth=2)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    formatter = plt.gca().get_yaxis().get_major_formatter()
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.gca().yaxis.get_offset_text().set_fontsize(14)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    formatter = plt.gca().get_xaxis().get_major_formatter()
    plt.gca().xaxis.set_major_formatter(formatter)
    plt.gca().xaxis.get_offset_text().set_fontsize(14)
    plt.xlim(0, np.percentile(data, 95) * 2)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("Contacts", fontsize=24)
    plt.ylabel("Density", fontsize=24)


    try:
        x_min = 0.0
        x_max = float(np.percentile(data, 99.5) * 1.5)
        x_grid = np.linspace(x_min, max(x_max, 1e-9), 2048)
        kde = gaussian_kde(data, bw_method='scott')
        y_kde = kde.evaluate(x_grid)
    except Exception:
        bins = min(256, max(32, int(np.sqrt(len(data)))))
        counts, edges = np.histogram(data, bins=bins, range=(0, x_max if 'x_max' in locals() else None), density=True)
        x_grid = (edges[:-1] + edges[1:]) / 2.0
        y_kde = counts


    try:
        y_smooth = gaussian_filter1d(y_kde, sigma=2)
    except Exception:
        y_smooth = y_kde


    q10 = np.quantile(data, 0.10)
    q99 = np.quantile(data, 0.99)
    mask = (x_grid >= max(q10 * 0.5, 0.0)) & (x_grid <= q99 * 1.1)
    x = x_grid[mask]
    y = y_smooth[mask]


    if len(x) == 0 or len(y) == 0:
        max_idx = 0
    else:
        height_thresh = 0.05 * float(y.max())  
        prom_thresh = 0.10 * float(y.max())    
        distance = max(5, len(y) // 200)  
        peaks, props = find_peaks(y, height=height_thresh, prominence=prom_thresh, distance=distance)

        if len(peaks) == 0:

            max_idx = int(np.argmax(y))
        else:
            prominences = props.get('prominences', np.zeros_like(peaks, dtype=float))
            max_idx = int(peaks[int(np.argmax(prominences))])

    x_peak = float(x[max_idx])

    ax.fill_between((x_peak * lower_value, x_peak * upper_value), 
                    0, ax.get_ylim()[1], alpha=0.5 , color='#bcbcbc')
    ax.axvline(x_peak * lower_value, linestyle='--', color='k')
    ax.axvline(x_peak * upper_value, linestyle='--', color='k')
    ax.axvline(x_peak, linestyle='--', color='#cb6e7f')
    ax.text(int(x_peak), ax.get_ylim()[1] / 4, str(int(x_peak)), fontsize=10, color='#cb6e7f')    

    plt.savefig(f'{output}.kde.plot.png', dpi=600, bbox_inches='tight')
    plt.savefig(f'{output}.kde.plot.pdf', dpi=600, bbox_inches='tight')
    logger.info(f"Output kde plot of contacts distribution in `{output}.kde.plot.png`")
    
    return int(x_peak), x_peak * lower_value, x_peak * upper_value


def hcr_by_contacts_cool(cool_file, output, lower_value=0.1, upper_value=1.75,
                    min_remove_whole_collapsed_contigs_rate=0.9):

    cool = cooler.Cooler(cool_file)
    binsize = cool.binsize
    bins = cool.bins()[:]
    
    contigsizes = cool.chromsizes
    matrix = cool.matrix(balance=False, sparse=True)[:]
    sum_values = np.array(matrix.sum(axis=1).T[0])
    logger.debug("Adjusting small bins value ...")
    small_bins = bins[bins['end'] - bins['start'] < cool.binsize]
    small_bins_sum_values = sum_values.T[small_bins.index]
    adjust_small_bins_sum_values = small_bins_sum_values.T / \
        ((small_bins['end'] - small_bins['start']) / binsize).values
    sum_values[:, small_bins.index] = adjust_small_bins_sum_values

    sum_values_nonzero = sum_values[sum_values > 0]

    min_value, max_value = plot(sum_values_nonzero, lower_value, 
                                upper_value, output=output.replace(".bed", ""))
    #median = np.median(sum_values)
    # max_value = np.percentile(sum_values_nonzero, percent)
    # logger.debug(f"Percent{percent} value is {max_value}")
    res = np.where((sum_values <= max_value) & (sum_values >= min_value))

    hcr_regions = bins.loc[res[1]]
    
    res = np.where((sum_values > max_value))
    collapsed_regions = bins.loc[res[1]]
    collapsed_regions = collapsed_regions[['chrom', 'start', 'end']]
    collapsed_regions.columns = ['Chromosome', 'Start', 'End']
    collapsed_regions['Chromosome'] = collapsed_regions['Chromosome'].astype(str)

    total_length = contigsizes.sum()
    num_hcr_regions = len(hcr_regions)
    logger.debug(f"Identify {num_hcr_regions} regions")
    hcr_length = sum(hcr_regions["end"] - hcr_regions["start"])
    logger.info(f"Identified {hcr_length/total_length:.2%} high-confidence regions")
    
    hcr_regions = hcr_regions[['chrom', 'start', 'end']]
    hcr_regions.columns = ['Chromosome', 'Start', 'End']
    hcr_regions_pr = pr.PyRanges(hcr_regions)
    hcr_regions_pr = hcr_regions_pr.merge()

    if min_remove_whole_collapsed_contigs_rate:
        contigsizes_db = contigsizes.to_dict()
        
        collapsed_regions = collapsed_regions.eval('Length = End - Start')
        
        collapsed_regions['Total_length'] = collapsed_regions['Chromosome'].map(contigsizes_db.get)
        collapsed = collapsed_regions.groupby('Chromosome')['Length'].sum().reset_index()
        collapsed['Total_length'] = collapsed['Chromosome'].map(contigsizes_db.get)
        
        collapsed = collapsed.eval('Rate = 1 - (Total_length - Length) / Total_length')
        collapsed = collapsed[collapsed['Rate'] > min_remove_whole_collapsed_contigs_rate]
        collapsed_df = contigsizes.loc[collapsed.Chromosome].reset_index()
        collapsed_df[1] = 0
        collapsed_df.columns = ['Chromosome', 'End', 'Start']
        collapsed_df = collapsed_df[['Chromosome', 'Start', 'End']]
        collapsed_pr = pr.PyRanges(collapsed_df)
        hcr_regions_pr = hcr_regions_pr.intersect(collapsed_pr, invert=True)


    hcr_regions_pr.df.to_csv(output, sep='\t', index=None, header=None)
    logger.info(f"Successful output HCRs into `{output}`.")

def hcr_by_contacts(depth_file, output, lower_value=0.1, upper_value=1.75,
                    min_remove_whole_collapsed_contigs_rate=0.9,
                    edge_length=None):

    
    depth = pl.read_csv(depth_file, separator='\t', has_header=False,
                        new_columns=['chrom', 'start', 'end', 'count']).to_pandas()

    binsize = (depth['end'] - depth['start']).max()
    bins = depth[['chrom', 'start', 'end']]
    contigsizes = depth[['chrom', 'end']].groupby('chrom')['end'].max()


  
    logger.debug("Adjusting small bins value ...")
    sum_values = depth['count'].values

    small_bins = depth[depth['end'] - depth['start'] < binsize]
    small_bins_sum_values = sum_values[small_bins.index]

    adjust_small_bins_sum_values = small_bins_sum_values.T / \
        ((small_bins['end'] - small_bins['start']) / binsize).values
    sum_values[small_bins.index] = adjust_small_bins_sum_values
  
    sum_values_nonzero = sum_values[sum_values > 0]

    peak_value, min_value, max_value = plot(sum_values_nonzero, lower_value, 
                                upper_value, output=output.replace(".bed", ""))
   
    res = np.where((sum_values <= max_value) & (sum_values >= min_value))
    
    bins['count'] = sum_values
    contig_counts = bins.groupby('chrom')['count'].mean().to_frame()
    contig_counts['CN'] = contig_counts['count'] / peak_value
    collapsed_contigs = contig_counts.query('count > @max_value')
    collapsed_contigs.to_csv(output.replace(".bed", ".high_coverage.contigs.txt"), sep='\t',
                             index=True, header=None)
    contigsizes.loc[contig_counts.query('count > @max_value').index].sum()

    hcr_regions = bins.loc[res[0]]
    
    res = np.where((sum_values > max_value))

    collapsed_regions = bins.loc[res[0]]
    collapsed_regions = collapsed_regions[['chrom', 'start', 'end']]
    collapsed_regions.columns = ['Chromosome', 'Start', 'End']
    collapsed_regions['Chromosome'] = collapsed_regions['Chromosome'].astype(str)

    total_length = contigsizes.sum()
    
    num_hcr_regions = len(hcr_regions)
    logger.debug(f"Identify {num_hcr_regions} regions")
    hcr_length = sum(hcr_regions["end"] - hcr_regions["start"])
    logger.info(f"Identified {hcr_length/total_length:.2%} high-confidence regions")
    
    hcr_regions['start'] = hcr_regions['start'].astype(np.int64)
    hcr_regions['end'] = hcr_regions['end'].astype(np.int64)
    hcr_regions.to_csv(f"tmp.{output}", sep='\t', index=None, header=None)

    cmd = f"bedtools merge -i tmp.{output} 2>/dev/null > tmp.merge.{output}"
    flag = os.system(cmd)
    assert flag == 0, f"Failed to merge HCRs by `{cmd}`"

    hcr_regions = pd.read_csv(f"tmp.merge.{output}", sep='\t', header=None)
    if Path(f"tmp.{output}").exists():
        os.remove(f"tmp.{output}")
    if Path(f"tmp.merge.{output}").exists():
        os.remove(f"tmp.merge.{output}")

    hcr_regions_df = hcr_regions
    # hcr_regions = hcr_regions[['chrom', 'start', 'end']]
    # hcr_regions.columns = ['Chromosome', 'Start', 'End']
    # hcr_regions_pr = pr.PyRanges(hcr_regions)
    # hcr_regions_pr = hcr_regions_pr.merge()

    if min_remove_whole_collapsed_contigs_rate:
        hcr_regions.columns = ['Chromosome', 'Start', 'End']
        hcr_regions_pr = pr.PyRanges(hcr_regions)
        contigsizes_db = contigsizes.to_dict()
        
        collapsed_regions = collapsed_regions.eval('Length = End - Start')
        
        collapsed_regions['Total_length'] = collapsed_regions['Chromosome'].map(contigsizes_db.get)
        collapsed = collapsed_regions.groupby('Chromosome')['Length'].sum().reset_index()
        collapsed['Total_length'] = collapsed['Chromosome'].map(contigsizes_db.get)
        
        collapsed = collapsed.eval('Rate = 1 - (Total_length - Length) / Total_length')
        collapsed = collapsed[collapsed['Rate'] > min_remove_whole_collapsed_contigs_rate]
        collapsed_df = contigsizes.loc[collapsed.Chromosome].reset_index()
        collapsed_df[1] = 0
        collapsed_df.columns = ['Chromosome', 'End', 'Start']
        collapsed_df = collapsed_df[['Chromosome', 'Start', 'End']]
        collapsed_pr = pr.PyRanges(collapsed_df)
        hcr_regions_pr = hcr_regions_pr.intersect(collapsed_pr, invert=True)
        hcr_regions_df = hcr_regions_pr.df

    if edge_length:
        contigsizes = contigsizes.rename("length").to_frame()
        large_contigs = contigsizes[contigsizes['length'] >= edge_length * 2]
        small_contigs = contigsizes[contigsizes['length'] < edge_length * 2]
        def func1(row):
            row['start'] = 0
            row['end'] = int(2e6)

            return row
    
        def func2(row2):
            row2['start'] = row2['length'] - 2e6
            row2['end'] = row2['length']

            return row2
        
        res_df = pd.concat([large_contigs.reset_index().apply(func1, axis=1),
              large_contigs.reset_index().apply(func2, axis=1),
              ], axis=0).drop(['length'], axis=1).astype({"start": np.int64, "end": np.int64})
        small_contigs['start'] = 0
        small_contigs = small_contigs.reset_index().rename(columns={"length": "end"})
        
        res_df = pd.concat([res_df, small_contigs], axis=0)
        res_df.columns = ["Chromosome", "Start", "End"]

        edge_gr = pr.PyRanges(res_df)
        report_regions_pr = edge_gr.join(hcr_regions_pr).new_position("intersection")
        report_regions_pr.df[['Chromosome', 'Start_b', 'End_b']].to_csv(
            output, sep='\t', index=None, header=None
        )
    
    else:
        hcr_regions_df.to_csv(output, sep='\t', index=None, header=None)
    logger.info(f"Successful output HCRs into `{output}`.")

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('cool_file', 
            help='')
   
    pOpt.add_argument('-o', '--output', type=str,
            default="output", help='output file [default: output]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    hcr_by_contacts(args.cool_file, args.output)

if __name__ == "__main__":
    main(sys.argv[1:])