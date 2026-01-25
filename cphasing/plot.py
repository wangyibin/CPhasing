#!/usr/bin/env python

import logging
import os
import os.path as op
import re
import sys
import time
import tempfile
import warnings
import gc

try:
    from pandas.core.common import SettingWithCopyWarning
    warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
except ImportError:
    pass

import cooler
import h5py 
import numpy as np
import pandas as pd
import polars as pl
import pyranges as pr
import scipy.sparse as sp

from collections import OrderedDict, defaultdict
from math import ceil
from itertools import product
from intervaltree import Interval, IntervalTree
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator, ScalarFormatter
from multiprocessing import Lock, Pool
from joblib import Parallel, delayed
from pandarallel import pandarallel
from pathlib import Path
from scipy.sparse import triu, coo_matrix, csr_matrix

from .agp import import_agp
from .utilities import to_humanized, to_humanized2, chrom_ticks_convert

# from line_profiler import profile 

logger = logging.getLogger(__name__)

HIC_METADATA = {}
try:
    HIC_METADATA['matrix-generated-by'] = np.string_(
        'CPhasing'
    )
    HIC_METADATA['matrix-generated-by-url'] = np.string_(
        'https://github.com/wangyibin/CPhasing'
    )
except AttributeError:
    HIC_METADATA['matrix-generated-by'] = np.bytes_(
        'CPhasing'
    )
    HIC_METADATA['matrix-generated-by-url'] = np.bytes_(
        'https://github.com/wangyibin/CPhasing'
    )

_BALANCE_BLOCKS = None
_BALANCE_BLOCK_RANGES = None

RED12 = ["#FFFFFF",
"#FEE8E8",
"#FED1D1",
"#FDB9B9",
"#FCA2A2",
"#FB8B8B",
"#FB7474",
"#FA5D5D",
"#F94646",
"#F82E2E",
"#F81717",
"#F70000",]


REDD12 = [
"#FFFFFF",
"#F6ECEB",
"#EED9D7",
"#E5C6C3",
"#DCB3AF",
"#D4A09B",
"#CB8E87",
"#C37B73",
"#BA685F",
"#B1554B",
"#A94237",
"#A02F23",
]

## https://github.com/XiaoTaoWang/NeoLoopFinder/blob/master/neoloop/visualize/core.py
RED5 = ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000']
whitered_cmap = LinearSegmentedColormap.from_list('whitered', ['#FFFFFF', '#FF0000'])

lock = Lock()

def half_colormaps(cmap):
    """
    Create a half colormap from the given colormap name.

    Parameters:
    ----------
    cmap_name : str
        The name of the colormap to create a half colormap from.

    Returns:
    -------
    half_cmap : LinearSegmentedColormap
        The half colormap.
    """
    import matplotlib.pyplot as plt
    from matplotlib import colors
    import colormaps as cmaps
    
    if cmap == 'whitered':
        whitered_cmap = LinearSegmentedColormap.from_list('interaction', ['#FFFFFF','#FFDFDF','#FF7575'])
        return whitered_cmap
    elif cmap == 'reds':
        colormap = getattr(cmaps, cmap)
        color_list = colormap.colors
        half_cmap = LinearSegmentedColormap.from_list(
                f"{cmap}_half", ['white', color_list[5]])

    else:
        try:
            """
            https://pratiman-91.github.io/colormaps/
            """
            colormap = getattr(cmaps, cmap)
            color_list = colormap.colors[:128]
            half_cmap = LinearSegmentedColormap.from_list(
                f"{cmap}_half", color_list)
        except:
            try:
                colormap = cmap 
                colormap = getattr(plt.cm, cmap)
                half_cmap = colors.LinearSegmentedColormap.from_list(
                    f"{cmap}_half", colormap(np.linspace(0, 0.5, 128)))
        
            except AttributeError:
                logger.warning(f"Colormap `{cmap}` not found, use `redp1_r` instead.")
                colormap = getattr(cmaps, "redp1_r")
                color_list = colormap.colors[:128]
                
                half_cmap = LinearSegmentedColormap.from_list(
                    f"{cmap}_half", color_list)


    return half_cmap


def format_float_dynamic(x, sig=4, min_dec=0, max_dec=6, sci_small=1e-4, sci_large=1e6):
    try:
        if not np.isfinite(x):
            return str(x)
        ax = abs(float(x))
        if ax != 0 and (ax < sci_small or ax >= sci_large):
            return f"{x:.{sig}g}"  
        if ax == 0:
            return f"{0:.{min_dec}f}"
        digits = int(np.floor(np.log10(ax))) + 1 
        dec = int(np.clip(sig - digits, min_dec, max_dec)) 
        return f"{x:.{dec}f}"
    except Exception:
        return str(x)
 
def get_bins(bin_size, chrom_size, start=0, orientation="+", 
                reverse=False, region=None):
    """
    Split the chromosomes into bins of length by specify bin_size

    Params:
    -------
    bin_size: `int` size of bins
    chrom_sizes: `list-like` list of chromosome size
    
    Returns:
    -------
    return a list of tuples 

    Examples:
    -------
    >>> chrom_sizes
    [('Chr1', 29999), ('Chr2', 23333)]
    >>> get_bins(50000, chrom_size, region='contig2’)
    [('Chr2', 0, 23333)]
    """
    bin_intervals = []
    # chrom_size_dict = dict(chrom_size)
    # if region:
    #     tmp = []
    #     for chroms in chrom_size:
    #         chrom, size = chroms
    #         if chrom == region:
    #             tmp.append(chroms)
    #     chrom_size = tmp 

    if region:
        chrom_size = [(chrom, size) for chrom, size in chrom_size if chrom == region]
    
    for chrom, size in chrom_size:
        length = size - start
        if orientation == "-" and length % bin_size > 0:
                old_start = start 
                start = start + (length % bin_size)
                bin_intervals.append((chrom, old_start, start))

        bin_intervals.extend(
            (chrom, interval, min(size, interval + bin_size))
            for interval in range(start, size, bin_size)
        )

    return bin_intervals if not reverse else bin_intervals[::-1]


def get_chrom_index_data(bin_interval):

    chrom_index_data = OrderedDict()
    previous_chrom = ''
    for i, (chrom, _, _) in enumerate(bin_interval):
        if chrom != previous_chrom:
            if previous_chrom:
                chrom_index_data[previous_chrom].append(i - 1)
            chrom_index_data[chrom] = [i]

        previous_chrom = chrom 
    
    return chrom_index_data   


def bedListToIntervalTree(bed_list):
    r"""
    convert a bed list to interval tree

    >>> bed_list = [('Chr1', 0, 10000), ('Chr1', 10000, 20000)]
    >>> res = bedListToIntervalTree(bed_list)
    >>> res['Chr1']
    [Interval(0, 10000), Interval(10000, 20000)]
    """
    bed_interval_tree = {}
    for i, _interval in enumerate(bed_list):
        chrom, start, end = _interval
        if chrom not in bed_interval_tree:
            bed_interval_tree[chrom] = IntervalTree()
        bed_interval_tree[chrom].add(Interval(start, end, i))

    return bed_interval_tree

def getContigIntervalTreeFromAGP(agp_df):
    r"""
    get contig intervaltree from agp
    
    Return:
        {'Chr1': IntervalTree((0, 100, (contig1, "+")))

    """
    contig_bin_interval_tree = {}
    for i, row in agp_df.iterrows():
        
        chrom = row.name
        contig = row.id
        start = row.start
        end = row.end
        orientation = row.orientation
        if chrom not in contig_bin_interval_tree:
            contig_bin_interval_tree[chrom] = IntervalTree()
        
        contig_bin_interval_tree[chrom].add(
            Interval(start, end, (contig, orientation)))

    return contig_bin_interval_tree

def splitAGPByBinSize(agp_df, bin_size=1000):
    r"""
    split chromosome regions into specifial intervals for agp file

    """
    logger.debug("Binnify the agp ...")
    agp_df['start'] -= 1
    agp_df['tig_start'] = agp_df['tig_start'].astype(int) - 1
    
    def process_row(row):
        if int(row.end - row.start + 1) <= bin_size:
            return pd.DataFrame([row])
        else:
            tmp_chrom_bins = get_bins(bin_size,  [(row.name, int(row.end))],
                                    start=row.start, orientation=row.orientation)
            tig_start = row.tig_start - row.tig_start % bin_size
            tmp_contig_bins = get_bins(
                bin_size, [(row.id, int(row.tig_end))], start=tig_start,
                reverse=False if row.orientation == "+" else True)

            rows = []
            for (_, start, end), (_, tig_start, tig_end) in \
                    zip(tmp_chrom_bins, tmp_contig_bins):
                # tmp_row = row.copy()
                # tmp_row[['start', 'end', 'tig_start', 'tig_end']] = start, \
                #         end, tig_start, tig_end
               
                # rows.append(tmp_row)
                rows.append({
                'chrom': row.name,
                'start': start,
                'end': end,
                'number': row.number,
                'type': row.type,
                'id': row.id,
                'tig_start': tig_start,
                'tig_end': tig_end,
                'orientation': row.orientation
                })


        return pd.DataFrame(rows)

    res_df = pd.concat(agp_df.parallel_apply(process_row, axis=1).tolist())
    res_df.set_index('chrom', inplace=True)

    return res_df
          

def findIntersectRegion(a, b, fraction=0.51):
    r"""
    find the region by minimum fraction
    >>> a
    IntervalTree([0, 100], [995, 1200])
    >>> b
    (0, 1000)
    >>> findIntersectRegion(a, b)
    [(0, 100)]
    """
    overlaps = sorted(a.overlap(*b))
    tree_begin, tree_end = b[0], b[1]
    
    if not overlaps:
        return overlaps
    start, end = overlaps[0], overlaps[-1]
    if (start.begin - tree_begin) < (start.length() * 0.5):
        overlaps.pop(0)
    if (end.end - tree_end) > (end.length() * 0.5):
        overlaps.pop(-1)
    
    return overlaps

def getContigOnChromBins(chrom_bins, contig_on_chrom_bins, dir="."):
    log_dir = "logs"
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    logger.debug("Get the coordinates of contig on chrom.")
    _file1 = tempfile.NamedTemporaryFile(suffix=".a", delete=True, dir=dir)
    _file2 = tempfile.NamedTemporaryFile(suffix=".b", delete=True, dir=dir)
    _file3 = tempfile.NamedTemporaryFile(suffix=".out", delete=True, dir=dir)

    chrom_bins.to_csv(_file1.name, sep='\t', header=None, index=None)

    contig_on_chrom_bins.to_csv(_file2.name,
                                sep='\t', index=True, header=None)

    cmd = "bedtools intersect -a {} -b {} -F 0.5 -wo > {} ".format(
                        _file1.name, _file2.name, _file3.name)
    # flag = os.system('bedtools intersect -a {} -b {} -F 0.5 -wo > {} 2>/dev/null'.format(
                        # _file1.name, _file2.name, _file3.name))
    flag = os.system(cmd + f" 2>{log_dir}/plot.bedtools.log")
    assert flag == 0 , f"Failed to execute command '{cmd}', please check log `plot.bedtools.log`."

    df = pd.read_csv(f"{_file3.name}", sep='\t',
                     header=None, index_col=None,
                     dtype={9: 'category'})
    
    _file1.close()
    _file2.close()
    _file3.close()

    df = df[~df.where(df == '.').any(axis=1)]

    # df_duplicated = df[df.duplicated([6, 7, 8], keep=False)]
    # df = df.drop(df_duplicated.groupby([6, 7, 8])[10].idxmin().values, axis=0)
    # df_duplicated = df[df.duplicated([0, 1, 2], keep=False)]
    # df = df.drop(df_duplicated.groupby([0, 1, 2])[10].idxmin().values, axis=0)
    df.drop([3, 4, 5, 10], axis=1, inplace=True)
    df.columns = ['chrom', 'start', 'end', 'contig',
                  'tig_start', 'tig_end', 'orientation']

    

    return df 

    _contig_on_chrom_bins = contig_on_chrom_bins.reset_index()
    _chrom_bins = chrom_bins.copy()
    _chrom_bins.reset_index(inplace=True)
    _chrom_bins.columns = ['Index', 'Chromosome', 'Start', 'End']
    _contig_on_chrom_bins.columns = ['Chromosome', 'Start', 'End', 'Contig',
                                    'Tig_start', 'Tig_end', 'Orientation']
    
    
    gr1 = pr.PyRanges(_chrom_bins)
    gr2 = pr.PyRanges(_contig_on_chrom_bins)
    res_df = gr1.join(gr2, report_overlap=True, how='right', preserve_order=True ).df 
    res_df.sort_values(by=['Index'], inplace=True)
    res_df['Length'] = res_df['End_b'] - res_df['Start_b']
    res_df['Fraction'] = res_df['Overlap'] / res_df['Length']
    res_df = res_df[res_df['Fraction'] >= 0.5]
    res_df.drop(['Overlap', 'Length', 'Fraction', 'Index', 'Start_b', 'End_b'], axis=1, inplace=True)
   
    res_df.columns = ['chrom', 'start', 'end', 'contig',
                      'tig_start', 'tig_end', 'orientation']
    
    res_df.reset_index(drop=True, inplace=True)
    return res_df

    def process_chromosome(args):
        gr1 = pr.PyRanges(args[0])
        gr2 = pr.PyRanges(args[1])
        res_df = gr1.join(gr2, report_overlap=True, how='right', preserve_order=True).df 
        res_df.sort_values(by=['Index'], inplace=True)
        res_df['Length'] = res_df['End_b'] - res_df['Start_b']
        res_df['Fraction'] = res_df['Overlap'] / res_df['Length']
        res_df = res_df[res_df['Fraction'] >= 0.5]
        res_df.drop(['Overlap', 'Length', 'Fraction', 'Index', 'Start_b', 'End_b'], axis=1, inplace=True)
        res_df.columns = ['chrom', 'start', 'end', 'contig',
                        'tig_start', 'tig_end', 'orientation']
        res_df.reset_index(drop=True, inplace=True)
        return res_df

    chromosomes = _chrom_bins['Chromosome'].unique()
    
    args = [(_chrom_bins[_chrom_bins['Chromosome'] == chrom], 
             _contig_on_chrom_bins[_contig_on_chrom_bins['Chromosome'] == chrom]) 
             for chrom in chromosomes]


    results = Parallel(n_jobs=4)(delayed(process_chromosome)(arg) for arg in args)
    final_res_df = pd.concat(results, ignore_index=True)

    return final_res_df

def chrRangeID(args, axis=0):
    """
    Chrom range transformation.
    Examples:
    --------
    >>> args = ["Chr1", 100, 200]
    >>> chrRangeID(args)
    "Chr1:100-200"
    >>> args = "Chr1:100-200"
    >>> chrRangeID(args, axis=1)
    ("Chr1", "100", "200")
    """
    if axis == 0:
        chrom, start, end = map(str, args)
        return "{}:{}-{}".format(chrom, start, end)
    elif axis == 1:
        chrom, ranges = args.split(':')
        start, end = ranges.split('-')
        return chrom, start, end
    else:
        return 

class OrderContigMatrix(object):
    """
    reorder the contig matrix to chromosome order
    """
    def __init__(self):
        pass
    

def aggregate_chunk_static(args):
    chunk, contig2chrom_index, agg_dict, index_columns = args
    
    new_bin1_id = np.searchsorted(contig2chrom_index, chunk['bin1_id'].to_numpy(), side='right') - 1
    new_bin2_id = np.searchsorted(contig2chrom_index, chunk['bin2_id'].to_numpy(), side='right') - 1
    
    chunk = chunk.with_columns([
        pl.Series('bin1_id', new_bin1_id),
        pl.Series('bin2_id', new_bin2_id)
    ])

    chunk = chunk.filter((pl.col('bin1_id') >= 0) & (pl.col('bin2_id') >= 0))
    
    if agg_dict['count'] == 'sum':
        result = (chunk.group_by(index_columns, maintain_order=True)
                    .agg(pl.sum('count')).to_pandas().reset_index())
    else:
        result = (chunk.to_pandas().groupby(index_columns, sort=True)
                        .agg(agg_dict)
                        .reset_index())
    
    result = result[(result['bin1_id'] >= 0) & (result['bin2_id'] >= 0)]

    return result

class sumSmallContig(object):
    """
    Sum small conitg count into one chromosome bin
    """
    def __init__(self, chrom_pixels, contig2chrom, edges,
                    columns, map, batchsize=10000):

        self._map = map
        # self.cool_path = cool_path
        self.contig2chrom = contig2chrom
        self.chromidx = contig2chrom.index.values
        self.contigidx = contig2chrom['contigidx'].values
        self.contig2chrom = self.contig2chrom.reset_index()

       

        contig2chrom_index = self.contig2chrom.groupby('chromidx')['contigidx'].first()
        self.contig2chrom_index = contig2chrom_index.values
        # self.contig_edges = list(zip(contig2chrom_index[:-1], 
        #                             contig2chrom_index[1:]))
        
        self.batchsize = batchsize 

        self.newbins_index = self.contig2chrom.index.values
        self.index_columns = ['bin1_id', 'bin2_id']
        self.value_coumns = list(columns)
        self.agg = {'count': 'sum'}
        
        edges = self.contig2chrom_index[np.arange(0, len(self.contig2chrom_index), batchsize)]
        edges = np.append(edges, self.contig2chrom_index[-1])
        self.edges = list(zip(edges[:-1], edges[1:]))

        self.pixels = chrom_pixels

        bin1 = self.pixels['bin1_id'].to_numpy()
        self.pixel_boundaries = np.searchsorted(bin1, edges, side='left')

    # # @profile
    # def _aggregate(self, span):

    #     # cool = cooler.Cooler(self.cool_path)
    
    #     # pixels = self.pixels
    #     # pixels = cool.matrix(balance=False, sparse=True, as_pixels=True)
    #     contig2chrom_index = self.contig2chrom_index
    #     lo, hi = span
        
    #     chunk = self.pixels.filter(pl.col('bin1_id').is_between(lo, hi))
        
    #     new_bin1_id = np.searchsorted(contig2chrom_index, chunk['bin1_id'].to_numpy(), side='right') - 1
    #     new_bin2_id = np.searchsorted(contig2chrom_index, chunk['bin2_id'].to_numpy(), side='right') - 1
        
    #     chunk = chunk.with_columns([
    #         pl.Series('bin1_id', new_bin1_id),
    #         pl.Series('bin2_id', new_bin2_id)
    #     ])

    #     chunk = chunk.filter((pl.col('bin1_id') >= 0) & (pl.col('bin2_id') >= 0))
    #     if self.agg['count'] == 'sum':
    #         result = (chunk.group_by(self.index_columns, maintain_order=True)
    #                     .agg(pl.sum('count')).to_pandas().reset_index())
    #     else:
    #         result = (chunk.to_pandas().groupby(self.index_columns, sort=True)
    #                         .agg(self.agg)
    #                         .reset_index())
        
    #     result = result[(result['bin1_id'] >= 0) & (result['bin2_id'] >= 0)]

    #     return result

    # def aggregate(self, span):
    #     try:
    #         chunk = self._aggregate(span)
    
    #     except MemoryError as e:
    #         raise RuntimeError(str(e))
    #     return chunk
    

    def __iter__(self):
        
        batchsize = self.batchsize
        boundaries = self.pixel_boundaries
        spans = list(zip(boundaries[:-1], boundaries[1:]))
        spans_length = len(spans)
        
        for i in range(0, spans_length, batchsize):
            current_spans = spans[i: i+batchsize]
            tasks = []
            for start, end in current_spans:
                if start >= end:
                    continue
                chunk = self.pixels[start:end]
                tasks.append((chunk, self.contig2chrom_index, self.agg, self.index_columns))
            
            if not tasks:
                continue

            try:
                if batchsize > 1:
                    lock.acquire()
                results = self._map(aggregate_chunk_static, tasks)
                
            finally:
                if batchsize > 1:
                    lock.release()
        
            for df in results:
                yield {k: v.values for k, v in df.items()}

    # def __iter__(self):
        
    #     batchsize = self.batchsize
    #     spans = self.edges
    #     spans_length = len(spans)
    #     for i in range(0, spans_length, batchsize):
    #         try:
    #             if batchsize > 1:
    #                 lock.acquire()
    #             results = self._map(self.aggregate, spans[i: i+batchsize])
                
    #         finally:
    #             if batchsize > 1:
    #                 lock.release()
        
    #         for df in results:
    #             yield {k: v.values for k, v in df.items()}
                
def sum_small_contig(chrom_pixels, contig2chrom, new_bins, output, 
                dtypes=None, columns=['count'], threads=1, **kwargs):
    from cooler.create import create

    # cool = cooler.Cooler(cool_path)

    chrom_counts = new_bins['chrom'].value_counts().sort_index()
    edges = np.r_[0, np.cumsum(chrom_counts)]
    # edges = np.r_[0, np.cumsum(new_bins.groupby(
    #     'chrom').count().reset_index()['start'].tolist())]
    edges = list(zip(edges[:-1], edges[1:]))
    
    if dtypes is None:
        dtypes = {'count': 'float64'}
    # input_dtypes = cool.pixels().dtypes
    input_dtypes = dtypes
    for col in columns:
        if col in input_dtypes:
            dtypes.setdefault(col, input_dtypes[col])
    
    try:
        if threads > 1:
            pool = Pool(threads)
            kwargs.setdefault('lock', lock)
               
            iterator = sumSmallContig(chrom_pixels, 
                contig2chrom,
                edges,
                columns,
                map=pool.map if threads > 1 else map, 
            )
        else:
            iterator = sumSmallContig(chrom_pixels, 
                    contig2chrom,
                    edges,
                    columns,
                    map=map, 
                )

        #kwargs.setdefault("append", True)
    
        create(output, new_bins, iterator, dtypes=dtypes,
                triucheck=False, dupcheck=False, boundscheck=False,
                **kwargs)
    
    finally:
        if threads > 1:
            pool.close()

def adjust_matrix_slow(matrix, agp, outprefix=None, chromSize=None, threads=4):

    start_time = time.time()
   
    cool = cooler.Cooler(matrix)
    agp_df, _ = import_agp(agp)
    bin_size = int(cool.binsize)

    if outprefix is None:
        agpprefix = op.basename(agp).rsplit(".", 1)[0]
        outprefix = op.basename(cool.filename).rsplit(".", 1)[0]
        outprefix = agpprefix + "." + outprefix
    ## get chromosome size database from arguments or agp file
    if not chromSize:
        chrom_sizes = agp_df.groupby(agp_df.index)['end'].max()
        chrom_sizes = pd.DataFrame(chrom_sizes)
        chrom_sizes.reset_index(inplace=True)
        chrom_sizes.columns = ['chrom', 'length']
    else: 
        chrom_sizes = pd.read_csv(chromSize, sep='\t', header=None,
                                index_col=0, names=['chrom', 'length'])

    chrom_sizes = [i for _, i in chrom_sizes.iterrows()]

    logger.info("Converting contig-level coordinate to chromosome-level")
    chrom_bin_interval_df = pd.DataFrame(get_bins(bin_size, chrom_sizes), 
                                    columns=['chrom', 'start', 'end'])
    
    pandarallel.initialize(nb_workers=threads, verbose=0)
    chrom_regions = chrom_bin_interval_df.parallel_apply(lambda x: chrRangeID(x.values), axis=1)
    contig_bins = cool.bins()
    # contig_idx_db = dict(zip(cool.chromnames, range(len(cool.chromnames))))
    new_agp = splitAGPByBinSize(agp_df, bin_size=bin_size)

    split_contig_on_chrom_df = getContigOnChromBins(
        chrom_bin_interval_df, new_agp.drop(['number', 'type'], axis=1))
    split_contig_on_chrom_df.drop(['orientation'], inplace=True, axis=1)

    
    contig_bins_index = contig_bins[:][['chrom', 'start' ,'end']].parallel_apply(
        lambda x: chrRangeID(x.values), axis=1)
    contig_bins_index = contig_bins_index.to_frame().reset_index()
    contig_bins_index.rename(
        columns={'index': 'contigidx', 0: 'contig_region'}, inplace=True)
    contig_bins_index.set_index('contig_region', inplace=True)

    chrom_bins_index = chrom_bin_interval_df.parallel_apply(
        lambda x: chrRangeID(x.values), axis=1)
    chrom_bins_index = chrom_bins_index.to_frame().reset_index()
    chrom_bins_index.rename(
        columns={'index': 'chromidx', 0: 'chrom_region'}, inplace=True)
    chrom_bins_index.set_index('chrom_region', inplace=True)

    def func(row): return chrRangeID(row.values)
    split_contig_on_chrom_df['chrom_region'] = split_contig_on_chrom_df[[
        'chrom', 'start', 'end']].parallel_apply(func, axis=1)
    # split_contig_on_chrom_df['contig'] = split_contig_on_chrom_df[
    #     'contig'].parallel_apply(lambda x: contig_idx_db.get(x, np.nan)).dropna()
    split_contig_on_chrom_df['contig_region'] = split_contig_on_chrom_df[[
        'contig', 'tig_start', 'tig_end']].parallel_apply(func, axis=1)

    cat_dtype = pd.CategoricalDtype(categories=chrom_regions,    
                                        ordered=True)
    split_contig_on_chrom_df['chrom_region'] = \
        split_contig_on_chrom_df['chrom_region'].astype(cat_dtype)
    split_contig_on_chrom_df['chromidx'] = \
        split_contig_on_chrom_df['chrom_region'].cat.codes.values
    
    contig_bins_index_set = set(contig_bins_index.index)
    split_contig_on_chrom_df = split_contig_on_chrom_df[
        split_contig_on_chrom_df['contig_region'].isin(contig_bins_index_set)]

    split_contig_on_chrom_df['contigidx'] = contig_bins_index.loc[
        split_contig_on_chrom_df['contig_region']]['contigidx'].values
    split_contig_on_chrom_df.sort_values(by=['chromidx'], inplace=True)

    split_contig_on_chrom_df.drop(['chrom', 'start', 'end',
                                   'contig', 'tig_start', 'tig_end',
                                   'chrom_region', 'contig_region'],
                                  inplace=True, axis=1)
    raw_chrom_interval_df = chrom_bin_interval_df.copy()
    ## Solve the problem that the dotted line of the heat map will drift when the bin is relatively large
    chrom_bin_interval_df = chrom_bin_interval_df.loc[split_contig_on_chrom_df['chromidx'].drop_duplicates()]
    chrom_bin_interval_df.reset_index(drop=True, inplace=True)
    shift_index = False
    if chrom_bin_interval_df.loc[0, 'start'] != 0:
        ## if the first bin is not start with 0, then we need to insert a new row before the first bin
        # logger.warning('The first bin of chromosome is not start with 0, '
        #                   'adjust the first bin to start with 0.')
        shift_index = True
        new_row = pd.DataFrame({
            'chrom': raw_chrom_interval_df.iloc[0]['chrom'],
            'start': 0,
            'end': raw_chrom_interval_df.iloc[0]['end']
        }, index=[0])
        chrom_bin_interval_df = pd.concat([new_row, chrom_bin_interval_df], ignore_index=True)
        
    
    logger.info('Starting to reorder matrix ...')

    matrix = cool.matrix(balance=False, sparse=True)
  
    dtypes = dict(cool.pixels().dtypes)
   
    contig2chrom = split_contig_on_chrom_df[['chromidx', 'contigidx']]
    contig2chrom.set_index('chromidx', inplace=True)
    # reorder matrix 
    # reordered_contigidx = contig2chrom['contigidx'].values
    # reordered_matrix = matrix[:].tocsr(
    #                      )[:, reordered_contigidx][reordered_contigidx, :]
    # reordered_matrix = triu(reordered_matrix).tocoo()
    n_bins = matrix.shape[0]
    reordered_contigidx = contig2chrom['contigidx'].values
    
    map_array = np.full(n_bins, -1, dtype=np.int32)
    map_array[reordered_contigidx] = np.arange(len(reordered_contigidx), dtype=np.int32)

    pixels = cool.pixels()[:]
    old_bin1 = pixels['bin1_id'].to_numpy()
    old_bin2 = pixels['bin2_id'].to_numpy()
    counts = pixels['count'].to_numpy()

    new_bin1 = map_array[old_bin1]
    new_bin2 = map_array[old_bin2]

    valid_mask = (new_bin1 != -1) & (new_bin2 != -1)
    
    new_bin1 = new_bin1[valid_mask]
    new_bin2 = new_bin2[valid_mask]
    counts = counts[valid_mask]

    swap_mask = new_bin1 > new_bin2
    final_bin1 = np.where(swap_mask, new_bin2, new_bin1)
    final_bin2 = np.where(swap_mask, new_bin1, new_bin2)
    
    os.environ["POLARS_MAX_THREADS"]  = str(threads)

    chrom_pixels = pl.DataFrame({
        'bin1_id': final_bin1,
        'bin2_id': final_bin2,
        'count': counts,
    })
    chrom_pixels = chrom_pixels.sort('bin1_id')

    # os.environ["POLARS_MAX_THREADS"]  = str(threads)

    # chrom_pixels = pl.DataFrame({
    #     'bin1_id': reordered_matrix.row,
    #     'bin2_id': reordered_matrix.col,
    #     'count': reordered_matrix.data,
    # })
    # chrom_pixels = chrom_pixels.sort('bin1_id')

    # del reordered_matrix
    # gc.collect()
    # order_cool_path = f"{outprefix}.ordered.cool"
    # cooler.create_cooler(order_cool_path, reordered_contig_bins,
    #                      chrom_pixels, dtypes=dtypes, triucheck=False,
    #                      dupcheck=False, boundscheck=False)#, metadata=HIC_METADATA)
    # logger.info('Successful, reorder the contig-level matrix, '
    #             f' and output into `{outprefix}.ordered.cool`')
    
    logger.info('Starting to convert contig bin to chromosome bin ...')
    contig2chrom['contigidx'] = range(len(contig2chrom))
    # contig2chrom = contig2chrom.reset_index().set_index('chromidx')

    if shift_index:
        contig2chrom['chromidx'] += 1
        dummy = pd.DataFrame({'chromidx': [0], 'contigidx': [0]})
        contig2chrom = pd.concat([dummy, contig2chrom], ignore_index=True)

    contig2chrom.reset_index().set_index('chromidx', inplace=True)

    sum_small_contig(chrom_pixels, contig2chrom, chrom_bin_interval_df, 
                     f'{outprefix}.chrom.cool', dtypes=dtypes)#, metadata=HIC_METADATA)
    logger.info('Successful, converted the contact into chromosome-level'
                f' and output into `{outprefix}.chrom.cool`')
    

    logger.info('Successful, adjusted matrix to chromosome-level, elasped time {:.2f}s'.format(time.time() - start_time))
    
    return f'{outprefix}.chrom.cool'

# def build_bin_mapping(cool, agp_df, bin_size, coarsen_factor=1):
#     target_bin_size = bin_size * coarsen_factor
    
#     contig_bins = cool.bins()[['chrom', 'start', 'end']][:]
#     n_old_bins = len(contig_bins)
#     dtype_id = np.int32 if n_old_bins < 2**31 - 1 else np.int64

#     agp_df = agp_df.copy()

#     chrom_lengths = agp_df.groupby(agp_df.index)['end'].max().to_dict()

#     new_bins_list = []
#     chrom_offset = {}
#     current_max_bin = 0

#     sorted_chroms = list(dict.fromkeys(agp_df.index))
    
#     for chrom in sorted_chroms:
#         length = chrom_lengths[chrom]

#         n_bins = int(np.ceil(length / target_bin_size))
#         chrom_offset[chrom] = current_max_bin

#         starts = np.arange(0, length, target_bin_size, dtype=np.int64) 
#         ends = np.clip(starts + target_bin_size, 0, length).astype(np.int64)
        
#         chrom_df = pd.DataFrame({
#             'chrom': chrom,
#             'start': starts,
#             'end': ends
#         })
#         new_bins_list.append(chrom_df)
#         current_max_bin += n_bins
        
#     new_bins_df = pd.concat(new_bins_list, ignore_index=True)
    
#     map_array = np.full(n_old_bins, -1, dtype=dtype_id)

#     if contig_bins['chrom'].dtype == 'object':
#         contig_bins['chrom'] = contig_bins['chrom'].astype('category')

#     contig_groups = contig_bins.reset_index().groupby('chrom')['index'].agg(['min', 'max'])
#     contig_id_map = contig_groups.to_dict('index')
    
#     contig_starts = contig_bins['start'].values
#     contig_ends = contig_bins['end'].values
    
#     for chrom_name, row in agp_df.iterrows():
#         contig_name = row.id
#         orientation = row.orientation
#         agp_start = row.start - 1 # 0-based
        
#         if contig_name not in contig_id_map:
#             continue
            
#         c_info = contig_id_map[contig_name]
#         old_start_id, old_end_id = c_info['min'], c_info['max']

#         tig_starts = contig_starts[old_start_id : old_end_id+1]
        
#         if orientation == '+':
#             chrom_coords = agp_start + tig_starts
#         else:
#             tig_len = row.tig_end
#             tig_ends = contig_ends[old_start_id : old_end_id+1]
#             chrom_coords = agp_start + (tig_len - tig_ends)
  
#         rel_bin_indices = (chrom_coords // target_bin_size).astype(dtype_id)
#         base_offset = chrom_offset[chrom_name]
#         new_bin_ids = base_offset + rel_bin_indices
        
#         map_array[old_start_id : old_end_id+1] = new_bin_ids
    
#     return map_array, new_bins_df

def build_bin_mapping(cool, agp_df, bin_size, coarsen_factor=1):
    target_bin_size = bin_size * coarsen_factor
    
    contig_bins = cool.bins()[['chrom', 'start', 'end']][:]
    n_old_bins = len(contig_bins)
    
    agp_df = agp_df.copy()
    chrom_lengths = agp_df.groupby(agp_df.index)['end'].max().to_dict()

    new_bins_list = []
    chrom_offset = {}
    current_max_bin = 0
    sorted_chroms = list(dict.fromkeys(agp_df.index))
    
    for chrom in sorted_chroms:
        length = chrom_lengths[chrom]
        n_bins = int(np.ceil(length / target_bin_size))
        chrom_offset[chrom] = current_max_bin
        starts = np.arange(0, length, target_bin_size, dtype=np.int64) 
        ends = np.clip(starts + target_bin_size, 0, length).astype(np.int64)
        chrom_df = pd.DataFrame({'chrom': chrom, 'start': starts, 'end': ends})
        new_bins_list.append(chrom_df)
        current_max_bin += n_bins
        
    new_bins_df = pd.concat(new_bins_list, ignore_index=True)
    
    if contig_bins['chrom'].dtype == 'object':
        contig_bins['chrom'] = contig_bins['chrom'].astype('category')

    contig_groups = contig_bins.reset_index().groupby('chrom')['index'].agg(['min', 'max'])
    contig_id_map = contig_groups.to_dict('index')
    contig_starts = contig_bins['start'].values
    contig_ends = contig_bins['end'].values
    
    mapping_dict = defaultdict(list)

    for chrom_name, row in agp_df.iterrows():
        contig_name = row.id
        orientation = row.orientation
        agp_start = row.start - 1 
        
        if contig_name not in contig_id_map:
            continue
            
        c_info = contig_id_map[contig_name]
        old_start_id, old_end_id = c_info['min'], c_info['max']

        num_bins = old_end_id - old_start_id + 1
        old_bin_ids = np.arange(old_start_id, old_end_id + 1)
        
        tig_starts = contig_starts[old_start_id : old_end_id+1]
        
        if orientation == '+':
            chrom_coords = agp_start + tig_starts
        else:
            tig_len = row.tig_end
            tig_ends = contig_ends[old_start_id : old_end_id+1]
            chrom_coords = agp_start + (tig_len - tig_ends)
  
        rel_bin_indices = (chrom_coords // target_bin_size).astype(np.int64)
        base_offset = chrom_offset[chrom_name]
        new_bin_ids = base_offset + rel_bin_indices

        for ob, nb in zip(old_bin_ids, new_bin_ids):
            mapping_dict[ob].append(nb)

    map_data = []
    for ob, nbs in mapping_dict.items():
        for nb in nbs:
            map_data.append((ob, nb))

    map_pl = pl.DataFrame(map_data, schema={'old': pl.Int32, 'new': pl.Int32})
    
    return map_pl, new_bins_df


# def process_chunk(chunk_df, map_array, temp_path, chunk_idx, bin_dtype, count_dtype):
#     chunk_df = chunk_df[(chunk_df['bin1_id'] > 0) & (chunk_df['bin2_id'] > 0)]
#     old_b1 = chunk_df['bin1_id'].values
#     old_b2 = chunk_df['bin2_id'].values
#     counts = chunk_df['count'].values
    
    
#     new_b1 = map_array[old_b1]
#     new_b2 = map_array[old_b2]
    
#     valid_mask = (new_b1 != -1) & (new_b2 != -1)
#     if not np.all(valid_mask):
#         new_b1 = new_b1[valid_mask]
#         new_b2 = new_b2[valid_mask]
#         counts = counts[valid_mask]
    
#     if len(new_b1) == 0:
#         return False

#     swap_mask = new_b1 > new_b2
#     final_b1 = np.where(swap_mask, new_b2, new_b1)
#     final_b2 = np.where(swap_mask, new_b1, new_b2)

#     chunk_pl = pl.DataFrame({
#         'bin1_id': final_b1,
#         'bin2_id': final_b2,
#         'count': counts
#     }).select([
#         pl.col('bin1_id').cast(bin_dtype),
#         pl.col('bin2_id').cast(bin_dtype),
#         pl.col('count').cast(count_dtype)
#     ])
    
#     chunk_pl.write_parquet(temp_path / f"chunk_{chunk_idx}.parquet")

#     return True

def process_chunk(chunk_df, map_pl, temp_path, chunk_idx, bin_dtype, count_dtype):
    chunk_pl = pl.from_pandas(chunk_df)

    chunk_pl = chunk_pl.filter((pl.col('bin1_id') >= 0) & (pl.col('bin2_id') >= 0))
    if chunk_pl.height == 0:
        return False

    chunk_pl = chunk_pl.with_columns([
        pl.col('bin1_id').cast(pl.Int32),
        pl.col('bin2_id').cast(pl.Int32)
    ])

    chunk_mapped = chunk_pl.join(
        map_pl, left_on='bin1_id', right_on='old', how='inner'
    ).rename({'new': 'new_bin1_id'})

    if chunk_mapped.height == 0:
        return False

    chunk_mapped = chunk_mapped.join(
        map_pl, left_on='bin2_id', right_on='old', how='inner'
    ).rename({'new': 'new_bin2_id'})

    if chunk_mapped.height == 0:
        return False
    
    final_pl = chunk_mapped.select([
        pl.when(pl.col('new_bin1_id') > pl.col('new_bin2_id'))
          .then(pl.col('new_bin2_id'))
          .otherwise(pl.col('new_bin1_id'))
          .cast(bin_dtype)
          .alias('bin1_id'),
          
        pl.when(pl.col('new_bin1_id') > pl.col('new_bin2_id'))
          .then(pl.col('new_bin1_id'))
          .otherwise(pl.col('new_bin2_id'))
          .cast(bin_dtype)
          .alias('bin2_id'),
          
        pl.col('count').cast(count_dtype)
    ])
    
    final_pl.write_parquet(temp_path / f"chunk_{chunk_idx}.parquet")

    return True

def adjust_matrix(matrix, agp, outprefix=None, chromSize=None, threads=4, coarsen_factor=1):

    from concurrent.futures import ThreadPoolExecutor

    start_time = time.time()
    
    cool = cooler.Cooler(matrix)
    agp_df, _ = import_agp(agp)
    bin_size = int(cool.binsize)

    target_bin_size = bin_size * coarsen_factor

    if outprefix is None:
        agpprefix = op.basename(agp).rsplit(".", 1)[0]
        outprefix = op.basename(cool.filename).rsplit(".", 1)[0]
        
        input_res_str = to_humanized2(bin_size)
        target_res_str = to_humanized2(target_bin_size)
        
        if input_res_str in outprefix:
            outprefix = outprefix.replace(input_res_str, target_res_str)
        
        outprefix = agpprefix + "." + outprefix

    logger.info(f"Generating bin mapping (Contig -> Chromosome) with coarsen factor {coarsen_factor}...")
    
    map_array, new_bins_df = build_bin_mapping(cool, agp_df, bin_size, coarsen_factor=coarsen_factor)
    new_bins_df['chrom'] = new_bins_df['chrom'].astype(str)
    new_bins_df['start'] = new_bins_df['start'].astype(np.int64)
    new_bins_df['end'] = new_bins_df['end'].astype(np.int64)
    # if (new_bins_df['end'] <= new_bins_df['start']).any():
    #     new_bins_df = new_bins_df[new_bins_df['end'] > new_bins_df['start']].reset_index(drop=True)

    max_bin_id = new_bins_df['end'].max()
    bin_dtype = pl.Int32 if max_bin_id < 2**32 / 2 else pl.Int64

    logger.info("Loading and aggregating pixels (Chunked Processing)...")
    
    os.environ["POLARS_MAX_THREADS"] = str(threads)
    
    chunksize = 10_000_000 
    nnz = cool.info['nnz']

    pixel_selector = cool.pixels()

    temp_dir = tempfile.TemporaryDirectory(
        dir=op.dirname(op.abspath(outprefix)) if outprefix else None, prefix="plot_tmp_"
    )
    temp_path = Path(temp_dir.name)
    logger.info(f"Using temporary directory: {temp_path}")

    chunk_idx = 0
    has_data = False

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = []
        
        for i in range(0, nnz, chunksize):
            lo, hi = i, min(i + chunksize, nnz)
            chunk_df = pixel_selector[lo:hi]
            
            if np.issubdtype(chunk_df['count'].values.dtype, np.integer):
                count_dtype = pl.Int32
            else:
                count_dtype = pl.Float32 

            future = executor.submit(
                process_chunk, 
                chunk_df, 
                map_array, 
                temp_path, 
                chunk_idx, 
                bin_dtype, 
                count_dtype
            )
            futures.append(future)
            chunk_idx += 1
        
        for future in futures:
            if future.result():
                has_data = True

    if not has_data or chunk_idx == 0:
        logger.warning("No valid pixels found after mapping.")
        agg_df = pl.DataFrame({'bin1_id': [], 'bin2_id': [], 'count': []}, 
                              schema={'bin1_id': bin_dtype, 'bin2_id': bin_dtype, 'count': count_dtype})
    else:
        logger.info(f"Processed {chunk_idx} chunks. Starting global aggregation...")
        
        try:
            agg_df = pl.scan_parquet(temp_path / "*.parquet") \
                .group_by(['bin1_id', 'bin2_id']) \
                .agg(pl.col('count').sum()) \
                .sort(['bin1_id', 'bin2_id']) \
                .collect(streaming=True)
        except Exception as e:
            logger.warning(f"Streaming aggregation failed ({e}), falling back to standard collection.")
            agg_df = pl.scan_parquet(temp_path / "*.parquet") \
                .group_by(['bin1_id', 'bin2_id']) \
                .agg(pl.col('count').sum()) \
                .sort(['bin1_id', 'bin2_id']) \
                .collect()

    temp_dir.cleanup()
    
    output_cool = f'{outprefix}.chrom.cool'
    logger.info(f"Writing to cooler file `{output_cool}`")
    
    b1 = agg_df['bin1_id'].to_numpy()
    b2 = agg_df['bin2_id'].to_numpy()
    c = agg_df['count'].to_numpy()
    
    del agg_df
    gc.collect()

    def pixel_iterator(chunksize=10_000_000):
        total_len = len(b1)
        for i in range(0, total_len, chunksize):
            yield {
                'bin1_id': b1[i:i+chunksize],
                'bin2_id': b2[i:i+chunksize],
                'count': c[i:i+chunksize]
            }

    Path(output_cool).parent.mkdir(parents=True, exist_ok=True)
    cooler.create.create(
        output_cool,
        new_bins_df,
        pixel_iterator(), 
        dtypes={'count': 'float64' if count_dtype == pl.Float32 else 'int64',
                'bin1_id': 'int32' if bin_dtype == pl.Int32 else 'int64',
                'bin2_id': 'int32' if bin_dtype == pl.Int32 else 'int64',
                },
        ordered=True,
        triucheck=False,
        dupcheck=False,
        boundscheck=False,
    )
    if coarsen_factor > 1:
        logger.info(f'Successful, adjusted and coarsen matrix to chromosome-level, elapsed time {time.time() - start_time:.2f}s')
    
    else:
        logger.info(f'Successful, adjusted matrix to chromosome-level, elapsed time {time.time() - start_time:.2f}s')
    
    return output_cool

def coarsen_matrix(cool, cool_binsize, k, out, threads):
    """
    coarsen a matrix

    Params:
    --------
    cool: str
         fine resolution cool path
    k: int
        factor 
    out: str
        output path of coarsened cool
    threads: int
        number of threads
    
    Returns:
    --------
    out: str
    
    Examples:
    --------
    >>> coarsen_matrix("sample.10k.cool", 50, "sample.500k.cool', 10)
    "sample.500k.cool"
    """
    from cooler.cli.coarsen import coarsen
    if not out:
        input_resolution = cool_binsize
        out_resolution = f"{k * input_resolution}"
        
        out = cool.replace(str(input_resolution), to_humanized2(out_resolution))
        out = out.replace(to_humanized2(input_resolution), to_humanized2(out_resolution))

        if to_humanized2(out_resolution) not in out:
            if "chrom.cool" in out:
                out = out.replace("chrom.cool", f"{to_humanized2(out_resolution)}.chrom.cool")
            else:
                out = out.replace("cool", f"{to_humanized2(out_resolution)}.chrom.cool")
            

    try:
        coarsen.main(args=[cool, 
                            '--out', out, 
                            '--factor', k, 
                            '--nproc', threads], 
                            prog_name='cooler')
        logger.info(f'Coarsen new cool into {out}.')
    except SystemExit as e:
        exc_info = sys.exc_info()
        exit_code = e.code
        if exit_code is None:
            exit_code = 0
        
        if exit_code != 0:
            raise e

    logger.info(f'Coarsen new cool into `{out}`.')

    return out

def _dot_block(args):
    idx, w_seg = args
    blk = _BALANCE_BLOCKS[idx]
    return idx, blk.dot(w_seg)

def balance_matrix(cool_path,
                        out_weights=None,
                        black_bed=None,
                        max_iter=200,
                        tol=1e-5,
                        min_count=0.0,
                        max_percentile=99.9,
                        threads=4,
                        column_name='weight',
                        report_every=50):

    logger.info(f"[balance_fast] load `{cool_path}`")
    c = cooler.Cooler(cool_path)
    n = c.shape[0]

    chrom_offset = c._load_dset('indexes/chrom_offset')
    bins_df = c.bins()[:][['chrom','start','end']]
    bins_df['idx'] = np.arange(len(bins_df))

    M_full = c.matrix(balance=False, sparse=True)[:].tocsr()
    logger.info(f"[balance_fast] nnz={M_full.nnz}, shape={M_full.shape}")

    blocks, ranges = [], []
    for i in range(len(chrom_offset) - 1):
        s, e = int(chrom_offset[i]), int(chrom_offset[i+1])
        if e > s:
            blocks.append(M_full[s:e, s:e].tocsr())
            ranges.append((s, e))
    del M_full
    logger.info(f"[balance_fast] blocks={len(blocks)}")

    raw_row_sum_parts = [np.array(blk.sum(axis=1)).ravel() for blk in blocks]
    raw_row_sum = np.zeros(n, dtype=float)
    for (s, e), part in zip(ranges, raw_row_sum_parts):
        raw_row_sum[s:e] = part

    mask = np.ones(n, dtype=bool)
    if black_bed and os.path.exists(black_bed):
        try:
            bed = pd.read_csv(black_bed, sep='\t', header=None, usecols=[0,1,2],
                              names=['chrom','start','end'])
            for chrom, start, end in bed.itertuples(index=False):
                sub = bins_df[(bins_df.chrom == chrom) &
                              (bins_df.end > start) &
                              (bins_df.start < end)]
                mask[sub.idx.values] = False
            logger.info(f"[balance_fast] blacklisted bins: {(~mask).sum()}")
        except Exception as e:
            logger.warning(f"[balance_fast] blacklist parse failed: {e}")

    low = raw_row_sum < float(min_count)
    if max_percentile < 100:
        nz = raw_row_sum[raw_row_sum > 0]
        hi_thresh = np.percentile(nz, max_percentile) if nz.size else np.inf
    else:
        hi_thresh = np.inf
    high = raw_row_sum > hi_thresh
    mask &= ~low & ~high
    valid = mask & (raw_row_sum > 0)
    logger.info(f"[balance_fast] valid bins: {valid.sum()} / {n}")

    w = np.ones(n, dtype=np.float64)
    w[~valid] = 0.0

    global _BALANCE_BLOCKS, _BALANCE_BLOCK_RANGES
    _BALANCE_BLOCKS = blocks
    _BALANCE_BLOCK_RANGES = ranges


    prev_change = None
    pool = None
    try:
        if threads and threads > 1:

            os.environ.setdefault("OMP_NUM_THREADS", "1")
            os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
            os.environ.setdefault("MKL_NUM_THREADS", "1")
            pool = Pool(processes=int(threads))

        for it in range(1, max_iter + 1):
            temp = np.zeros(n, dtype=np.float64)
            if pool:
                tasks = []
                for bi, (s, e) in enumerate(ranges):
                    tasks.append((bi, w[s:e].copy()))
                for bi, vec in pool.map(_dot_block, tasks, chunksize=max(1, len(tasks)//(threads or 1))):
                    s, e = ranges[bi]
                    temp[s:e] = vec
            else:
                for bi, (s, e) in enumerate(ranges):
                    temp[s:e] = blocks[bi].dot(w[s:e])

            scaled_row_sum = w * temp
            denom = scaled_row_sum.copy()
            denom[denom == 0] = 1.0
            upd = 1.0 / denom
            w[valid] *= upd[valid]

            med = np.median(w[valid][w[valid] > 0]) if np.any(w[valid] > 0) else 1.0
            if not np.isfinite(med) or med == 0:
                med = 1.0
            w[valid] /= med

            if pool:
                tasks = [(bi, w[s:e].copy()) for bi, (s, e) in enumerate(ranges)]
                temp[:] = 0
                for bi, vec in pool.map(_dot_block, tasks, chunksize=max(1, len(tasks)//(threads or 1))):
                    s, e = ranges[bi]
                    temp[s:e] = vec
            else:
                for bi, (s, e) in enumerate(ranges):
                    temp[s:e] = blocks[bi].dot(w[s:e])
            scaled_row_sum = w * temp
            cur = scaled_row_sum[valid]
            diff = np.abs(cur - 1.0)
            rel_change = diff.mean() if cur.size else 0.0

            if it % report_every == 0 or it == 1 or it == max_iter:
                logger.debug(f"[balance_fast] iter {it:02d}: mean(|r-1|)={rel_change:.3e}, valid={valid.sum()}")
            if prev_change is not None and rel_change < tol:
                logger.debug(f"[balance_fast] converged at iter {it}, mean(|r-1|)={rel_change:.3e}")
                break
            prev_change = rel_change
    finally:
        if pool:
            pool.close()
            pool.join()

    w[~valid] = np.nan

    if out_weights:
        pd.DataFrame({column_name: w}).to_csv(out_weights, sep='\t', index=False)
        logger.info(f"[balance_fast] write weights → {out_weights}")

    try:
        with h5py.File(cool_path, 'r+') as h5:
            bins_grp = h5['bins']
            if column_name in bins_grp:
                del bins_grp[column_name]
            bins_grp.create_dataset(column_name, data=w, compression='gzip')
        logger.info(f"[balance_fast] appended `{column_name}` to bins.")
    except Exception as e:
        logger.warning(f"[balance_fast] append weights failed: {e}")

    return w

def balance_matrix_cooler(cool, black_bed=None, threads=4):
    from cooler.cli.balance import balance 
    args = [cool, "-p", threads, '--cis-only', '--force']

    if black_bed:
        args.extend(['--blacklist', black_bed])
    
    try:
        logger.info(f'Balancing matrix ...')
        balance.main(args=args, 
                        prog_name='cooler')
       
    except SystemExit as e:
        exc_info = sys.exc_info()
        exit_code = e.code
        if exit_code is None:
            exit_code = 0
        
        if exit_code != 0:
            raise e
    
    if black_bed:
        if Path(black_bed).exists():
            Path(black_bed).unlink()
    
def cluster_chrom_by_inter_contacts(matrix, bins,
                                    hap_pattern=r'(Chr\d+)g(\d+)',
                                    chrom_offset=None,
                                    chromnames=None,
                                    include_cis=True,
                                    agg='sum',  
                                    method='euclidean',
                                    normalize='row',
                                    linkage_method='average',
                                    optimal_order=True):

    import re
    import scipy.sparse as sp
    from scipy.cluster.hierarchy import (
        linkage,
        optimal_leaf_ordering,
        leaves_list,
        fcluster,
    )
    from scipy.spatial.distance import pdist, squareform

    if chromnames is None or chrom_offset is None:
        if isinstance(bins, pd.DataFrame):
            if 'chrom' in bins.columns:
                bdf = bins[['chrom', 'index']].copy()
            else:
                bdf = bins.reset_index()[['chrom', 'index']].copy()
        else:
            raise ValueError("bins need contain 'chrom' and 'index'")
        grp = bdf.groupby('chrom', sort=False)['index'].count()
        chromnames = grp.index.tolist()
        counts = grp.values.astype(int)
        chrom_offset = np.r_[0, np.cumsum(counts)].tolist()
    n_chr = len(chromnames)


    group_members = {}
    group_order = [] 
    pat = re.compile(hap_pattern)
    for i, name in enumerate(chromnames):
        m = pat.match(str(name))
        key = m.group(1) if m else str(name)  
        if key not in group_members:
            group_members[key] = []
            group_order.append(key)
        group_members[key].append(i)

    def block_sum(i, j):
        i0, i1 = chrom_offset[i], chrom_offset[i+1]
        j0, j1 = chrom_offset[j], chrom_offset[j+1]
        blk = matrix[i0:i1, j0:j1]
        s = blk.sum() if sp.issparse(blk) else float(np.nansum(blk))
        if agg == 'mean':
            denom = max((i1 - i0) * (j1 - j0), 1)
            s = s / denom
        return float(s)

    new_order = []
    for key in group_order:
        members = group_members[key]
        k = len(members)
        if k <= 1:
            new_order.extend(members)
            continue

        S = np.zeros((k, k), dtype=float)
        for a in range(k):
            ia = members[a]
            for b in range(k):
                ib = members[b]
                if (not include_cis) and ia == ib:
                    S[a, b] = 0.0
                else:
                    S[a, b] = block_sum(ia, ib)

        if k == 2:
            new_order.extend(members)
            continue
        S = np.nan_to_num(S, nan=0.0, posinf=0.0, neginf=0.0)
        S_max = float(np.nanmax(S)) if np.isfinite(S).any() else 0.0
        D = (S_max - S).astype(float)
        np.fill_diagonal(D, 0.0)
        inter = D
        if normalize == 'row':
            rs = inter.sum(axis=1, keepdims=True)
            rs[rs == 0] = 1.0
            inter = inter / rs
        elif normalize == 'sqrt':
            rs = inter.sum(axis=1, keepdims=True)
            cs = inter.sum(axis=0, keepdims=True)
            rs[rs == 0], cs[cs == 0] = 1.0, 1.0
            inter = inter / np.sqrt(rs @ cs)

        inter = np.nan_to_num(inter, nan=0.0, posinf=0.0, neginf=0.0)

        if method == 'corr':
            with np.errstate(invalid='ignore'):
                C = np.corrcoef(inter + 1e-12) 
            C = np.nan_to_num(C, nan=0.0)
            np.fill_diagonal(C, 1.0)
            D = 1.0 - C
            D[D < 0] = 0.0
            dist_condensed = squareform(D, checks=False)
        elif method == 'cosine':
            dist_condensed = pdist(inter, metric='cosine')
        else:
            dist_condensed = pdist(inter, metric='euclidean')

        Z = linkage(dist_condensed, method=linkage_method)
        if optimal_order:
            Z = optimal_leaf_ordering(Z, dist_condensed)
    
        leaves = leaves_list(Z).astype(int)

        # labels = fcluster(Z, t=0.7, criterion='distance')
        # groups = OrderedDict()
        # for leaf in leaves.tolist():
        #     lab = int(labels[leaf])
        #     groups.setdefault(lab, []).append(leaf)
        # order_labs = sorted(groups.keys(),
        #                         key=lambda lab: (-len(groups[lab]), leaves.tolist().index(groups[lab][0])))
        # leaves_ordered = [leaf for lab in order_labs for leaf in groups[lab]]
        # leaves = leaves_ordered

        ordered_members = [members[i] for i in leaves]
        new_order.extend(ordered_members)

    order_idx = np.asarray(new_order, dtype=int)
    order_labels = [chromnames[i] for i in order_idx]
    return order_idx, order_labels

def apply_chrom_order_to_matrix(matrix, chrom_offset, order_idx):
    import scipy.sparse as sp
    ranges = [(chrom_offset[i], chrom_offset[i+1]) for i in range(len(chrom_offset)-1)]
    bin_order = np.concatenate([np.arange(*ranges[i], dtype=int) for i in order_idx])
    if sp.issparse(matrix):
        newM = matrix[bin_order, :][:, bin_order]
    else:
        newM = matrix[np.ix_(bin_order, bin_order)]
    lengths = [ranges[i][1] - ranges[i][0] for i in order_idx]
    new_offset = np.r_[0, np.cumsum(lengths)].tolist()
    return newM, new_offset


def plot_heatmap(matrix, output, 
                    vmin=None, vmax=None,
                    scale="log1p", 
                    triangle=False,
                    xlabel=None, ylabel=None, 
                    xticks=True, yticks=True,
                    rotate_xticks=False, rotate_yticks=False,
                    ytick_min_dist=20,
                    avoid_overlap_yticks=True,
                    remove_short_bin=True,
                    add_hap_border=False,
                    hap_border_color='black',
                    hap_border_width=1.0,
                    add_hap_shadow=False,
                    hap_shadow_color='grey',
                    hap_shadow_alpha=0.1,
                    add_lines=False,
                    line_color='black',
                    line_width=0.5,
                    line_style='--',
                    chromosomes=None, 
                    plot_cis_only=False,
                    plot_hap_only=True,
                    sort_chromosome_by_inter_contacts=False,
                    hap_pattern=r'(Chr\d+)g(\d+)',
                    per_chromosomes=False,
                    group_by_homologs=False,
                    chrom_per_row=4,
                    factor=None,
                    factor_formula="genome_log",
                    size_scale=(3.5, 4.0),
                    tri_size_scale=(7.0, 3.2), 
                    figwidth=None,
                    figheight=None,
                    fontsize=None,
                    dpi=1200, 
                    cmap="redp1_r",
                    cbar_bottom=False,
                    balanced=False,
                    remove_empty=True,
                    threads=1):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.use('Agg')
    from matplotlib.colors import LogNorm

    if hap_pattern:
        hap_pattern = r"{}".format(hap_pattern)

    logger.info(f"Load contact matrix `{matrix}`.")
    cool = cooler.Cooler(matrix)
    if balanced:
        log1p = False 
        log = False 
        scale = ""
        
    if scale == "log1p":
        log1p = True
        log = False
    elif scale == "log":
        log = True 
        log1p = False
    else:
        log = False
        log1p =False

    logger.info("Plotting heatmap ...")

    if not per_chromosomes:
        chrom_offset = cool._load_dset('indexes/chrom_offset')
        chromsizes = cool.chromsizes
        if chromsizes.dtype == np.int32:
            chromsizes = chromsizes.astype(np.uint32).astype(np.int64)

        binsize = cool.binsize
        bins = cool.bins()[:].reset_index(drop=False)
        if bins['start'].dtype == np.int32:
            bins['start'] = bins['start'].astype(np.uint32).astype(np.int64)
        if bins['end'].dtype == np.int32:
            bins['end'] = bins['end'].astype(np.uint32).astype(np.int64)
        bins['chrom'] = bins['chrom'].astype('str')
        bins = bins.set_index('chrom')
        if binsize is None:
            binsize = np.argmax(np.bincount(bins['end'] - bins['start']))

            
        chromnames = cool.chromnames

        if chromosomes:
            if len(chromosomes) == 1:
                matrix = cool.matrix(balance=balanced, sparse=True).fetch(chromosomes[0]).tocsr()
                
            else:
                matrix = cool.matrix(balance=balanced, sparse=True)[:].tocsr()
        else:
            matrix = cool.matrix(balance=balanced, sparse=True)[:].tocsr()
 
        if chromosomes:
            if len(chromosomes) == 1:
                bins = cool.bins().fetch(chromosomes[0])
                new_idx = bins.index.values.tolist()
               
                chrom_offset = [0, len(new_idx)]
                chromnames = chromosomes
                try:
                    chromsizes = chromsizes.loc[chromnames]
                except KeyError:
                    size =  bins['end'].max() - bins['start'].min()
                    chromsizes = pd.DataFrame({"chrom": chromnames, "length": size}).set_index('chrom')
                   
            else:
                bins = bins.loc[chromosomes]
                new_idx = bins['index']
                bins.drop('index', axis=1, inplace=True)
                bins.reset_index(drop=False, inplace=True)
                bins.reset_index(drop=False, inplace=True)
                bins.set_index('chrom', inplace=True)
                matrix = matrix[new_idx, :][:, new_idx]
                chromnames = chromosomes
                chromsizes = chromsizes.loc[chromnames]
            
                chrom_offset = np.r_[0, 
                                    np.cumsum(new_idx
                                                .reset_index()
                                                .groupby('chrom', sort=False)
                                                .count()
                                                .values
                                                .flatten())
                                                    ].tolist()
        # else:
        if remove_short_bin:
            logger.debug("Removing the short bin ...")
            retain_chroms = chromsizes[chromsizes >= binsize].index.values
            if len(retain_chroms) < len(chromsizes): 
                bins = bins.loc[retain_chroms]

                new_idx = bins['index']
                matrix = matrix[new_idx, :][:, new_idx]
                chromnames = retain_chroms.tolist()
                chrom_offset = np.r_[0, 
                                np.cumsum(new_idx
                                            .reset_index()
                                            .groupby('chrom', sort=False)
                                            .count()
                                            .values
                                            .flatten())
                                                ].tolist()
                chromsizes = chromsizes.loc[chromnames]
            else:
                chrom_offset = chrom_offset

        if balanced and remove_empty:
             if 'weight' in bins.columns:
                valid_mask = bins['weight'].notna().values
                if not np.all(valid_mask):
                    logger.info(f"Removing {np.sum(~valid_mask)} bins masked by balancing.")
                    matrix = matrix[valid_mask, :][:, valid_mask]
                    bins = bins.iloc[valid_mask]
                    
                    new_counts = bins.groupby(level=0, sort=False).size()
                    new_counts = new_counts.reindex(chromnames).fillna(0).astype(int)
                    
                    keep_chroms = new_counts[new_counts > 0].index.tolist()
                    
                    if len(keep_chroms) < len(chromnames):
                        chromnames = keep_chroms
                        new_counts = new_counts.loc[chromnames]
                        chromsizes = chromsizes.loc[chromnames]
                        
                    chrom_offset = np.r_[0, np.cumsum(new_counts.values)].tolist()


        if sort_chromosome_by_inter_contacts and len(chromnames) > 1:
            order_idx, order_labels = cluster_chrom_by_inter_contacts(
                                            matrix, bins, hap_pattern=hap_pattern,
                                      )
            logger.info(
                "Reorder homologous chromosomes by calculating the distance of inter-chromosomal contacts: \n    "
                + ",".join(order_labels)
            )
            matrix, chrom_offset = apply_chrom_order_to_matrix(matrix, chrom_offset, order_idx)
            chromnames = np.array(chromnames)[order_idx].tolist()

        chromnames_df = pd.DataFrame(chromnames, columns=['chrom'])
    
        try:
            hap_name_df = chromnames_df['chrom'].str.extract(hap_pattern, expand=True).dropna()
            hap_name_df['source_chrom'] = chromnames_df['chrom']
            hap_name_df.reset_index(drop=False, inplace=True)
 
            hap_name_df.columns = ['index', 'chrom', 'hap', 'source_chrom' ]
            source_chrom_db = dict(zip(hap_name_df['source_chrom'], hap_name_df['chrom']))

            hap_name_df = hap_name_df.groupby('chrom', sort=False)['index'].agg(lambda x: list(x)).to_frame()
            
            hap_name_df['source_chrom'] = hap_name_df.index.map(
                lambda x: [chrom for chrom, hap in source_chrom_db.items() if hap == x]
            )
            
        
            if not hap_name_df.empty:
                ## find chromosomes not in hap_name_df

                missing_chroms = set(chromnames_df['chrom'].values) - set(hap_name_df['source_chrom'].values.sum())

                for mc in missing_chroms:
                    new_idx = chromnames_df[chromnames_df['chrom'] == mc].index.values[0]
                    hap_name_df.loc[mc] = [[new_idx], mc]

                ## sort by first index
                hap_name_df = hap_name_df.sort_values(
                    by='index', 
                    key=lambda x: x.apply(lambda y: y[0])
                )

            hap_name_df = hap_name_df['index']
                
    
        except ValueError:
            hap_name_df = pd.Series()
        
        if len(hap_name_df) <= 1:
            hap_name_df = chromnames_df['chrom'].to_frame()
            hap_name_df.reset_index(drop=False, inplace=True)
            hap_name_df.columns = ['index', 'chrom']
            hap_name_df['hap'] = hap_name_df['chrom']
            hap_name_df = hap_name_df.groupby('chrom', sort=False)['index'].agg(lambda x: list(x))
            hap_names =hap_names = hap_name_df.index.values
            median_hap_count = 1
        else:
            median_hap_count = hap_name_df.apply(lambda x: len(x)).values.min()
            hap_names = hap_name_df.index.values

        matrix = matrix.todense()

        if log1p or log:
            logger.debug("Masking the zero ...")
            if sp.issparse(matrix):
                mask_nan = np.isnan(matrix.data)
                mask_inf = np.isinf(matrix.data)
                
                if mask_nan.any() or mask_inf.any():
                    try:
                        min_value = np.nanmin(matrix.data[~(mask_nan | mask_inf)])
                        matrix.data[mask_nan] = min_value
                        matrix.data[mask_inf] = min_value
                    except Exception:
                        pass
            else:
                mask = matrix == 0
                mask_nan = np.isnan(matrix)
                mask_inf = np.isinf(matrix)

                try:
                    min_value = np.nanmin(matrix[~(mask | mask_nan | mask_inf)])
                    matrix[mask] = min_value
                    matrix[mask_nan] = min_value
                    matrix[mask_inf] = min_value
                except Exception:
                    pass 

        if log1p:
            if sp.issparse(matrix):
                matrix.data = np.log10(matrix.data + 1)
            else:
                matrix += 1 
                # norm = LogNorm(vmax=vmax, vmin=vmin)
                matrix = np.log10(matrix)
            norm = "log1p"
        elif log:
            norm = "log"
            if sp.issparse(matrix):
                matrix.data = np.log(matrix.data)
            else:
                matrix = np.log(matrix)
            # norm = LogNorm(vmax=vmax, vmin=vmin)
        else:
            if balanced:
                norm = "balanced"
            else:
                norm = None

        # factor = round(matrix.shape[0] / 5000 + 1, 2)  
        # factor = np.log2(round((matrix.shape[0] * binsize) / 1e9  + 1, 2)) + 1
    
        # logger.debug(f"Figure factor is `{factor}`")
        # if triangle:
        #     fig_width = round(7 * factor, 2) 
        #     fig_height = round(3.2 * factor, 2)
        # else:
        #     if figwidth is None:
        #         fig_width = round(7 * factor, 2) 
        #     else:
        #         fig_width = figwidth 
        #     if figheight is None:
        #         fig_height = round(8 * factor, 2)
        #     else:
        #         fig_height = figheight
        if factor is None:
            if factor_formula == "bins_linear":
                factor = round(matrix.shape[0] / 5000.0 + 1, 2)
            elif factor_formula == "genome_log":
                genome_len_gb = (matrix.shape[0] * float(binsize)) / 1e9
                factor = float(np.log2(max(round(genome_len_gb + 1, 2), 1e-9)) + 1)
            else:
                factor = 1.0
        else:
            try:
                factor = float(factor)
            except Exception:
                logger.warning(f"Invalid factor `{factor}`, fallback to 1.0.")
                factor = 1.0
            if factor <= 0:
                logger.warning(f"Non-positive factor `{factor}`, fallback to 1.0.")
                factor = 1.0

        logger.debug(f"Figure factor is `{factor}` (formula={factor_formula})")

        if triangle:
            base_w, base_h = tri_size_scale
        else:
            base_w, base_h = size_scale

        scaled_w = round(base_w * factor, 2)
        scaled_h = round(base_h * factor, 2)

        fig_width = scaled_w if figwidth is None else figwidth
        fig_height = scaled_h if figheight is None else figheight

        logger.info(f"  Set figsize: (W: {fig_width}, H: {fig_height})")
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['svg.fonttype'] = 'none'
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)

        # if hap_names is None or not add_hap_border:
        #     font_factor = np.log10(len(chromnames) ) if len(chromnames) >= 10 else 1
        #     tick_fontsize = int(15 / font_factor) if not fontsize else fontsize
           
        # else:
        #     font_factor = np.log10(len(hap_names) ) if len(hap_names) > 5 else 1
        #     tick_fontsize = int(15 / font_factor) if not fontsize else fontsize
        # if tick_fontsize == np.inf:
        #     tick_fontsize = 16
        # logger.info(f"  Set fontsize to `{tick_fontsize}`.")

        if fontsize is None:
            if hap_names is None or (not add_hap_border and not add_hap_shadow):
                n_labels = len(chromnames) if chromnames is not None else 1
            else:
                n_labels = len(hap_names) if hap_names is not None else 1
            density = np.log10(n_labels) if n_labels >= 10 else 1.0

            f = float(factor) if factor is not None else 1.0
            f_scale = float(np.sqrt(max(f, 1e-6))) 
            base_font = 6.0
            tick_fontsize = int(np.clip(base_font * f_scale / max(density, 1e-6), 6, 18))
        else:
            tick_fontsize = int(fontsize)
        logger.info(f"  Set fontsize to `{tick_fontsize}`.")

        ax = plot_heatmap_core(matrix, ax, bins=bins, 
                               chromnames=chromnames, 
                               hap_name_df=hap_name_df,
                               chromsizes=chromsizes,
                            chrom_offset=chrom_offset, norm=norm,
                            triangle=triangle,
                            plot_cis_only=plot_cis_only,
                            plot_hap_only=plot_hap_only,
                            xticks=xticks, yticks=yticks, 
                            rotate_xticks=rotate_xticks, rotate_yticks=rotate_yticks,
                            vmin=vmin, vmax=vmax, 
                            cmap=cmap, cbar_bottom=cbar_bottom, add_lines=add_lines,
                            add_hap_border=add_hap_border,
                            hap_border_color=hap_border_color,
                            hap_border_width=hap_border_width,
                            add_hap_shadow=add_hap_shadow,
                            hap_shadow_color=hap_shadow_color,
                            hap_shadow_alpha=hap_shadow_alpha,
                            line_color=line_color,
                            line_width=line_width,
                            line_style=line_style,
                            tick_fontsize=tick_fontsize,
                            ytick_min_dist=ytick_min_dist,
                            avoid_overlap_yticks=avoid_overlap_yticks)
    
    else: 
        if fontsize is None:
            tick_fontsize = 12
        else:
            tick_fontsize = fontsize

        ax = plot_per_chromosome_heatmap(cool, chromosomes, fontsize=tick_fontsize,
                                         chrom_per_row=chrom_per_row, triangle=triangle,
                                        cmap=cmap, balanced=balanced,
                                        vmin=vmin, vmax=vmax, group_by_homologs=group_by_homologs,
                                        hap_pattern=hap_pattern, threads=threads)

    logger.info(f"  Set dpi to `{dpi}`.")
    plt.savefig(output, dpi=dpi, bbox_inches='tight')
    
    logger.info(f'Successful, plotted the heatmap into `{output}`')
    if output.endswith('.png'):
        plt.savefig(output.replace('.png', '.pdf'), dpi=dpi, bbox_inches='tight')
        logger.info(f"Successful, also saved the heatmap into `{output.replace('.png', '.pdf')}`")

    return ax


def plot_per_chromosome_heatmap(cool, chromosomes, log1p=True, fontsize=10,
                                    chrom_per_row=4, remove_short_bin=True, triangle=False,
                                    cmap='redp1_r', balanced=False, vmin=None, vmax=None,
                                    group_by_homologs=False, hap_pattern=r'(Chr\d+)g(\d+)',
                                    threads=1):
    """
    modified from hicPlotMatrix plotPerChr
    """
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from matplotlib.colors import LogNorm
    

    if balanced:
        log1p = False 
        log = True

    chrom_offset = cool._load_dset('indexes/chrom_offset')
    chromsizes = cool.chromsizes
    if chromsizes.dtype == np.int32:
        chromsizes = chromsizes.astype(np.uint32).astype(np.int64)

    binsize = cool.binsize

    bins = cool.bins()[:].reset_index(drop=False)
    if bins['start'].dtype == np.int32:
        bins['start'] = bins['start'].astype(np.uint32).astype(np.int64)
        bins['end'] = bins['end'].astype(np.uint32).astype(np.int64)
        
    bins['chrom'] = bins['chrom'].astype('str')
    bins = bins.set_index('chrom')
    if binsize is None:
        binsize = np.argmax(np.bincount(bins['end'] - bins['start']))

    matrix = cool.matrix(balance=balanced, sparse=True)[:].tocsr()

    if chromosomes:
        chromsizes = chromsizes.loc[chromosomes]
        bins = bins.loc[chromosomes]
        new_idx = bins['index']
        matrix = matrix[new_idx, :][:, new_idx]
        chrom_offset = np.r_[0, 
                            np.cumsum(new_idx
                                        .reset_index()
                                        .groupby('chrom', sort=False)
                                        .count()
                                        .values
                                        .flatten())
                                        ].tolist()
    else:
        if remove_short_bin:
            retain_chroms = chromsizes[chromsizes >= binsize].index.values
            bins = bins.loc[retain_chroms]
            new_idx = bins.loc[retain_chroms]['index']
            matrix = matrix[new_idx, :][:, new_idx]
            chromosomes = retain_chroms.tolist()
            grouped_counts = new_idx.reset_index().groupby('chrom', sort=False).count().values.flatten()
            chrom_offset = np.r_[0, np.cumsum(grouped_counts)].tolist()
            chromsizes = chromsizes.loc[chromosomes]

    matrix = matrix.todense()
    chromosomes = chromosomes if chromosomes else cool.chromnames


    if group_by_homologs:

        def _compute_chrom_offset(bins_df, chrom_list):
            b = bins_df.reset_index() if 'chrom' not in bins_df.columns else bins_df.reset_index()
            if 'chrom' not in b.columns:
                b.columns.values[0] = 'chrom'
            counts = b.groupby('chrom', sort=False)['index'].count() if 'index' in b.columns else b.groupby('chrom', sort=False).size()
            counts = counts.reindex(chrom_list).fillna(0).astype(int)
            return np.r_[0, np.cumsum(counts.values)].tolist()
        
        pat = re.compile(hap_pattern)
        groups = OrderedDict()
        order_seen = []
        for ch in chromosomes:
            m = pat.match(str(ch))
            key = m.group(1) if m else str(ch)
            if key not in groups:
                groups[key] = []
                order_seen.append(key)
            groups[key].append(ch)

        ordered_chroms = [c for key in order_seen for c in groups[key]]

        try:
            chromsizes = chromsizes.loc[ordered_chroms]
            bins = bins.loc[ordered_chroms]
            chromosomes = ordered_chroms
            chrom_offset = _compute_chrom_offset(bins, chromosomes)

        except Exception:
            groups = OrderedDict()
            groups_list = [[c] for c in chromosomes]
            groups = OrderedDict((str(i), g) for i, g in enumerate(groups_list))

        num_rows = len(groups)
        num_cols = max(len(g) for g in groups.values()) if groups else 0
        width_ratios = [1.0] * num_cols + [0.05]
    
        grid_chrom = []
        for g in groups.values():
            row = list(g) + [None] * (num_cols - len(g))
            grid_chrom.append(row)
    else:
        num_rows = int(ceil(float(len(chromosomes)) / chrom_per_row))
        num_cols = min(chrom_per_row, len(chromosomes))
        width_ratios = [1.0] * num_cols + [0.05]


    if triangle:
        grids = gridspec.GridSpec(num_rows, num_cols + 1,
                                width_ratios=width_ratios,
                                height_ratios=[1] * num_rows,
                                hspace=0.4, wspace=0.2)
     
        fig_height = 1.9 * num_rows
        fig_width = 3.5 * num_cols
        fig = plt.figure(figsize=(fig_width, fig_height))

    else:
        grids = gridspec.GridSpec(num_rows, num_cols + 1,
                                width_ratios=width_ratios,
                                height_ratios=[1] * num_rows,
                                hspace=0.4, wspace=0.2)
        fig_height = 3 * num_rows
        fig_width = sum((np.array(width_ratios) + 0.05) * 3)
        
        fig = plt.figure(figsize=(fig_width, fig_height))

    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype'] = 42

    def _plot(ax, chrom, matrix, chrom_range, plot_cbar=False, fontsize=12):
        chrom_matrix = matrix[chrom_range[0]:chrom_range[1], :][:, chrom_range[0]:chrom_range[1]]
        chrom_matrix = chrom_matrix
        if log1p or log:
            mask = chrom_matrix == 0
            mask_nan = np.isnan(chrom_matrix)
            mask_inf = np.isinf(chrom_matrix)

            try:
                min_value = np.nanmin(matrix[~(mask | mask_nan | mask_inf)])
                matrix[mask] = min_value
                matrix[mask_nan] = min_value
                matrix[mask_inf] = min_value
            except Exception:
                pass 

        if log1p:
            chrom_matrix += 1
            chrom_matrix = np.log10(chrom_matrix)
            norm = "log1p"
            # norm = LogNorm()
        elif log:
            norm = "log"
            chrom_matrix = np.log(chrom_matrix)
            # norm = LogNorm()
        else:
            norm = None


        plot_heatmap_core(chrom_matrix, ax, tick_fontsize=fontsize,
                          chromsizes=chromsizes.loc[chrom], bins=bins, chrom_offset=chrom_offset,
                          norm=norm, xlabel=chrom, triangle=triangle,
                            cmap=cmap, xticks=True, yticks=False,
                            vmin=vmin, vmax=vmax, plot_cbar=plot_cbar)
    
    if vmax is None:
        plot_cbar = True
    else:
        plot_cbar = False

    args = []
    if group_by_homologs:
        index_map = {ch: idx for idx, ch in enumerate(chromosomes)}
        for r in range(num_rows):
            for c in range(num_cols):
                chrom = grid_chrom[r][c]
              
                if chrom is None:
                    continue
                i = index_map[chrom]
                ax = plt.subplot(grids[r, c])
                
                chrom_range = (chrom_offset[i], chrom_offset[i+1])
                args.append((ax, chrom, matrix, chrom_range))
                logger.debug(f"Plotting the heatmap of `{chrom}` (group row {r}, col {c}) ...")
                _plot(ax, chrom, matrix, chrom_range, plot_cbar=plot_cbar, fontsize=fontsize)
    else:
        for i, chrom in enumerate(chromosomes):
            row = i // chrom_per_row
            col = i % chrom_per_row
            ax = plt.subplot(grids[row, col])
            chrom_range = (chrom_offset[i], chrom_offset[i+1])
            args.append((ax, chrom, matrix, chrom_range))
            logger.debug(f"Plotting the heatmap of `{chrom}` ...")
            _plot(ax, chrom, matrix, chrom_range, plot_cbar=plot_cbar, fontsize=fontsize)

    if plot_cbar is False:
        try:
            import colormaps as cmaps
            from matplotlib import colors as mcolors
        except Exception:
            from matplotlib import colors as mcolors
            cmaps = None

        if cmap.endswith('_half'):
            colormap = half_colormaps(cmap.replace('_half', ''))    
        else:
            if cmap == 'whitered':
                colormap = whitered_cmap
            else:
                try:
                    colormap = getattr(cmaps, cmap)
                except Exception:
                    try:
                        colormap = getattr(plt.cm, cmap)
                    except Exception:
                        logger.warning(f"Colormap `{cmap}` not found, use `redp1_r` instead.")
                        colormap = getattr(cmaps, "redp1_r") if cmaps is not None else plt.cm.get_cmap('Reds')

        img = None
        for ax in fig.axes:
            ims = ax.get_images()
            if ims:
                img = ims[0]
                break

        from matplotlib import colors as mcolors

        if img is not None:
            cmap_getter = getattr(img, 'get_cmap', None)
            try:
                cmap_used = cmap_getter() if callable(cmap_getter) else (cmap_getter if cmap_getter is not None else None)
            except Exception:
                cmap_used = None
            if cmap_used is None:
                try:
                    cmap_used = getattr(plt.cm, cmap)
                except Exception:
                    cmap_used = plt.cm.viridis

            norm_used = None
            get_norm = getattr(img, 'get_norm', None)
            if callable(get_norm):
                try:
                    norm_used = get_norm()
                except Exception:
                    norm_used = None

            if norm_used is None and hasattr(img, 'norm'):
                try:
                    norm_used = getattr(img, 'norm')
                except Exception:
                    norm_used = None

            if norm_used is None:
                try:
                    vmin_img, vmax_img = img.get_clim()
                except Exception:
                    vmin_img, vmax_img = (None, None)
                vmin_final = vmin if vmin is not None else (vmin_img if vmin_img is not None else 0.0)
                vmax_final = vmax if vmax is not None else (vmax_img if vmax_img is not None else 1.0)
                norm_used = mcolors.Normalize(vmin=vmin_final, vmax=vmax_final)
            else:
                try:
                    if vmin is not None:
                        setattr(norm_used, 'vmin', vmin)
                    if vmax is not None: 
                        setattr(norm_used, 'vmax', vmax)
                except Exception:
                    pass

            sm = plt.cm.ScalarMappable(cmap=cmap_used, norm=norm_used)
            v0 = getattr(norm_used, 'vmin', (vmin if vmin is not None else 0.0))
            v1 = getattr(norm_used, 'vmax', (vmax if vmax is not None else 1.0))
            sm.set_array(np.linspace(v0, v1, 256))
        else:
            norm_obj = mcolors.Normalize(vmin=vmin if vmin is not None else 0.0,
                                        vmax=vmax if vmax is not None else 1.0)
            sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm_obj)
            sm.set_array(np.linspace(norm_obj.vmin, norm_obj.vmax, 256))
       
        from matplotlib.transforms import Bbox
        try:
            ax_positions = [a.get_position() for a in fig.axes if a is not None and a.get_axes_locator() is None]
        except Exception:
            ax_positions = [a.get_position() for a in fig.axes]
        if ax_positions:
            bbox = Bbox.union(ax_positions)
            pad = 0.01  
            cbar_width = 0.015  
            cbar_height = 0.18  
            x = bbox.x1 + pad
            y = bbox.y1 - cbar_height - pad
           
            if x + cbar_width > 0.98:
                x = bbox.x1 - cbar_width - pad
            if y < 0.02:
                y = bbox.y0 + pad
            cax = fig.add_axes([x, y, cbar_width, cbar_height])
            cbar = fig.colorbar(sm, cax=cax, orientation='vertical')
            cbar.ax.tick_params(labelsize=fontsize)
            try:
                cbar.set_label("Contacts", fontsize=fontsize)
            except Exception:
                pass
        else:
            cax = plt.subplot(grids[:, -1])
            cbar = fig.colorbar(sm, cax=cax, shrink=0.4, pad=0.3)



        if log1p:
            cbar.set_label("Log$_{10}$(Contact + 1)")
        elif log:
            cbar.set_label("Log(Contact)")
        elif balanced:
            cbar.set_label("Normalized Contacts")
        else:
            cbar.set_label("Contacts")

    return fig

def plot_heatmap_core(matrix, 
                        ax,
                        bins=None,
                        chromnames=None,
                        hap_name_df=None,
                        chromsizes=None,
                        chrom_offset=None,
                        norm=None, vmin=None, vmax=None,
                        triangle=False,
                        plot_cis_only=False,
                        plot_hap_only=False,
                        xlabel=None, ylabel=None, 
                        xticks=True, yticks=True,
                        rotate_xticks=False, rotate_yticks=False,
                        avoid_overlap_yticks = True,
                        ytick_min_dist=20,
                        tick_fontsize=16,
                        cmap="redp1_r",
                        add_hap_border=False,
                        hap_border_color='black',
                        hap_border_width=1.0,
                        add_hap_shadow=False,
                        hap_shadow_color='grey',
                        hap_shadow_alpha=0.1,
                        add_lines=False,
                        line_color='black',
                        line_width=0.5,
                        line_style='--',
                        plot_cbar=True,
                        cbar_bottom=False):
    import colormaps as cmaps
    import matplotlib.pyplot as plt 
    import matplotlib
    matplotlib.use('Agg')       
    from matplotlib import colors
    from matplotlib.patches import Rectangle
    import seaborn as sns 

    if cmap.endswith('_half'):
        colormap = half_colormaps(cmap.replace('_half', ''))    
    else:
        if cmap == 'whitered':
            colormap = whitered_cmap
        else:
            try:
                """
                https://pratiman-91.github.io/colormaps/
                """
                colormap = getattr(cmaps, cmap)
            except:
                try:
                    colormap = cmap 
                    colormap = getattr(plt.cm, cmap)
                except AttributeError:
                    logger.warning(f"Colormap `{cmap}` not found, use `redp1_r` instead.")
                    colormap = getattr(cmaps, "redp1_r")

    
    if plot_cis_only:
        new_matrix = np.zeros(matrix.shape)
        for i in range(len(chrom_offset) - 1):
            new_matrix[chrom_offset[i]:chrom_offset[i+1], 
                                            chrom_offset[i]:chrom_offset[i+1]] = \
                                                matrix[chrom_offset[i]:chrom_offset[i+1], 
                                            chrom_offset[i]:chrom_offset[i+1]]
        matrix = new_matrix

    if plot_hap_only and hap_name_df is not None:
        new_matrix = np.zeros(matrix.shape)
        for chrom, idx_list in hap_name_df.items():
            idx_list = list(idx_list)
            if not idx_list:
                continue
            start = chrom_offset[idx_list[0]]
            
            end = chrom_offset[idx_list[-1] + 1]
            new_matrix[start:end, start:end] = matrix[start:end, start:end]
        matrix = new_matrix
        
    
    if norm is None or norm == "balanced":
        if vmax is None:
            cis_matrix = []
            if sp.issparse(matrix):
                for i in range(len(chrom_offset) - 1):
                    sub_mat = matrix[chrom_offset[i]:chrom_offset[i+1], 
                                    chrom_offset[i]:chrom_offset[i+1]]
                    sub_triu = sp.triu(sub_mat, k=1)
                    cis_matrix.append(sub_triu.data)

                cis_matrix = np.concatenate(cis_matrix)
                cis_matrix = cis_matrix[cis_matrix > 0] 
            else:
                for i in range(len(chrom_offset) - 1):
                    cis_matrix.append(
                                    np.array(matrix[chrom_offset[i]:chrom_offset[i+1], 
                                            chrom_offset[i]:chrom_offset[i+1]]))
                
                cis_matrix = [arr[np.triu_indices(arr.shape[0], k=1)] for arr in cis_matrix]
                cis_matrix = np.concatenate([arr.flatten().ravel() for arr in cis_matrix])
                cis_matrix = cis_matrix[~np.isnan(cis_matrix)]
                cis_matrix = cis_matrix[cis_matrix > 0] 
            
            def vmax_by_iqr(vals, k):
                q1, q3 = np.percentile(vals, [25, 75])
                iqr = q3 - q1
                return float(q3 + k * iqr)

            def vmax_by_mad(vals, k):
                med = np.median(vals)
                mad = np.median(np.abs(vals - med)) * 1.4826
                return float(med + k * mad)
         
            # vmax = vmax_by_mad(cis_matrix, 5)
            vmax = vmax_by_iqr(cis_matrix, 1.5)
          
            if vmax == 0:
                try:
                    vmax = np.min(cis_matrix[cis_matrix != 0])
                
                    tiny_vmax_thr = 1e-20
                    if vmax < 1e-20:
                        if (norm is None or norm == "balanced") and (vmax is not None) and np.isfinite(vmax) and (vmax > 0) and (vmax < tiny_vmax_thr):
                            with np.errstate(divide='ignore', invalid='ignore'):
                                matrix = np.asarray(matrix, dtype=float) / max(vmax, 1e-12)
                            vmin = 0.0
                            vmax = 1.0
                            logger.info("  Rescaled by tiny vmax; set vmin=0, vmax=1 for colorbar.")
                except:
                    vmax = 0
       
            try:
                if isinstance(vmax, (int, np.integer)):
                    s_vmax = str(int(vmax))
                else:
                    s_vmax = format_float_dynamic(vmax, sig=4, min_dec=0, max_dec=6)
                logger.info(f"  Set vmax to `{s_vmax}`.")
            except Exception:
                logger.info(f"  Set vmax to `{vmax}`.")
           
        

            if norm is None and vmin is None:
                vmin = 0
            del cis_matrix
    else:
        # if norm == 'log1p':
        #     if vmax is None:
        #         vmax = np.ceil(np.max(matrix))
        #         logger.info(f"  Set vmax to `{vmax}`.")

        # elif norm == 'log':
        #     if vmax is None:
        #         vmax = np.ceil(np.max(matrix))
        #         logger.info(f"  Set vmax to `{vmax}`.")

        if norm == 'log1p' or norm == 'log':
            if vmax is None:
                cis_matrix = []
                if sp.issparse(matrix):
                    for i in range(len(chrom_offset) - 1):
                        sub_mat = matrix[chrom_offset[i]:chrom_offset[i+1], 
                                        chrom_offset[i]:chrom_offset[i+1]]
                        sub_triu = sp.triu(sub_mat, k=1)
                        cis_matrix.append(sub_triu.data)

                    cis_matrix = np.concatenate(cis_matrix)
                
                    vmax = np.percentile(cis_matrix, 99) if matrix.nnz > 0 else 1.0
                else:
                    for i in range(len(chrom_offset) - 1):
                       
                        sub_mat = np.array(matrix[chrom_offset[i]:chrom_offset[i+1], 
                                            chrom_offset[i]:chrom_offset[i+1]])
                      
                        rows, cols = np.triu_indices(sub_mat.shape[0], k=1)
                        cis_matrix.append(sub_mat[rows, cols])
                        
                    cis_matrix = np.concatenate(cis_matrix)
                    vmax = np.percentile(np.asarray(cis_matrix), 99)
                
                vmax = np.ceil(vmax)
                logger.info(f"  Set vmax to `{vmax}`.")
    
    if triangle:
        ## https://github.com/XiaoTaoWang/NeoLoopFinder/blob/master/neoloop/visualize/core.py
        # fig = plt.figure(figsize=(7, 4.2))
        # grid = GridSpec(1, 1, figure=fig, left=0.1, right=0.9,
        #             bottom=0.1, top=0.9, hspace=0.04, height_ratios=[5])
        # ax = fig.add_subplot(grid[0])
        if sp.issparse(matrix):
            matrix = matrix.todense()

        add_lines = False
        # fig, ax = plt.subplots(figsize=(7, 3.2))
        n = matrix.shape[0]
        t = np.array([[1, 0.5], [-1, 0.5]])
        A = np.dot(np.array([(i[1], i[0]) for i in product(range(n, -1, -1),range(0, n+1, 1))]), t)

        x = A[:, 1].reshape(n+1, n+1)
        y = A[:, 0].reshape(n+1, n+1)
        y[y<0] = -y[y<0]
    
        cax = ax.pcolormesh(x, y, np.flipud(np.array(matrix)), cmap=colormap, edgecolor='none',
                            snap=True, linewidth=.001, rasterized=True)
        fig = plt.gcf()
        # ax.plot(100, 3, marker='v', markersize=2, linestyle=None, color='blue')
        # ax.axis('off')
        # ax.xaxis.set_visible(False)
        # ax.yaxis.set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
    else:
        if sp.issparse(matrix):
            matrix = matrix.todense()
        fig = plt.gcf()
        cax = ax.imshow(matrix, cmap=colormap, aspect='equal',
                            interpolation=None, rasterized=True)
      

        
    binsize = np.argmax(np.bincount(bins['end'] - bins['start']))
    
    if plot_cbar:
        cax.set_norm(colors.Normalize(vmin=vmin, vmax=vmax))

        ax.set_xlim(0, matrix.shape[0])
        ax.set_ylim(0, matrix.shape[1])

        if cbar_bottom:
            cbar = fig.colorbar(cax, ax=ax, orientation='horizontal', shrink=.3, pad=0.06, location='bottom')
        else:
            cbar = fig.colorbar(cax, ax=ax, shrink=.4, pad=0.03)
        if tick_fontsize is not None:
            cbar.ax.tick_params(labelsize=tick_fontsize*0.8)
        else:
            cbar.ax.tick_params(labelsize=12)

        cbar.locator = plt.MaxNLocator(3)
        fmt = ScalarFormatter(useMathText=True)
        fmt.set_scientific(True)
        fmt.set_powerlimits((-3, 3))
        fmt.set_useOffset(False)
        cbar.formatter = fmt 
        if cbar.orientation == 'vertical':
            cbar.ax.yaxis.set_major_formatter(fmt)
            cbar.ax.ticklabel_format(style='sci', axis='y', scilimits=(-3, 3))
            try:
                cbar.ax.yaxis.set_offset_position('right')
                cbar.ax.yaxis.get_offset_text().set_fontsize(int(tick_fontsize * 0.75))
            except Exception:
                pass
        else:
            cbar.ax.xaxis.set_major_formatter(fmt)
            cbar.ax.ticklabel_format(style='sci', axis='x', scilimits=(-3, 3))
            try:
                cbar.ax.xaxis.set_offset_position('top')
                cbar.ax.xaxis.get_offset_text().set_fontsize(int(tick_fontsize * 0.75))
            except Exception:
                pass
        cbar.update_ticks()
        if norm == 'log1p':
            cbar.set_label("Log$_{10}$(Contact + 1)", fontsize=tick_fontsize)
        elif norm == 'log':
            cbar.set_label("Log(Contact)", fontsize=tick_fontsize)
        elif norm == 'balanced':
            cbar.set_label("Normalized Contacts", fontsize=tick_fontsize)
        else:
            cbar.set_label("Contacts", fontsize=tick_fontsize)
    else:
        cbar = None

    if add_lines and chrom_offset is not None:
        ax.hlines(np.array(chrom_offset[1:-1]) - 0.5, *ax.get_xlim(), 
                    linewidth=line_width, color=line_color, linestyles=line_style)
        ax.vlines(np.array(chrom_offset[1:-1]) - 0.5, *ax.get_ylim(), 
                    linewidth=line_width, color=line_color, linestyles=line_style)
    
    if chrom_offset is not None:
        mid_tick_pos = list((np.array(chrom_offset)[:-1] + np.array(chrom_offset)[1:]) / 2)
    else:
        mid_tick_pos = [0]

    ax.tick_params(width=0)
    if xticks and chromnames:
        rotation = "horizontal" if rotate_xticks else "vertical" 
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
        ax.tick_params(axis='x', length=5, width=1)
        x_limit = ax.get_xlim()[1]
        _xticks = ax.get_xticks()
        _xticks = _xticks[_xticks < x_limit]
        
        if len(_xticks) > 1:
            dist = _xticks[1] - _xticks[0]
            if len(_xticks) > 0 and (x_limit - _xticks[-1]) <= (dist / 2):
                _xticks = _xticks[:-1]

        _xticks = np.r_[_xticks, x_limit]

        dist = _xticks[1] - _xticks[0]

        ax.set_xticks(_xticks)
        _xticks = _xticks * binsize 
        total_size = chromsizes.sum()
        _xticks[-1] = total_size
        xticklabels = chrom_ticks_convert(_xticks)

        ax.set_xticklabels(xticklabels, fontsize=tick_fontsize, rotation=rotation)

    else:
        ax.set_xticks([])
    if chromnames:
        if hap_name_df is not None and (add_hap_border or add_hap_shadow):
            add_hap_internal_lines = True 
            hap_internal_mode = 'grid'
            new_mid_tick_pos = []
            hap_names = []
            new_chrom_offset = []
            for chrom, idx in hap_name_df.items():
                hap_names.append(chrom)
                start, end = chrom_offset[idx[0]], chrom_offset[idx[-1] + 1]
                new_mid_tick_pos.append((end - start) / 2 + start - 0.5)
                new_chrom_offset.append(start)
            else:
                new_chrom_offset.append(chrom_offset[-1])
            
            if yticks:
                ax.set_yticks(new_mid_tick_pos)
                rotation = "vertical" if rotate_yticks else "horizontal"
                ax.set_yticklabels(hap_names, fontsize=tick_fontsize, rotation=rotation)
            else:
                ax.set_yticks([])
            

            for start, end in zip(new_chrom_offset[:-1], new_chrom_offset[1:]):
                if add_hap_border:
                    plt.plot([start, start], [start, end], color=hap_border_color, linestyle='-', linewidth=hap_border_width)
                    plt.plot([start, end], [start, start], color=hap_border_color, linestyle='-', linewidth=hap_border_width)
                    plt.plot([end, end], [start, end], color=hap_border_color, linestyle='-', linewidth=hap_border_width)
                    plt.plot([start, end], [end, end], color=hap_border_color, linestyle='-', linewidth=hap_border_width)

                if add_hap_shadow:
                    rect = Rectangle(
                            xy=(start, start), 
                            width=end - start, 
                            height=end - start, 
                            facecolor=hap_shadow_color, 
                            alpha=hap_shadow_alpha,        
                            edgecolor=None,
                        )
                    ax.add_patch(rect)
            if add_hap_internal_lines and not add_lines:
                small_lw = line_width
                small_color = line_color
                if hap_internal_mode == "box":
                    for chrom, idx in hap_name_df.items():
                        local_offsets = [chrom_offset[i] for i in idx] + [chrom_offset[idx[-1] + 1]]
                        for s, e in zip(local_offsets[:-1], local_offsets[1:]):
                            plt.plot([s, s], [s, e], color=small_color, linestyle='-', linewidth=small_lw)
                            plt.plot([s, e], [s, s], color=small_color, linestyle='-', linewidth=small_lw)
                            plt.plot([e, e], [s, e], color=small_color, linestyle='-', linewidth=small_lw)
                            plt.plot([s, e], [e, e], color=small_color, linestyle='-', linewidth=small_lw)
                elif hap_internal_mode == "grid":
                    for chrom, idx in hap_name_df.items():
                        group_start = chrom_offset[idx[0]]
                        group_end = chrom_offset[idx[-1] + 1]
                        boundaries = [chrom_offset[i] for i in idx] + [group_end]
                        for b in boundaries[1:-1]:
                            plt.plot([b, b], [group_start, group_end], color=small_color, linestyle=line_style, linewidth=small_lw)
                            plt.plot([group_start, group_end], [b, b], color=small_color, linestyle=line_style, linewidth=small_lw)

        else:
            if yticks:
                ax.set_yticks(mid_tick_pos)
                rotation = "vertical" if rotate_yticks else "horizontal"
                ax.set_yticklabels(chromnames, fontsize=tick_fontsize, rotation=rotation)
            else:
                ax.set_yticks([])
    else:
        ax.set_yticks([])


    if yticks and avoid_overlap_yticks:
        fig = ax.figure
        try:
            fig.canvas.draw()
        except Exception:
            pass
        kept_last_y = None
        for lab in ax.get_yticklabels():
            if not lab.get_visible():
                continue
            if not lab.get_text():
                continue
            y_data = lab.get_position()[1]
            pixel = ax.transData.transform((0, y_data))
            y_px = pixel[1]
            if kept_last_y is None or abs(y_px - kept_last_y) >= ytick_min_dist:
                kept_last_y = y_px
            else:
                lab.set_visible(False)

    if xlabel:
        ax.set_xlabel(xlabel, fontsize=18, labelpad=15)
    
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=18, labelpad=15)

    if not triangle:
        sns.despine(top=False, right=False)


    if triangle:
        ax.axis([0, matrix.shape[0], 0, matrix.shape[1]])
    else:
        ax.axis([0, matrix.shape[0], matrix.shape[1], 0])

    return ax


def sort_chromosome_by_interaction(cool, size_penalty_coef: float = 0.5, collapse_small_threshold_factor: float = 0.5):
    import numpy as np
    from scipy.cluster.hierarchy import linkage, optimal_leaf_ordering, leaves_list
    from scipy.spatial.distance import squareform

    c = cooler.Cooler(cool)
    chromnames = list(c.chromnames)
    chrom_offset = np.asarray(c._load_dset('indexes/chrom_offset'), dtype=np.int64)
    n_chr = len(chromnames)
    if n_chr <= 1:
        return np.arange(n_chr, dtype=int), chromnames

    pixels = c.pixels()
    nnz = len(pixels)
    if nnz == 0:
        return np.arange(n_chr, dtype=int), chromnames

    try:
        chrom_sizes = np.asarray([float(c.chromsizes[name]) for name in chromnames], dtype=float)
    except Exception:
        chrom_sizes = None

    collapsed_map = None
    keep_idx = None
    new_n_chr = n_chr
    if chrom_sizes is not None and n_chr > 200:
        thr = max(np.median(chrom_sizes) * collapse_small_threshold_factor, 1e4)
        small_idx = np.where(chrom_sizes < thr)[0]
        if 0 < small_idx.size < n_chr:
            small_set = set(small_idx.tolist())
            keep_idx = [i for i in range(n_chr) if i not in small_set]
            other_bucket = len(keep_idx)
            new_n_chr = other_bucket + 1
            collapsed_map = np.full(n_chr, other_bucket, dtype=np.int32)
            for newi, oldi in enumerate(keep_idx):
                collapsed_map[oldi] = newi
        else:
            collapsed_map = None
            new_n_chr = n_chr
    else:
        collapsed_map = None
        new_n_chr = n_chr

    S = np.zeros((new_n_chr, new_n_chr), dtype=float)

    chunksize = max(1_000_000, nnz // 50)  
    for i in range(0, nnz, chunksize):
        chunk = pixels[i: i + chunksize]
        b1 = np.asarray(chunk['bin1_id'], dtype=np.int64)
        b2 = np.asarray(chunk['bin2_id'], dtype=np.int64)
        cnt = np.asarray(chunk['count'], dtype=float)

        chr1 = np.searchsorted(chrom_offset, b1, side='right') - 1
        chr2 = np.searchsorted(chrom_offset, b2, side='right') - 1

        mask = chr1 != chr2
        if not np.any(mask):
            continue
        i1 = chr1[mask].astype(np.int64)
        i2 = chr2[mask].astype(np.int64)
        w = cnt[mask].astype(float)

        if collapsed_map is not None:
            i1 = collapsed_map[i1]
            i2 = collapsed_map[i2]
            valid = (i1 >= 0) & (i2 >= 0)
            if not np.any(valid):
                continue
            i1 = i1[valid].astype(np.int64)
            i2 = i2[valid].astype(np.int64)
            w = w[valid]

        flat_idx = i1 * new_n_chr + i2
        flat = np.bincount(flat_idx, weights=w, minlength=new_n_chr * new_n_chr)
        S += flat.reshape((new_n_chr, new_n_chr))

    S = S + S.T

    if np.all(S == 0):
        if collapsed_map is not None and keep_idx is not None:
            labels = [chromnames[i] for i in keep_idx] + ["small_contigs"]
            return np.arange(new_n_chr, dtype=int), labels
        return np.arange(n_chr, dtype=int), chromnames

    if collapsed_map is not None and keep_idx is not None:
        chromnames_new = [chromnames[i] for i in keep_idx] + ["small_contigs"]
        chromnames = chromnames_new
        n_chr = new_n_chr
    else:
        n_chr = new_n_chr

    S_max = float(np.nanmax(S))
    D = (S_max - S).astype(float)
    np.fill_diagonal(D, 0.0)

    try:
        if chrom_sizes is not None and size_penalty_coef > 0.0:
            if collapsed_map is not None and keep_idx is not None:
                keep_sizes = [chrom_sizes[i] for i in keep_idx]
                small_sizes = chrom_sizes[np.where(collapsed_map == collapsed_map[keep_idx[0]])] if False else chrom_sizes[np.isin(np.arange(len(chrom_sizes)), np.where(collapsed_map == (len(keep_idx)))[0])]
                small_sizes = chrom_sizes[[i for i in range(len(chrom_sizes)) if i not in keep_idx]] if len(keep_idx) < len(chrom_sizes) else np.array([1.0])
                new_sizes = np.concatenate([np.array(keep_sizes, dtype=float), [float(np.median(small_sizes)) if small_sizes.size else 1.0]])
                log_sizes = np.log1p(new_sizes)
            else:
                log_sizes = np.log1p(chrom_sizes)
            size_diff = np.abs(log_sizes[:, None] - log_sizes[None, :])
            max_diff = float(np.nanmax(size_diff)) if size_diff.size else 0.0
            if max_diff > 0:
                size_penalty_norm = size_diff / max_diff
                D = D + (size_penalty_coef * S_max) * size_penalty_norm
    except Exception:
        pass

    try:
        dist_vec = squareform(D, checks=False)
    except Exception:
        jitter = 1e-12 * np.random.RandomState(0).rand(D.shape[0], D.shape[1])
        D = D + jitter
        dist_vec = squareform(D, checks=False)

    Z = linkage(dist_vec, method='average')
    try:
        Z_opt = optimal_leaf_ordering(Z, dist_vec)
        leaves = leaves_list(Z_opt).astype(int)
    except Exception:
        leaves = leaves_list(Z).astype(int)

    order_idx = np.asarray(leaves, dtype=int)
    order_labels = [chromnames[i] for i in order_idx]

    return order_idx, order_labels


def sort_chromosome_by_interaction_fast(
    cool,
    method: str = "greedy",
    size_penalty_coef: float = 0.5,
    collapse_small_threshold_factor: float = 0.5,
    collapse_min_chrs: int = 200,
):
    import numpy as np
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import minimum_spanning_tree

    c = cooler.Cooler(cool)
    chromnames = list(c.chromnames)
    chrom_offset = np.asarray(c._load_dset("indexes/chrom_offset"), dtype=np.int64)
    n_chr = len(chromnames)
    if n_chr <= 1:
        return np.arange(n_chr, dtype=int), chromnames

    pixels = c.pixels()
    nnz = len(pixels)
    if nnz == 0:
        return np.arange(n_chr, dtype=int), chromnames

    # optional: collapse many small contigs -> one bucket
    try:
        chrom_sizes = np.asarray([float(c.chromsizes[name]) for name in chromnames], dtype=float)
    except Exception:
        chrom_sizes = None

    collapsed_map = None
    keep_idx = None
    new_n_chr = n_chr
    if chrom_sizes is not None and n_chr >= collapse_min_chrs:
        thr = max(np.median(chrom_sizes) * collapse_small_threshold_factor, 1e4)
        small_idx = np.where(chrom_sizes < thr)[0]
        if 0 < small_idx.size < n_chr:
            small_set = set(small_idx.tolist())
            keep_idx = [i for i in range(n_chr) if i not in small_set]
            other_bucket = len(keep_idx)
            new_n_chr = other_bucket + 1
            collapsed_map = np.full(n_chr, other_bucket, dtype=np.int32)
            for newi, oldi in enumerate(keep_idx):
                collapsed_map[oldi] = newi

    S = np.zeros((new_n_chr, new_n_chr), dtype=float)
    chunksize = max(1_000_000, nnz // 50)
    for i in range(0, nnz, chunksize):
        chunk = pixels[i : i + chunksize]
        b1 = np.asarray(chunk["bin1_id"], dtype=np.int64)
        b2 = np.asarray(chunk["bin2_id"], dtype=np.int64)
        cnt = np.asarray(chunk["count"], dtype=float)

        chr1 = np.searchsorted(chrom_offset, b1, side="right") - 1
        chr2 = np.searchsorted(chrom_offset, b2, side="right") - 1
        mask = chr1 != chr2
        if not np.any(mask):
            continue
        i1 = chr1[mask].astype(np.int64)
        i2 = chr2[mask].astype(np.int64)
        w = cnt[mask].astype(float)

        if collapsed_map is not None:
            i1 = collapsed_map[i1]
            i2 = collapsed_map[i2]
            valid = (i1 >= 0) & (i2 >= 0)
            if not np.any(valid):
                continue
            i1 = i1[valid].astype(np.int64)
            i2 = i2[valid].astype(np.int64)
            w = w[valid]

        flat_idx = i1 * new_n_chr + i2
        flat = np.bincount(flat_idx, weights=w, minlength=new_n_chr * new_n_chr)
        S += flat.reshape((new_n_chr, new_n_chr))

    S = S + S.T
    if np.all(S == 0):
        # no inter-chromosomal signal
        if collapsed_map is not None and keep_idx is not None:
            labels = [chromnames[i] for i in keep_idx] + ["small_contigs"]
            return np.arange(new_n_chr, dtype=int), labels
        return np.arange(n_chr, dtype=int), chromnames

    # optionally incorporate size penalty into a modified score matrix (only used by greedy/mst)
    if chrom_sizes is not None and size_penalty_coef > 0.0:
        if collapsed_map is not None and keep_idx is not None:
            keep_sizes = np.array([chrom_sizes[i] for i in keep_idx], dtype=float)
            small_sizes = chrom_sizes[[i for i in range(len(chrom_sizes)) if i not in keep_idx]]
            small_val = float(np.median(small_sizes)) if small_sizes.size else 1.0
            sizes = np.concatenate([keep_sizes, [small_val]])
        else:
            sizes = chrom_sizes.astype(float)
        log_sizes = np.log1p(sizes)
        size_diff = np.abs(log_sizes[:, None] - log_sizes[None, :])
        # normalize and scale relative to S max
        max_diff = float(np.nanmax(size_diff)) if size_diff.size else 0.0
        if max_diff > 0:
            size_penalty = (size_diff / max_diff) * (size_penalty_coef * float(np.nanmax(S)))
        else:
            size_penalty = np.zeros_like(S)
    else:
        size_penalty = np.zeros_like(S)

    # fast strategies
    if method == "sum":
        scores = S.sum(axis=1)
        order_idx = np.argsort(-scores)
    elif method == "size_then_sum" and chrom_sizes is not None:
        if collapsed_map is not None and keep_idx is not None:
            sizes = np.array([float(chrom_sizes[i]) for i in keep_idx] + [np.median(chrom_sizes[[i for i in range(len(chrom_sizes)) if i not in keep_idx]])])
        else:
            sizes = chrom_sizes.astype(float)
        sums = S.sum(axis=1)
        # sort by size desc, tie-breaker sums desc
        order_idx = np.lexsort(( -sums, -sizes ))
    elif method == "greedy":
        # seed by highest total inter-contact (or largest size)
        totals = S.sum(axis=1)
        seed = int(np.argmax(totals))
        visited = [seed]
        unvisited = set(range(new_n_chr))
        unvisited.remove(seed)
        while unvisited:
            last = visited[-1]
            # score = S[last] - size_penalty[last]  (apply penalty)
            scores = S[last].copy() - size_penalty[last]
            # ignore visited
            for v in visited:
                scores[v] = -np.inf
            nxt = int(np.argmax(scores[list(unvisited)]))
            # map back index in list(unvisited)
            nxt = list(unvisited)[nxt]
            visited.append(nxt)
            unvisited.remove(nxt)
        order_idx = np.array(visited, dtype=int)
    elif method == "mst":
        # maximum spanning tree -> DFS order
        # use minimum_spanning_tree on negative weights to get maximum spanning tree
        W = -S
        W[np.isinf(W)] = 0.0
        csr = csr_matrix(W)
        T = minimum_spanning_tree(csr)  # returns sparse tree for min; since W=-S, this is max tree
        # build adjacency
        T = T.toarray()
        T = T + T.T  # undirected
        # adjacency weights positive
        adj = (T != 0).astype(int)
        # DFS
        root = int(np.argmax(S.sum(axis=1)))
        stack = [root]
        visited = set()
        order = []
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            order.append(node)
            nbrs = np.where(adj[node] != 0)[0]
            # sort neighbors by descending weight to make DFS greedy
            nbrs_sorted = sorted(nbrs.tolist(), key=lambda x: -S[node, x])
            for nb in nbrs_sorted[::-1]:
                if nb not in visited:
                    stack.append(nb)
        # append any disconnected nodes
        for i in range(new_n_chr):
            if i not in visited:
                order.append(i)
        order_idx = np.array(order, dtype=int)
    else:
        raise ValueError(f"Unknown method `{method}`")

    # map back collapsed labels if needed
    if collapsed_map is not None and keep_idx is not None:
        labels = [chromnames[i] for i in keep_idx] + ["small_contigs"]
        order_labels = [labels[i] for i in order_idx]
    else:
        order_labels = [chromnames[i] for i in order_idx]
    return order_idx, order_labels


# def compute_chrom_interaction_matrix_from_pixels(cool, chrom_offset, collapse_map=None, chunksize=5_000_000):
#     import numpy as np
#     pixels = cool.pixels()
#     nnz = len(pixels)
#     if nnz == 0:
#         n_new = (np.max(collapse_map) + 1) if collapse_map is not None else (len(chrom_offset) - 1)
#         return np.zeros((n_new, n_new), dtype=float)

#     if collapse_map is None:
#         n_new = len(chrom_offset) - 1
#     else:
#         n_new = int(collapse_map.max() + 1)
#     S = np.zeros((n_new, n_new), dtype=float)
#     for i in range(0, nnz, chunksize):
#         chunk = pixels[i: i + chunksize]
#         b1 = np.asarray(chunk['bin1_id'], dtype=np.int64)
#         b2 = np.asarray(chunk['bin2_id'], dtype=np.int64)
#         w = np.asarray(chunk['count'], dtype=float)
#         chr1 = np.searchsorted(chrom_offset, b1, side='right') - 1
#         chr2 = np.searchsorted(chrom_offset, b2, side='right') - 1
#         mask = chr1 != chr2
#         if not np.any(mask):
#             continue
#         i1 = chr1[mask].astype(np.int64)
#         i2 = chr2[mask].astype(np.int64)
#         weights = w[mask]
#         if collapse_map is not None:
#             i1 = collapse_map[i1]
#             i2 = collapse_map[i2]
#             valid = (i1 >= 0) & (i2 >= 0)
#             if not np.any(valid):
#                 continue
#             i1 = i1[valid]; i2 = i2[valid]; weights = weights[valid]
#         flat = i1 * n_new + i2
#         binc = np.bincount(flat, weights=weights, minlength=n_new * n_new)
#         S += binc.reshape((n_new, n_new))

#     S = S + S.T
#     return S

# def sort_chromosomes_grouped(cool, bins,
#                              hap_pattern=r'(Chr\d+)g(\d+)',
#                              intra_method='sum',
#                              inter_method='sum',
#                              contig_method='size',
#                              compute_chunksize=5_000_000):

#     import re
#     import numpy as np
#     from collections import OrderedDict
#     from scipy.sparse import csr_matrix, issparse
#     from scipy.sparse.csgraph import minimum_spanning_tree

#     if isinstance(cool, str):
#         c = cooler.Cooler(cool)
#     else:
#         c = cool


#     if isinstance(bins, np.ndarray) or isinstance(bins, list):
#         raise ValueError("bins must be a DataFrame-like with 'chrom' and 'index'")
#     if 'chrom' in bins.columns:
#         bdf = bins[['chrom', 'index']].copy()
#     else:
#         bdf = bins.reset_index()[['chrom', 'index']].copy()
#     grp = bdf.groupby('chrom', sort=False)['index'].count()
#     chromnames = grp.index.tolist()
#     counts = grp.values.astype(int)
#     chrom_offset = np.r_[0, np.cumsum(counts)].astype(np.int64)

#     n_chr = len(chromnames)
#     if n_chr <= 1:
#         return np.arange(n_chr, dtype=int), chromnames

#     S = compute_chrom_interaction_matrix_from_pixels(c, chrom_offset, collapse_map=None, chunksize=compute_chunksize)

#     S = np.nan_to_num(np.asarray(S, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
#     if S.shape[0] != n_chr:

#         S = np.zeros((n_chr, n_chr), dtype=float)

#     def sort_members(members, method='sum'):
#         k = len(members)
#         if k <= 1:
#             return members[:]
#         sub = S[np.ix_(members, members)]
#         if method == 'sum':
#             scores = np.asarray(sub.sum(axis=1)).ravel()
#             order = np.argsort(-scores)
#             return [members[i] for i in order]
#         elif method == 'greedy':
#             totals = np.asarray(sub.sum(axis=1)).ravel()
#             seed = int(np.argmax(totals))
#             visited = [members[seed]]
#             unvisited = set(members) - {visited[0]}
#             # map value->pos helper
#             idx_map = {m: i for i, m in enumerate(members)}
#             while unvisited:
#                 last = visited[-1]
#                 li = idx_map[last]
#                 # scores to unvisited
#                 cand_list = list(unvisited)
#                 cand_idx = [idx_map[c] for c in cand_list]
#                 cand_scores = sub[li, cand_idx]
#                 best_rel = int(np.argmax(cand_scores))
#                 best = cand_list[best_rel]
#                 visited.append(best)
#                 unvisited.remove(best)
#             return visited
#         elif method == 'mst':
#             # build complete weight matrix W for these members
#             W = np.asarray(sub, dtype=float)
#             # negative for minimum spanning tree
#             M = -W
#             csr = csr_matrix(M)
#             T = minimum_spanning_tree(csr).toarray()
#             T = T + T.T
#             adj = (T != 0).astype(int)
#             # root = node with largest total weight
#             root = int(np.argmax(W.sum(axis=1)))
#             stack = [root]
#             visited = set()
#             order = []
#             while stack:
#                 node = stack.pop()
#                 if node in visited:
#                     continue
#                 visited.add(node)
#                 order.append(members[node])
#                 nbrs = np.where(adj[node] != 0)[0]
#                 nbrs_sorted = sorted(nbrs.tolist(), key=lambda x: -W[node, x])
#                 for nb in nbrs_sorted[::-1]:
#                     if nb not in visited:
#                         stack.append(nb)
#             # append any missing
#             for i in range(k):
#                 if members[i] not in order:
#                     order.append(members[i])
#             return order
#         else:
#             raise ValueError(f"unknown intra_method `{method}`")


#     pat = re.compile(hap_pattern)
#     group_members = OrderedDict()
#     group_order = []
#     for i, name in enumerate(chromnames):
#         m = pat.match(str(name))
#         key = m.group(1) if m else str(name)
#         if key not in group_members:
#             group_members[key] = []
#             group_order.append(key)
#         group_members[key].append(i)

   
#     ordered_group_members = {}
#     for key, members in group_members.items():
#         if len(members) <= 1:
#             ordered_group_members[key] = members[:]
#             continue
#         ordered_group_members[key] = sort_members(members, method=intra_method)

#     group_keys = list(group_order)
#     Gk = len(group_keys)
#     G = np.zeros((Gk, Gk), dtype=float)

#     for a, ka in enumerate(group_keys):
#         ma = ordered_group_members[ka]
#         for b, kb in enumerate(group_keys[a+1:], start=a+1):
#             mb = ordered_group_members[kb]
#             # sum over cross-block using S and numpy indexing
#             sa = S[np.ix_(ma, mb)].sum()
#             G[a, b] = sa
#             G[b, a] = sa

#     def sort_groups(group_idx_list, method='sum'):
#         m = len(group_idx_list)
#         if m <= 1:
#             return group_idx_list[:]
#         if method == 'sum':
#             scores = G.sum(axis=1)
#             order = np.argsort(-scores)
#             return [group_idx_list[i] for i in order]
#         elif method == 'greedy':
#             totals = G.sum(axis=1)
#             seed = int(np.argmax(totals))
#             visited = [group_idx_list[seed]]
#             unvisited = set(group_idx_list) - {visited[0]}
#             while unvisited:
#                 last = visited[-1]
#                 last_idx = group_idx_list.index(last)
#                 best = max(unvisited, key=lambda x: G[last_idx, group_idx_list.index(x)])
#                 visited.append(best)
#                 unvisited.remove(best)
#             return visited
#         elif method == 'mst':
#             M = -G
#             csr = csr_matrix(M)
#             T = minimum_spanning_tree(csr).toarray()
#             T = T + T.T
#             adj = (T != 0).astype(int)
#             root = int(np.argmax(G.sum(axis=1)))
#             stack = [root]
#             visited = set()
#             order = []
#             while stack:
#                 node = stack.pop()
#                 if node in visited:
#                     continue
#                 visited.add(node)
#                 order.append(group_idx_list[node])
#                 nbrs = np.where(adj[node] != 0)[0]
#                 nbrs_sorted = sorted(nbrs.tolist(), key=lambda x: -G[node, x])
#                 for nb in nbrs_sorted[::-1]:
#                     if nb not in visited:
#                         stack.append(nb)
#             for i in range(m):
#                 if group_idx_list[i] not in order:
#                     order.append(group_idx_list[i])
#             return order
#         else:
#             raise ValueError(f"unknown inter_method `{method}`")

#     group_idx_list = list(range(Gk))
#     ordered_group_idx = sort_groups(group_idx_list, method=inter_method)

#     final_order = []
#     final_labels = []
#     for gi in ordered_group_idx:
#         key = group_keys[gi]
#         members = ordered_group_members[key]
#         final_order.extend(members)
#         final_labels.extend([chromnames[i] for i in members])

#     return np.asarray(final_order, dtype=int), final_labels

# def compute_chrom_interaction_matrix_from_pixels(cool, chrom_offset, collapse_map=None, chunksize=5_000_000):
#     import numpy as np
#     from scipy.sparse import coo_matrix, csr_matrix

#     pixels = cool.pixels()
#     nnz = len(pixels)
#     if nnz == 0:
#         n_new = (int(collapse_map.max()) + 1) if collapse_map is not None else (len(chrom_offset) - 1)
#         return csr_matrix((n_new, n_new), dtype=float)

#     if collapse_map is None:
#         n_new = len(chrom_offset) - 1
#     else:
#         n_new = int(collapse_map.max() + 1)

#     rows = []
#     cols = []
#     data = []

#     for i in range(0, nnz, chunksize):
#         chunk = pixels[i: i + chunksize]
#         b1 = np.asarray(chunk['bin1_id'], dtype=np.int64)
#         b2 = np.asarray(chunk['bin2_id'], dtype=np.int64)
#         w = np.asarray(chunk['count'], dtype=float)

#         chr1 = np.searchsorted(chrom_offset, b1, side='right') - 1
#         chr2 = np.searchsorted(chrom_offset, b2, side='right') - 1

#         mask = chr1 != chr2
#         if not np.any(mask):
#             continue

#         i1 = chr1[mask].astype(np.int64)
#         i2 = chr2[mask].astype(np.int64)
#         weights = w[mask]

#         if collapse_map is not None:
#             i1 = collapse_map[i1]
#             i2 = collapse_map[i2]
#             valid = (i1 >= 0) & (i2 >= 0)
#             if not np.any(valid):
#                 continue
#             i1 = i1[valid].astype(np.int64)
#             i2 = i2[valid].astype(np.int64)
#             weights = weights[valid]

#         flat = i1 * n_new + i2

#         uniq, inv = np.unique(flat, return_inverse=True)
#         sums = np.bincount(inv, weights=weights)
#         if sums.size == 0:
#             continue
#         rows_chunk = (uniq // n_new).astype(np.int32)
#         cols_chunk = (uniq % n_new).astype(np.int32)

#         rows.extend(rows_chunk.tolist())
#         cols.extend(cols_chunk.tolist())
#         data.extend(sums.tolist())

#     if len(data) == 0:
#         return csr_matrix((n_new, n_new), dtype=float)

#     S = coo_matrix((data, (rows, cols)), shape=(n_new, n_new), dtype=float)

#     S = (S + S.T).tocsr()
#     return S


def compute_chrom_interaction_matrix_from_pixels(
    cool, chrom_offset, collapse_map=None, chunksize=5_000_000):

    # pixels = cool.pixels()
    nnz = cool.info['nnz']
    
    if collapse_map is None:
        n_new = len(chrom_offset) - 1
    else:
        n_new = int(collapse_map.max() + 1)

    if nnz == 0:
        return csr_matrix((n_new, n_new), dtype=float)

    n_bins = cool.info['nbins']
    bin_to_chrom = np.full(n_bins, -1, dtype=np.int32)
    for i, (start, end) in enumerate(zip(chrom_offset[:-1], chrom_offset[1:])):
        if start < n_bins:
            # Ensure we don't go out of bounds if chrom_offset is weird
            end = min(end, n_bins)
            bin_to_chrom[start:end] = i

    use_dense = (n_new <= 5000)
    
    if use_dense:
        S_dense = np.zeros((n_new, n_new), dtype=np.float64)
    else:
        rows_list = []
        cols_list = []
        data_list = []

    pixels = cool.pixels()
    
    for i in range(0, nnz, chunksize):
        chunk = pixels[i: i + chunksize]
        b1 = np.asarray(chunk['bin1_id'], dtype=np.int64)
        b2 = np.asarray(chunk['bin2_id'], dtype=np.int64)
        w = np.asarray(chunk['count'], dtype=float)

        chr1 = bin_to_chrom[b1]
        chr2 = bin_to_chrom[b2]

        mask = (chr1 != -1) & (chr2 != -1) & (chr1 != chr2)
        if not np.any(mask):
            continue

        i1 = chr1[mask]
        i2 = chr2[mask]
        weights = w[mask]

        if collapse_map is not None:
            i1 = collapse_map[i1]
            i2 = collapse_map[i2]
            valid = (i1 >= 0) & (i2 >= 0)
            if not np.any(valid):
                continue
            i1 = i1[valid]
            i2 = i2[valid]
            weights = weights[valid]

        flat = i1.astype(np.int64) * n_new + i2.astype(np.int64)

        if use_dense:
            if flat.size == 0:
                continue
            counts = np.bincount(flat, weights=weights, minlength=n_new*n_new)
            if counts.size > n_new*n_new:
                counts = counts[:n_new*n_new]
            S_dense += counts.reshape((n_new, n_new))
        else:
            uniq, inv = np.unique(flat, return_inverse=True)
            sums = np.bincount(inv, weights=weights)
            
            rows_chunk = (uniq // n_new).astype(np.int32)
            cols_chunk = (uniq % n_new).astype(np.int32)

            rows_list.append(rows_chunk)
            cols_list.append(cols_chunk)
            data_list.append(sums)

    if use_dense:
        S_dense += S_dense.T
        return csr_matrix(S_dense)
    else:
        if not data_list:
            return csr_matrix((n_new, n_new), dtype=float)
        
        rows = np.concatenate(rows_list)
        cols = np.concatenate(cols_list)
        data = np.concatenate(data_list)

        S = coo_matrix((data, (rows, cols)), shape=(n_new, n_new), dtype=float)
        S = (S + S.T).tocsr()
        return S


def sort_chromosomes_grouped(cool, 
                             hap_pattern=r'(Chr\d+)g(\d+)',
                             intra_method='greedy',
                             inter_method='greedy',
                             contig_method='size',
                             compute_chunksize=5_000_000,
                             min_inter_contacts=25):
    import re
    import numpy as np
    from collections import OrderedDict
    from scipy.sparse import csr_matrix, issparse
    from scipy.sparse.csgraph import minimum_spanning_tree

    if isinstance(cool, str):
        c = cooler.Cooler(cool)
    else:
        c = cool

    bins = c.bins()[:]['chrom'].reset_index()

    if isinstance(bins, np.ndarray) or isinstance(bins, list):
        raise ValueError("bins must be a DataFrame-like with 'chrom' and 'index'")
    if 'chrom' in bins.columns:
        bdf = bins[['chrom', 'index']].copy()
    else:
        bdf = bins.reset_index()[['chrom', 'index']].copy()
    grp = bdf.groupby('chrom', sort=False)['index'].count()
    chromnames = grp.index.tolist()
    counts = grp.values.astype(int)
    chrom_offset = np.r_[0, np.cumsum(counts)].astype(np.int64)

    n_chr = len(chromnames)
    if n_chr <= 1:
        return np.arange(n_chr, dtype=int), chromnames


    S = compute_chrom_interaction_matrix_from_pixels(
        c, chrom_offset, collapse_map=None, chunksize=compute_chunksize)
    if not issparse(S):
        S = csr_matrix(S)
    

    inter_sum = np.asarray(S.sum(axis=1)).ravel()
    keep_mask = inter_sum >= min_inter_contacts
    removed_names = []
    removed_indices = []
    if not np.all(keep_mask):
        removed_indices = np.where(~keep_mask)[0]
        removed_names = [chromnames[i] for i in removed_indices]
        valid_indices = np.where(keep_mask)[0]
        S = S[valid_indices, :][:, valid_indices]
        chromnames = [chromnames[i] for i in valid_indices]
        counts = counts[valid_indices]
        chrom_offset = np.r_[0, np.cumsum(counts)].astype(np.int64)
        n_chr = len(chromnames)
        logger.info(f"Filtered out {len(keep_mask) - n_chr} contigs with < {min_inter_contacts} inter contacts.")

    S = S.tocsr()  # ensure csr

 
    def sort_members(members, method='sum'):
        k = len(members)
        if k <= 1:
            return members[:]
        sub = S[members, :][:, members]  # sparse k x k
        sub_arr = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub, dtype=float)
        if method == 'sum':
            scores = sub_arr.sum(axis=1)
            order = np.argsort(-scores)
            return [members[i] for i in order]
        elif method == 'greedy':
            n = sub_arr.shape[0]
            if n == 0:
                return []
            
            # Vectorized greedy sort
            # totals = sub_arr.sum(axis=1)
            # current = np.argmax(totals)
            
            # path = np.zeros(n, dtype=int)
            # path[0] = current
            
            # mask = np.ones(n, dtype=bool)
            # mask[current] = False
            
            # for i in range(1, n):
            #     # Get interactions of current node
            #     scores = sub_arr[current].copy()
            #     # Mask visited nodes with -1.0 (assuming interactions >= 0)
            #     scores[~mask] = -1.0
                
            #     nxt = np.argmax(scores)
            #     path[i] = nxt
            #     mask[nxt] = False
            #     current = nxt
                
            # return [members[i] for i in path]
            path = []
            mask = np.ones(n, dtype=bool)
            
            while len(path) < n:
                remaining_indices = np.where(mask)[0]
                if len(path) == 0:
                    current = remaining_indices[np.argmax(sub_arr[remaining_indices].sum(axis=1))]
                else:
                    scores = sub[path[-1]].toarray().flatten()
                    scores[~mask] = -1.0
                    if np.max(scores) <= 0: 
                        current = remaining_indices[np.argmax(np.asarray(sub[remaining_indices].sum(axis=1)).ravel())]

                    else:
                        current = np.argmax(scores)
                
                path.append(current)
                mask[current] = False
                
            return [members[i] for i in path]
        elif method == 'mst':
            W = sub_arr.copy()
            M = -W
            csr = csr_matrix(M)
            T = minimum_spanning_tree(csr).toarray()
            T = T + T.T
            adj = (T != 0).astype(int)
            root = int(np.argmax(W.sum(axis=1)))
            stack = [root]
            visited_nodes = set()
            order = []
            while stack:
                node = stack.pop()
                if node in visited_nodes:
                    continue
                visited_nodes.add(node)
                order.append(members[node])
                nbrs = np.where(adj[node] != 0)[0]
                nbrs_sorted = sorted(nbrs.tolist(), key=lambda x: -W[node, x])
                for nb in nbrs_sorted[::-1]:
                    if nb not in visited_nodes:
                        stack.append(nb)
            for i in range(k):
                if members[i] not in order:
                    order.append(members[i])
            return order
        else:
            raise ValueError(f"unknown intra_method `{method}`")


    pat = re.compile(hap_pattern)
    group_members = OrderedDict()
    group_order = []
    for i, name in enumerate(chromnames):
        m = pat.match(str(name))
        key = m.group(1) if m else str(name)
        if key not in group_members:
            group_members[key] = []
            group_order.append(key)
        group_members[key].append(i)

    ordered_group_members = {}
    for key, members in group_members.items():
        if len(members) <= 1:
            ordered_group_members[key] = members[:]
        else:
            ordered_group_members[key] = sort_members(members, method=intra_method)

    group_keys = list(group_order)
    Gk = len(group_keys)
    
    group_lengths = []
    for key in group_keys:
        mems = group_members[key]
        group_lengths.append(sum(counts[m] for m in mems))
    group_lengths = np.array(group_lengths, dtype=float)

    row_ind = []
    col_ind = []
    for g_idx, key in enumerate(group_keys):
        mems = ordered_group_members[key]
        row_ind.extend(mems)
        col_ind.extend([g_idx] * len(mems))
    
    data = np.ones(len(row_ind), dtype=float)
    P = csr_matrix((data, (row_ind, col_ind)), shape=(n_chr, Gk))
    
    # G = P.T @ S @ P
    G_mat = P.T @ S @ P
    
    use_dense_G = (Gk <= 5000)
    
    if issparse(G_mat):
        G_coo = G_mat.tocoo()
       
        norm_data = G_coo.data / (group_lengths[G_coo.row] * group_lengths[G_coo.col])
        G_coo.data = np.log1p(norm_data * 1e6)
        G = G_coo.tocsr()
    else:
        G = np.asarray(G_mat)
        L = np.outer(group_lengths, group_lengths)
        L[L == 0] = 1.0
        G = np.log1p((G / L) * 1e6)
    
    if use_dense_G and issparse(G):
        G = G.toarray()
        
    if not issparse(G):
        np.fill_diagonal(G, 0.0)


    L = np.outer(group_lengths, group_lengths)
    L[L == 0] = 1.0 

    G = G / L
    G = np.log1p(G * 1e6)

    def sort_groups(group_idx_list, method='sum'):
        m = len(group_idx_list)
        if m <= 1:
            return group_idx_list[:]
        nonlocal G
        if issparse(G) and not hasattr(G, "__getitem__"):
            G = G.tocsr()
        if method == 'sum':
            scores = G.sum(axis=1)
            order = np.argsort(-scores)
            return [group_idx_list[i] for i in order]
        elif method == 'greedy':
            n = G.shape[0]
            if n == 0:
                return []
            
            path = []
            mask = np.ones(n, dtype=bool)
            if issparse(G):
                threshold = np.percentile(G.data, 50) if G.nnz > 0 else 0
            else:
                threshold = np.percentile(G, 50)

            while len(path) < n:
                remaining_indices = np.where(mask)[0]
                if len(path) == 0:
                    sums = np.asarray(G[remaining_indices].sum(axis=1)).ravel()
                    current = remaining_indices[np.argmax(sums)]
                else:
                  
                    if issparse(G):
                        scores = G[path[-1]].toarray().flatten()
                    else:
                        scores = G[path[-1]].copy()
                        
                    scores[~mask] = -1.0
                    if np.max(scores) <= threshold: 
                        sums = np.asarray(G[remaining_indices].sum(axis=1)).ravel()
                        current = remaining_indices[np.argmax(sums)]
                    else:
                        current = np.argmax(scores)
                
                path.append(current)
                mask[current] = False
                
            return [group_idx_list[i] for i in path]
        elif method == 'mst':
            M = -G
            csr = csr_matrix(M)
            T = minimum_spanning_tree(csr).toarray()
            T = T + T.T
            adj = (T != 0).astype(int)
            root = int(np.argmax(G.sum(axis=1)))
            stack = [root]
            visited = set()
            order = []
            while stack:
                node = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                order.append(group_idx_list[node])
                nbrs = np.where(adj[node] != 0)[0]
                nbrs_sorted = sorted(nbrs.tolist(), key=lambda x: -G[node, x])
                for nb in nbrs_sorted[::-1]:
                    if nb not in visited:
                        stack.append(nb)
            for i in range(m):
                if group_idx_list[i] not in order:
                    order.append(group_idx_list[i])
            return order
        else:
            raise ValueError(f"unknown inter_method `{method}`")

    group_idx_list = list(range(Gk))
    ordered_group_idx = sort_groups(group_idx_list, method=inter_method)
 
    final_order = []
    final_labels = []
    for gi in ordered_group_idx:
        key = group_keys[gi]
        members = ordered_group_members[key]
        final_order.extend(members)
        final_labels.extend([chromnames[i] for i in members])

    if removed_names:
        start_idx = len(chromnames)
        for i, name in enumerate(removed_names):
            final_labels.append(name)

    final_order = np.concatenate([final_order, removed_indices])

    return np.asarray(final_order, dtype=int), final_labels
