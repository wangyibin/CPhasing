#!/usr/bin/env python

import logging
import os
import os.path as op
import sys
import time
import tempfile
import warnings

try:
    from pandas.core.common import SettingWithCopyWarning
    warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
except ImportError:
    pass


import cooler
import numpy as np
import pandas as pd

from collections import OrderedDict
from math import ceil
from itertools import product
from intervaltree import Interval, IntervalTree
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator
from multiprocessing import Lock, Pool
from joblib import Parallel, delayed
from pandarallel import pandarallel
from scipy.sparse import triu

from .agp import import_agp
from .utilities import to_humanized, chrom_ticks_convert

logger = logging.getLogger(__name__)

HIC_METADATA = {}
HIC_METADATA['matrix-generated-by'] = np.string_(
    'CPhasing'
)
HIC_METADATA['matrix-generated-by-url'] = np.string_(
    'https://github.com/wangyibin/CPhasing'
)

## https://github.com/XiaoTaoWang/NeoLoopFinder/blob/master/neoloop/visualize/core.py
whitered_cmap = LinearSegmentedColormap.from_list('interaction', ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000'])

lock = Lock()
 
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
    chrom_size_dict = dict(chrom_size)
    if region:
        # chrom_size = [(region, chrom_size_dict[region])]
        tmp = []
        for chroms in chrom_size:
            chrom, size = chroms
            if chrom == region:
                tmp.append(chroms)
        chrom_size = tmp 
    
    for chrom, size in chrom_size:
        length = size - start
        if orientation == "-":
      
            if length % bin_size > 0:
                old_start = start 
                start = start + (length % bin_size)
                bin_intervals.append((chrom, old_start, start))
        for interval in range(start, size, bin_size):
                bin_intervals.append((chrom, interval, 
                                    min(size, interval + bin_size)))

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
    db = []
    agp_df['start'] -= 1
    agp_df['tig_start'] = agp_df['tig_start'].astype(int) - 1
    def process_row(row):
        if int(row.end - row.start + 1) <= bin_size:
            return pd.DataFrame([row])
        else:
            tmp_chrom_bins = get_bins(bin_size,  [(row.name, int(row.end))],
                                    start=row.start, orientation=row.orientation)
            tmp_contig_bins = get_bins(
                bin_size, [(row.id, int(row.tig_end))], start=0,
                reverse=False if row.orientation == "+" else True)

            rows = []
            for (_, start, end), (_, tig_start, tig_end) in \
                    zip(tmp_chrom_bins, tmp_contig_bins):
                tmp_row = row.copy()
                tmp_row[['start', 'end', 'tig_start', 'tig_end']] = start, \
                        end, tig_start, tig_end
               
                rows.append(tmp_row)

        return pd.DataFrame(rows)

    res_df = pd.concat(agp_df.parallel_apply(process_row, axis=1).tolist())

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
    logger.debug("Get the coordinates of contig on chrom.")
    _file1 = tempfile.NamedTemporaryFile(delete=True, dir=dir)
    _file2 = tempfile.NamedTemporaryFile(delete=True, dir=dir)
    _file3 = tempfile.NamedTemporaryFile(delete=True, dir=dir)

    chrom_bins.to_csv(_file1.name, sep='\t', header=None, index=None)

    contig_on_chrom_bins.to_csv(_file2.name,
                                sep='\t', index=True, header=None)

    os.system('bedtools intersect -a {} -b {} -F 0.5 -wo > {} 2>/dev/null'.format(
        _file1.name, _file2.name, _file3.name))

    df = pd.read_csv(_file3.name, sep='\t',
                     header=None, index_col=None,
                     dtype={9: 'category'})
    
    _file1.close()
    _file2.close()
    _file3.close()

    df = df[~df.where(df == '.').any(axis=1)]
    df.drop([3, 4, 5, 10], axis=1, inplace=True)
    df.columns = ['chrom', 'start', 'end', 'contig',
                  'tig_start', 'tig_end', 'orientation']
    
    return df

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
    

class sumSmallContig(object):
    """
    Sum small conitg count into one chromosome bin
    """
    def __init__(self, cool_path, contig2chrom, edges,
                    columns, map, batchsize=1000):

        self._map = map
        self.cool_path = cool_path
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
        

    def _aggregate(self, span):

        cool = cooler.Cooler(self.cool_path)
        pixels = cool.matrix(balance=False, sparse=True, as_pixels=True)
        contig2chrom_index = self.contig2chrom_index

        lo, hi = span
        chunk = pixels[lo: hi+1]
        
        old_bin1_id = chunk['bin1_id'].values
        old_bin2_id = chunk['bin2_id'].values
     
        chunk['bin1_id'] = np.searchsorted(contig2chrom_index, old_bin1_id, 
                                            side='right') - 1
        chunk['bin2_id'] = np.searchsorted(contig2chrom_index, old_bin2_id, 
                                            side='right') - 1

        return (chunk.groupby(self.index_columns, sort=True)
                    .aggregate(self.agg)
                    .reset_index())
      

    def aggregate(self, span):
        try:
            chunk = self._aggregate(span)
    
        except MemoryError as e:
            raise RuntimeError(str(e))
        return chunk
    
    def __iter__(self):
        
        batchsize = self.batchsize
        spans = self.edges
        spans_length = len(spans)
        for i in range(0, spans_length, batchsize):
            try:
                if batchsize > 1:
                    lock.acquire()
                results = self._map(self.aggregate, spans[i: i+batchsize])
                
            finally:
                if batchsize > 1:
                    lock.release()
        
            for df in results:
                yield {k: v.values for k, v in df.items()}
                

def sum_small_contig(cool_path, contig2chrom, new_bins, output, 
                dtypes=None, columns=['count'], threads=1, **kwargs):
    from cooler.create import create

    cool = cooler.Cooler(cool_path)

    edges = np.r_[0, np.cumsum(new_bins.groupby(
        'chrom').count().reset_index()['start'].tolist())]
    edges = list(zip(edges[:-1], edges[1:]))
    
    if dtypes is None:
        dtypes = {}
    input_dtypes = cool.pixels().dtypes
 
    for col in columns:
        if col in input_dtypes:
            dtypes.setdefault(col, input_dtypes[col])
    
    try:
        if threads > 1:
            pool = Pool(threads)
            kwargs.setdefault('lock', lock)
        
        iterator = sumSmallContig(cool_path, 
            contig2chrom,
            edges,
            columns,
            map=pool.map if threads > 1 else map, 
        )

        #kwargs.setdefault("append", True)
    
        create(output, new_bins, iterator, dtypes=dtypes, **kwargs)
    
    finally:
        if threads > 1:
            pool.close()


def adjust_matrix(matrix, agp, outprefix=None, chromSize=None, threads=4):

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

    split_contig_on_chrom_df.drop(['chrom', 'start', 'end',
                                   'contig', 'tig_start', 'tig_end',
                                   'chrom_region', 'contig_region'],
                                  inplace=True, axis=1)

    
    logger.info('Starting to reorder matrix ...')

    matrix = cool.matrix(balance=False, sparse=True)
    dtypes = dict(cool.pixels().dtypes)
    
    contig2chrom = split_contig_on_chrom_df[['chromidx', 'contigidx']]
    contig2chrom.set_index('chromidx', inplace=True)
    grouped_contig_idx = split_contig_on_chrom_df.groupby('chromidx')[
                                        'contigidx'].aggregate(tuple)
    grouped_contig_idx = grouped_contig_idx.tolist()
    
    # reorder matrix 
    reordered_contigidx = contig2chrom['contigidx'].values

    reordered_matrix = matrix[:].tocsr(
                         )[:, reordered_contigidx][reordered_contigidx, :]
    reordered_contig_bins = contig_bins[:].loc[reordered_contigidx].reset_index(
        drop=True)
    
    reordered_matrix = triu(reordered_matrix).tocoo()

    chrom_pixels = {
        'bin1_id': reordered_matrix.row,
        'bin2_id': reordered_matrix.col,
        'count': reordered_matrix.data
    }

    order_cool_path = f"{outprefix}.ordered.cool"
    cooler.create_cooler(order_cool_path, reordered_contig_bins,
                         chrom_pixels, dtypes=dtypes)#, metadata=HIC_METADATA)
    logger.info('Successful, reorder the contig-level matrix, '
                f' and output into `{outprefix}.ordered.cool`')
    
    logger.info('Starting to convert contig bin to chromosome bin ...')
    contig2chrom['contigidx'] = range(len(contig2chrom))
    contig2chrom = contig2chrom.reset_index().set_index('chromidx')
    
    sum_small_contig(order_cool_path, contig2chrom, chrom_bin_interval_df, 
                     f'{outprefix}.chrom.cool', dtypes=dtypes)#, metadata=HIC_METADATA)
    logger.info('Successful, converted the contact into chromosome-level'
                f' and output into `{outprefix}.chrom.cool`')
    

    logger.info('Successful, adjusted matrix, elasped time {:.2f}s'.format(time.time() - start_time))

    logger.info(f'Removed `{order_cool_path}`')
    os.remove(order_cool_path)
    
    return f'{outprefix}.chrom.cool'

def coarsen_matrix(cool, k, out, threads):
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
        input_resolution = cooler.Cooler(cool).binsize
        out_resolution = f"{k * input_resolution}"
        
        out = cool.replace(str(input_resolution), to_humanized(out_resolution))
        out = out.replace(to_humanized(input_resolution), to_humanized(out_resolution))

        if to_humanized(out_resolution) not in out:
            if "chrom.cool" in out:
                out = out.replace("chrom.cool", f"{to_humanized(out_resolution)}.chrom.cool")
            else:
                out = out.replace("cool", f"{to_humanized(out_resolution)}.chrom.cool")
            

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


def balance_matrix(cool, force=False, threads=4):
    from cooler.cli.balance import balance 
    
    try:
        logger.info(f'Balancing matrix ...')
        if force:
            balance.main(args=[cool, 
                            '-p', threads,
                             '--force' ], 
                            prog_name='cooler')
        else:
            balance.main(args=[cool, 
                            '-p', threads, ], 
                            prog_name='cooler')
        
    except SystemExit as e:
        exc_info = sys.exc_info()
        exit_code = e.code
        if exit_code is None:
            exit_code = 0
        
        if exit_code != 0:
            raise e


# def plot_matrix(matrix, output, chroms, 
#                 per_chromosomes, cmap, dpi):
#     from hicexplorer import hicPlotMatrix

#     addition_options = []
#     if per_chromosomes:
#         addition_options.append('--perChromosome')

#     if chroms and chroms[0] != "":
#         addition_options.append('--chromosomeOrder')
#         addition_options.extend(list(chroms))
#         # cool = cooler.Cooler(matrix)
        
  
    
#     hicPlotMatrix.main(args=['-m', matrix,
#                             '--dpi', str(dpi),
#                             '--outFileName', output,
#                             '--colorMap', cmap,
#                             '--log1p'] 
#                             + 
#                             addition_options)
    
#     logger.info(f'Successful, plotted the heatmap into `{output}`')

def plot_heatmap(matrix, output, 
                    vmin=None, vmax=None,
                    scale="log1p", 
                    triangle=False,
                    xlabel=None, ylabel=None, 
                    xticks=True, yticks=True,
                    rotate_xticks=False, rotate_yticks=False,
                    remove_short_bin=True,
                    add_lines=False,
                    chromosomes=None, 
                    per_chromosomes=False,
                    chrom_per_row=4,
                    fontsize=None,
                    dpi=1200, 
                    cmap="redp1_r",
                    balanced=False,
                    threads=1):
    import matplotlib.pyplot as plt

    from matplotlib.colors import LogNorm

    logger.info(f"Load contact matrix `{matrix}`.")
    cool = cooler.Cooler(matrix)
    if balanced:
        log1p = False 
        log = False 
        scale = ""
        
    if scale == "log1p":
        log1p = True
    elif scale == "log":
        log = True 
    else:
        log = False
        log1p =False



    if not per_chromosomes:
        chrom_offset = cool._load_dset('indexes/chrom_offset')
        chromsizes = cool.chromsizes
        binsize = cool.binsize
        bins = cool.bins()[:].reset_index(drop=False)
        bins['chrom'] = bins['chrom'].astype('str')
        bins = bins.set_index('chrom')

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
            else:
                bins = bins.loc[chromosomes]
                new_idx = bins['index']
                matrix = matrix[new_idx, :][:, new_idx]
                chromnames = chromosomes
            
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
                logger.debug("Removing the short bin ...")
                retain_chroms = chromsizes[chromsizes >= binsize].index.values
                if len(retain_chroms) < len(chromsizes): 
                    bins = bins.loc[retain_chroms]
                    new_idx = bins['index']
                    matrix = matrix[new_idx, :][:, new_idx]
                    chromnames = retain_chroms.tolist()
                    grouped_counts = new_idx.reset_index().groupby('chrom', sort=False).count().values.flatten()
                    chrom_offset = np.r_[0, np.cumsum(grouped_counts)].tolist()
                else:
                    chrom_offset = chrom_offset.tolist()

        matrix = matrix.todense()

        if log1p or log:
            logger.debug("Masking the zero ...")
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
            matrix += 1 
            # norm = LogNorm(vmax=vmax, vmin=vmin)
            norm = "log1p"
            matrix = np.log10(matrix)
        elif log:
            norm = "log"
            matrix = np.log(matrix)
            # norm = LogNorm(vmax=vmax, vmin=vmin)
        else:
            if balanced:
                norm = "balanced"
            else:
                norm = None

        factor = round(matrix.shape[0] / 5000 + 1, 2)   
        
        fig_width = round(7 * factor, 2) 
        fig_height = round(8 * factor, 2)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)
        
        tick_fontsize = 15 / np.log10(len(chromnames)) if not fontsize else fontsize
        
        logger.info("Plotting heatmap ...")
        ax = plot_heatmap_core(matrix, ax, bins=bins, chromnames=chromnames, 
                            chrom_offset=chrom_offset, norm=norm,
                            triangle=triangle,
                            xticks=xticks, yticks=yticks, 
                            rotate_xticks=rotate_xticks, rotate_yticks=rotate_yticks,
                            vmin=vmin, vmax=vmax, 
                            cmap=cmap, add_lines=add_lines,
                            tick_fontsize=tick_fontsize)
    
    else: 
        ax = plot_per_chromosome_heatmap(cool, chromosomes,
                                         chrom_per_row=chrom_per_row,
                                        cmap=cmap, balanced=balanced, threads=threads)

    plt.savefig(output, dpi=dpi, bbox_inches='tight')

    logger.info(f'Successful, plotted the heatmap into `{output}`')

    return ax


def plot_per_chromosome_heatmap(cool, chromosomes, log1p=True, 
                                    chrom_per_row=4, remove_short_bin=True,
                                    cmap='redp1_r', balanced=False, threads=1):
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
    binsize = cool.binsize
    bins = cool.bins()[:].reset_index(drop=False)
    bins['chrom'] = bins['chrom'].astype('str')
    bins = bins.set_index('chrom')

    matrix = cool.matrix(balance=balanced, sparse=True)[:].tocsr()

    if chromosomes:
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

    matrix = matrix.todense()
    chromosomes = chromosomes if chromosomes else cool.chromnames
    
    num_rows = int(ceil(float(len(chromosomes)) / chrom_per_row))
    num_cols = min(chrom_per_row, len(chromosomes))
    width_ratios = [1.0] * num_cols + [0.05]
    grids = gridspec.GridSpec(num_rows, num_cols + 1,
                                width_ratios=width_ratios,
                                height_ratios=[1] * num_rows)

    fig_height = 6 * num_rows
    fig_width = sum((np.array(width_ratios) + 0.05) * 6)
    

    fig = plt.figure(figsize=(fig_width, fig_height))

    def _plot(ax, chrom, matrix, chrom_range):
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


        plot_heatmap_core(chrom_matrix, ax, bins=bins, chrom_offset=chrom_offset,
                          norm=norm, xlabel=chrom, 
                            cmap=cmap, xticks=False, yticks=False)
    args = []
    for i, chrom in enumerate(chromosomes):
        row = i // chrom_per_row
        col = i % chrom_per_row
    
        ax = plt.subplot(grids[row, col])
        chrom_range = (chrom_offset[i], chrom_offset[i+1])
        args.append((ax, chrom, matrix, chrom_range))
        logger.debug(f"Plotting the heatmap of `{chrom}` ...")
        _plot(ax, chrom, matrix, chrom_range)
    # Parallel(n_jobs=threads)(
    #     delayed(_plot)(i, j, k,l) for i, j, k, l in args)

    return fig

def plot_heatmap_core(matrix, 
                        ax,
                        bins=None,
                        chromnames=None,
                        chrom_offset=None,
                        norm=None, vmin=None, vmax=None,
                        triangle=True,
                        xlabel=None, ylabel=None, 
                        xticks=True, yticks=True,
                        rotate_xticks=False, rotate_yticks=False,
                        tick_fontsize=16,
                        cmap="redp1_r",
                        add_lines=False):
    import colormaps as cmaps
    import matplotlib.pyplot as plt 
    import matplotlib
    matplotlib.use('Agg')       
    from matplotlib import colors
    import seaborn as sns 

    if cmap == 'whitered':
        colormap = whitered_cmap
    else:
        try:
            """
            https://pratiman-91.github.io/colormaps/
            """
            colormap = getattr(cmaps, cmap)
        except:
            colormap = cmap 
    
    if norm is None:
        if vmax is None:
            vmax = np.percentile(np.array(matrix), 95)

   
    # cax = make_axes_locatable(ax).append_axes("right", size="2%", pad=0.09)
    # sns.heatmap(matrix, ax=ax, cmap=colormap, square=True, 
    #                 norm=norm,
    #                 cbar=True, 
    #                 cbar_kws=dict(shrink=.4, pad=0.03))
    
    # cbar = ax.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=10)
    
    if triangle:
        ## https://github.com/XiaoTaoWang/NeoLoopFinder/blob/master/neoloop/visualize/core.py
        # fig = plt.figure(figsize=(7, 4.2))
        # grid = GridSpec(1, 1, figure=fig, left=0.1, right=0.9,
        #             bottom=0.1, top=0.9, hspace=0.04, height_ratios=[5])
        # ax = fig.add_subplot(grid[0])
        add_lines = False
        fig, ax = plt.subplots(figsize=(7, 3.2))
        n = matrix.shape[0]
        t = np.array([[1, 0.5], [-1, 0.5]])
        A = np.dot(np.array([(i[1], i[0]) for i in product(range(n, -1, -1),range(0, n+1, 1))]), t)

        x = A[:, 1].reshape(n+1, n+1)
        y = A[:, 0].reshape(n+1, n+1)
        y[y<0] = -y[y<0]
    
        cax = ax.pcolormesh(x, y, np.flipud(np.array(matrix)), cmap=colormap, edgecolor='none',
                            snap=True, linewidth=.001, rasterized=False)
        fig = plt.gcf()
        # ax.plot(100, 3, marker='v', markersize=2, linestyle=None, color='blue')
        # ax.axis('off')
        
    else:
        fig = plt.gcf()
        cax = ax.imshow(matrix, cmap=colormap, aspect='equal',
                        interpolation=None)
    
    binsize = (bins['end'] - bins['start']).median()
    tmp_chrom_offset = np.array(chrom_offset.copy())
    tmp_chrom_range = list(zip(bins.iloc[tmp_chrom_offset[:-1]]['start'].values, bins.iloc[(tmp_chrom_offset - 1)[1:]]['end'].values))
    tmp_chrom_start = tmp_chrom_range[0][0]
    tmp_chrom_size = tmp_chrom_range[0][0] + sum(list(map(lambda x:(x[1] - x[0]), tmp_chrom_range)))
        
    cax.set_norm(colors.Normalize(vmin=vmin, vmax=vmax))

    ax.set_xlim(0, matrix.shape[0])
    ax.set_ylim(0, matrix.shape[1])


    cbar = fig.colorbar(cax, ax=ax, shrink=.4, pad=0.03)
    cbar.ax.tick_params(labelsize=10)
    cbar.locator = plt.MaxNLocator(5)
    
    if norm == 'log1p':
        cbar.set_label("Log$_{10}$(Contact + 1)", fontsize=12)
    elif norm == 'log':
        cbar.set_label("Log(Contact)", fontsize=12)
    elif norm == 'balanced':
        cbar.set_label("Normalized Contacts", fontsize=12)
    else:
        cbar.set_label("Contacts", fontsize=12)

    if add_lines and chrom_offset:
        ax.hlines(chrom_offset[1:], *ax.get_xlim(), 
                    linewidth=0.5, color='black', linestyles="--")
        ax.vlines(chrom_offset[1:], *ax.get_ylim(), 
                    linewidth=0.5, color='black', linestyles="--")
    
    if chrom_offset:
        mid_tick_pos = list((np.array(chrom_offset)[:-1] + np.array(chrom_offset)[1:]) / 2)
    else:
        mid_tick_pos = [0]

    ax.tick_params(width=0)
    if xticks and chromnames:
        # ax.set_xticks(mid_tick_pos)
        rotation = "horizontal" if rotate_xticks else "vertical" 
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
        ax.tick_params(axis='x', length=5, width=2)
        xticklabels = chrom_ticks_convert(ax.get_xticks()[:-1] * binsize)
        ax.set_xticklabels(xticklabels, fontsize=tick_fontsize, rotation=rotation)

    else:
        ax.set_xticks([])
    if yticks and chromnames:
        ax.set_yticks(mid_tick_pos)
        rotation = "vertical" if rotate_yticks else "horizontal"
        ax.set_yticklabels(chromnames, fontsize=tick_fontsize, rotation=rotation)
    else:
        ax.set_yticks([])

    if xlabel:
        ax.set_xlabel(xlabel, fontsize=18, labelpad=15)
    
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=18, labelpad=15)

    sns.despine(top=False, right=False)


    if triangle:
        ax.axis([0, matrix.shape[0], 0, matrix.shape[1]])
    else:
        ax.axis([0, matrix.shape[0], matrix.shape[1], 0])

    return ax