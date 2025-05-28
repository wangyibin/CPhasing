#!/usr/bin/env python

import logging
import os
import os.path as op
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
import numpy as np
import pandas as pd
import polars as pl
import pyranges as pr

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
    flag = os.system(cmd + " 2>plot.bedtools.log")
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
    

class sumSmallContig(object):
    """
    Sum small conitg count into one chromosome bin
    """
    def __init__(self, chrom_pixels, contig2chrom, edges,
                    columns, map, batchsize=5000):

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

    # @profile
    def _aggregate(self, span):

        # cool = cooler.Cooler(self.cool_path)
    
        # pixels = self.pixels
        # pixels = cool.matrix(balance=False, sparse=True, as_pixels=True)
        contig2chrom_index = self.contig2chrom_index
        lo, hi = span
        
        chunk = self.pixels.filter(pl.col('bin1_id').is_between(lo, hi))
        
        new_bin1_id = np.searchsorted(contig2chrom_index, chunk['bin1_id'].to_numpy(), side='right') - 1
        new_bin2_id = np.searchsorted(contig2chrom_index, chunk['bin2_id'].to_numpy(), side='right') - 1
        
        chunk = chunk.with_columns([
            pl.Series('bin1_id', new_bin1_id),
            pl.Series('bin2_id', new_bin2_id)
        ])

        chunk = chunk.filter((pl.col('bin1_id') >= 0) & (pl.col('bin2_id') >= 0))
        if self.agg['count'] == 'sum':
            result = (chunk.group_by(self.index_columns, maintain_order=True)
                        .agg(pl.sum('count')).to_pandas().reset_index())
        else:
            result = (chunk.to_pandas().groupby(self.index_columns, sort=True)
                            .agg(self.agg)
                            .reset_index())
        
        result = result[(result['bin1_id'] >= 0) & (result['bin2_id'] >= 0)]

        return result

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
        
            # from mpire import WorkerPool
            # with WorkerPool(threads) as pool:
        
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
    raw_chrom_interval_df = chrom_bin_interval_df.copy()
    ## Solve the problem that the dotted line of the heat map will drift when the bin is relatively large
    chrom_bin_interval_df = chrom_bin_interval_df.loc[split_contig_on_chrom_df['chromidx'].drop_duplicates()]
    chrom_bin_interval_df.reset_index(drop=True, inplace=True)
    if chrom_bin_interval_df.loc[0, 'start'] != 0:
        ## if the first bin is not start with 0, then we need to insert a new row before the first bin
        # logger.warning('The first bin of chromosome is not start with 0, '
        #                   'adjust the first bin to start with 0.')
        new_row = pd.DataFrame({
            'chrom': raw_chrom_interval_df.iloc[0]['chrom'],
            'start': 0,
            'end': raw_chrom_interval_df.iloc[0]['start']
        }, index=[0])
        chrom_bin_interval_df = pd.concat([new_row, chrom_bin_interval_df], ignore_index=True)
        
    
    logger.info('Starting to reorder matrix ...')

    matrix = cool.matrix(balance=False, sparse=True)
  
    dtypes = dict(cool.pixels().dtypes)
   
    contig2chrom = split_contig_on_chrom_df[['chromidx', 'contigidx']]
    contig2chrom.set_index('chromidx', inplace=True)
    # reorder matrix 
    reordered_contigidx = contig2chrom['contigidx'].values
    reordered_matrix = matrix[:].tocsr(
                         )[:, reordered_contigidx][reordered_contigidx, :]
    reordered_matrix = triu(reordered_matrix).tocoo()
    os.environ["POLARS_MAX_THREADS"]  = str(threads)

    chrom_pixels = pl.DataFrame({
        'bin1_id': reordered_matrix.row,
        'bin2_id': reordered_matrix.col,
        'count': reordered_matrix.data,
    })
    
    del reordered_matrix
    gc.collect()
    # order_cool_path = f"{outprefix}.ordered.cool"
    # cooler.create_cooler(order_cool_path, reordered_contig_bins,
    #                      chrom_pixels, dtypes=dtypes, triucheck=False,
    #                      dupcheck=False, boundscheck=False)#, metadata=HIC_METADATA)
    # logger.info('Successful, reorder the contig-level matrix, '
    #             f' and output into `{outprefix}.ordered.cool`')
    
    logger.info('Starting to convert contig bin to chromosome bin ...')
    contig2chrom['contigidx'] = range(len(contig2chrom))
    contig2chrom = contig2chrom.reset_index().set_index('chromidx')

    sum_small_contig(chrom_pixels, contig2chrom, chrom_bin_interval_df, 
                     f'{outprefix}.chrom.cool', dtypes=dtypes)#, metadata=HIC_METADATA)
    logger.info('Successful, converted the contact into chromosome-level'
                f' and output into `{outprefix}.chrom.cool`')
    

    logger.info('Successful, adjusted matrix to chromosome-level, elasped time {:.2f}s'.format(time.time() - start_time))
    
    return f'{outprefix}.chrom.cool'

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
                    add_hap_border=False,
                    add_lines=False,
                    chromosomes=None, 
                    hap_pattern=r'(Chr\d+)g(\d+)',
                    per_chromosomes=False,
                    chrom_per_row=4,
                    fontsize=None,
                    dpi=1200, 
                    cmap="redp1_r",
                    balanced=False,
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
    elif scale == "log":
        log = True 
    else:
        log = False
        log1p =False

    logger.info("Plotting heatmap ...")

    if not per_chromosomes:
        chrom_offset = cool._load_dset('indexes/chrom_offset')
        chromsizes = cool.chromsizes
        binsize = cool.binsize
        bins = cool.bins()[:].reset_index(drop=False)
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
                    chromsizes = chromsizes.loc[chromnames]
                else:
                    chrom_offset = chrom_offset.tolist()

        chromnames_df = pd.DataFrame(chromnames, columns=['chrom'])
        hap_name_df = chromnames_df['chrom'].str.extract(hap_pattern, expand=True).dropna()
        hap_name_df.reset_index(drop=False, inplace=True)
        hap_name_df.columns = ['index', 'chrom', 'hap', ]
        hap_name_df = hap_name_df.groupby('chrom')['index'].agg(lambda x: list(x))
        
        if len(hap_name_df) <= 1:
            hap_name_df = None
            hap_names = None
            median_hap_count = 1
        else:
            median_hap_count = hap_name_df.apply(lambda x: len(x)).values.min()
            hap_names = hap_name_df.index.values

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
        factor = np.log2(round((matrix.shape[0] * binsize) / 1e9  + 1, 2)) + 1
    
        logger.debug(f"Figure factor is `{factor}`")
        if triangle:
            fig_width = round(7 * factor, 2) 
            fig_height = round(3.2 * factor, 2)
        else:
            fig_width = round(7 * factor, 2) 
            fig_height = round(8 * factor, 2)

        logger.info(f"  Set figsize: (W: {fig_width}, H: {fig_height})")
        plt.rcParams['font.family'] = 'Arial'
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)

        if hap_names is None or not add_hap_border:
            font_factor = np.log10(len(chromnames) ) if len(chromnames) >= 10 else 1
            tick_fontsize = int(15 / font_factor) if not fontsize else fontsize
           
        else:
            font_factor = np.log10(len(hap_names) ) if len(hap_names) > 5 else 1
            tick_fontsize = int(15 / font_factor) if not fontsize else fontsize
        if tick_fontsize == np.inf:
            tick_fontsize = 16
        logger.info(f"  Set fontsize to `{tick_fontsize}`.")

        
        ax = plot_heatmap_core(matrix, ax, bins=bins, 
                               chromnames=chromnames, 
                               hap_name_df=hap_name_df,
                               chromsizes=chromsizes,
                            chrom_offset=chrom_offset, norm=norm,
                            triangle=triangle,
                            xticks=xticks, yticks=yticks, 
                            rotate_xticks=rotate_xticks, rotate_yticks=rotate_yticks,
                            vmin=vmin, vmax=vmax, 
                            cmap=cmap, add_lines=add_lines,
                            add_hap_border=add_hap_border,
                            tick_fontsize=tick_fontsize)
    
    else: 
        ax = plot_per_chromosome_heatmap(cool, chromosomes,
                                         chrom_per_row=chrom_per_row, triangle=triangle,
                                        cmap=cmap, balanced=balanced, threads=threads)

    logger.info(f"  Set dpi to `{dpi}`.")
    plt.savefig(output, dpi=dpi, bbox_inches='tight')

    logger.info(f'Successful, plotted the heatmap into `{output}`')

    return ax


def plot_per_chromosome_heatmap(cool, chromosomes, log1p=True, 
                                    chrom_per_row=4, remove_short_bin=True, triangle=False,
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
    
    num_rows = int(ceil(float(len(chromosomes)) / chrom_per_row))
    num_cols = min(chrom_per_row, len(chromosomes))
    width_ratios = [1.0] * num_cols + [0.05]
    grids = gridspec.GridSpec(num_rows, num_cols + 1,
                                width_ratios=width_ratios,
                                height_ratios=[1] * num_rows)

    if triangle:
        grids = gridspec.GridSpec(num_rows, num_cols + 1,
                                width_ratios=width_ratios,
                                height_ratios=[1] * num_rows,
                                hspace=0.4, wspace=0.1)
     
        fig_height = 3.8 * num_rows
        fig_width = 7 * num_cols
        fig = plt.figure(figsize=(fig_width, fig_height))

    else:
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


        plot_heatmap_core(chrom_matrix, ax, chromsizes=chromsizes.loc[chrom], bins=bins, chrom_offset=chrom_offset,
                          norm=norm, xlabel=chrom, triangle=triangle,
                            cmap=cmap, xticks=True, yticks=False)
    args = []
    for i, chrom in enumerate(chromosomes):
        row = i // chrom_per_row
        col = i % chrom_per_row
    
        ax = plt.subplot(grids[row, col])
        chrom_range = (chrom_offset[i], chrom_offset[i+1])
        args.append((ax, chrom, matrix, chrom_range))
        logger.debug(f"Plotting the heatmap of `{chrom}` ...")
        _plot(ax, chrom, matrix, chrom_range)

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
                        xlabel=None, ylabel=None, 
                        xticks=True, yticks=True,
                        rotate_xticks=False, rotate_yticks=False,
                        tick_fontsize=16,
                        cmap="redp1_r",
                        add_hap_border=False,
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
            try:
                colormap = cmap 
                colormap = getattr(plt.cm, cmap)
            except AttributeError:
                logger.warning(f"Colormap `{cmap}` not found, use `redp1_r` instead.")
                colormap = getattr(cmaps, "redp1_r")

            
    
    if norm is None or norm == "balanced":
        if vmax is None:
            cis_matrix = []
            for i in range(len(chrom_offset) - 1):
                cis_matrix.append(
                                np.array(matrix[chrom_offset[i]:chrom_offset[i+1], 
                                        chrom_offset[i]:chrom_offset[i+1]]))

            # cis_matrix_median = np.median(np.concatenate([arr.ravel() for arr in cis_matrix]))
            # cis_matrix_std = np.std(np.concatenate([arr.ravel() for arr in cis_matrix]))
            ## remove diagonal
            cis_matrix = [arr[np.triu_indices(arr.shape[0], k=1)] for arr in cis_matrix]
            cis_matrix = np.concatenate([arr.flatten().ravel() for arr in cis_matrix])
            cis_matrix = cis_matrix[~np.isnan(cis_matrix)]
            vmax = np.percentile(cis_matrix, 98)
            
            if vmax == 0:
                vmax = 0.1
            if isinstance(vmax, int):
                logger.info(f"  Set vmax to `{vmax}`.")
            else:
                logger.info(f"  Set vmax to `{vmax:.6f}`.")

            if norm is None and vmin is None:
                vmin = 0
            del cis_matrix
            

    #     if vmax is None:
    #         vmax = np.percentile(np.array(matrix), 95)

   
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
        # fig, ax = plt.subplots(figsize=(7, 3.2))
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
        # ax.xaxis.set_visible(False)
        # ax.yaxis.set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
    else:
        fig = plt.gcf()
        cax = ax.imshow(matrix, cmap=colormap, aspect='equal',
                            interpolation=None)
        
    binsize = np.argmax(np.bincount(bins['end'] - bins['start']))

    # tmp_chrom_offset = np.array(chrom_offset.copy())
    # tmp_chrom_range = list(zip(bins.iloc[tmp_chrom_offset[:-1]]['start'].values, bins.iloc[(tmp_chrom_offset - 1)[1:]]['end'].values))
    # tmp_chrom_start = tmp_chrom_range[0][0]
    # tmp_chrom_size = tmp_chrom_range[0][0] + sum(list(map(lambda x:(x[1] - x[0]), tmp_chrom_range)))
    
    cax.set_norm(colors.Normalize(vmin=vmin, vmax=vmax))

    ax.set_xlim(0, matrix.shape[0])
    ax.set_ylim(0, matrix.shape[1])


    cbar = fig.colorbar(cax, ax=ax, shrink=.4, pad=0.03)
    cbar.ax.tick_params(labelsize=tick_fontsize*0.8)
    cbar.locator = plt.MaxNLocator(5)
    
    if norm == 'log1p':
        cbar.set_label("Log$_{10}$(Contact + 1)", fontsize=tick_fontsize)
    elif norm == 'log':
        cbar.set_label("Log(Contact)", fontsize=tick_fontsize)
    elif norm == 'balanced':
        cbar.set_label("Normalized Contacts", fontsize=tick_fontsize)
    else:
        cbar.set_label("Contacts", fontsize=tick_fontsize)

    if add_lines and chrom_offset:
        ax.hlines(np.array(chrom_offset[1:-1]) - 0.5, *ax.get_xlim(), 
                    linewidth=0.5, color='black', linestyles="--")
        ax.vlines(np.array(chrom_offset[1:-1]) - 0.5, *ax.get_ylim(), 
                    linewidth=0.5, color='black', linestyles="--")
    
    if chrom_offset:
        mid_tick_pos = list((np.array(chrom_offset)[:-1] + np.array(chrom_offset)[1:]) / 2)
    else:
        mid_tick_pos = [0]

    ax.tick_params(width=0)
    if xticks and chromnames:
        rotation = "horizontal" if rotate_xticks else "vertical" 
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
        ax.tick_params(axis='x', length=5, width=2)
        _xticks = ax.get_xticks()[:-1]
        dist = _xticks[1] - _xticks[0]
 
        if (ax.get_xlim()[1] - _xticks[-1]) <= (dist / 2):
            _xticks = _xticks[:-1]

        _xticks = np.r_[_xticks, ax.get_xlim()[1]]
        ax.set_xticks(_xticks)
        _xticks = _xticks * binsize 
        total_size = chromsizes.sum()
        _xticks[-1] = total_size
        xticklabels = chrom_ticks_convert(_xticks)

        ax.set_xticklabels(xticklabels, fontsize=tick_fontsize, rotation=rotation)

    else:
        ax.set_xticks([])
    if yticks and chromnames:
        if hap_name_df is not None and add_hap_border:
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
            ax.set_yticks(new_mid_tick_pos)
            rotation = "vertical" if rotate_yticks else "horizontal"
            ax.set_yticklabels(hap_names, fontsize=tick_fontsize, rotation=rotation)

            ## add lines by start
            for start, end in zip(new_chrom_offset[:-1], new_chrom_offset[1:]):
                plt.plot([start, start], [start, end], color='black', linestyle='-', linewidth=1.0)
                plt.plot([start, end], [start, start], color='black', linestyle='-', linewidth=1.0)
                plt.plot([end, end], [start, end], color='black', linestyle='-', linewidth=1.0)
                plt.plot([start, end], [end, end], color='black', linestyle='-', linewidth=1.0)
            
            # ax.hlines(np.array(new_chrom_offset[:-1]), *ax.get_xlim(), 
            #         linewidth=1.0, color='black', linestyles="-")
            # ax.vlines(np.array(new_chrom_offset[:-1]), *ax.get_ylim(),
            #         linewidth=1.0, color='black', linestyles="-")
                
        else:
            ax.set_yticks(mid_tick_pos)
            rotation = "vertical" if rotate_yticks else "horizontal"
            ax.set_yticklabels(chromnames, fontsize=tick_fontsize, rotation=rotation)
    else:
        ax.set_yticks([])

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