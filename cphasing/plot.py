#!/usr/bin/env python

import logging
import os
import os.path as op
import sys
import time
import warnings


from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
import gc
import tempfile
from collections import OrderedDict
from itertools import combinations

import cooler
import numpy as np
import pandas as pd
from hicexplorer.reduceMatrix import reduce_matrix
from intervaltree import Interval, IntervalTree
from multiprocessing import Lock, Pool
from pandarallel import pandarallel
from scipy.sparse import coo_matrix, csr_matrix, triu

logger = logging.getLogger(__name__)

HIC_METADATA = {}
HIC_METADATA['matrix-generated-by'] = np.string_(
    'CPhasing_plot'
)
HIC_METADATA['matrix-generated-by-url'] = np.string_(
    'https://github.com/wangyibin/CPhasing'
)


lock = Lock()

# def parse_arguments(args=None):
#     p = argparse.ArgumentParser(prog=__file__,
#                         description=__doc__,
#                         formatter_class=argparse.RawTextHelpFormatter,
#                         conflict_handler='resolve')
#     pReq = p.add_argument_group('Required arguments')
#     pOpt = p.add_argument_group('Optional arguments')
#     pReq.add_argument('-m', '--matrix', 
#             help='cool matrix from ALLHiC_buildMatrix',
#             metavar='cool',
#             required=True)
#     pReq.add_argument('-a', '--agp', 
#             help='agp file from ALLHiC.',
#             metavar='agp',
#             required=True)
#     pOpt.add_argument('--binSize', '-bs',
#             default=500000, type=int,
#             help='Size in bp for the plotting heatmap bins. [default: %(default)s]')
#     pOpt.add_argument('--chromSize', default=None,
#             help='chromosome size and chromosome order of heatmap, [default: from agp]')
#     pOpt.add_argument('-t', '--threads', type=int, default=8,
#             help='number of program threads[default:%(default)s]')
#     pOpt.add_argument('--cmap', default='YlOrRd',
#             help='colormap of heatmap [default: %(default)s]')
#     pOpt.add_argument('-o', '--outprefix', 
#             help='prefix of output heatmap [default: prefix.cool]')
#     pOpt.add_argument('-om', '--outmatrix', default=None,
#             help='output the chromosomal-level matrix [default: %(default)s]')
#     pOpt.add_argument('-h', '--help', action='help',
#             help='show help message and exit.')
    
#     return p  

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
    logger.info('Load agp file: `{}`'.format(agpfile))
    
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
    if region:
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
    db = []
    for i, row in agp_df.iterrows():
        row.start = row.start - 1
        row.tig_start = int(row.tig_start) - 1
        #row.tig_end = int(row.tig_end) + 1
        if int(row.end - row.start + 1) <= bin_size:
            tmp_row = row.copy()
            db.append(tmp_row)
        else:
            tmp_chrom_bins = get_bins(bin_size,  [(row.name, int(row.end))],
                                      start=row.start, orientation=row.orientation)
            tmp_contig_bins = get_bins(
                bin_size, [(row.id, int(row.tig_end))], start=0,
                reverse=False if row.orientation == "+" else True)

            for (_, start, end), (_, tig_start, tig_end) in \
                    zip(tmp_chrom_bins, tmp_contig_bins):

                tmp_row = row.copy()
                tmp_row['start', 'end', 'tig_start', 'tig_end'] = start, \
                    end, tig_start, tig_end
                db.append(tmp_row)

    res_df = pd.concat(db, axis=1).T

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

def getContigOnChromBins(chrom_bins, contig_on_chrom_bins):

    _file1 = tempfile.NamedTemporaryFile(delete=False)
    _file2 = tempfile.NamedTemporaryFile(delete=False)
    _file3 = tempfile.NamedTemporaryFile(delete=False)

    chrom_bins.to_csv(_file1.name, sep='\t', header=None, index=None)
    contig_on_chrom_bins.to_csv(_file2.name,
                                sep='\t', index=True, header=None)

    os.system('bedtools intersect -a {} -b {} -F 0.5 -wo > {} 2>/dev/null'.format(
        _file1.name, _file2.name, _file3.name))

    df = pd.read_csv(_file3.name, sep='\t',
                     header=None, index_col=None,
                     dtype={0: 'category',
                               9: 'category'})
    _file1.close()
    _file2.close()
    _file3.close()

    df = df.drop([3, 4, 5, 10], axis=1)
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

       
        contig2chrom_index = self.contig2chrom.groupby('chromidx')['contigidx'].aggregate(list)

        contig2chrom_index = contig2chrom_index.apply(lambda x: x[0])

        self.contig2chrom_index = contig2chrom_index.values
        # self.contig_edges = list(zip(contig2chrom_index[:-1], 
        #                             contig2chrom_index[1:]))
        
        self.batchsize = batchsize 

        self.newbins_index = self.contig2chrom.index.values
        self.index_columns = ['bin1_id', 'bin2_id']
        self.value_coumns = list(columns)
        self.agg = {'count': 'sum'}
        
        edges = []
        for i in range(0, len(self.contig2chrom_index), batchsize):
            tmp_list = self.contig2chrom_index[i: i+batchsize]
            edges.append(tmp_list[-1])
        edges = np.r_[0, edges]
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
    
        for i in range(0, len(spans), batchsize):
            try:
                if batchsize > 1:
                    lock.acquire()
                results = self._map(self.aggregate, spans[i: i+batchsize])
                
            finally:
                if batchsize > 1:
                    lock.release()
        
            for df in results:
                yield {k: v.values for k, v in df.iteritems()}
                

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
        outprefix = op.basename(cool.filename).rsplit(".", 1)[0]
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
    chrom_bin_interval_df = pd.DataFrame(get_bins(bin_size, chrom_sizes), 
                                    columns=['chrom', 'start', 'end'])
    chrom_regions = chrom_bin_interval_df.apply(lambda x: chrRangeID(x.values), axis=1)
    contig_bins = cool.bins()

    new_agp = splitAGPByBinSize(agp_df, bin_size=bin_size)
    
    split_contig_on_chrom_df = getContigOnChromBins(
        chrom_bin_interval_df, new_agp.drop(['number', 'type'], axis=1))
    split_contig_on_chrom_df.drop(['orientation'], inplace=True, axis=1)

    pandarallel.initialize(nb_workers=threads, verbose=0)
    contig_bins_index = contig_bins[:].parallel_apply(
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
    split_contig_on_chrom_df['contig_region'] = split_contig_on_chrom_df[[
        'contig', 'tig_start', 'tig_end']].parallel_apply(func, axis=1)

    cat_dtype = pd.CategoricalDtype(categories=chrom_regions,    
                                        ordered=True)
    split_contig_on_chrom_df['chrom_region'] = \
        split_contig_on_chrom_df['chrom_region'].astype(cat_dtype)
    split_contig_on_chrom_df['chromidx'] = \
        split_contig_on_chrom_df['chrom_region'].cat.codes.values
   
    
    split_contig_on_chrom_df['contigidx'] = contig_bins_index.loc[
        split_contig_on_chrom_df['contig_region']]['contigidx'].values

    split_contig_on_chrom_df.drop(['chrom', 'start', 'end',
                                   'contig', 'tig_start', 'tig_end',
                                   'chrom_region', 'contig_region'],
                                  inplace=True, axis=1)

    
    logger.info('Starting to reorder matrix ...')

    matrix = cool.matrix(balance=False, sparse=True)
    
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

    chrom_pixels = dict(zip(['bin1_id', 'bin2_id', 'count'],
                            [reordered_matrix.row,
                            reordered_matrix.col,
                            reordered_matrix.data]))
    order_cool_path = f"{outprefix}.ordered.cool"
    cooler.create_cooler(order_cool_path, reordered_contig_bins,
                         chrom_pixels, metadata=HIC_METADATA)
    logger.info('Successful, reorder the contig-level matrix, '
                f'and output into `{outprefix}.ordered.cool`')
    
    logger.info('Staring to collaspe chromosome bin ...')
    contig2chrom['contigidx'] = range(len(contig2chrom))
    contig2chrom = contig2chrom.reset_index().set_index('chromidx')
    
    sum_small_contig(order_cool_path, contig2chrom, chrom_bin_interval_df, 
                     f'{outprefix}.chrom.cool', metadata=HIC_METADATA)
    logger.info('Successful, collasped the contact into chromosome-level'
                f'and output into {outprefix}.chrom.cool')
    

    logger.info('Done, elasped time {:.2f}s'.format(time.time() - start_time))

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
        out = cool.replace(str(input_resolution), str(out_resolution))

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

def plot_matrix(matrix, output, chroms, 
                per_chromosomes, cmap, dpi):
    from hicexplorer import hicPlotMatrix

    addition_options = []
    
    if chroms and chroms[0] != "":
        addition_options.append('--chromosomeOrder')
        addition_options.extend(list(chroms))

    if per_chromosomes:
        addition_options.append('--perChromosome')
    
    hicPlotMatrix.main(args=['-m', matrix,
                            '--dpi', str(dpi),
                            '-o', output,
                            '--colorMap', cmap,
                            '--log1p'] 
                            + 
                            addition_options)
    
    logger.info(f'Successful, plotted the heatmap into `{output}`')