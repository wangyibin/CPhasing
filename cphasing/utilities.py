#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
utility libraries
"""

import logging
import os
import os.path as op
import sys
import re

import cooler
import numpy as np
import pandas as pd

from Bio import SeqIO 
from collections import OrderedDict
from pathlib import Path, PosixPath


logger = logging.getLogger(__name__)

def listify(item):
    """
    To return a list or tuple value.
    From https://github.com/tanghaibao/jcvi/blob/master/jcvi/apps/base.py listify
    """
    return item if (isinstance(item, list) or 
                   isinstance(item, tuple)) else [item]


def tail(infile, n, offset=0):
    """
    output n lines in file tail

    Params:
    --------
    infile: str
        input file
    n: int
        number of lines
    offset: int
        offset number [0]

    Returns:
    --------
    lines: list
        list of per line

    Examples:
    --------
    >>> tail("sample.txt", 1)
    ["hello"]
    """

    from subprocess import Popen, PIPE
    p = Popen(['tail', '-n', f'{n + offset}', infile], stdout=PIPE)
    lines = p.stdout.readlines()
    lines = list(map(lambda x: str(x, "utf-8"), lines))

    return lines

def cmd_exists(program):
    """
    Check program is exists.
    """
    from shutil import which
    return which(program) is not None

def run_cmd(command, log=sys.stderr, out2err=False):
    """
    run command on shell

    Params:
    --------
    command: list or tuple
        command list
    log: str, default sys.stderr
        the file of log
    out2err: bool, default False
    Returns:
    --------
    None

    Examples:
    --------
    >>> command = ["ls", "-l"]
    >>> run_cmd(command)
    example.txt result.txt
    0
    """
    from subprocess import Popen, PIPE
    import _io

    logger.info('Running command:')
    logger.info('\t' + ' '.join(command))

    if not isinstance(log, _io.TextIOWrapper):
        errout = open(log, 'w')
    else:
        errout = log
    
    pipelines = []
    try:
        if out2err:
            pipelines.append(
                Popen(command, 
                        stderr=errout,
                        stdout=errout,
                        bufsize=-1)
            )
        else:
            pipelines.append(
                Popen(command, 
                        stderr=errout,
                        bufsize=-1)
            )
        pipelines[-1].wait()
    except:
        return 2
    finally:
        for p in pipelines:
            if p.poll() is None:
                p.terminate()
    
    return 0

def is_empty(_file):
    if os.path.getsize(_file) == 0:
        return True 
    else:
        return False

def is_compressed_table_empty(_file):
    try:
        pd.read_csv(_file, chunksize=1, comment="#")
    except EOFError:
        return True 
    except pd.errors.EmptyDataError:
        return True 

    return False


def xopen(infile, mode='r'):
    """
    open file 

    Params:
    --------
    infile: `str`, input file
    mode: `str`, mode of open ["r"]

    Returns:
    --------
    handle: `_io.TextIOWrapper`

    Examples:
    --------
    >>> xopen('input.fastq.gz', 'r')
    
    """
    import gzip

    if not isinstance(infile, PosixPath):
        infile = Path(infile)

    if infile.suffix == ".gz":
        handle = gzip.open(infile, mode + 't')
    else:
        handle = open(infile, mode)
    return handle

def gz_reader():
    from shutil import which 
    if which('pigz'):
        return 'pigz'
    else:
        return 'gzip'

def list_flatten(list_2d):
    """
    convert 2d list into 1d list

    Params:
    --------
    list_2d: `list` or `array-like`
            2d list [[1, 2, 3], [4, 5, 6]]
    
    Returns:
    --------
    1d list 

    Examples:
    --------
    >>> l = [[1, 2, 3], [4, 5, 6]]
    >>> list_flatten(l)
    [1, 2, 3, 4, 5, 6]
    """

    return [i for item in list_2d for i in item]

def rm_orientation(contig):
    """
    remove orientation in the suffix with contig

    Params:
    --------
    contig: str
        contig with orientation
    
    Returns:
    --------
    str:
        contig without orientation

    Examples:
    --------
    >>> contig = 'utg001+'
    >>> rm_orientation(contig)
    'utg001'

    """
    import re
    return re.sub('[+-]$', '', contig)

def digest(fasta_records, enzyme):
    """
    Divide a genome into restriction fragments. 

    Params:
    --------
    fasta_records: dict
        Dictionary of chromosome names to sequence records.
    enzyme: str
        Name of restriction enzyme. i.e. HindIII, MboI, Arima
    
    Returns:
    --------
    Dataframe with columns: "chrom", "start", "end".

    Examples:
    --------
    >>> fasta_records = {'Chr1': 'AAGCAAAGCGGGATCGATC', 
                        'Chr2': 'AAGCTTGATCGATC'}
    >>> digest(fasta_records, "MboI")
      chrom  start  end
    0  Chr1      0   13
    1  Chr1     13   17
    2  Chr1     17   19
    3  Chr2      0    8
    4  Chr2      8   12
    5  Chr2     12   14
    """
    from Bio import Restriction, Seq


    chroms = fasta_records.keys()
    try:
        if enzyme.lower() == 'arima':
            cut_finder = (
                Restriction.RestrictionBatch(['MboI', 'HinfI'])
                .search
                )
        else:
            cut_finder = getattr(Restriction, enzyme).search
    
    except AttributeError:
        raise ValueError(f'Unknown enzyme name: {enzyme}')
    
    def _each(chrom):
        seq = Seq.Seq(str(fasta_records[chrom]))
        cut_res = cut_finder(seq)
        if isinstance(cut_res, list):
            cut_sites = cut_res
        elif isinstance(cut_res, dict):
            cut_sites = []
            for enzyme in cut_res:
                cut_sites.extend(cut_sites[enzyme])
            
            cut_sites.sort()

        cuts = np.r_[0, np.array(cut_sites) + 1, len(seq)].astype(int)
        n_frags = len(cuts) - 1
        
        fragments = pd.DataFrame(
            {'chrom': [chrom] * n_frags, 
            'start': cuts[:-1],
            'end': cuts[1:]},
            columns=['chrom', 'start', 'end']
        )

        return fragments

    fragments_res_df = pd.concat(map(_each, chroms), axis=0, ignore_index=True)
    
    return fragments_res_df

def restriction_site(enzyme):
    """
    Get the restriction site pattern.

    enzyme: str
        restriction enzyme, i.e. MboI, HindIII.
    
    Returns:
    --------
    list:
        list of ligation sites

    Examples:
    --------
    >>> ligation_site("MboI")
    ['GATC']
    >>> ligation_site("Arima")
    ['GATCGATC', 'GATCGANT', 'GANTGATC', 'GANTGANT']
    """
    from Bio import Restriction, Seq
    
    sites = []
    if isinstance(enzyme, str):
        if enzyme.lower() == 'arima':
            site_pattern = ['GATCGATC', 'GATCGANT', 
                            'GANTGATC', 'GANTGANT']
            return site_pattern
        elif enzyme.lower() == 'dnpii':
            enzyme = 'MboI'
            
        enzyme = getattr(Restriction, enzyme)
        site_pattern = enzyme.site
        sites.append(site_pattern)
    else:
        site_pattern = enzyme.site
        sites.append(site_pattern)

    return sites
    
def ligation_site(enzyme):
    """
    Get the ligation site.

    Params:
    --------
    enzyme: str
        restriction enzyme, i.e. MboI, HindIII.
    
    Returns:
    --------
    list:
        list of ligation sites

    Examples:
    --------
    >>> ligation_site("MboI")
    ['GATCGATC']
    >>> ligation_site("Arima")
    ['GATCGATC', 'GATCGANT', 'GANTGATC', 'GANTGANT']
    """
    from Bio import Restriction

    sites = []
    if isinstance(enzyme, str):
        if enzyme.lower() == 'arima':
            ligated_site = ['GATCGATC', 'GATCGANT', 
                            'GANTGATC', 'GANTGANT']
            return ligated_site
        elif enzyme.lower() == 'dnpii':
            enzyme = 'MboI'
            
        enzyme = getattr(Restriction, enzyme)
        cut_pattern = enzyme.elucidate()
        sites.append(cut_pattern)
    else:
        cut_pattern = enzyme.elucidate()
        sites.append(cut_pattern)

    ligated_sites = []
    for pattern in sites:
        left_side = []
        right_side = []
        for character in pattern:
            if not(len(left_side) > 0 and left_side[-1] == '_'):
                if not character == '^':
                    left_side.append(character)
            
            if character == '^' or len(right_side) > 0:
                if not character == '_':
                    right_side.append(character)
        
        left_side_re = ''.join(left_side[:-1]).strip("N")
        right_side_re = ''.join(right_side[1:]).strip("N")

        ligated_site = left_side_re + right_side_re
        ligated_sites.append(ligated_site)
    
    return ligated_sites

def get_genome_size(fasta):
    """
    Get the size of genome.

    Params:
    --------
    fasta: str
        genome in fasta format.
    
    Returns:
    --------
    int:
        size of genome
    
    Examples:
    --------
    >>> get_genome_size("sample.fasta")
    100000
    """
    from pyfaidx import Fasta

    if str(fasta)[-3:] == ".gz":
        handle = xopen(fasta)
        fasta = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        sizes = [len(record) for record in fasta]
    else:
        fasta = Fasta(fasta)
        sizes = [record.unpadded_len for record in fasta]
    
    return sum(sizes)

def get_contigs(fasta):
    """
    get contigs from fasta
    """
    length_db = get_contig_length(fasta)
    contigs = list(length_db.keys())

    return contigs

def get_contig_length(fasta):
    """
    get contig length database from fasta
    """
    from pyfaidx import Fasta
    fasta = Fasta(fasta)
    contigs = list(map(lambda x: x.name, fasta))
    lengths = list(map(len, list(fasta)))
    length_db = OrderedDict(zip(contigs, lengths))

    return length_db 

def get_contig_idx(fasta):
    """
    get contigs and it's idx
    """
    contigs = get_contigs(fasta)
    return dict(zip(contigs, range(len(contigs))))

def read_chrom_sizes(chrom_size):
    """
    read chrom size file

    Params:
    --------
    chrom_size: str
        path of chrom size (two columns)
    
    Returns:
    --------
    chrom_sizes: pd.DataFrame
        dataframe of chromsizes
    
    Examples:
    --------
    >>> df = read_chrom_sizes("sample.chromsizes")
    """
    # logger.info(f"Loading contig sizes `{chrom_size}`...")
    df = pd.read_csv(chrom_size, sep='\t', header=None, index_col=0, names=['chrom', 'length'],
                     usecols=[0, 1])

    return df 

def get_contig_size_from_fasta(fasta, force=False):
    fasta_path = Path(fasta)
    contigsizes = Path(f"{fasta_path.stem}.contigsizes")

    if contigsizes.exists() and force is False:
        logger.info(f"Load exists contigsizes: `{str(contigsizes)}`")
        return contigsizes
        
    cmd = ["cphasing-rs", "chromsizes", str(fasta_path), 
                "-o", str(contigsizes)]
    
    flag = run_cmd(cmd, log=os.devnull)
    assert flag == 0, "Failed to execute command, please check log."           

    return contigsizes  

def calculate_Nx_from_contigsizes(contigsizes: pd.DataFrame, Nx: int) -> int:
    """
    calculate N10-N90 from a contigsizes dataframe

    Params:
    --------
    contigsizes: pd.DataFrame

    Nx: int
        value of Nx
    
    Returns:
    --------
    int
    """
    
    total_length = contigsizes['length'].sum()
    
    length_array = contigsizes['length'].values
    length_array.sort()
    length_array = length_array[::-1]
    
    N_length = int(total_length * (Nx / 100))

    sum_length = 0
    for length in length_array:
        sum_length += length 
        if sum_length >= N_length:
            return length

def stat_contigsizes(contigsizes: pd.DataFrame) -> pd.DataFrame:
    """
    stat the summary of contigs fasta from a contigsizes
        include the Nx, min, max, total length.
    
    Params:
    ---------
    contigsizes: `pd.DataFrame`

    Returns:
    --------
    pd.DataFrame

    """

    Nx_items = list(range(10, 100, 10))
    Nx_results = []
    for Nx_item in Nx_items:
        Nx_results.append(calculate_Nx_from_contigsizes(contigsizes, Nx_item))
    
    min_length, max_length = min(contigsizes['length']), max(contigsizes['length'])
    total_length = sum(contigsizes['length'])

    Nx_items = [f"N{item}" for item in Nx_items]
    
    result_dict = OrderedDict(zip(Nx_items, Nx_results))
   
    result_dict['Max'] = max_length
    result_dict['Min'] = min_length
    result_dict['Total'] = total_length
    

    df = pd.DataFrame.from_dict(result_dict, orient='index')

    return df 
    
def parse_corrected_contig(contig):
    try:
        contig, ranges = contig.split(":")
        start, end = ranges.split("-")
        start, end = int(start), int(end)

        return contig, start, end
    
    except:
        return None 

def check_platform():
    import platform
    machine = platform.uname().machine

    return machine.split("_")[0]

def choose_software_by_platform(software):
    if check_platform() == 'arm':
        return f'{software}'
    else:
        return software

def check_allhic_version():
    from shutil import which
    allhic = choose_software_by_platform('allhic')
    assert which(allhic) is not None, "No such command of `allhic`" 
    
    for i in os.popen(f'{allhic} --version'):
        version = i.strip().split()[-1]
    
    def version_check(version):
        first_version, second_version, thrid_version =\
                list(map(int, version.split(".")))

        if first_version == 0:
            if second_version == 9:
                if thrid_version >= 14:
                    return True
                else:
                    return False 
            elif second_version < 9:
                return False
            else:
                return True
        else:
            return True
    
    assert version_check(version), "the version of `allhic` must be >= 0.9.14."


def read_fasta(fasta: str) ->OrderedDict:
    """
    Create a directory of fasta record.

    Params:
    --------
    fasta: str
        Path of fasta file.
    
    Returns:
    --------
    OrderedDict:
        {"seq_id", Seq()}
    
    Examples:
    --------
    >>> read_fasta('sample.fasta')
    OrderedDict(('ctg1', 1000))
    """
    from Bio import SeqIO
    
    logger.info(f'Load fasta file: `{fasta}`.')
    db = OrderedDict()
    
    fasta = SeqIO.parse(xopen(fasta), 'fasta')
    for record in fasta:
        if record.id not in db:
            db[record.id] = record.seq
        
    return db


def _zero_diags(chunk, n_diags):
    """
    zero diag in cool pixels table
    """
    if n_diags > 0:
        if n_diags > 0:
            mask = np.abs(chunk['pixels']['bin1_id'] 
                            - chunk['pixels']['bin2_id']) < n_diags
            chunk['pixels']['count'][mask] = 0

    return chunk

def delete_row_lil(mat, i):
    """
    https://stackoverflow.com/questions/13077527/is-there-a-numpy-delete-equivalent-for-sparse-matrices
    """
    from scipy.sparse import lil_matrix
    
    if not isinstance(mat, lil_matrix):
        raise ValueError("works only for LIL format -- use .tolil() first")
    mat.rows = np.delete(mat.rows, i)
    mat.data = np.delete(mat.data, i)
    mat._shape = (mat._shape[0] - 1, mat._shape[1])


def humanized2numeric(size):
    """
    convert the humanized to chromosome size
    >>> humanized2numeric("10k")
    10000
    """
    if not isinstance(size, str):
        logger.error("Please input correct string")
        raise ValueError("Value error")
    
    size = size.lower()
    size = size.replace(" ", "")

    try:
        if size[-1] == "k":
            res_size = int(size[:-1]) * 1000
        elif size[-1] == "m":
            res_size = int(size[:-1]) * 1000000
        elif size[-1] == "g":
            res_size = int(size[:-1]) * 1000000000
        else:
            res_size = int(size) 
    
    except ValueError:
        logger.error("Please input correct string")
        raise ValueError("Value error")
    
    return res_size

def to_humanized(size):
    """
    Convert the unit of chromosome size to suitable unit.
    >>> chrom_size_convert(100000)
    100 Kbp
    >>> chrom_size_convert(1000000)
    1 Mbp
    """
    size = int(size)
    if size <= 1e3:
        label = "{:,.0f}".format((size)) + ""
    elif size < 1e6:
        label = "{:,.0f}".format((size / 1e3)) + "k"
    elif size < 1e9:
        label = "{:,.1f}".format((size / 1e6)) + "m"
    
    else:
        label = "{:,.0f}".format((size / 1e9)) + "g"

    return label

def to_humanized2(size):
    """
    Convert the unit of chromosome size to suitable unit.
    >>> chrom_size_convert(100000)
    100 Kbp
    >>> chrom_size_convert(1000000)
    1 Mbp
    """
    size = int(size)
    if size <= 1e3:
        label = "{:,.0f}".format((size)) + ""
    elif size < 1e6:
        label = "{:,.0f}".format((size / 1e3)) + "k"
    elif size < 1e9:
        label = "{:,.0f}".format((size / 1e6)) + "m"
    
    else:
        label = "{:,.0f}".format((size / 1e9)) + "g"

    return label

def chrom_ticks_convert(ticks, add_suffix=True):
    """
    Convert a list of  chromosome size to suitable unit.
    >>> ticks = [10000, 20000, 30000]
    >>> chrom_ticks_convert(ticks)
    ['10', '20', '30Kbp']
    """
    if ticks[-1]  - ticks[1] <= 1e3:
        labels = ["{:,.0f}".format((x)) 
                  for x in ticks] 
        if add_suffix:
            labels[-1] += " bp"
    elif ticks[-1]  - ticks[1] <= 4e5:
        labels = ["{:,.0f}".format((x / 1e3)) 
                  for x in ticks]
        if add_suffix:
            labels[-1] += 'Kbp'
    else:
        labels = ["{:,.1f}".format((x / 1e6)) 
                  for x in ticks]
        if add_suffix:
            labels[-1] += " Mbp"
    
    return labels

def cool2depth(coolfile, output):
    cool = cooler.Cooler(coolfile)
    binsize = cool.binsize
    bins = cool.bins()[:]
    
    contigsizes = cool.chromsizes
    matrix = cool.matrix(balance=False, sparse=True)[:]
    sum_values = np.array(matrix.sum(axis=1).T[0])

    bins['values'] = sum_values[0]
    bins.to_csv(output, sep='\t', header=None, index=None)
    

def merge_matrix(coolfile, 
                outcool, 
                min_contacts=3,
                no_mask_nan=False,
                symmetric_upper=True,
                balance=False,
                no_dia=True):
    """
    merge slidewindows matrix into whole contig matrix.
    
    INPUT_COOL_PATH : Path to COOL file.

    OUTPUT_COOL_PATH : Path to output COOL file.

    """
    from scipy.sparse import coo_matrix, dia_matrix, triu

    logger.info(f'Load `{coolfile}` ...')
    cool = cooler.Cooler(coolfile)

    pixels = cool.matrix(balance=balance, sparse=True, 
                            as_pixels=True, join=True)[:]
    if balance:
        pixels['count'] = pixels['count'] * pixels['balanced']
    if no_mask_nan:
        matrix = cool.matrix(balance=balance, sparse=True)

    logger.info('Merging matrix ...')
    pix_counts = pixels.groupby(by=['chrom1', 'chrom2'], 
                            observed=True)['count'].sum().dropna()
    
    pix_counts = pix_counts.to_frame().reset_index()
    
    pix_counts = pix_counts[pix_counts['count'] >= min_contacts]
    
    ## reset chromosome into bin id
    pix_counts['chrom1'] = pix_counts['chrom1'].cat.codes.values
    pix_counts['chrom2'] = pix_counts['chrom2'].cat.codes.values
    pix_counts.columns = ['bin1_id', 'bin2_id', 'count']
    
    ## create sparse matrix 
    row = pix_counts['bin1_id']
    col = pix_counts['bin2_id']
    data = pix_counts['count']
    
    m, n = len(cool.chromnames), len(cool.chromnames)
    matrix = coo_matrix((data, (row, col)), shape=(m, n))
    dia = dia_matrix(([matrix.diagonal()], [0]), shape=matrix.shape)
    
    if no_dia:
        matrix = matrix + matrix.T - (2 * dia)
    else:
        matrix = matrix + matrix.T - dia
    
    if symmetric_upper:
        matrix = triu(matrix).tocoo()
    
        new_pixels = dict(zip(['bin1_id', 'bin2_id', 'count'],
                            [matrix.row, matrix.col, matrix.data]))
    
    ## new bins 
    new_bins = pd.DataFrame(cool.chromsizes)
    new_bins = new_bins.reset_index()
    new_bins.columns = ['chrom', 'end']
    new_bins['start'] = 0
    new_bins = new_bins[['chrom', 'start', 'end']]
    
    cooler.create_cooler(outcool, new_bins, new_pixels, 
                            dtypes=dict(count='float64'),
                            symmetric_upper=symmetric_upper)
    logger.info(f'Successful merge matrix into whole contig: `{outcool}`')


def extract_matrix(coolfile, chrom_list, out):
    """
    Extract matrix from a cool by chrom list.
    """
    from hicmatrix import HiCMatrix as hm
    num = len(chrom_list)
    hic_matrix = hm.hiCMatrix(coolfile)
    hic_matrix.reorderChromosomes(chrom_list)
    hic_matrix.save(out)
    logger.info(f'Extract {num} chromosomes into `{out}`.') 

def prune_matrix(coolfile, prunepairs, outcool):
    """
    Prune matrix by a chrom/contig pairs.
    """
    cool = cooler.Cooler(coolfile)
   
    chrom_index = cool.bins()[:]['chrom'].reset_index()
    chrom_index.set_index('chrom', inplace=True)
    chrom_index = chrom_index.to_dict()['index']
    bins = cool.bins()[:]
    bins.set_index('chrom')
    pairs = list(set(prunepairs))

    new_pairs = []
    for i, pair in enumerate(pairs):
        try:
            if chrom_index[pair[0]] > chrom_index[pair[1]]:
                pairs = pairs[::-1] 
        except KeyError:
            continue

        new_pairs.append(pair)

    new_pairs = list(set(new_pairs))

    pixels = cool.matrix(balance=False, sparse=True, as_pixels=True, join=True)
    pixels = pixels[:]

    pixels = pixels.reset_index(drop=False).set_index(['chrom1', 'chrom2'])

    pruned_idx = []
    for pair in new_pairs:
        try:
            pruned_idx.extend(pixels.loc[pair]['index'].values.tolist())
        except KeyError:
            continue

    pixels = cool.matrix(balance=False, sparse=True, as_pixels=True)[:]

    pixels.loc[pruned_idx] = 0
    pixels = pixels[pixels['count'] != 0]
  
    cooler.create_cooler(outcool, cool.bins()[:], pixels)


def trim_axes(axes, N):
    """
    little helper to message the axs list to have correct length...

    Params:
    --------
    axes: `list`
            list of axes
    N: `int`
            number of axes to return
    
    Returns:
    --------
    axes: `list`
            list of trimed axes
    
    Examples:
    --------
    >>> fig, axes = plt.subplots(5, 2)
    >>> axes = trim_axes(axes, 7)
    array([<matplotlib.axes._subplots.AxesSubplot object at 0x7f5f48365198>,
       <matplotlib.axes._subplots.AxesSubplot object at 0x7f5f49275438>,
       <matplotlib.axes._subplots.AxesSubplot object at 0x7f5f4712e320>, ...]
    """

    axes = axes.flat
    for ax in axes[N:]:
        ax.remove()
    
    return axes[:N]

def generate_to_hic_cmd(agp, pairs, n=0, _3ddna_path="~/software/3d-dna", output="to_hic.cmd.sh"):

    pairs_prefix = str(Path(pairs).name).replace(".gz", "").replace(".pairs", "")
    agp_prefix = str(Path(agp).name).replace(".agp", "")
    cmd = f"""
## Please submit the following command yourself

_3ddna_path={_3ddna_path}
cphasing-rs pairs2mnd {pairs} -o {pairs_prefix}.mnd.txt
cphasing utils agp2assembly {agp} -o {agp_prefix}.assembly
bash $_3ddna_path/visualize/run-assembly-visualizer.sh -p true {agp_prefix}.assembly {pairs_prefix}.mnd.txt

## After curation, convert review.assembly to new agp 
# cphasing utils assembly2agp groups.review.assembly -o groups.review -n {n}
# cphasing utils agp2fasta groups.review.agp -o groups.review.fasta
    """

    with open(output, 'w') as out:
        out.write(cmd)

