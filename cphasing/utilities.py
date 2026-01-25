#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
utility libraries
"""

import logging
import os
import os.path as op
import shutil
import sys
import re
import time

import cooler
import numpy as np
import pandas as pd

from Bio import SeqIO 
from collections import OrderedDict
from pathlib import Path, PosixPath

from psutil import Process, NoSuchProcess
from threading import Thread

logger = logging.getLogger(__name__)

## https://joblib.readthedocs.io/en/latest/auto_examples/parallel_generator.html#sphx-glr-auto-examples-parallel-generator-py
class MemoryMonitor(Thread):
    """Monitor the memory usage in MB in a separate thread.

    Note that this class is good enough to highlight the memory profile of
    Parallel in this example, but is not a general purpose profiler fit for
    all cases.
    """
    def __init__(self):
        super().__init__()
        self.stop = False
        self.memory_buffer = []
        self.daemon = True
        self.start()

    def get_memory(self):
        "Get memory of a process and its children."
        p = Process()
        memory = p.memory_info().rss
        for c in p.children():
            memory += c.memory_info().rss
        return memory

    def run(self):
        memory_start = self.get_memory()
        while not self.stop:
            try:
                self.memory_buffer.append(self.get_memory() - memory_start)
            except NoSuchProcess:
                continue 
            time.sleep(0.2)

    def join(self):
        self.stop = True
        super().join()


def export_log(logfile):
    from . import console, console_html
    from rich.terminal_theme import MONOKAI
    
    try:
        console.save_text(logfile.replace(".txt", ""), styles=True)
        console_html.save_text(logfile)
        # out_html = logfile.replace(".txt", ".html")
        # console_html.save_html(out_html, theme=MONOKAI)

    except AssertionError:
        pass 

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

def choose_compressor():
    if not cmd_exists('crabz') and not cmd_exists('pigz'):
            raise ValueError(f"pigz or crabz: command not found")
    
    if cmd_exists('crabz'):
        return 'crabz'
    else:
        if cmd_exists('pigz'):
            return 'pigz'
        else:
            raise ValueError(f"pigz: command not found")

def compress_cmd(infile, threads=10):
    compressor = choose_compressor()

    if not Path(infile).exists():
        raise FileNotFoundError(f"File `{infile}` not found")

    if not Path(infile).is_file():
        raise ValueError(f"File `{infile}` is not a file")

    if compressor == 'crabz':
        cmd = ["crabz", "-p", str(threads), infile]
    else:
        cmd = ["pigz", "-c", "-p", str(threads), infile]

    return cmd 

def decompress_cmd(infile, threads=10):
    compressor = choose_compressor()

    if infile.endswith(".gz"):
        if compressor == "crabz":
            cmd = ["crabz", "-d", "-p", str(threads), infile]
        else:
            cmd = ["pigz", "-c", "-d", "-p", str(threads), infile]
    else:
        cmd = ["cat", infile]

    return cmd


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
    from .pqs import PQS
    if Path(_file).is_dir():
        p = PQS()
        if p.is_pqs(_file):
            return False
        else:
            return True
    else:
        try:
            pd.read_csv(_file, chunksize=1, comment="#")
        except EOFError:
            return True 
        except pd.errors.EmptyDataError:
            return True 
        except UnicodeDecodeError:
            raise ValueError(f"File `{_file}` is compressed file, you should add `.gz` suffix.")

    return False

def is_file_changed(input_file):
    from .pqs import PQS
    if input_file is None:
        return False
    
    if not Path(input_file).exists():
        return True 

    if Path(input_file).is_dir():
        p = PQS()
        if p.is_pqs(input_file):
            return p.is_changed(input_file)
        else:
            return True

    input_file_path = Path(input_file).absolute().parent
    file_name = Path(input_file).name
    prefix = Path(input_file).stem

    if Path(f"{input_file_path}/.{prefix}.md5.txt").exists():
        text = os.popen(f"md5sum -c {input_file_path}/.{prefix}.md5.txt 2>/dev/null").read()
        if not text.strip():
            return True
        if text.strip().split()[-1] == "OK":
            return False
        else:
            os.system(f"md5sum {input_file_path}/{file_name} > {input_file_path}/.{prefix}.md5.txt")
            return True
    else:
        os.system(f"md5sum {input_file_path}/{file_name} > {input_file_path}/.{prefix}.md5.txt")
        return True


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
        # try:
        #     import mgzip 
        # except ImportError:
        handle = gzip.open(infile, mode + 't')
        # finally:
            # handle = mgzip.open(str(infile), mode + 't', thread=8)
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

    try:
        from needletail import parse_fastx_file, NeedletailError
    except ImportError:
        handle = xopen(fasta)
        fasta = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        sizes = [len(fasta[record]) for record in fasta]
    else:
        try:
            sizes = []
            for record in parse_fastx_file(fasta):
                sizes.append(len(record.seq))
        except NeedletailError:
            handle = xopen(fasta)
            fasta = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            sizes = [len(fasta[record]) for record in fasta]
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


def read_fasta(fasta: str) -> OrderedDict:
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
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    from Bio.Seq import Seq

    logger.info(f'Load fasta file: `{fasta}`.')
    db = OrderedDict()
    
    try:
        from needletail import parse_fastx_file, NeedletailError
    except ImportError:
        with xopen(fasta) as handle:
            for title, seq in SimpleFastaParser(handle):
                seq_id = title.split(None, 1)[0]
                if seq_id not in db:
                    db[seq_id] = Seq(seq)
    else:
        try:
            for record in parse_fastx_file(fasta):
                seq_id = record.name
                if seq_id not in db:
                    db[seq_id] = Seq(record.seq)
        except NeedletailError:
            logger.error(f"NeedletailError: File `{fasta}` may be corrupted.")
            sys.exit(-1)
 
    return db

def read_fasta_yield(fasta: str):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    from Bio.Seq import Seq
    logger.info(f'Load fasta file: `{fasta}`.')
    try:
        from needletail import parse_fastx_file, NeedletailError
    except ImportError:
        with xopen(fasta) as handle:
            for title, seq in SimpleFastaParser(handle):
                seq_id = title.split(None, 1)[0]
                yield seq_id, Seq(seq)
    else:
        try:
            for record in parse_fastx_file(fasta):
                seq_id = record.name
                yield seq_id, Seq(record.seq)
        except NeedletailError:
            logger.error(f"NeedletailError: File `{fasta}` may be corrupted.")
            sys.exit(-1)


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
    if isinstance(size, int):
        return size
    elif isinstance(size, float):
        return int(size)

    elif not isinstance(size, str):
        logger.error("Please input correct string")
        raise ValueError("Value error")
    
    size = size.lower()
    size = size.replace(" ", "")

    try:
        if size[-1] == "k":
            res_size = float(size[:-1]) * 1000
        elif size[-1] == "m":
            res_size = float(size[:-1]) * 1000000
        elif size[-1] == "g":
            res_size = float(size[:-1]) * 1000000000
        else:
            res_size = float(size) 
    
    except ValueError:
        logger.error("Please input correct string")
        raise ValueError("Value error")
    
    res_size = int(res_size)

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
    if size < 1e3:
        label = "{:,.0f}".format((size)) + ""
    elif size < 1e6:
        label = "{:,.0f}".format((size / 1e3)) + "k"
    elif size < 1e9:
        label = "{:,.0f}".format((size / 1e6)) + "m"
    
    else:
        label = "{:,.0f}".format((size / 1e9)) + "g"

    return label

def to_humanized3(size):
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
        label = "{:,.2f}".format((size / 1e3)) + "k"
    elif size < 1e9:
        label = "{:,.2f}".format((size / 1e6)) + "m"
    
    else:
        label = "{:,.2f}".format((size / 1e9)) + "g"

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
    elif ticks[-1]  - ticks[1] <= 1e10:
        labels = ["{:,.1f}".format((x / 1e6)) 
                  for x in ticks]
        if add_suffix:
            labels[-1] = f"      {labels[-1]} Mbp"
    else:
        labels = ["{:,.2f}".format((x / 1e9)) 
                  for x in ticks]
        if add_suffix:
            labels[-1] = f"      {labels[-1]} Gbp"
    
    return labels

def cool2depth(coolfile, output, remove_cis=False):
    cool = cooler.Cooler(coolfile)
    binsize = cool.binsize
    
    bins = cool.bins()[:]
    
    contigsizes = cool.chromsizes
    matrix = cool.matrix(balance=False, sparse=True)[:]
    
    matrix = matrix + matrix.T
    if remove_cis:
        cis_bins = bins.reset_index().groupby('chrom').apply(lambda x: x.index.tolist())

        ## calculate cis sum

        matrix = matrix.tocsr()
        for chrom, idxes in cis_bins.items():
            res = matrix[idxes, :][:, idxes].sum(axis=1)
            bins.loc[idxes, 'values'] = res

        matrix = matrix.tocoo()

        sum_values = np.array(matrix.sum(axis=1).T[0])
        bins['values2'] = sum_values[0]
        bins['values'] = bins['values2'] - bins['values']
        bins = bins.drop(columns=['values2'])

    else:
        sum_values = np.array(matrix.sum(axis=1).T[0])
        sum_values = sum_values - matrix.diagonal() / 2
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

def pretty_cmd(cmd_list, n=4):
    new_cmd = []

    for i, item in enumerate(cmd_list):
        if i % n == 0 and i != 0:
            new_cmd.append("\\\n    ")
            
        new_cmd.append(item)
    
    return list(map(str, new_cmd ))


def generate_to_hic_cmd(agp, fasta, pairs, mapq=1, n=0, _3ddna_path="~/software/3d-dna", 
                        output="to_hic.cmd.sh"):

    pairs_prefix = str(Path(pairs).name).replace(".gz", "").replace(".pairs", "").replace(".pqs", "")
    agp_prefix = str(Path(agp).name).replace(".agp", "")
    cmd = f"""#!/usr/bin/bash
## Please submit the following command in yourself

_3ddna_path={_3ddna_path}
min_quality={mapq}

cphasing-rs pairs2mnd -q ${{min_quality}} {pairs} -o {pairs_prefix}.mnd.txt
cphasing utils agp2assembly {agp} -o {agp_prefix}.assembly
bash $_3ddna_path/visualize/run-assembly-visualizer.sh -p true -c -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000 {agp_prefix}.assembly {pairs_prefix}.mnd.txt

## After curation, convert review.assembly to new agp 
# cphasing utils assembly2agp groups.review.assembly -o groups.review -n {n}
# cphasing utils agp2fasta groups.review.agp {fasta} -o groups.review.fasta
    """


    if output is not None:
        with open(output, 'w') as out:
            out.write(cmd)
    else:
        return cmd


def generate_curation_cmd(agp, fasta,
                        pairs, binsize, cool_binsize,
                        mapq=1, scaffolding_dir="4.scaffolding", 
                        plot_dir="5.plot",
                        n=0, _3ddna_path="~/software/3d-dna", 
                        output="curation.cmd.sh"):

    pairs_prefix = str(Path(pairs).name).replace(".gz", "").replace(".pairs", "").replace(".pqs", "")
    agp_prefix = str(Path(agp).name).replace(".agp", "")
    cmd = f"""#!/usr/bin/bash
## Please submit the following command in yourself

_3ddna_path={_3ddna_path}
min_quality={mapq}

## Step 1 [Optional]: sort chromosomes based on interaction map, which may help the curation
cphasing sort-chromosomes -a ../{str(scaffolding_dir)}/{agp} -m ../{str(plot_dir)}/{agp_prefix}.{pairs_prefix}.q{mapq}.{to_humanized2(binsize)}.chrom.cool -o {agp_prefix}.sorted.agp
cphasing plot -a {agp_prefix}.sorted.agp -m ../{str(plot_dir)}/{pairs_prefix}.q{mapq}.{to_humanized2(cool_binsize)}.cool -o {agp_prefix}.sorted.{pairs_prefix}.q${{min_quality}}.{to_humanized2(binsize)}.wg.png -bs {to_humanized2(binsize)} --add-hap-border --no-lines --disable-natural-sort -oc 

## Step 2: visualize the interaction map and curate the assembly
cphasing-rs pairs2mnd -q ${{min_quality}} {pairs} -o {pairs_prefix}.mnd.txt
cphasing utils agp2assembly {agp_prefix}.sorted.agp -o {agp_prefix}.sorted.assembly

### Step2.1: you can seperately curation each homologous group
cphasing utils split-agp {agp_prefix}.sorted.agp -o separate_groups

cd separate_groups
for group_agp in *.agp; do
    cphasing utils agp2assembly $group_agp -o $(basename $group_agp .agp).assembly
    echo "bash $_3ddna_path/visualize/run-assembly-visualizer.sh -p true -c -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000 $(basename $group_agp .agp).assembly ../{pairs_prefix}.mnd.txt" 
done > to_hic_separate.cmd.sh
cd ..

### Step2.2: or directly curate the whole assembly

echo "bash $_3ddna_path/visualize/run-assembly-visualizer.sh -p true -c -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000 {agp_prefix}.sorted.assembly {pairs_prefix}.mnd.txt" > to_hic.cmd.sh
echo
echo -e "\e[1mPlease submit the command 'to_hic.cmd.sh' or 'seperate_groups/to_hic_separate.cmd.sh' for manually adjust assembly in Juicbox\e[0m"



## Step 3: After curation, convert review.assembly to new agp 
# cphasing utils assembly2agp groups.sorted.review.assembly -o groups.review -n {n}
# cphasing utils agp2fasta groups.review.agp {fasta} -o groups.review.fasta
    """


    if output is not None:
        with open(output, 'w') as out:
            out.write(cmd)
        
        return output
    else:
        return cmd


def generate_plot_cmd(pairs, pairs_prefix, contigsizes, agp, 
                        min_quality, init_binsize,
                      binsize, colormap, whitered, mode,
                      output="plot.cmd.sh"):

    out_heatmap =  f"groups.q${{min_quality}}.${{heatmap_binsize}}.wg.png"
    plot_another_args = ""
    if colormap != "redp1_r":
        plot_another_args += f" -cm {colormap}"
    if whitered:
        plot_another_args += " --whitered"

    if mode == "phasing":
        plot_another_args += " --add-hap-border --no-lines"
    
    cmd = f"""#!/usr/bin/bash

min_quality={min_quality}
cool_binsize={to_humanized2(init_binsize)}
heatmap_binsize={to_humanized2(binsize)}

if [ ! -f {pairs_prefix}.q${{min_quality}}.${{cool_binsize}}.cool ]; then
    cphasing pairs2cool {pairs} \\
        {contigsizes} {pairs_prefix}.q${{min_quality}}.${{cool_binsize}}.cool \\
        -q ${{min_quality}} -bs ${{cool_binsize}}
fi
cphasing plot -a {agp} \\
    -m {pairs_prefix}.q${{min_quality}}.${{cool_binsize}}.cool \\
        -o {out_heatmap} -bs ${{heatmap_binsize}} -oc {plot_another_args}
    """

    with open(output, 'w') as out:
        out.write(cmd)

def recommend_binsize_by_genomesize(genomesize):
    """
    Recommend the binsize by genome size.
    """
    
    if genomesize < 5e8:
        binsize = 100000 
        init_binsize = 10000

    elif genomesize < 5e9:
        binsize = 500000
        init_binsize = 10000
    
    elif genomesize < 2e10:
        binsize = 1000000
        init_binsize = 20000
    
    elif genomesize < 5e10:
        binsize = 5000000
        init_binsize = 100000

    else:
        binsize = 10000000
        init_binsize = 200000
    
    return init_binsize, binsize



def is_pairs_with_mapq(pairs):
    """
    Check if the pairs file contains mapq.
    """
    with xopen(pairs) as fp:
        for line in fp:
            if line.startswith("#"):
                continue
            line_list = line.strip().split("\t")
            if len(line_list) <= 7:
                return False 
            else:
                if line_list[7].isdigit():
                    return True
                else:
                    return False

def binnify(chromsizes, binsize, trim_length=0):
    """
    Divide a genome into evenly sized bins.

    Parameters
    ----------
    chromsizes : Series
        pandas Series indexed by chromosome name with chromosome lengths in bp.
    binsize : int
        size of bins in bp
    trim_length: int
        length to trim at the termini of the contig
    Returns
    -------
    bins : `pandas.DataFrame`
        Dataframe with columns: ``chrom``, ``start``, ``end``.

    """

    chrom_list = []
    start_list = []
    end_list = []
    binsize = int(binsize)
    for chrom, clen in sorted(chromsizes.items(), key=lambda x: x[0]):
        n_bins = int(np.ceil(clen / binsize))
        binedges = np.arange(0, (n_bins + 1)) * binsize
        binedges[-1] = clen
        chrom_list.extend([chrom] * n_bins)
        if trim_length:
            if clen > trim_length * 3:
                # binedges[0] = binedges[0] + trim_length
                # if (binedges[-1] - binedges[-2]) > trim_length:
                #     binedges[-2] = binedges[-2] - trim_length
                binedges = np.clip(binedges, trim_length, clen - trim_length)

        start_list.extend(binedges[:-1])
        end_list.extend(binedges[1:])
        

    bintable = pd.DataFrame({
        "chrom": chrom_list,
        "start": start_list,
        "end": end_list
    })

    # if trim_length:
    #     bintable = bintable[bintable['start'] < bintable['end']]

    # bintable = bintable.sort_values(["chrom", "start"]).reset_index(drop=True)
    bintable["chrom"] = pd.Categorical(
        bintable["chrom"], categories=list(chromsizes.index), ordered=True
    )

    return bintable


def porec_downsample(porec_table, n, threads=10):
    import polars as pl

    os.environ['POLARS_MAX_THREADS'] = str(threads)

    HEADER = ["read_idx", "read_length",
              "read_start", "read_end", 
              "strand", "chrom", "start", 
              "end", "mapping_quality", 
              "identity", "filter_reason"]
    df = pl.read_csv(porec_table, separator='\t', 
                        has_header=False, 
                        new_columns=HEADER)

    total_read_idxes = df.unique("read_idx")['read_idx']
   
    total_read_count = len(total_read_idxes)
    df_grouped_nuinque = df.groupby(['read_idx']).agg(pl.col("chrom").count().alias("n_unique"))

    df_grouped_nuinque_db = df_grouped_nuinque.to_pandas()
    df_grouped_nuinque_db.set_index('read_idx', inplace=True)
    df_grouped_nuinque_db = df_grouped_nuinque_db.to_dict()['n_unique']
    total_read_idxes = total_read_idxes.to_pandas().to_frame()

    total_read_idxes['count'] = total_read_idxes['read_idx'].map(df_grouped_nuinque_db.get)
    total_pairs = (total_read_idxes['count'] * (total_read_idxes['count'] - 1) / 2).sum()
    median_count = int(total_read_idxes['count'].median())
    random_state = np.random.RandomState(12345)

    for n in listify(n):
        if n == 0:
            continue
        logger.debug(f"Total read count: {total_read_count}")
        if (n // median_count) >= total_read_count:
            read_idxes = total_read_idxes
        else:
            read_idxes = total_read_idxes.sample(n // median_count, random_state=random_state)
        read_idxes = read_idxes.sort_values(by='read_idx')

        read_idxes.set_index('read_idx', inplace=True)
        

        m = int((read_idxes['count'] * (read_idxes['count'] - 1) // 2).sum())
    
        logger.info(f"Initially generated {m} pairs")
        while  m > n and not read_idxes.empty:
            old_m = m
            
            batch_size = (m - n) // 10
            logger.debug(f"Batch size {batch_size}")
            read_idx = read_idxes.sample(batch_size, random_state=random_state).index
            logger.debug(f"Remove {len(read_idx)}")
            read_idxes.drop(read_idx, inplace=True)
        
            m = int((read_idxes['count'] * (read_idxes['count'] - 1) / 2).sum())
            logger.debug(f"Generated {m} pairs")
            if m == old_m:
                break

        read_idxes = set(read_idxes.index.tolist())

        _df = df.filter(pl.col("read_idx").is_in(read_idxes))
        with xopen(f"random.{n}.{porec_table}", 'w') as out:
            _df.write_csv(out, separator='\t', include_header=False)


def pairs_pqs_downsample(pairs, n, min_mapq=0, threads=10):
    import polars as pl
    pl.enable_string_cache()
    from joblib import Parallel, delayed
    from pytools import natsorted
    from .pqs import PQS 
    from .core import Pairs2

    os.environ["POLARS_MAX_THREADS"] = str(threads)
    
    if Path(pairs).is_dir():
        p = PQS(pairs)
        p.init_read()
        chunksize = p._metadata['chunksize']
        chunks = natsorted(p.read(pairs, min_mapq=min_mapq, return_as="files"), key=str)

        if min_mapq == 0:
           
            last_file = pl.read_parquet(chunks[-1])
            last_file_length = len(last_file)

            length = chunksize * len(chunks[:-1]) + last_file_length

            np.random.seed(12345)
            for _n in listify(n):
                if length > _n:
                    
                    idxes = np.random.choice(length, size=_n, replace=False)
                    idxes.sort()

                    output_dir = f"random.{_n}.pairs.pqs"
                    Path(output_dir).mkdir(exist_ok=True)
                    Path(f"{output_dir}/q0").mkdir(parents=True, exist_ok=True)
                    Path(f"{output_dir}/q1").mkdir(parents=True, exist_ok=True)

                    shutil.copy(f"{pairs}/_metadata", output_dir)
                    shutil.copy(f"{pairs}/_readme", output_dir)
                    shutil.copy(f"{pairs}/_contigsizes", output_dir)
                    

                    def process_parquet(parquet_path, idxes, chunksize, output_dir):
                        parquet_idx = int(Path(parquet_path).stem)
                        local_idxes = idxes[(idxes >= parquet_idx * chunksize) & (idxes < (parquet_idx + 1) * chunksize)] - (parquet_idx * chunksize)
                        df = pl.read_parquet(parquet_path)
                        res_df = df[local_idxes]

                        res_df.write_parquet(f"{output_dir}/q0/{parquet_idx}.parquet")
                        res_df.filter(pl.col("mapq") > 0).write_parquet(f"{output_dir}/q1/{parquet_idx}.parquet")

                    args = [(parquet, idxes, chunksize, output_dir) for parquet in chunks]

                    Parallel(n_jobs=threads)(
                        delayed(process_parquet)(*i) for i in args
                    )
                                    
                else:
                    logger.warning(f"The total length of the input pairs is not enough to sample {_n} records, exit")
                    sys.exit(-1)
                    
        else:
            df = pl.scan_parquet(f"{pairs}/q1/*.parquet")
            length = len(df.collect())

            np.random.seed(12345)

            for _n in listify(n):
                if length > _n:
                    idxes = np.random.choice(length, size=_n, replace=False)
                    idxes.sort()

                    output_dir = f"random.{_n}.pairs.pqs"
                    Path(output_dir).mkdir(exist_ok=True)
                    Path(f"{output_dir}/q0").mkdir(parents=True, exist_ok=True)
                    Path(f"{output_dir}/q1").mkdir(parents=True, exist_ok=True)

                    shutil.copy(f"{pairs}/_metadata", output_dir)
                    shutil.copy(f"{pairs}/_readme", output_dir)
                    shutil.copy(f"{pairs}/_contigsizes", output_dir)

                    for parquet in chunks:
                        df = pl.read_parquet(parquet)
                        length = len(df)
                        res_df = df[idxes[idxes < length]]

                        res_df.write_parquet(f"{output_dir}/q0/{Path(parquet).stem}.parquet")
                        res_df.write_parquet(f"{output_dir}/q1/{Path(parquet).stem}.parquet")

                        idxes = idxes - length
                        idxes = idxes[idxes >= 0]
                else:
                    logger.warning(f"The total length of the input pairs is not enough to sample {_n} records, exit")
                    sys.exit(-1)
                
    else:

        p = Pairs2(pairs) 
        
        for n in listify(n):
            if n == 0:
                continue
        
            cmd = ["cphasing-rs", "pairs-downsample", str(pairs),
                    "-n", str(n), "-o", f"random.{n}.{pairs}"] 
            
            flag = run_cmd(cmd)



def parse_split_contigs(contig):
    contig, _range = contig.rsplit("|", 1)
    start, end = _range.split("_")
    start, end = int(start), int(end)

    return contig, start, end

def get_fasta_from_split_contig(fasta, contigsizes, output):
    # fasta_dict = read_fasta(fasta)
    
    # raw_contigs = contigsizes.index.tolist()
    # contigs = contigsizes.index.map(parse_split_contigs).tolist()
    # with open(output, 'w') as out:
    #    for raw_contig, (contig, start, end) in zip(raw_contigs, contigs):
    #        if start == end:
    #             continue
    #        out.write(f">{raw_contig}\n{fasta_dict[contig][start: end]}\n")

    from Bio.SeqIO.FastaIO import SimpleFastaParser
    from collections import defaultdict

    regions = defaultdict(list)

    for raw_contig in contigsizes.index:
        contig, start, end = parse_split_contigs(raw_contig)
        if start != end:
            regions[contig].append((raw_contig, start, end))

    with xopen(fasta) as handle, open(output, 'w') as out:
        for title, seq in SimpleFastaParser(handle):
            contig = title.split(None, 1)[0]
            
            if contig in regions:
                for raw_contig, start, end in regions[contig]:
                    sub_seq = seq[start:end]
                    out.write(f">{raw_contig}\n{sub_seq}\n")
                
    
    return output 


def determine_split_length(
    contigsizes_df: pd.DataFrame,
    max_split_len: int = 20_000_000,
    min_split_len: int = 1_000_000,
    outlier_factor: float = 5.0 # 3.0
):
    if contigsizes_df.empty:
        logger.warning("Contig sizes data is empty. Cannot determine split_length.")
        return None

    lengths = contigsizes_df['length']
    
    if len(lengths) < 10: 
        logger.info("Too few contigs to perform automatic splitting.")
        return None

    max_len = lengths.max()
    p95_len = lengths.quantile(0.95)
    median_len = lengths.median()

    logger.info(f"Contig length stats: Max={max_len/1e6:.2f}Mb, 95th-percentile={p95_len/1e6:.2f}Mb, Median={median_len/1e3:.2f}Kb.")

   
    if max_len > p95_len * outlier_factor and max_len > max_split_len:

        split_len = int(min(p95_len, max_split_len))
        split_len = int(max(split_len, min_split_len))
        split_len = (split_len // 1000_000 + 1) * 1000_000 
        logger.info(f"Outlier contig detected (max: {max_len/1e6:.2f}Mb). Activating splitting.")
        logger.info(f"Automatic split_length set to: {split_len/1e6:.2f}Mb.")
        return split_len
    else:
        logger.info("Contig length distribution is relatively uniform. Deactivating splitting.")
        return None
