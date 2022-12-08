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

import numpy as np
import pandas as pd

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
        else:
            return p.returncode
    return 0


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
    fasta = Fasta(fasta)
    sizes = [record.unpadded_len for record in fasta]
    
    return sum(sizes)

def check_allhic_version():
    from shutil import which
    assert which('allhic') is not None, "No such command of `allhic`" 
    
    for i in os.popen('allhic --version'):
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


def read_fasta(fasta: str) -> dict:
    """
    Create a directory of fasta record.

    Params:
    --------
    fasta: str
        Path of fasta file.
    
    Returns:
    --------
    dict:
        {"seq_id", Seq()}
    
    Examples:
    --------
    """
    from Bio import SeqIO
    
    logger.info(f'Load fasta file: `{fasta}`.')
    db = {}
    
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