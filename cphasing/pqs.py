#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Library of .pqs file format, which is dictory contain several parquet files and a _contigsizes file.
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

import polars as pl 

from pathlib import Path
from .utilities import read_chrom_sizes

_README = """
# .pqs file format
The _contigsizes file contains the size of each contig in the .pqs file.
The _metadata file contains the metadata of the .pqs file.
The q0, q1 directory contains the parquet files of the data.

q0 mean the mapping quality of data >= 0.
q1 mean the mapping quality of data >= 1.

.pqs/  
|-- _contigsizes 
|-- _metadata
|-- _readme 
|-- q0/
|   |-- q0_0.parquet
|   |-- q0_1.parquet
|   |-- ...
|-- q1/
|   |-- q1_0.parquet
|   |-- q1_1.parquet
|   |-- ...
"""

class PQS:
    """
    A class to represent a .pqs file.
    """
    def __init__(self):
        self._readme = _README 
        
        pass 
    
    @property
    def contigsizes(self):
        """
        Get the contig sizes of the .pqs file.
        """

        return self._contigsizes

    @property
    def metadata(self):
        """
        Get the metadata of the .pqs file.
        """

        return self._metadata
    
    @property
    def readme(self):
        """
        Get the readme of the .pqs file.
        """
        return self._readme
    
    @property
    def q0(self):
        """
        Get the q0 directory of the .pqs file.
        """
        return self._q0
    
    @property
    def q1(self):
        """
        Get the q1 directory of the .pqs file.
        """
        return self._q1
    
    def read(self, path):
        """
        Read the .pqs file from the given path.
        """
        pass
    
    def write(self, path):
        """
        Write the .pqs file to the given path.
        """
        pass


    def from_pairs(self, pairs, path):
        """
        Create a .pqs file from the given pairs.
        """
        pass

    def from_porec_table(self, porec, path):
        """
        Create a .pqs file from the given porec table.
        """
        pass
