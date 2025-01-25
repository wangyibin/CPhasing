#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
recluster the contigs on existing clusters
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

from ..core import AlleleTable, PruneTable, ClusterTable


class Recluster:
    """
    Developping...
    """
    def __init__(self, clustertable, contacts, countre, prunetable, min_length=50000):
        self.clustertable = clustertable
        self.contacts = contacts
        self.countre = countre
        self.prunetable = prunetable

        self.min_length = min_length
    
    


