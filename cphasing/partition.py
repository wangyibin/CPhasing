#!/usr/bin/env python
# -*- coding:utf-8 -*-


import logging
import os
import os.path as op
import sys

import glob

from .utilities import run_cmd

logger = logging.getLogger(__name__)

def get_partition_res(count_re, k, targetDir):

    prefix = count_re.replace('.txt', '')
    targets = glob.glob('{}/{}.{}g*txt'.format(targetDir, prefix, k))
    targets = sorted(targets)
    
    return targets

class Partitioner:
    def __init__(self, count_re, pairtable, k, 
                maxLinkDensity=2, minREs=10, 
                nonInformativeRatio=3):
                self.count_re = count_re
                self.pairtable = pairtable
                self.k = k
                self.maxLinkDensity = maxLinkDensity
                self.minREs = minREs
                self.nonInformativeRatio = nonInformativeRatio
    
    @staticmethod 
    def partition(count_re, pairtable, k, 
                maxLinkDensity=2, minREs=10, 
                nonInformativeRatio=3):

        cmd = ['allhic', 'partition', count_re, 
                pairtable, str(k), 
                '--maxLinkDensity', str(maxLinkDensity),
                '--minREs', str(minREs),
                '--nonInformativeRatio', str(nonInformativeRatio)]

        run_cmd(cmd)

        targets = get_partition_res(count_re, k, "./")

        return targets 



class AdativePartitioner(Partitioner):
    def __init__(self, *args, **kwargs):
        super(Partitioner, self).__init__(*args, **kwargs)

    

