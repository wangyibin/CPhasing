#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
cut and trim hic reads on ligation site
"""

import logging
import os
import os.path as op
import sys

import re

from Bio import SeqIO, Seq, SeqRecord
from .utilities import xopen

## deprecated
def findPattern(seq, RE):
    regex = re.compile(RE)
    for m in regex.finditer(seq):
        pos = m.start() + len(RE) // 2
        break
    else:
        pos = None
    return pos

def cutsite_trimming(fastqfile, RE, output):
    """
    HindIII: AAGCTAGCTT
    MboI: GATCGATC
    """
    enzyme_length = len(RE) // 2

    with xopen(fastqfile) as handle:
        out_handle = xopen(output, 'w')
        fastq = SeqIO.parse(handle, "fastq")
        for record in fastq:
           
            pos = record.seq.find(RE) 
            if pos != -1:
                record = record[:pos + enzyme_length]
            SeqIO.write(record, out_handle, "fastq")   

if __name__ == "__main__":
    cutsite_trimming(sys.argv[1], sys.argv[2], sys.argv[3])
    
