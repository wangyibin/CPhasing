#!/usr/bin/env python
# -*- coding:utf-8 -*-


import argparse
import logging
import os
import os.path as op
import sys

from pathlib import Path

from ..utilities import run_cmd




logger = logging.getLogger(__name__)

def refine(bam, fasta, bedgraph, 
           penalty, 
           ref_prob_cutoff, 
           prob_cutoff,
           designate_mapq, 
           output, threads=4):
    """
    Refine the alignments by methylation signal.
    """
    cmd = ["cphasing-rs", "methalign",
           bam, "-f", fasta,
           "-b", bedgraph,
           "--ref-penalty", str(penalty),
           "--read-penalty", str(penalty),
           "-r", str(ref_prob_cutoff),
           "-c", str(prob_cutoff),
           "--designate-mapq", str(designate_mapq),
           "-o", output,
           "-t", str(threads)
           ]
    log_dir = Path("logs").mkdir(exist_ok=True)
    logs = op.join("logs", f"{Path(output).stem}.methalign.log")

    flag = run_cmd(cmd, log=logs)
    assert flag == 0, "Failed to refine the alignments."
    logger.info(f"Output refined bam file is saved to {output}")

    return output