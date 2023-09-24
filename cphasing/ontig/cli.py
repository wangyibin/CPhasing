#!/usr/bin/env python

"""
cli for ultra-long ont pipelines
"""

import click
import logging
import sys


from ..cli import CommandGroup
from ..cli import cli 

@cli.group(cls=CommandGroup, short_help='Sub-command for the UL ONT pipeline.')
@click.pass_context
def ontig(ctx):
    pass

