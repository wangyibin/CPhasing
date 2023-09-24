#!/usr/bin/env python

"""
cli for methylation align pipelines
"""

import click
import logging
import sys

from ..cli import CommandGroup
from ..cli import cli 

@cli.group(cls=CommandGroup, short_help='Sub-command for the UL ONT pipeline.')
@click.pass_context
def methalign(ctx):
    pass
