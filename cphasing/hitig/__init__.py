import logging
import warnings
import sys

from rich.console import Console
from rich.logging import RichHandler

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=Console(stderr=True))]
)

__author__ = ("Jiaxin Yu", "Yibin Wang", "Xingtan Zhang")