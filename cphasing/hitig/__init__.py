import logging
import warnings
import sys

from rich.logging import Console, RichHandler

logging.basicConfig(
    # level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=Console(stderr=True))]
)
