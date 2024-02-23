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

logging.getLogger('matplotlib').setLevel(logging.ERROR)
logging.getLogger('cooler').setLevel(logging.ERROR)
logging.getLogger('hicmatrix').setLevel(logging.ERROR)
logging.getLogger('numexpr').setLevel(logging.ERROR)
if not sys.warnoptions:
    warnings.simplefilter("ignore")

warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
warnings.simplefilter(action="ignore", category=DeprecationWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

__author__ = "Yibin Wang"
__copyright__ = "Copyright (c) 2024, Yibin Wang"
__email__ = "yibinwang96@outlook.com"
__license__ = "BSD"
__status__ = "Development"
__version__ = "0.0.63.r161"
