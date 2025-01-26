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

__author__ = ("Yibin Wang", "Xingtan Zhang")
__copyright__ = "Copyright (c) 2025, tanger-lab"
__email__ = ("yibinwang96@outlook.com", "zhangxingtan@caas.cn")
__license__ = "BSD"
__status__ = "Development"
__version__ = "0.2.4.r290"
__epilog__ =  f"""
            \b
            Version: {__version__} | \n
            Author: tanger-lab | \n
            Please check out the docs at: [https://wangyibin.github.io/CPhasing](https://wangyibin.github.io/CPhasing)
            """