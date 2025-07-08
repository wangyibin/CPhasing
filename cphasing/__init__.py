import logging
import warnings
import io 
import sys
import click

from rich.logging import Console, RichHandler

console = Console(stderr=True, record=True, file=io.StringIO())
console_html = Console(stderr=True, record=True)
logging.basicConfig(
    # level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=console, enable_link_path=False, rich_tracebacks=True), #]
            RichHandler(console=console_html, rich_tracebacks=True)]
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
__version__ = "0.2.6.r303"
__url__ = "https://github.com/wangyibin/CPhasing"
__doc_url__ = "https://wangyibin.github.io/CPhasing"
__epilog__ =  f"""
            \b
            Version: {__version__} | \n
            Author: tanger-lab | \n
            Please check out the docs at: [https://wangyibin.github.io/CPhasing](https://wangyibin.github.io/CPhasing)
            """