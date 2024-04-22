import pandas as pd

PQ_ENGINE = "pyarrow"
PQ_VERSION = "2.6"

CSV_ENGINE = "pyarrow"

## Datatype for PAFTable or PoreCTable
CHROM_COORD_DTYPE = "uint32"
READ_COORD_DTYPE = "uint32"
READ_LENGTH_DTYPE = "uint32"
MAPPING_QUALITY_DTYPE = "uint32"
PERCENTAGE_DTYPE = "float32"

## Datatype of HyperGraph order
HYPERGRAPH_ORDER_DTYPE = "int8"
HYPERGRAPH_COL_DTYPE = "int32"