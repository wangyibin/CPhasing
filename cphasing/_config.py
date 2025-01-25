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


## Default parameters for CLI 
ALLELES_TRIM_LENGTH = 25_000

## HyperGraph construction
EDGE_LENGTH = "500k"


## HyperPartition
HCR_LOWER = 0.01
HCR_UPPER = 1.5
MIN_CONTACTS = 5
MIN_CIS_WEIGHT = 5.0
MIN_WEIGHT = 0.1

MIN_SCAFFOLD_LENGTH = 2e6 
MIN_LENGTH = 10_000


