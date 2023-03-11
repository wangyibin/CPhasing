# Release notes #

## [v0.0.12]
Date: 2023-03-11
### New features
- `HyperGraph`, instead of hypernetx to speed up 
- `prepare`, a cli function for some preparation for C-Phasing
- `optimize.cpp`, optimize algorithm v1
### Enhancement
- `HyperEdges`, restruct to {idx, row, col}
- `Extractor` and `HyperExtractor`, change to export `HyperEdges` 
- `process_pore_c_table`, speed up
- `HyperPartition`, use custom `HyperGraph` to speed up
- `merge_matrix`, merge binning matrix into whole contigs
- `optimize`, cli for optimize score cacluation
### Bug fixes
- `alleles`, fix the bug of empty allele table

## [v0.0.11]
Date: 2023-03-08
### Enhancement
- `Extractor` and `HyperExtractor`, change to export dict to list
- `HyperPartition`, add parallel in hypergraph generation 
## [v0.0.10]
Date: 2023-03-08
### New features
- `Extractor`, extract edges from pairs file
- `HyperExtractor`, extract hyperedges from pore-c table
- `HyperEdges`, msgspec Struct for hyperedges serialization

### Enhancement
- `HyperPartition`, divesting edges extract function to improve usability

## [v0.0.9]
Date: 2023-01-06
`refactor the frame of cphasing that seperating the hic and pore-c pipelines`
### New features
- `agp2assembly`, convert agp to assembly file
- `pairs2mnd`, convert pairs to mnd file
- `paf2table`, convert paf to table
- `pore_c_cchrom2contig`, convert chromosome-level to contig-level
- `pairs_chrom2contig`, convert chromosome-level to contig-level
- `optimzie`, init add optimize
### Enhancement
- `PAFTable`, speed up
- `PoreCTable`, speed up and add `chrom2contig`
- `HyperPartition`, speed up the pore-c reads process
- `paf2pairs`, speed up
### Bug fixes
- `plot`, to humanized


## [v0.0.8]
Date: 2022-11-17
### New features
- `hypergraph`, contain a hypergraph cluster algorithms
- `HyperPartition`, partition contigs in diploid or allopolyploid
- `PartigAllele`, enable using sequences similarity to generate allele table
- `PAFTable`, filter alignment
- `PoreCTable`, stat alignment information
### Enhancement
- `agp2fasta`, `build`, using dictory instead of faidx to speed up
### Bug fixes
- `PairHeader`, save header without "\n".
- `rescue`, fix CountRE only default load >=3.
- `plot`, fix bug of hicexplorer outputFilename error.

## [v0.0.7]
Date: 2022-10-31
### New features
- `plot`, provide adjust matrix and plot matrix function.
- `AdapativePartition`, partition by adapative pipeline for find best results
- `PoreCMapper`, Pore-C reads mapping pipeline.
- `PAFTable`, a function for processing paf file. (include `to_pairs`)

### Enhancement
- `correct`, increase recall and precision to 70%
- `pairs`, add `chrom2contig` to convert a chrom-level pairs to contig-level
## [v0.0.6]
Date: 2022-10-07

### New features
- Change packages name to CPhasing.
- New function of 
    - `correct`
    - `alleles`
    - `rescue`
    - `partition`
- New API of `Pairs` in `core`.

### Bug fixes
- Using normalized signal in `prune`.
