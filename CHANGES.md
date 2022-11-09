# Release notes #
## [v0.0.8]
Date: 2022-11-09
### New features
- `hypergraph`, contain a hypergraph cluster algorithms
- `HyperPartition`, partition contigs in diploid or allopolyploid
### Bug fixes
- `PairHeader`, save header without "\n".


## [v0.0.7]
Date: 2022-10-31
### New features
- `plot`, provide adjust matrix and plot matrix function.
- `AdapativePartition`, partition by adapative pipeline for find best results
- `PoreCMapper`, Pore-C reads mapping pipeline.
- `PAFTable`, a function for processing paf file. (include `to_pairs`)
### Bug fixes

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
