# Release notes #
## [v0.1.8]
Date: 2024-05-15
## Enhancement
- `hypergraph`, `Extractor`, add mapq to hypergraph
- `remove_misassembly`, add remove_misassembly of error between homologous chromosome
## Bug fixes
- `plot`, fixed bug of resolution not in filename

## [v0.1.7]
Date: 2024-05-11
## Enhancement
- `hypergraph`, speed up data loading
- `scaffolding`, speed up clm loadding 

## [v0.1.6]
Date: 2024-05-07
## Enhancement
- `hcr_by_contacts`, remove whole collapsed contigs 
- `stat_porec_table`, reduce memory usage
## Bug fixes
- `scaffolding`, length db not load contig that RE count < 3
- `simulate_collapse`, fixed bug of contig position in collapsed contigs

## [v0.1.5]
Date: 2024-04-29
## Enhancement 
- environment, add pigz 
- `hic mapper`, add "remove pcr duplicates" as default parameters

## [v0.1.4]
Date: 2024-04-23
## Enhancement
- colorful the help text
- `hcr`, use kde peaks 
- `pipeline`, use 


## [v0.1.3]
Date: 2024-04-22
## New features
- `_chromap`, add a modified `chromap`
## Enhancement
- `hypergraph`, reduce memory usage


## [v0.1.2]
Date: 2024-04-18
## Enhancement 
- `scaffolding`, optimized
- `cphasing-rs pairs2clm`, speed up

## [v0.1.1]
Date: 2024-04-11
## New features
- `hitig`
    - `scaffolding`, which enable to scaffolding by ultra-long reads.

## Enhancement
- `plot`, speed up
- `alleles`, remove filter function
- `hypergraph`, speed up by porec table load.


## [v0.1.0]
Date: 2024-04-07
## New features
- `alleles2`, a method for allelic contig identification by self mapping
- `cphasing-rs prune`, a method to prune by raw allele table (ALLHiC allele table)

## Enhancement
- `scaffolding`, increase the performance of long contigs

## Bug fixes
- `plot`, fixed bug that it can not be used in pandas v2.0

## [v0.0.64]
Date: 2024-03-11
## Bug fixes
- `hyperpartition`, fixed the inconsistent of phasing results, when specified alleletable and prunetable
- `alleles`, filterring by valid kmer length of contigs

## [v0.0.62]
Date: 2024-02-18
## New features
- `cphasing-rs` `kprune`, increase the recall of cross-allelic

## [v0.0.61]
Date: 2024-01-19
## New features
- `plot_high_order_distribution.py`, plot the contact order distribution of pore-c alignments

## Enhancement
- `scaffolding`, adopt the HapHiC_sort to scaffolding, and add HaplotypeAlign
- `partig`, change `d` to `0.1`
- `AlleleTable`, add `strand` column to the format of `allele2`
- `Tour`, add `to_dict`, `backup`, `save`

## [v0.0.60]
## Enhancement
- `hyperpartition`, add prune based on hypergraph
- `hyperpartition`, add different n setting
- `pipeline`, remove single kprune, change it to hypergraph prune


## [v0.0.59]
## Enhancement
- `hyperpartition`, add automatic search resolution
- `kprune`, update the formulate of normalization

## [v0.0.58]

## Enhancement
- `mapper`, add restriction site filter and realign
- `hyperpartition`, new mode by phasing
- `hyperpartition`, support import pore-c table or pairs
- `HyperEdges`, add contigsizes, mapq

## [v0.0.57]
## Enhancement
- `pipeline`, add `hcr`
- `hcr`, optimized

## [v0.0.56]
## Enhancement 
- `cphasing-rs`, update `pairs-intersect` and `porec-intersect`
- `kprune`, remove `count_re` argument

## [v0.0.55]
## New features
- `hcr_from_contacts`, add new function of hcr
- `plot_lines`, plot lines of evaluation

## [v0.0.54]
## Bug fixes
- `pipeline`, fixed bugs


## [v0.0.53]
## New features
- `pipeline`, add a pipeline of C-Phasing
## Enhancement
- `prepare`, restruct 
- `pairs2cool`, move to second command
## [v0.0.52]
## Enhancement
- `kprune`, load data into tempfile to decrease the memory usage.

##[v0.0.51]
Date: 2023-11-11
## New features
- `agp`, add `pseudo-agp`, create a pseudo agp from simulation contigs

## Bug fixed
- `hyperpartition`, add `cross_allelic_factor`


##[v0.0.50]
Date: 2023-11-10
## Enhancement
- `kprune`, add multiprocessing

##[v0.0.49]
Date: 2023-11-10
## Enhancement 
- `cphasing-rs`, update it into v0.0.10

## New features
- `porec2csv`, import porec table into pao csv
- `PoreCTable.binnify`, binnify the contig by binsize
- `PoreCTable.divide_contig_by_nparts`, divide contig by the number of parts
## Bug fixed 
- `hyperpartition`, fixed bug that phasing mode report error when the prune table not apply

## [v0.0.48]
Date: 2023-11-05
## Enhancement
- `kprune`, advanced it to identity more cross-allelic 
## New features
- `evaluate_prune.py`, evaluate the result of `kprune`
- `methalign/pipe`, init the pipeline of methalign


## [v0.0.47]
Date: 2023-10-15
## Enhancement
- `higig`, mv `ontig` to `hitig`.

## [v0.0.46]
## Enhancement
- `hcr`, add `break_pos` to correct hcr by break positions
## Bug fixed
- `mapper`, fixed min_quality can not be adopt

## [v0.0.45]
## Enhancement
- `hcr`, intergated `bed2depth` into `hcr`
## Bug fixed
- `bed2depth`, fixed bug
## [v0.0.44]
## Enhancement
- `ontig`, restruct the framework and add `split-reads`, `find-chimeric`, `hcr` three functions.

## [v0.0.43]
Date: 2023-09-24
## New features
- `methalign`, init
## Enhancement
- rename `ultra_long` to `ontig`


## [v0.0.42]
Date: 2023-09-22
## Enhancement
- `ultra_long`, expose ultra_long tutorial to C-Phasing README.md

## [v0.0.41]
Date: 2023-08-25
## Bug fixed
- `scaffolding`, fixed cp groups.agp error when change the name of agp
## [v0.0.40]
Date: 2023-08-19
## Enhancement
- add `pyproject.toml`
## Bug fixed
- `setup.py`, fixed the setup_helper error.

## [v0.0.39]
## Enhancement
- before this version `cphasing-rs` can not found
## [v0.0.38]
Date: 2023-08-18
## Bug fixed
- `build`, fixed bug that can not specify the output of agp
- `scaffolding`, fixed bug that can not specify the output of agp
## [v0.0.37]
Date: 2023-08-16
## Enhancement
- `README.md`, add detail describtions
- `scaffolding`, directly output `groups.agp`
## Bug fixed
- fixed bug that resulted some misassembly from v0.0.17

## [v0.0.36]
Date: 2023-08-14
## Bug fixed
- `hyperpartition`  
    - fixed bug from v0.0.35 that default run ultra-complex mode
    - fixed bug that whitelist can not used in single partition
- `recluster`
    - fixed bug of `KeyError`

## [v0.0.35]
Date: 2023-08-12
## Enhancement
- `pairs2cool`, add `--fofn` parameter
- `hyperpartition`, mask `ultra-complex` parameter
## [v0.0.34]
Date: 2023-08-11
## Enhancement
- `mapper`, compress hic mapper result
- `hyperpartition`, add contig filter when import first cluster result
Date: 2023-08-04
## Bug fixed
- `hypergraph`, fixed the column index error
## [v0.0.32]
Date: 2023-08-04
## Bug fixed
- `hypergraph`, fixed it can not parse single porec table 
## [v0.0.31]
Date: 2023-08-02
## New features
- `ultra_long`, add v1 pipeline
## [v0.0.30]
Date: 2023-08-02
## Enhancement
- `hypergraph`, speed up
## Bug fixed
- `porec-intersection`, fixed it can not output result
## [v0.0.29]
Date: 2023-07-31
## Bug fixed
- `hyperpartition`, fixed key word of `edges` not yet change to `hypergraph`
## [v0.0.28]
Date: 2023-07-29
## New features
- `collapse`, create a collapse rescued contact matrix
## Enhancement
- update `environment.yml`
- update `cphasing-rs`, which reduce the size of binary and update to `v0.0.4`
## Bug fixed
- `scaffolding`, fixed path error
## [v0.0.27]
Date: 2023-07-14
## Enhancement
- remove `check_allhic_version`
## [v0.0.26]
Date: 2023-07-12
## New features
- `statcluster`, stat the clustertable
- `plot_hist`, plot the distribution of somethings
## Enhancement
- `hypergraph`, add prune allelic hyperedges
- `statagp`, change `agpstat` to `statagp`
- change `extract` to `hypergraph`
- change `optimize` to `scaffolding`
- `hyperpartition`, out "Chr??g?" as chromosome name

## Bug fixed
- `hyprpartition`, not add `-pt` or `at`, the `-inc` not work

## [v0.0.25]
Date: 2023-07-05
## New features
- `evalutate_known_assembly.py`, evalutate the simulated human genome
- `cis_trans_by_contigs.py`, calculate the cis and trans contacts by pairs
- `logo`, add logo
## Enhancement
- `hyperpartition`, automatic using the exists first clusters results

## [v0.0.24]
Date: 2023-06-20
## New features
- `ultra_long`, init ultra_long method
## Enhancement 
- `mapper`, add `outprefix` prameter
## [v0.0.23]
Date: 2023-06-17
## New features
- `pictures`, add the heatmap of M1
## Enhancement
- `hyperpartition`, expose `--min-allelic-overlap`
- `kprune`, add sort and normalize function
## [v0.0.22]
Date: 2023-06-16
## New features
- `tmpoptimze`, a temp function for order and orientation
- `PoreCMapper`, based on minimap2 and cphasing-rs
## Enhancement
- `plot`, add `--no-ticks` to remove ticks from picture

## [v0.0.21]
Date: 2023-06-12
## Enhancement
- `hyperpartition`
    - change `zero-allelic` to `allelic-factor`
    - add `merge_cluster`
    - init `ultra_complex`
## Bug fixed 
- `plot`, enable `factor` can set to 1

## [v0.0.20]
Date: 2023-06-07
## Enhancement 
- `AlleleTable`, load info
- `PartigRecord`, output allele table with additional information
## New features
- `AlleleInfo`, object function of allele table information
## Bug fixed
- `HyperPartition`, fixed bug of prunt -> prunetable

## [v0.0.19]
Date: 2023-06-06
## Enhancement 
- `KPruneHyperGraph`, speed up
- `HyperPartition`, add the get_prunetable_KPrune

## [v0.0.18]
Data: 2023-06-06
## Enhancement
- `hyperpartition`
    - `min_weight`, add a parameter of the graph minimum weight 
- `IRMM`, add a function of the graph filter by minimum weight
## New features
- `KPruneHyperGraph`, add a function that implement pruning on hypergraph

## [v0.0.17]
Date: 2023-06-01
## Enhancement
- `environment.yml`, update
- `hyperpartition`
    - add `--first-cluster` to load exists first cluster results
    - expose `--zero-allelic` parameter
    - expose `--allelic-similarity` parameter
- `PoreCTable`, change to cphasing-rs porec table
## New features
- `docs`, add docs
- `utils`
    - `extract_matrix`, extract matrix by th contig list
    - `prune_matrix`, prune matrix by prune table
## Bug fixes
- `merge-matrix`, fixed the empty contact in result

## [v0.0.16]
Date: 2023-05-05
## Enhancement
- `hyperpartition`, add mutual exclusion merge
## Bug fixes
- `hyperpartition`, fixed the bug that incorrect scaffold length in`_incremental_partition` 
## [v0.0.15]
Date: 2023-05-01
### New features
- `PruneTable`, Object for prune table.
### Enhancement
- `kprune`, output prune table
- `hyperpartition`
    - add k parameter to output the specified number groups
## [v0.0.14]
Date: 2023-03-22
### New features
- `plot_heatmap`, custom function to plot heatmap 
### Enhancement
- `HyperPartition`
    - add post check to increase the accuracy of partition
    - add min-contacts parameters to remove low contact contig pairs
    - add negative allelic algorithm
    - change multi to incremental
    - add `whitelist` and `blacklist` parameters
## Bug fixes
- `KPrune`, dropna in pixels which will let igraph fall into loops

## [v0.0.13]
Date: 2023-03-15
### New features
- `KPruner`, prune contig contacts according similarity allele table
### Enhancement
- `extract_incidence_matrix`, restruct to speed up 
- `extract`, `hyperpartition`, import contig sizes instead of fasta
- `HyperPartition`, add filter function of minimum scaffolding length  
- `build`, add only-agp option
### Bug fixes
- `HyperExtractor`, fixed the problem of read index for multi dataframe concat
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
