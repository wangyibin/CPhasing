# Change logs #
## [v0.2.6] - 2025-03-11
#### Enhancement
- Parameter optimization.
- Supported pairs.pqs processing.
- `pipeline`, if input pairs or pairs.gz, it will first convert it to pairs.pqs to speed up subsequence pairs data load. About a 15 percent increase in speed.
- `plot`, enable plot the border or haplotypes (`--add-hap-border`)
#### Bug fixes
- `pipeline`, fixed bug that when user input uncommon suffix of Hi-C data (e.g. "_R1.f.fastq.gz") the pairs file could not be found. #issue17
- `hypergraph`  
    - fixed bug that `cphasing-rs pairs-intersect` can not filter contacts by mapq  
    - fixed bug that report error when edge_length set to 0
- `pairs2cool`, fixed bug that loss contacts when the chromsizes in `pairs.gz` is unsorted. 
- `rename`
    - fixed bug that when the chromname in agp is named group??, the `rename` cannot be executed
    - generate the unrenamed chromosome, #issue16

## [v0.2.5] - 2025-02-07
#### New features
- `PQS`, custom format of contacts to speed up the load and parse, do not used in this version.
#### Enhancement
- `cphasing-rs`, slightly speed up `pairs2clm`
- `rename`, use g1, g2, g3 ... to rename duplicates chromosomes, when `--unphased` specified.
- `hcr_by_contacts`, speed up it
#### Bug fixes
- `hypergraph`, fixed bug that report error of polars prompt of set "skip_nulls=False"


## [v0.2.4] - 2025-01-25
#### Enhancement
- `pairs2depth`, speed up by pipeline
- `alleles`, add `--trim-length` to trim the both end of contigs to remove the 
                effect of overlapping from hifiasm assembly graph.
- `pipeline`
    - `--preset precision`: Optimize parameters to improve accuracy at the expense of anchor rate
    - `--preset sensitive`: Using in some complex genome, which contain many fragmented contigs and low signals contigs

#### Bug fixes
- `activate_cphasing`, fixed bug that exit shell window when pixi failed to install
- `pipeline`
    - fixed bug that program can not exit when error occurred 
    - fixed bug that alleles parameters cannot affect in phasing mode
    - fixed bug that reported in issue #14, which `--chimeric-correct` mode cannot use correct fasta in 3.hyperpartition
- `plot`, fixed bug that program can not coarsen adjusted matrix to plot another binsize matrix

## [v0.2.3] - 2025-01-14
#### Enhancement
- `pipeline`, report peak memory usage
- `hyperpartition`
    - add `--min-cis-weight=5.0` to remove low contacts contigs
    - increase `--min-contacts` to `25.0` to remove low contacts contigs

#### Bug fixes
- `pairs2cool`, fixed bug that "bin1_id > bin2_id" in some cases

## [v0.2.2] - 2025-01-10
#### Enhancement
- `pipeline`, When mode=phasing, pipeline integrating 1.alleles into 3.hyperpartition to speed up
#### Bug fixes
- `prepare`, fixed bug of "got unexpected keyword argument `has_header`'


## [v0.2.1] - 2025-01-07
#### Enhancement
- `hypergraph`, add hcr_bed filter step.
- `hyperpartition`, speed up the `phasing` mode
- `activate_cphasing`, use pixi to activate the environment of CPhasing
- `pairs2cool`, speed up it at the expense of memory consumption
- `plot`, speed up it
#### Bug fixes    
- `alleles`, a cheat method to fix bug of partig can not parse contig > 130 Mb
- `rename`, bug of read agp
- `pipeline`
    - input hic with `_1` or `_2` suffix can not load successful
    - input hic can not find fastq path in `output dir`


## [v0.2.0] - 2024-10-1
#### Enhancement
- `pipeline`
    - change the normalization method of the hic pipeline
    - change pipeline output directories
- `hyperpartition`, improve the performance of cluster merge
#### Bug fixes
- `plot`, fixed bug of split contigs

## [v0.1.9] - 2024-07-09
#### New features
- `collapsed-rescue`, init a collapsed contigs rescue function
#### Enhancement
- `higig`
    - support for hifi data. Moreover, support for junk and collapsed identification
    - `correct-alignments`, increase performance by filtering low quality LIS
- `pairs2cool`, add min_mapq filtering
- `plot`, add triangle plot 
- `scaffolding`, speed up split clm by `cphasing-rs splitclm`
- `alleles`, speed up
#### Bug fixes
- `pipeline` 
    - fixed bug of contacts generating
    - fixed bug of the `corrected.agp` not output when chimeric correct mode
- `scaffolding`, fixed bug of countre extract not found contigs 
- `hcr_from_contacts`, find error peak in low coverage regions

## [v0.1.8] - 2024-05-15
#### Enhancement
- `hypergraph`, `Extractor`, add mapq to hypergraph
- `remove_misassembly`, add remove_misassembly of error between homologous chromosome
#### Bug fixes
- `plot`, fixed bug of resolution not in filename

## [v0.1.7] - 2024-05-11
#### Enhancement
- `hypergraph`, speed up data loading
- `scaffolding`, speed up clm loadding 

## [v0.1.6] - 2024-05-07
#### Enhancement
- `hcr_by_contacts`, remove whole collapsed contigs 
- `stat_porec_table`, reduce memory usage
#### Bug fixes
- `scaffolding`, length db not load contig that RE count < 3
- `simulate_collapse`, fixed bug of contig position in collapsed contigs

## [v0.1.5] - 2024-04-29
#### Enhancement 
- environment, add pigz 
- `hic mapper`, add "remove pcr duplicates" as default parameters

## [v0.1.4] - 2024-04-23
#### Enhancement
- colorful the help text
- `hcr`, use kde peaks 
- `pipeline`, use 


## [v0.1.3] - 2024-04-22
#### New features
- `_chromap`, add a modified `chromap`
#### Enhancement
- `hypergraph`, reduce memory usage


## [v0.1.2] - 2024-04-18
#### Enhancement 
- `scaffolding`, optimized
- `cphasing-rs pairs2clm`, speed up

## [v0.1.1] - 2024-04-11
#### New features
- `hitig`
    - `scaffolding`, which enable to scaffolding by ultra-long reads.

#### Enhancement
- `plot`, speed up
- `alleles`, remove filter function
- `hypergraph`, speed up by porec table load.


## [v0.1.0] - 2024-04-07
#### New features
- `alleles2`, a method for allelic contig identification by self mapping
- `cphasing-rs prune`, a method to prune by raw allele table (ALLHiC allele table)

#### Enhancement
- `scaffolding`, increase the performance of long contigs

#### Bug fixes
- `plot`, fixed bug that it can not be used in pandas v2.0


#### [v0.0.64] - 2024-03-11
#### Bug fixes
- `hyperpartition`, fixed the inconsistent of phasing results, when specified alleletable and prunetable
- `alleles`, filterring by valid kmer length of contigs

#### [v0.0.62] - 2024-02-18
#### New features
- `cphasing-rs` `kprune`, increase the recall of cross-allelic

#### [v0.0.61] - 2024-01-19
#### New features
- `plot_high_order_distribution.py`, plot the contact order distribution of pore-c alignments

#### Enhancement
- `scaffolding`, adopt the HapHiC_sort to scaffolding, and add HaplotypeAlign
- `partig`, change `d` to `0.1`
- `AlleleTable`, add `strand` column to the format of `allele2`
- `Tour`, add `to_dict`, `backup`, `save`

## [v0.0.60]
#### Enhancement
- `hyperpartition`, add prune based on hypergraph
- `hyperpartition`, add different n setting
- `pipeline`, remove single kprune, change it to hypergraph prune


## [v0.0.59]
#### Enhancement
- `hyperpartition`, add automatic search resolution
- `kprune`, update the formulate of normalization

## [v0.0.58]

#### Enhancement
- `mapper`, add restriction site filter and realign
- `hyperpartition`, new mode by phasing
- `hyperpartition`, support import pore-c table or pairs
- `HyperEdges`, add contigsizes, mapq

## [v0.0.57]
#### Enhancement
- `pipeline`, add `hcr`
- `hcr`, optimized

## [v0.0.56]
#### Enhancement 
- `cphasing-rs`, update `pairs-intersect` and `porec-intersect`
- `kprune`, remove `count_re` argument

## [v0.0.55]
#### New features
- `hcr_from_contacts`, add new function of hcr
- `plot_lines`, plot lines of evaluation

## [v0.0.54]
#### Bug fixes
- `pipeline`, fixed bugs


## [v0.0.53]
#### New features
- `pipeline`, add a pipeline of C-Phasing
#### Enhancement
- `prepare`, restruct 
- `pairs2cool`, move to second command
## [v0.0.52]
#### Enhancement
- `kprune`, load data into tempfile to decrease the memory usage.

##[v0.0.51]
- 2023-11-11
#### New features
- `agp`, add `pseudo-agp`, create a pseudo agp from simulation contigs

#### Bug fixed
- `hyperpartition`, add `cross_allelic_factor`


##[v0.0.50]
- 2023-11-10
#### Enhancement
- `kprune`, add multiprocessing

##[v0.0.49]
- 2023-11-10
#### Enhancement 
- `cphasing-rs`, update it into v0.0.10

#### New features
- `porec2csv`, import porec table into pao csv
- `PoreCTable.binnify`, binnify the contig by binsize
- `PoreCTable.divide_contig_by_nparts`, divide contig by the number of parts
#### Bug fixed
- `hyperpartition`, fixed bug that phasing mode report error when the prune table not apply

## [v0.0.48]
- 2023-11-05
#### Enhancement
- `kprune`, advanced it to identity more cross-allelic 
#### New features
- `evaluate_prune.py`, evaluate the result of `kprune`
- `methalign/pipe`, init the pipeline of methalign


## [v0.0.47]
- 2023-10-15
#### Enhancement
- `higig`, mv `ontig` to `hitig`.

## [v0.0.46]
#### Enhancement
- `hcr`, add `break_pos` to correct hcr by break positions
#### Bug fixed
- `mapper`, fixed min_quality can not be adopt

## [v0.0.45]
#### Enhancement
- `hcr`, intergated `bed2depth` into `hcr`
#### Bug fixed
- `bed2depth`, fixed bug
## [v0.0.44]
#### Enhancement
- `ontig`, restruct the framework and add `split-reads`, `find-chimeric`, `hcr` three functions.

## [v0.0.43] - 2023-09-24
#### New features
- `methalign`, init
#### Enhancement
- rename `ultra_long` to `ontig`


## [v0.0.42] - 2023-09-22
#### Enhancement
- `ultra_long`, expose ultra_long tutorial to C-Phasing README.md

## [v0.0.41] - 2023-08-25
#### Bug fixed
- `scaffolding`, fixed cp groups.agp error when change the name of agp
## [v0.0.40]
- 2023-08-19
#### Enhancement
- add `pyproject.toml`
#### Bug fixed
- `setup.py`, fixed the setup_helper error.

## [v0.0.39]
#### Enhancement
- before this version `cphasing-rs` can not found
## [v0.0.38]
- 2023-08-18
#### Bug fixed
- `build`, fixed bug that can not specify the output of agp
- `scaffolding`, fixed bug that can not specify the output of agp
## [v0.0.37]
- 2023-08-16
#### Enhancement
- `README.md`, add detail describtions
- `scaffolding`, directly output `groups.agp`
#### Bug fixed
- fixed bug that resulted some misassembly from v0.0.17

## [v0.0.36]
- 2023-08-14
#### Bug fixed
- `hyperpartition`  
    - fixed bug from v0.0.35 that default run ultra-complex mode
    - fixed bug that whitelist can not used in single partition
- `recluster`
    - fixed bug of `KeyError`

## [v0.0.35]
- 2023-08-12
#### Enhancement
- `pairs2cool`, add `--fofn` parameter
- `hyperpartition`, mask `ultra-complex` parameter
## [v0.0.34]
- 2023-08-11
#### Enhancement
- `mapper`, compress hic mapper result
- `hyperpartition`, add contig filter when import first cluster result
- 2023-08-04
#### Bug fixed
- `hypergraph`, fixed the column index error
## [v0.0.32]
- 2023-08-04
#### Bug fixed
- `hypergraph`, fixed it can not parse single porec table 
## [v0.0.31]
- 2023-08-02
#### New features
- `ultra_long`, add v1 pipeline
## [v0.0.30]
- 2023-08-02
#### Enhancement
- `hypergraph`, speed up
#### Bug fixed
- `porec-intersection`, fixed it can not output result
## [v0.0.29]
- 2023-07-31
#### Bug fixed
- `hyperpartition`, fixed key word of `edges` not yet change to `hypergraph`
## [v0.0.28]
- 2023-07-29
#### New features
- `collapse`, create a collapse rescued contact matrix
#### Enhancement
- update `environment.yml`
- update `cphasing-rs`, which reduce the size of binary and update to `v0.0.4`
#### Bug fixed
- `scaffolding`, fixed path error
## [v0.0.27]
- 2023-07-14
#### Enhancement
- remove `check_allhic_version`
## [v0.0.26]
- 2023-07-12
#### New features
- `statcluster`, stat the clustertable
- `plot_hist`, plot the distribution of somethings
#### Enhancement
- `hypergraph`, add prune allelic hyperedges
- `statagp`, change `agpstat` to `statagp`
- change `extract` to `hypergraph`
- change `optimize` to `scaffolding`
- `hyperpartition`, out "Chr??g?" as chromosome name

#### Bug fixed
- `hyprpartition`, not add `-pt` or `at`, the `-inc` not work

## [v0.0.25]
- 2023-07-05
#### New features
- `evalutate_known_assembly.py`, evalutate the simulated human genome
- `cis_trans_by_contigs.py`, calculate the cis and trans contacts by pairs
- `logo`, add logo
#### Enhancement
- `hyperpartition`, automatic using the exists first clusters results

## [v0.0.24]
- 2023-06-20
#### New features
- `ultra_long`, init ultra_long method
#### Enhancement 
- `mapper`, add `outprefix` prameter
## [v0.0.23]
- 2023-06-17
#### New features
- `pictures`, add the heatmap of M1
#### Enhancement
- `hyperpartition`, expose `--min-allelic-overlap`
- `kprune`, add sort and normalize function
## [v0.0.22]
- 2023-06-16
#### New features
- `tmpoptimze`, a temp function for order and orientation
- `PoreCMapper`, based on minimap2 and cphasing-rs
#### Enhancement
- `plot`, add `--no-ticks` to remove ticks from picture

## [v0.0.21]
- 2023-06-12
#### Enhancement
- `hyperpartition`
    - change `zero-allelic` to `allelic-factor`
    - add `merge_cluster`
    - init `ultra_complex`
###### Bug fixed
- `plot`, enable `factor` can set to 1

## [v0.0.20]
- 2023-06-07
#### Enhancement 
- `AlleleTable`, load info
- `PartigRecord`, output allele table with additional information
#### New features
- `AlleleInfo`, object function of allele table information
#### Bug fixed
- `HyperPartition`, fixed bug of prunt -> prunetable

## [v0.0.19]
- 2023-06-06
#### Enhancement 
- `KPruneHyperGraph`, speed up
- `HyperPartition`, add the get_prunetable_KPrune

## [v0.0.18]
- 2023-06-06
#### Enhancement
- `hyperpartition`
    - `min_weight`, add a parameter of the graph minimum weight 
- `IRMM`, add a function of the graph filter by minimum weight
#### New features
- `KPruneHyperGraph`, add a function that implement pruning on hypergraph

## [v0.0.17]
- 2023-06-01
#### Enhancement
- `environment.yml`, update
- `hyperpartition`
    - add `--first-cluster` to load exists first cluster results
    - expose `--zero-allelic` parameter
    - expose `--allelic-similarity` parameter
- `PoreCTable`, change to cphasing-rs porec table
#### New features
- `docs`, add docs
- `utils`
    - `extract_matrix`, extract matrix by th contig list
    - `prune_matrix`, prune matrix by prune table
#### Bug fixes
- `merge-matrix`, fixed the empty contact in result

## [v0.0.16]
- 2023-05-05
#### Enhancement
- `hyperpartition`, add mutual exclusion merge
#### Bug fixes
- `hyperpartition`, fixed the bug that incorrect scaffold length in`_incremental_partition` 
## [v0.0.15]
- 2023-05-01
#### New features
- `PruneTable`, Object for prune table.
#### Enhancement
- `kprune`, output prune table
- `hyperpartition`
    - add k parameter to output the specified number groups
## [v0.0.14]
- 2023-03-22
#### New features
- `plot_heatmap`, custom function to plot heatmap 
#### Enhancement
- `HyperPartition`
    - add post check to increase the accuracy of partition
    - add min-contacts parameters to remove low contact contig pairs
    - add negative allelic algorithm
    - change multi to incremental
    - add `whitelist` and `blacklist` parameters
#### Bug fixes
- `KPrune`, dropna in pixels which will let igraph fall into loops

## [v0.0.13]
- 2023-03-15
#### New features
- `KPruner`, prune contig contacts according similarity allele table
#### Enhancement
- `extract_incidence_matrix`, restruct to speed up 
- `extract`, `hyperpartition`, import contig sizes instead of fasta
- `HyperPartition`, add filter function of minimum scaffolding length  
- `build`, add only-agp option
##### Bug fixes
- `HyperExtractor`, fixed the problem of read index for multi dataframe concat
## [v0.0.12]
- 2023-03-11
#### New features
- `HyperGraph`, instead of hypernetx to speed up 
- `prepare`, a cli function for some preparation for C-Phasing
- `optimize.cpp`, optimize algorithm v1
#### Enhancement
- `HyperEdges`, restruct to {idx, row, col}
- `Extractor` and `HyperExtractor`, change to export `HyperEdges` 
- `process_pore_c_table`, speed up
- `HyperPartition`, use custom `HyperGraph` to speed up
- `merge_matrix`, merge binning matrix into whole contigs
- `optimize`, cli for optimize score cacluation
#### Bug fixes
- `alleles`, fix the bug of empty allele table

## [v0.0.11]
 - 2023-03-08
#### Enhancement
- `Extractor` and `HyperExtractor`, change to export dict to list
- `HyperPartition`, add parallel in hypergraph generation 
## [v0.0.10]
 - 2023-03-08
#### New features
- `Extractor`, extract edges from pairs file
- `HyperExtractor`, extract hyperedges from pore-c table
- `HyperEdges`, msgspec Struct for hyperedges serialization

#### Enhancement
- `HyperPartition`, divesting edges extract function to improve usability

## [v0.0.9]
- 2023-01-06
`refactor the frame of cphasing that seperating the hic and pore-c pipelines`
#### New features
- `agp2assembly`, convert agp to assembly file
- `pairs2mnd`, convert pairs to mnd file
- `paf2table`, convert paf to table
- `pore_c_cchrom2contig`, convert chromosome-level to contig-level
- `pairs_chrom2contig`, convert chromosome-level to contig-level
- `optimzie`, init add optimize
#### Enhancement
- `PAFTable`, speed up
- `PoreCTable`, speed up and add `chrom2contig`
- `HyperPartition`, speed up the pore-c reads process
- `paf2pairs`, speed up
#### Bug fixes
- `plot`, to humanized


## [v0.0.8]
- 2022-11-17
#### New features
- `hypergraph`, contain a hypergraph cluster algorithms
- `HyperPartition`, partition contigs in diploid or allopolyploid
- `PartigAllele`, enable using sequences similarity to generate allele table
- `PAFTable`, filter alignment
- `PoreCTable`, stat alignment information
#### Enhancement
- `agp2fasta`, `build`, using dictory instead of faidx to speed up
#### Bug fixes
- `PairHeader`, save header without "\n".
- `rescue`, fix CountRE only default load >=3.
- `plot`, fix bug of hicexplorer outputFilename error.

## [v0.0.7]
- 2022-10-31
#### New features
- `plot`, provide adjust matrix and plot matrix function.
- `AdapativePartition`, partition by adapative pipeline for find best results
- `PoreCMapper`, Pore-C reads mapping pipeline.
- `PAFTable`, a function for processing paf file. (include `to_pairs`)

##### Enhancement
- `correct`, increase recall and precision to 70%
- `pairs`, add `chrom2contig` to convert a chrom-level pairs to contig-level
## [v0.0.6]
- 2022-10-07

#### New features
- Change packages name to CPhasing.
- New function of 
    - `correct`
    - `alleles`
    - `rescue`
    - `partition`
- New API of `Pairs` in `core`.

##### Bug fixes
- Using normalized signal in `prune`.
