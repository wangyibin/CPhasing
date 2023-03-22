# **C**-Phasing: **Phasing** and scaffolding polyploid genomes based on Pore-**C** or Hi-**C** data.
```
   ____      ____  _               _             
  / ___|    |  _ \| |__   __ _ ___(_)_ __   __ _ 
 | |   _____| |_) | '_ \ / _` / __| | '_ \ / _` |
 | |__|_____|  __/| | | | (_| \__ \ | | | | (_| |
  \____|    |_|   |_| |_|\__,_|___/_|_| |_|\__, |
                                           |___/ 
```
## Introduction
The major problem of scaffolding polyploid genome by Hi-C is that the lower unique mapping. Now, the long reads based chromosome conformation capture technology, called Pore-C, provide a fine way to solving this problem. So, we developed a new pipeline, called `C-Phasing`, specificially tailored to the polyploid phasing and scaffolding by Pore-C data. Also it can be used to scaffolding by Hi-C data, but will be slowly. 
  
The advantages of `C-Phasing`:   
- High speed.   
- High anchor rate of genome. 
- High accuracy of polyploid phasing. 
## Dependencies
### For core function.
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
### For Hi-C pipeline.
- [chromap](https://github.com/haowenz/chromap)


## Installation
```
git clone https://github.com/wangyibin/CPhasing.git
cd CPhasing
conda env create -f environment.yml

export PATH=/path/to/CPhasing/bin:$PATH
export PYTHONPATH=/path/to/CPhasing:$PYTHONPATH
```

## Pore-C 
### Polyploid
1. **mapping**
    - pore-c-snakemake 
        > Time consuming and memory consuming.   
    - minimap2  
        > This will lose a lot of signals and lead to unsatisfactory results.

2. **prepare**
    ```bash
    cphasing prepare genome draft.asm.fasta HindIII
    cphasing prepare cool sample.pairs draft.asm.chromsizes 
    ```

3. **prunning**
    - `alleles`
    ```bash
    cphasing alleles -f draft.asm.fasta
    ```
    - `kprune`
    > for polyploid or diploid phasing
    ```bash
    cphasing kprune allele.table contigs.whole.cool 
    ```
    - `extract`
    ```bash
    cphasing extract sample.merge.pq draft.asm.contigs sample.hyperedges
    ```
    - `hyperpartition`
    ```bash
    ## for haploid scaffolding 
    cphasing hyperpartition sample.hyperedges draft.asm.contigsizes output.clusters.txt 
    ## for polyploid or diploid phasing must add prune information and use the multi partition mode
    cphasing hyperpartition sample.hyperedges draft.asm.contigsizes output.clusters.txt --prune prune.contig.list --multi  
    ```
    - `ordering and orientation`
    ```
    cphasing optimize group1.count_HindIII.txt sample.10000.cool 
    ```
    - `build`
    ```bash
    cphasing build draft.asm.fasta
    ```
    - `plot`
    ```bash
    cphasing plot -a groups.agp -m sample.10000.cool -o groups.wg.png
    ```
    




