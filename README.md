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
### For Pore-C pipeline.
- [minimap2](https://github.com/lh3/minimap2)
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


1. **prepare**
    ```bash
    cphasing-rs chromsizes draft.asm.fasta > draft.asm.contigsizes
    ```
2. **mapping**
    ```bash
    minimap2 -x map-ont --secondary=no -c draft.asm.fasta sample.fastq.gz \
         | cphasing-rs paf2table - -o sample.porec.gz
    
    cphasing-rs porec2pairs sample.porec.gz draft.asm.contigsizes -o sample.pairs.gz 
    cphasing prepare pairs2cool sample.pairs.gz draft.asm.contigsizes
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
4. **partition**
    - `extract`
    ```bash
    cphasing extract sample.porec.gz draft.asm.contigsizes sample.hyperedges
    ```
    - `hyperpartition`
    ```bash
    ## for haploid scaffolding 
    cphasing hyperpartition sample.hyperedges draft.asm.contigsizes output.clusters.txt 
    ## for polyploid or diploid phasing must add prune information and use the multi partition mode
    ### auto generate groups
    cphasing hyperpartition sample.hyperedges draft.asm.contigsizes output.clusters.txt --prune prune.contig.list -inc
    ### k aware, 8:4 indicate that this polyploid is a tetraploid with 8 chromosome in each haplotype
    cphasing hyperpartition sample.hyperedges draft.asm.contigsizes output.clusters.txt --prune prune.contig.list -inc -k 8:4
    ```
5. **optimize**
    - `ordering and orientation` **incomplete**
    ```bash
    cphasing optimize group1.count_HindIII.txt sample.10000.cool 
    ```
    - `build`
    ```bash
    cphasing build draft.asm.fasta
    ```
6. **plot**
    - `plot`
    ```bash
    cphasing plot -a groups.agp -m sample.10000.cool -o groups.wg.png
    ```

### Evalutation 
- A diploid that simulated from two rice genome

- Simulation data from haplotype resolved human genome

- Autotetraploid genome
![AP heatmap](pictures/AP/groups.wg.png)
- Hexaploid genome
![SP heatmap](pictures/SP/groups.wg.png)

## ToDo list
- [ ] mapper
- [x] kprune
- [x] hyperpartition
- [ ] optimize
- [x] plot
- [ ] pipeline