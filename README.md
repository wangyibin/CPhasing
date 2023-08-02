<img src="pictures/logo/C-Phasing_logo3.jpg" alt="C-Phasing logo" width="140px" align="left" />
<h1 align="center"><b>C</b>-Phasing</h1>
<p align="center"> <b>Phasing</b> and scaffolding polyploid genomes based on Pore-<b>C</b> or Hi-<b>C</b> data</p>

***  

## Introduction
The main problem with Hi-C scaffolding of polyploid genomes is the lower unambiguous mapping. Now, the long reads based chromosome conformation capture technology, called Pore-C, provides a fine way to solve this problem. So we have developed a new pipeline, called `C-Phasing`, which is specifically tailored to polyploid phasing and scaffolding by Pore-C data. It can also be used to scaffold assembly by Hi-C data, but will be slow.  
  
The advantages of `C-Phasing`:   
- High speed.   
- High anchor rate of genome. 
- High accuracy of polyploid phasing. 

## Installation
```
git clone https://github.com/wangyibin/CPhasing.git
cd CPhasing
conda env create -f environment.yml

export PATH=/path/to/CPhasing/bin:$PATH
export PYTHONPATH=/path/to/CPhasing:$PYTHONPATH
```
### Dependencies
> All dependencies already add into `environment.yml`.
#### For core function
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
#### For Pore-C pipeline
- [minimap2](https://github.com/lh3/minimap2)
- [pigz](https://github.com/madler/pigz)
#### For Hi-C pipeline
- [chromap](https://github.com/haowenz/chromap)


## Pipeline of Pore-C data
1. **mapping** 
    ```bash
    cphasing mapper draft.asm.fasta sample.fastq.gz
    ```
    > For Hi-C data please use `cphasing hic mapper` 
2. **pairs2cool**
    - `pairs2cool`
    > Convert 4DN pairs to cool and calculate contacts between the whole contigs. 
    > The `sample.whole.cool` will be generated at the same time.
    ```
    cphasing prepare pairs2cool sample.pairs.gz draft.asm.contigsizes sample.10000.cool
    ```
3. **alleles** (Optional for phasing mode)
   
    - `alleles`
    ```bash
    cphasing alleles -f draft.asm.fasta
    ```
    - `kprune`
    > for polyploid or diploid phasing
    ```bash
    cphasing kprune allele.table contigs.whole.cool -c sample.counts_HindIII.txt
    ```
4. **partition**
    - `hypergraph`
    > Construct a hypergraph of pore-c alignments and store as hypergraph format
    ```bash
    cphasing hypergraph sample.porec.gz draft.asm.contigsizes sample.hg -t 4
    ```
    - `hyperpartition`
    ```bash
    ## for haploid scaffolding 
    cphasing hyperpartition sample.hg draft.asm.contigsizes output.clusters.txt
    ## for polyploid or diploid phasing must add prune information and use the incremental partition mode
    ### auto generate groups
    cphasing hyperpartition sample.hg draft.asm.contigsizes output.clusters.txt --prune prune.contig.list -inc -t 4
    ### k aware, 8:4 indicate that this polyploid is a tetraploid with 8 chromosome in each haplotype
    cphasing hyperpartition sample.hg draft.asm.contigsizes output.clusters.txt --prune prune.contig.list -inc -k 8:4 -t 4
    ```
5. **scaffolding**
    ```bash
    cphasing scaffolding output.clusters.txt sample.counts_AAGCTT.txt sample.clm -t 4
    ```
6. **build**
    ```bash
    cphasing build draft.asm.fasta
    ```
7. **plot**
    ```bash
    cphasing plot -a groups.agp -m sample.10000.cool -o groups.wg.png
    ```

### Evalutation 
- A diploid that simulated from two rice genome

- Simulation data from haplotype resolved human genome

- Autotetraploid genome  
    AP (left) and alfalfa (right)
    <p float="center">
        <img src="pictures/AP/groups.wg.png" width="250" />
        <img src="pictures/M1/groups.wg.png" width="250" />
    </p>
    

- Hexaploid genome  
    sweet potato
    <p float="center">
        <img src="pictures/SP/groups.wg.png" width="400" />
    </p>
- Ultra-complex polyploid

## ToDo list
- [x] mapper
- [x] kprune
- [x] hyperpartition
- [x] optimize
- [x] plot
- [ ] pipeline