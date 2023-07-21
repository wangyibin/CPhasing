<img src="pictures/logo/C-Phasing_logo2.jpg" alt="C-Phasing logo" width="160px" align="right" />
<h1 align="center">C-Phasing</h1>
<p align="center">Phasing and scaffolding polyploid genomes based on Pore-C or Hi-C data</p>

***  

**C**-Phasing: **Phasing** and scaffolding polyploid genomes based on Pore-**C** or Hi-**C** data.
## Introduction
The major problem of scaffolding polyploid genome by Hi-C is that the lower unique mapping. Now, the long reads based chromosome conformation capture technology, called Pore-C, provide a fine way to solve this problem. So, we developed a new pipeline, called `C-Phasing`, specificially tailored to the polyploid phasing and scaffolding by Pore-C data. Also it can be used to scaffold by Hi-C data, but will be slowly. 
  
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
#### For core function.
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
#### For Pore-C pipeline.
- [minimap2](https://github.com/lh3/minimap2)
- [pigz](https://github.com/madler/pigz)
#### For Hi-C pipeline.
- [chromap](https://github.com/haowenz/chromap)


## Pipeline of Pore-C data
### Polyploid
1. **mapping**
    ```bash
    cphasing mapper draft.asm.fasta sample.fastq.gz
    ```
2. **alleles**
    - `alleles`
    ```bash
    cphasing alleles -f draft.asm.fasta
    ```
    - `kprune`
    > for polyploid or diploid phasing
    ```bash
    cphasing kprune allele.table contigs.whole.cool -c sample.counts_HindIII.txt
    ```
3. **partition**
    - `hypergraph`
    ```bash
    cphasing hypergraph sample.porec.gz draft.asm.contigsizes sample.hyperedges
    ```
    - `hyperpartition`
    ```bash
    ## for haploid scaffolding 
    cphasing hyperpartition sample.hyperedges draft.asm.contigsizes output.clusters.txt 
    ## for polyploid or diploid phasing must add prune information and use the incremental partition mode
    ### auto generate groups
    cphasing hyperpartition sample.hyperedges draft.asm.contigsizes output.clusters.txt --prune prune.contig.list -inc
    ### k aware, 8:4 indicate that this polyploid is a tetraploid with 8 chromosome in each haplotype
    cphasing hyperpartition sample.hyperedges draft.asm.contigsizes output.clusters.txt --prune prune.contig.list -inc -k 8:4
    ```
4. **scaffolding**
    ```bash
    cphasing scaffolding output.clusters.txt sample.counts_AAGCTT.txt sample.clm -t 10
    ```
5. **build**
    ```bash
    cphasing build draft.asm.fasta
    ```
6. **plot**
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