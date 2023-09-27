<img src="pictures/logo/C-Phasing_logo3.jpg" alt="C-Phasing logo" width="140px" align="left" />
<h1 align="center"><b>C</b>-Phasing</h1>
<p align="center"> <b>Phasing</b> and scaffolding polyploid genomes based on Pore-<b>C</b>, Ultra-long, or Hi-<b>C</b> data</p>

***  

## Introduction
The main problem with Hi-C scaffolding of polyploid genomes is the lower level of unambiguous mapping. Now, the long reads based chromosome conformation capture technology, called Pore-C, provides an effective way to solve this problem. So we have developed a new pipeline, called `C-Phasing`, which is specifically tailored for polyploid phasing and scaffolding using Pore-C data. It can also be used to scaffold assembly using Hi-C data, but it will be slow.  
  
The advantages of `C-Phasing`:   
- High speed.   
- High anchor rate of genome. 
- High accuracy of polyploid phasing. 

## Installation
### Via anaconda (recommend)
```bash
## Download C-Phasing and install all dependencies
git clone https://github.com/wangyibin/CPhasing.git
cd CPhasing
conda env create -f environment.yml
conda activate cphasing

## Add these command into .bash_profile or .bashrc
export PATH=/path/to/CPhasing/bin:$PATH
export PYTHONPATH=/path/to/CPhasing:$PYTHONPATH
```
### Custom install the dependencies
#### Install C-Phasing
```bash
## Download C-Phasing and install python dependencies
git clone https://github.com/wangyibin/CPhasing.git
cd CPhasing
pip install .

## Add these into .bash_profile or .bashrc
export PATH=/path/to/CPhasing/bin:$PATH
```
#### Dependencies
1. For core function
    - [bedtools](https://bedtools.readthedocs.io/en/latest/)
    - [seqkit](https://bioinf.shenwei.me/seqkit/)
    - [pigz](https://github.com/madler/pigz)
2. For Pore-C pipeline
    - [minimap2](https://github.com/lh3/minimap2)
3. For Hi-C pipeline
    - [chromap](https://github.com/haowenz/chromap)

## Pipeline of Ultra-long data [Optional]
C-Phasing enable to use ultra-long to correct chimeric and use ultra-long data to identify the high confidence regions (HCRs) to help assembly.  
### **[ontig tutorial](cphasing/ontig/README.md)**

## Pipeline of Pore-C data
1. **mapping**  
    > Mapping pore-c data into contig-level assembly and output the pore-c table and 4DN pairs.
    ```bash
    ## results are `sample.porec.gz` and `sample.pairs.gz`
    cphasing mapper draft.asm.fasta sample.fastq.gz -t 10
    ```  
    > Note: At first, only one data can be run until the index is successfully created.  

    > For Hi-C data please use `cphasing hic mapper`

    Note: If you mapping multiple pore-c datas, the multiple `pairs.gz` files should be merged by following steps:
    ```bash
    zgrep "^#" sample-1.pairs.gz > header
    cat header <(zcat sample-1.pairs.gz sample-2.pairs.gz | grep -v "^#") | pigz -p 4 -c > sample.pairs.gz 
    ```
2. **prepare**
    > Prepare some data for subsequence analysis
    - **Step1** `pairs2cool`
    ```bash
    ## results are `sample.10000.cool` and `sample.whole.cool`
    cphasing prepare pairs2cool sample.pairs.gz draft.asm.contigsizes sample.10000.cool
    ```
    - **Step2** `extract`  
    ```bash
    ## results are `sample.counts_AAGCTT.txt` and `sample.clm`
    cphasing prepare extract sample.pairs.gz draft.asm.fasta
    ```
3. **alleles** (Optional for phasing mode)
    > This step is specific to diploid and polyploid phasing. If you only want to scaffolding a haploid, ignore this step.
    - **Step1** `alleles`
    ```bash
    ## result is `draft.asm.allele.table`
    cphasing alleles -f draft.asm.fasta
    ```
    - **Step2** `kprune`
    ```bash
    ## result is `prune.contig.table`
    cphasing kprune draft.asm.allele.table sample.whole.cool -c sample.counts_AAGCTT.txt
    ```
4. **partition**  
    - **Step1** `hypergraph`
    > Construct a hypergraph of pore-c alignments and store as hypergraph format
    ```bash
    ## result is `sample.hg`
    ## for single file
    cphasing hypergraph sample.porec.gz draft.asm.contigsizes sample.hg -t 4
    ## for multiple files
    cphasing hypergraph <(ls *porec.gz) draft.asm.contigsizes sample.hg -t 4 --fofn 
    ```
    Note: The parameter of `--pairs` should be added when you want to import 4DN pairs format.
    - **Step2** `hyperpartition`
    ```bash
    ## result is `output.clusters.txt`
    ## for haploid scaffolding 
    cphasing hyperpartition sample.hg draft.asm.contigsizes output.clusters.txt
    ## for polyploid or diploid phasing must add prune information and use the incremental partition mode
    ### k aware, 8:4 indicate that this is a tetraploid with 8 chromosome in each haplotype
    cphasing hyperpartition sample.hg draft.asm.contigsizes output.clusters.txt -pt prune.contig.table -inc -k 8:4 -t 4
    ### auto generate groups
    cphasing hyperpartition sample.hg draft.asm.contigsizes output.clusters.txt -pt prune.contig.table -inc -t 4
    ```
    Notes: If you are not satisfied with the number of groups, you can adjust the resolution parameters (`-r1` or `-r2`)  to make difference groups. `-r1` controls the first cluster while `-r2` controls the second cluster in phasing mode. Increasing the resolution can generate more groups while decreasing it will get fewer groups.  
    If you want to group some contigs, you can input a contig list with the `-wl` parameter.
5. **scaffolding**
    ```bash
    ## result is `groups.agp`
    cphasing scaffolding output.clusters.txt sample.counts_AAGCTT.txt sample.clm -f draft.asm.fasta -t 4
    ```
6. **plot**
    > Adjust the contig-level contacts matrix into chromosome-level and plot a heatmap.
    ```bash
    cphasing plot -a groups.agp -m sample.10000.cool -o groups.wg.png
    ```


