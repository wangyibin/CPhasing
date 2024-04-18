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
conda env create -f environment.py11.yml
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
    - [wfmash](https://github.com/waveygang/wfmash)(>=0.13.0)
    - [bedtools](https://bedtools.readthedocs.io/en/latest/)
    - [seqkit](https://bioinf.shenwei.me/seqkit/)
    - [pigz](https://github.com/madler/pigz)
2. For Pore-C pipeline
    - [minimap2](https://github.com/lh3/minimap2)(>= v2.24)
3. For Hi-C pipeline
    - [chromap](https://github.com/haowenz/chromap)

## Pipeline of Ultra-long data [Optional]
C-Phasing enable to use ultra-long to correct chimeric and identify the high confidence regions (HCRs) to help assembly.  
### **[hitig tutorial](cphasing/hitig)**

## One command pipeline of C-Phasing
- Start from a pore-c data  
Current only support single fastq file, if you have multi cells, you can use `cphasing mapper` and `cphasing-rs porec-merge` to generate the `porec.gz`.  

```bash
cphasing pipeline -f draft.asm.fasta -pcd sample.fastq.gz -t 10 -s all -n "8:4"
```
- Start from a pore-c table  
```bash
cphasing pipeline -f draft.asm.fasta -pct sample.porec.gz -t 10 -s all 
```  
- Start from a Hi-C 4DN pairs file  
```bash
cphasing pipeline -f draft.asm.fasta -prs sample.pairs.gz -t 10
```  
- Skip some steps  
```bash
## skip 1.alleles and 2.prepare steps 
cphasing pipeline -f draft.asm.fasta -pct sample.porec.gz -t 10 -ss 1,2
```

- Only run specified steps  
```bash
## run 3.hyperpartition 
cphasing pipeline -f draft.asm.fasta -pct sample.porec.gz -t 10 -s 3
```

## Step by step pipeline of Pore-C data  
0. **mapping**  
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
0.1. **hcr** (Optional)

    > Only retain high confidence regions (HCRs) to subsequencing analysis
    ```bash
    ## results are `sample.hcr.bed`, `sample.hcr.porec.gz` and `sample.hcr.pairs.gz`
    cphasing -pct sample.porec.gz -cs drfat.asm.contigsizes 
    ```
1. **alleles** (Optional for phasing mode)  
    > This step is specific to diploid and polyploid phasing. If you only want to scaffolding a haploid, ignore this step.
    - **Step1** `alleles`
    ```bash
    ## result is `draft.asm.allele.table`
    cphasing alleles -f draft.asm.fasta
    ```

2. **prepare**  
    > Prepare some data for subsequence analysis  
    ```bash  
    ## results are `sample.counts_AAGCTT.txt`, `sample.clm`, `sample.split.contacts`, `sample.contacts`
    cphasing prepare draft.asm.fasta sample.pairs.gz 
    ```  

3. **hyperpartition**  
    ```bash  
    ## result is `output.clusters.txt`
    ## for haploid scaffolding 
    cphasing hyperpartition sample.porec.gz draft.asm.contigsizes output.clusters.txt
    ## for polyploid or diploid phasing must add prune information and use the incremental partition mode
    ### chromsome number aware, 8:4 indicate that this is a tetraploid with 8 chromosome in each haplotype
    cphasing hyperpartition sample.porec.gz draft.asm.contigsizes output.clusters.txt -pt prune.contig.table -inc -n 8:4 -t 4
    ### auto generate groups
    cphasing hyperpartition sample.porec.gz draft.asm.contigsizes output.clusters.txt -pt prune.contig.table -inc -t 4
    ```  

    Notes: If you are not satisfied with the number of groups, you can adjust the resolution parameters (`-r1` or `-r2`)  to make difference groups. `-r1` controls the first cluster while `-r2` controls the second cluster in phasing mode. Increasing the resolution can generate more groups while decreasing it will get fewer groups.  
    If you want to group some contigs, you can input a contig list with the `-wl` parameter.

4. **scaffolding**  
    ```bash  
    ## result is `groups.agp`
    cphasing scaffolding output.clusters.txt sample.counts_AAGCTT.txt sample.clm -sc sample.split.contacts -f draft.asm.fasta -t 4

    ## for polyploid can specified allele table to adjust the orientation of different haplotypes
    cphasing scaffolding output.clusters.txt sample.counts_AAGCTT.txt sample.clm -at draft.asm.allele.table -sc sample.split.contacts -f draft.asm.fasta 
    ```  
5. **plot**  
    > Adjust the contig-level contacts matrix into chromosome-level and plot a heatmap.  

    ```bash  
    ## result is `sample.10000.cool` and `groups.wg.png`
    cphasing pairs2cool sample.pairs.gz draft.asm.contigsizes sample.10000.cool
    cphasing plot -a groups.agp -m sample.10000.cool -o groups.wg.png
    ```  


### Curation by Juicebox
```bash
cphasing alignments pairs2mnd sample.pairs.gz sample.mnd.txt
cphasing utils agp2assembly groups.agp > groups.assembly
bash ~/software/3d-dna/visualize/run-assembly-visualizer.sh sample.assembly sample.mnd.txt
```


