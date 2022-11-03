# **C**-Phasing
**Phasing** and scaffolding polyploid genomes based on Pore-**C** or Hi-**C** data.

## Dependencies
### For core function.
- [gmap](http://research-pub.gene.com/gmap/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [asmkit](https://github.com/wangyibin/asmkit)

### For ultra-fast Hi-C pipeline.
- [chromap](https://github.com/haowenz/chromap)

### For fine polyploid Hi-C pipeline (developing function).
- [hisat2](http://daehwankimlab.github.io/hisat2/)
- [samtools](http://www.htslib.org/)

## Installation
```
git clone https://github.com/wangyibin/CPhasing.git
cd CPhasing
pip install -r requirements.txt

export PATH=/path/to/CPhasing/bin:$PATH
export PYTHONPATH=/path/to/CPhasing:$PYTHONPATH
```

## Hi-C data
1. mapping
```
cphasing hic mapper -r draft.asm.fasta -1  Lib_R1.fastq.gz -2 Lib_R2.fastq.gz -t 20
```
2. correct
```
cphasing correct draft.asm.fasta Lib.pairs -t 20 -o corrected.fasta -ob corrected.bed -op corrected.pairs
```
3. Generate allele table
```
cphasing alleles -f corrected.fasta -c mono.cds -b mono.bed -k 4 -t 20
```
4. extract contact for partition
```
cphasing extract corrected.pairs.gz corrected.fasta -e HindIII
```
5. pregroup by homologous
```
cphasing pregroup Allele.ctg.table corrected.counts_AAGCTT.txt corrected.pairs.txt -f corrected.fasta
```
6. pruning in each homologous group
```
cphasing prune Chr01.allele.table Chr01.counts_AAGCTT.txt Chr01.pairs.txt
...
```
7. partition in each homologous group
```
cphasing partition Chr01.counts_GATC.txt Chr01.pairs.prune.txt 4 --adaptive --threads 10
...
```
4 is repersent the number of groups (ploidy in there)  

8. recluster in each homologous group
```

```

## Citiing  
**If you use the hic data to scaffold genomes please also cite ALLHiC.**  
Zhang, X. ,  Zhang, S. ,  Zhao, Q. ,  Ming, R. , &  Tang, H. . (2019). Assembly of allele-aware, chromosomal-scale autopolyploid genomes based on hi-c data. Nature Plants, 5(5). doi: [10.1038/s41477-019-0487-8](https://doi.org/10.1038/s41477-019-0487-8)

```bibtex
@article{2019Assembly,
        title={Assembly of allele-aware, chromosomal-scale autopolyploid genomes based on Hi-C data},
        author={ Zhang, Xingtan  and  Zhang, Shengcheng  and  Zhao, Qian  and  Ming, Ray  and  Tang, Haibao },
        journal={Nature Plants},
        volume={5},
        number={5},
        year={2019},
}
```