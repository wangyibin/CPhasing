## legacy from ALLHiC
    These functions were developing for ALLHiC2 but deprecated.


    - [gmap](http://research-pub.gene.com/gmap/)


## Hi-C scaffolding based on pregroup by homologous
### Autopolyploid
1. **mapping**
```
cphasing hic mapper -r draft.asm.fasta -1 Lib_R1.fastq.gz -2 Lib_R2.fastq.gz -t 20
```
2. **correct** [Optional]
```
cphasing correct draft.asm.fasta Lib.pairs -t 20 -o corrected.fasta -ob corrected.bed -op corrected.pairs
```

3. **alleles**  
Generate allele table.
```
cphasing alleles -f corrected.fasta -c mono.cds -b mono.bed -k 4 -t 20 --method gene
```
4. **extract**  
Extract contacts for partition.
```
cphasing extract corrected.pairs corrected.fasta -e HindIII
```
5. **pregroup**  
Pregroup by homologous
```
cphasing pregroup Allele.ctg.table corrected.counts_AAGCTT.txt corrected.pairs.txt -f corrected.fasta
```
6. **prune**  
Pruning in each homologous group.
```
cphasing prune Chr01.allele.table Chr01.counts_AAGCTT.txt Chr01.pairs.txt
...
```
7. **partition**  
Partition in each homologous group
```
cphasing partition Chr01.counts_GATC.txt Chr01.pairs.prune.txt 4 --adaptive --threads 10
...
```
`4` is repersent the number of groups (ploidy in there)  

8. **recluster**  
Allelic information was used to further optimize the clustering results.
```
cphasing recluster clusters.txt Chr01.allele.table Chr01.counts_GATC.txt Chr01.pairs.txt 4
...
```
9. **rescue**  
Assign unplaced contigs into partitioned clusters.
```
cphasing rescue clusters.txt corrected.counts_AAGCTT.txt corrected.pairs.txt -o rescued.clusters.txt
```
10. **optimize**  
Ordering and orientation for each group.
```
allhic optimze sample.counts_GATC.group1.txt corrected.clm
...
```
11. **build**  
```
cphasing build corrected.fasta
```
12. **plot**  
```
cooler cload pairs corrected.fasta.chromsizes:10000 corrected.pairs corrected.10k.cool -c1 2 -p1 3 -c2 4 -p2 5
cphasing plot -a groups.agp -m corrected.10k.cool -o groups.wg.png
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