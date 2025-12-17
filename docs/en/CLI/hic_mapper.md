# HiC Mapper
## one cell
```shell
cphasing hic mapper -f draft.asm.fasta -1 hic_R1.fastq.gz -2 hic_R2.fastq.gz -t 40
```
!!!note
    If the total length of your input genome is larger than 8 Gb, the `-k 27 -w 14` should be specified, to avoid the error of `chromap`. 



## mutiple cells
#### Submit by specified multiple times
```shell title="run_mapping.sh"
cphasing hic mapper -f draft.asm.fasta -1 hic-1_R1.fastq.gz -1 hic-2_R1.fastq.gz -2 hic-1_R2.fastq.gz -2 hic-2_R2.fastq.gz -t 40
```
!!! note
    The input Hi-C reads must paired by position

#### Submit each cell to the cluster 
- use multiple scripts to run mapper for each sample  
```shell  title="run_sample1.sh"
cphasing hic mapper -f draft.asm.fasta -1 hic-1_R1.fastq.gz -2 hic-1_R2.fastq.gz -t 40
```

    ```shell  title="run_sample2.sh"
    cphasing hic mapper -f draft.asm.fasta -1 hic-2_R1.fastq.gz -2 hic-2_R2.fastq.gz -t 40
    ```

    ```shell  title="run_sample3.sh"
    cphasing hic mapper -f draft.asm.fasta -1 hic-3_R1.fastq.gz -2 hic-3_R2.fastq.gz -t 40
    ```
!!!note
    Please submit one first; wait until the index is successfully created and the remaining jobs are enabled to submit to the cluster.

- merge results
```shell
cphasing pairs-merge hic-*.pairs.pqs -o hic.merge.pairs.pqs
```


