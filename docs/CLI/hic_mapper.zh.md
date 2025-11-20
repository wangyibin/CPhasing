# HiC Mapper
处理Hi-C数据
## one cell
```shell
cphasing hic mapper -f draft.asm.fasta -1 hic_R1.fastq.gz -2 hic_R2.fastq.gz -t 40
```
!!!note
    如果基因组大小大于8Gb，应指定更高的kmer和window大小，例如`-k 27 -w 14`，避免`chromap`报错。
   


## mutiple cells

#### 将不同cell的Hi-C 数据提交到集群
- 将多个脚本提交到不同的节点
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
    因为有一步创建索引的步骤，因此用户需要先提交一个任务，直到索引创建完成再提交剩下的任务，或者在不同的目录下分别提交。

- 合并结果
将多个`.pairs.gz`文件合并成一个。
```shell
cphasing pairs-merge hic-*.pairs.pqs -o hic.merge.pairs.pqs
```
