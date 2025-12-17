# 单倍体或单套基因组的挂载

C-Phasing也支持单倍体的挂载，只需设置`-n n`, 或者设置`--mode haploid`. 

例如一个Diploid (2n=2x=24)：

```shell
cphasing pipeline -f monoploid_ctg.fasta -hic1 hic_R1.fastq.gz -hic2 hic_R2.fastq.gz -t 40 -n 12
```
