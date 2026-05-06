# Prepartition (预分组)

`prepartition` 步骤通过将 contig 比对到单倍体参考基因组来实现第一轮聚类。该过程会在 Hi-C/Pore-C phasing（定相/分型）之前，将 contig 预先划分到不同的同源组中。

## 1. 运行 Prepartition

基于参考单倍体基因组对组装好的 contig 进行聚类。

```shell
cphasing prepartition \
  monoploid.fasta \
  contigs.fasta \
  -t 60 \
  -o prepartition.clusters.txt

```
## 2. 传递给 Pipeline（主流程） 
获得预划分结果后，使用 -fc (first cluster) 参数将其传递给主 cphasing pipeline 流程。


```shell
cphasing pipeline \
  -f contigs.fasta \
  -pcd porec_reads.fastq.gz \
  -fc prepartition.clusters.txt \
  -t 60 
```

- `-fc / --first-cluster`: 预划分步骤中生成的聚类文件。