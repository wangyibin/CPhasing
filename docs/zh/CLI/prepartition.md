# 同源染色体聚类 (Prepartition)

`prepartition` 命令旨在通过将 contig 比对到近缘物种的参考基因组来实现同源染色体聚类（Pre-cluster）。

当你拥有一个高质量的近缘物种参考基因组（或单倍体参考基因组）时，此步骤可以基于同源共线性为 contig 提供一个可靠的初始分组。随后，可将此该聚类结果传递给 CPhasing 流程。

## 基本用法

通过提供参考基因组和你的草图 contigs 来运行 `prepartition` 命令：

```bash
cphasing prepartition monoploid.fasta contigs.fasta -t 40 -o prepartition.clusters.txt
```

## 与 `cphasing pipeline`结合使用
生成 `prepartition.clusters.tx`t 文件后，你可以通过在主流程中指定 `-fc (first cluster)` 参数将其引入。这允许程序跳过第一轮的同源染色体聚类，直接使用基于参考基因组引导的分组来进行后续的同源染色体分组。


```bash
cphasing pipeline -f contigs.fasta \
                  -hic1 hic_R1.fastq.gz \
                  -hic2 hic_R2.fastq.gz \
                  -t 40 \
                  -n 12:4 \
                  -hcr \
                  -p GATC \
                  -fc prepartition.clusters.txt
```