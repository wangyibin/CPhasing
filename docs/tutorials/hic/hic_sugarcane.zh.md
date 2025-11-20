# 现代栽培甘蔗组装(Hi-C)


非整倍 (2n=114), 预估基因组大小 ~ 10 Gb 

## Hi-C数据处理
由于Hi-C数据量大，我们将每个cell的Hi-C数据分别提交到不同的节点上进行比对，然后再合并成一个文件。

!!!note 
    因为`hic mapper`里面有一步索引构建，用户需要先提交一个任务，等到索引创建完成再提交剩余的任务。

```shell
cphasing hic mapper -f sh_hifi.bp.p_utg.fasta -hic1 hic-1_R1.fastq.gz -hic2 hic-1_R2.fastq.gz -t 40 -k 27 -w 14 
```

```shell
cphasing hic mapper -f sh_hifi.bp.p_utg.fasta -hic1 hic-2_R1.fastq.gz -hic2 hic-2_R2.fastq.gz -t 40 -k 27 -w 14 
cphasing hic mapper -f sh_hifi.bp.p_utg.fasta -hic1 hic-3_R1.fastq.gz -hic2 hic-3_R2.fastq.gz -t 40 -k 27 -w 14 
cphasing hic mapper -f sh_hifi.bp.p_utg.fasta -hic1 hic-4_R1.fastq.gz -hic2 hic-4_R2.fastq.gz -t 40 -k 27 -w 14 
```

```shell
cphasing pairs-merge hic-*.pairs.pqs -o hic.merge.pairs.pqs
```


## 组装流程
- 现代栽培甘蔗属于非整倍体，不同同源染色体组内的染色体数量不同，因此我们倾向于先让程序自行分组看看（`-n 0:0`)。
```shell
cphasing pipeline -f sh_hifi.bp.p_utg.fasta -pct hic.mrege.pairs.pqs -t 40 -n 0:0 -hcr -p AAGCTT 
```

- 或者使用`--merge-use-allele`通过等位contig信息，辅助第一次分组
```shell
cphasing pipeline -f sh_hifi.bp.p_utg.fasta -pct porec.mrege.porec.gz -t 40 -n 10:0 -hcr -p AAGCTT --merge-use-allele
```