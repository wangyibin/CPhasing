# 现代栽培甘蔗组装(Pore-C)

非整倍 (2n=114), 预估基因组大小 ~ 10 Gb 

## Pore-C数据处理
由于Pore-C数据量大，我们将每个cell的Pore-C数据分别提交到不同的节点上进行比对，然后再合并成一个文件。
`cphasing porec-merge`
```shell
cphasing mapper sh_hifi.bp.p_utg.fasta porec-1.fastq.gz -t 40 
cphasing mapper sh_hifi.bp.p_utg.fasta porec-2.fastq.gz -t 40 
cphasing mapper sh_hifi.bp.p_utg.fasta porec-3.fastq.gz -t 40 
cphasing mapper sh_hifi.bp.p_utg.fasta porec-4.fastq.gz -t 40 
cphasing mapper sh_hifi.bp.p_utg.fasta porec-5.fastq.gz -t 40 
```

```shell
cphasing porec-merge porec-*.porec.gz -o porec.merge.porec.gz 
```


## 组装流程
- 现代栽培甘蔗属于非整倍体，不同同源染色体组内的染色体数量不同，因此我们倾向于先让程序自行分组看看（`-n 10:0`)。
```shell
cphasing pipeline -f sh_hifi.bp.p_utg.fasta -pct porec.mrege.porec.gz -t 40 -n 10:0 -hcr -p AAGCTT 
```
