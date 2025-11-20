
`cphasing mapper` 是用来处理**Pore-C**或者**HiFi-C**数据的，它输出两个文件:    
    （1）存储高阶互作信息的porec table(`.porec.gz`)  
    （2）存储虚拟成对（VPC）互作的4DN pairs (`.pairs.pqs`) 
!!!note
    HiFi-C数据目前也以`.porec.gz`后缀命名，内容是一样的。

## 示例

### 处理单个cell的Pore-C数据
```shell
cphasing mapper draft.contigs.fasta sample.porec.fastq.gz -t 40 
```

### 处理多个cell的Pore-C数据
#### 一个脚本完成比对
```shell
cphasing mapper draft.contigs.fasta sample1.porec.fastq.gz sample2.porec.fastq.gz -t 40 
```
or 
```shell
cphasing mapper draft.contigs.fasta sample*.porec.fastq.gz -t 40 
```
!!! note
    结果使用第一个文件的前缀作为输出前缀


#### 分别提交到集群

- 使用多个脚本提交任务到不同节点
```shell title="run_sample1.sh"
cphasing mapper draft.contigs.fasta sample1.porec.fastq.gz -t 40
```
```shell title="run_sample2.sh"
cphasing mapper draft.contigs.fasta sample2.porec.fastq.gz -t 40
```
```shell title="run_sample3.sh"
cphasing mapper draft.contigs.fasta sample3.porec.fastq.gz -t 40
```

- 合并结果
```shell
cphasing-rs porec-merge sample1.porec.porec.gz sample2.porec.porec.gz sample3.porec.porec.gz -o sample.merge.porec.gz
cphasing-rs pairs-merge sample1.porec.pairs.pqs sample2.porec.pairs.pqs sample3.porec.pairs.pqs -o sample.merge.pairs.pqs
```


### 处理单个cell的HiFi-C数据
```shell
cphasing mapper draft.contigs.fasta sample.porec.fastq.gz -t 40 --mm2-params "-x map-hifi"
```

### 处理多个cell的HiFi-C数据
与处理Pore-C数据的步骤一致，但需要加`--mm2-params "-x map-hifi"`参数。

