
`cphasing mapper` 是用来处理**Pore-C**或者**HiFi-C**数据的，它输出两个文件：
is designed for process Pore-C data, its output two files:  
    （1）存储高阶互作信息的porec table(`.porec.gz`)  
    （2）存储虚拟成对（VPC）互作的4DN pairs (`.pairs.gz`) 
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
cphasing-rs pairs-merge sample1.porec.pairs.gz sample2.porec.pairs.gz sample3.porec.pairs.gz -o sample.merge.pairs.gz
```


### 处理单个cell的HiFi-C数据
```shell
cphasing mapper draft.contigs.fasta sample.porec.fastq.gz -t 40 --mm2-params "-x map-hifi"
```

### 处理多个cell的HiFi-C数据
与处理Pore-C数据的步骤一致，但需要加`--mm2-params "-x map-hifi"`参数。

## Parameters
```shell title="cphasing mapper -h"
 Usage: cphasing mapper [OPTIONS] REFERENCE FASTQ...   
 
 Mapper for pore-c reads.                                                  
 REFERENCE: Path of reference                                              
                                                                           
 FASTQ: Path of pore-c reads, multiple file enabled, the prefix of output  
 default only use sample 1.                                                
                                                                           
╭─ Arguments ─────────────────────────────────────────────────────────────╮
│ *  REFERENCE    (PATH) [required]                                       │
│ *  FASTQ        (PATH) [required]                                       │
╰─────────────────────────────────────────────────────────────────────────╯
╭─ Options ───────────────────────────────────────────────────────────────╮
│ --enzyme        -e        Restrict site pattern, use comma to separate  │
│                           multiple patterns.                            │
│                           (STR)                                         │
│ --mm2-params              additional parameters for minimap2            │
│                           (STR)                                         │
│                           [default: -x map-ont]                         │
│ --mapq          -q        Minimum quality of mapping [0, 60].           │
│                           (INT)                                         │
│                           [default: 0; 0<=x<=60]                        │
│ --min-identity  -p        Minimum percentage identity of alignments [0, │
│                           1.0].                                         │
│                           (FLOAT)                                       │
│                           [default: 0.8; 0.0<=x<=1.0]                   │
│ --min-length    -l        Minimum length of fragments.                  │
│                           (INT)                                         │
│                           [default: 150]                                │
│ --max-edge      -me       Maximum length of fragment located in the     │
│                           edge of contigs.                              │
│                           (INT)                                         │
│                           [default: 2000]                               │
│ --force         -f        Force run all the command, ignore existing    │
│                           results. The index file also will be removed. │
│ --outprefix     -o        output prefix, if none use the prefix of      │
│                           fastq                                         │
│                           (TEXT)                                        │
│ --threads       -t        Number of threads.                            │
│                           (INT)                                         │
│                           [default: 4]                                  │
│ --help          -h,-help  Show this message and exit.                   │
╰─────────────────────────────────────────────────────────────────────────╯
```
