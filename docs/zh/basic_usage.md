
# 流程化的命令
## pipeline命令
以下案例中的参数`-n 8:4`表示组装一个染色体基数为**8**的四倍体（**4**），指定程序第一次分组输出8组，第二次再每一组里面输出4组，最终输出32条染色体。我们也支持自动分组，只需设置`-n 0:0`, 这表示两轮自动分组。也可以设置为第一分成8组第二次自动分组（`-n 8:0`）或者第一次自动分组第二次分4组（`-n 0:8`)
!!! note
    CPhasing也支持单套基因组的挂载，只需设置一个分组数，例如`-n 8`，染色体数为8，期望程序将contig分成8条染色体。`pipeline`会自动跳过`1.allele`这步，并且只执行一次contig分组.

!!! note
    如果用户输入的基因组为异源多倍体，亚基因组的相似性较低，可能无法实现理想的第一次分组，用户应该设置不同的`-n`参数，例如异源四倍体(2n=4x=32)，组装两套亚基因组的话，应该把它当成无需分型的二倍体组装，设置`-n 16`, 需要对两套亚基因组分别分型的话可以设置`-n 16:2`。  
    如果为异源六倍体(AAABBB, 2n=6x=48), 可以设置`-n 16:3`，另一种类型的异源六倍体（AABBCC，2n=6x=48），则设置`-n 24:2`。

!!! tip
    如果是hifiasm组装的contig，组装分型的基因组选择输入`p_utg.fasta`, 或者用户可以尝试合并hifiasm (hic)输出的多个hap*.p_utg.fasta。


### 输入 **Pore-C 测序reads**:
```bash
cphasing pipeline -f draft.asm.fasta -pcd sample.fastq.gz -t 10 -n 8:4
```

### 输入多组 **Pore-C data**: 
每个数据指定一次`-pcd`参数
```bash
cphasing pipeline -f draft.asm.fasta -pcd sample1.fastq.gz -pcd sample2.fastq.gz -t 10 -n 8:4
```  
    
!!! note
    如果你希望，将不同cell的数据分别投递到不同的节点上跑，可以使用`cphasing mapper`参数，并最后使用`cphasing-rs porec-merge`将输出的`porec.gz` 文件合并为一个`.merge.porec.gz`, 然后再通过`-pct`参数输入。详细文档[:octicons-arrow-right-24:Mapper](CLI/mapper.md)
   
### 输入 **Pore-C table (`.porec.gz`)**:
这种文件是由`cphasing mapper`产生的.
```bash
cphasing pipeline -f draft.asm.fasta -pct sample.porec.gz -t 10 -n 8:4
```


### 输入**HiFi-C**数据
仅需比对的时候更改为`pipeline`或者`mapper`的`--mm2-params "-x map-hifi" `参数. 剩下步骤和生成的文件都与Pore-C数据一致。
```bash
cphasing pipeline -f draft.asm.fasta -pcd hific.fastq.gz --mm2-params "-x map-hifi" -t 10 -n 8:4
```
!!! note 
    HiFi-C的比对结果以`.porec.gz`结尾，处理方式与Pore-C一致，均可用`porec-merge`、`porec-intersect`等命令处理。

### 输入**Hi-C data**数据
```bash
cphasing pipeline -f draft.asm.fasta -hic1 Lib_R1.fastq.gz -hic2 Lib_R2.fastq.gz -t 10 -n 8:4
```
!!! note
    - **1** | `cphasing pipeline` 不支持输入多组Hi-C数据，但是用户可以自行使用`cphasing hic mapper` 进行比对，然后再使用`cphasing-rs pairs-merge` 合并所有的`.pairs.gz`文件, 然后通过 `-prs` 参数输入.  
    - **2** | 如果基因组的大小大于8Gb，需要提高比对的`k` 和`w`的大小，避免`chromap`报错，例如 `-hic-mapper-k 27 -hic-mapper-w 14`.


### 输入4DN pairs (`.pairs.gz` or `.pairs.pqs`)文件

=== "`pairs.pqs`"
    ```bash
    cphasing pipeline -f draft.asm.fasta -prs sample.pairs.pqs -t 10 -n 8:4
    ```
=== "`pairs.gz`"
    ```bash
    cphasing pipeline -f draft.asm.fasta -prs sample.pairs.gz -t 10 -n 8:4
    ```

### 跳过一些步骤
```bash
## skip steps 1.alleles and 2.prepare steps 
cphasing pipeline -f draft.asm.fasta -pct sample.porec.gz -t 10 -ss 1,2
```
### 只跑某一步
```bash
## run 3.hyperpartition 
cphasing pipeline -f draft.asm.fasta -pct sample.porec.gz -t 10 -s 3
```
!!!tips
    `-s 3,4` to run step 3 and 4.

### 提升组装质量
大部分情况下，设置`-hcr`参数会使得组装质量更高，虽然这会让运行速度变慢，但是该参数的设置可以有效去除一些在全基因组频繁互作的区域对分型的影响。同时由于酶切位点分布的偏好性，建议指定`-p AAGCTT` 进行标准化。
```bash
cphasing pipeline -f draft.asm.fasta -pct sample.porec.gz -t 10 -hcr -p AAGCTT 
```  

## 通过Juicebox 人工调整基因组组装
- 首先需要生成两个文件：`.assembly` 和 `.hic`, 以下步骤依赖于[3d-dna](https://github.com/aidenlab/3d-dna) 软件，用户需自行安装。

```shell
cphasing pairs2mnd sample.pairs.gz -o sample.mnd.txt
cphasing utils agp2assembly groups.agp > groups.assembly
bash ~/software/3d-dna/visualize/run-assembly-visualizer.sh sample.assembly sample.mnd.txt
```
!!! note
    如果contig进行嵌合纠错，请使用`groups.corrected.agp`，并从原始的`pairs.gz`通过`cphasing-rs pairs-break`  生成一个新的 `corrected.pairs.gz`。
    

- 调整完后生成agp和fasta文件
```shell
## convert assembly to agp
cphasing utils assembly2agp groups.review.assembly -n 8:4 
## or haploid or a homologous group
cphasing utils assembly2agp groups.review.assembly -n 8
## extract contigs from agp 
cphasing agp2fasta groups.review.agp draft.asm.fasta --contigs > contigs.fasta
## extract chromosome-level fasta from agp
cphasing agp2fasta groups.review.agp draft.asm.fasta > groups.review.asm.fasta
```


## 根据参考基因组重新命名和调整整条染色体的方向
我们内置了`cphasing rename`命令，该命令可以根据一套单倍型或者近缘物种的染色体水平基因组去调整多倍体分型组装的染色体编号和整条染色体的方向方向，可以查看详细的文档[:octicons-arrow-right-24:Rename](CLI/rename.md)
```bash
cphasing rename -r mono.fasta -f draft.asm.fasta -a groups.review.agp -t 20
```
!!! note 
    为了减少时间消耗，我们只对每套同源染色体组内的第一条进行比对，然后根据此比对结果，对整个同源组进行重命名和调整方向。这是因为我们在`scaffolding`步骤的时候已经把同源染色体的方向调整为一致了。如果你需要对所有同源组内的染色体进行比对，可以设置`--unphased`参数实现。

## 绘制热图
查看详细文档[:octicons-arrow-right-24:Plot](CLI/plot.md)