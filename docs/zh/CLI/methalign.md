# Methalign 

为了提升比对的精度，特别是高度同源区域的互作数量，我们开发了一个模块**`Methalign`**。尝试使用甲基化信息(5mCG)，对模糊的比对进行二次分配。这个模块主要针对超复杂多倍体，它们存在大量的高度同源区域，这些区域，即使使用长度长的互作数据也难以区分。普通多倍体请使用`Mapper`处理Pore-C数据[:octicons-arrow-right-24:Mapper](mapper.md)。

!!!info 
    如果你不需要使用甲基化信息，请使用 [:octicons-arrow-right-24:Mapper](mapper.md)

!!! note
    输入的bam文件应该包含`MM/ML`tags

## 激活 methalign环境
```shell
source activate_cphasing methalign
```
!!! info
    第一次激活，请在有网络的情况下运行该命令

## 计算contig水平基因组的5mCG信息

=== "HiFi reads"
    
    Align the HiFi reads by pbmm2

    ```shell
    pbmm2 index --preset CCS contigs.fasta index.mmi 
    pbmm2 align --preset CCS index.mmi HiFi_reads.bam | samtools view - -b -o HiFi.align.bam
    samtools sort HiFi.align.bam -o HiFi.align.sorted.bam
    samtools index HiFi.align.sorted.bam
    ```

    Calculate the 5mC sites by pb-cpg-tools
    ```shell
    aligned_bam_to_cpg_scores --bam HiFi.align.sorted.bam --ref contigs.fasta --model pileup_calling_model.v1.tflite --modsites-mode reference -q 2 
    awk '{print $1,$2,$3,$9}' OFS='\t' HiFi_align_sorted_bam_to_cpg_scores.combined.bed > output.methyl.bg
    ```

=== "ONT reads" 
    Align the ONT reads by dorado

    ```shell
    dorado aligner contigs.fasta ont_reads.bam -t 20 | samtools view -q 2 -@ 10 -b | samtools sort -@ 10 > ont.align.sorted.bam 
    samtools index ont.align.sorted.bam 
    ```

    Calculate the 5mC sites by modkit

    ```shell
    modkit pileup ont.align.sorted.bam output.methyl.bed --log-filepath pileup.log --cpg --ref contigs.fasta -t 60 --combine-strands
    awk '{print $1,$2,$3,$11}' OFS='\t' output.methyl.bed > output.methyl.bg
    ```

## 将Pore-C或者HiFi-C比对到contig水平基因组上
```shell
dorado aligner contigs.fasta porec_reads.bam --secondary=yes > porec.align.bam
```
!!! note
    Replace `--secondary=yes` to `--mm2-opts "--secondary=yes"` when using dorado >= 0.8.0

## 二次分配比对片段

- 二次分配比对片段
```shell
methalign -t 20 contigs.fasta output.methyl.bg porec.reads*.bam -o porec.align.refined.paf.gz
```
!!! note
    这步会生成两个文件`methalign.refined.paf.gz`和`methalign.refined.porec.gz`。

- 之后使用`cphasing pipeline`完成分型组装

```shell
cphasing pipeline -f contigs.fasta -pct methalign.refined.porec.gz -t 10 -n 8:4 -hcr -p AAGCTT 
```

!!! tips
    - 将bam文件分割成多份以加快处理速度
    ```shell
    cphasing-rs split-bam porec.align.bam -o split_outdir/output.split
    cd split_outdir
    ```
