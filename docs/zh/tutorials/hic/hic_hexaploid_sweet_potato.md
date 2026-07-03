# 六倍体甘薯组装与挂载教程

本教程演示了六倍体甘薯 (*Ipomoea batatas*, 2n = 6x = 90) 的基因组组装与挂载（Scaffolding）工作流程。

---

## 1. 使用 Hifiasm 进行 Contig 组装（Hi-C 模式）

首先，利用 PacBio HiFi 测序数据和双端 Hi-C 数据运行 `hifiasm`（Hi-C 模式），以生成单倍型分辨的 unitigs：

```bash
~/software/hifiasm-0.25.0/hifiasm -t 190 \
    -o hifi.asm \
    --h1 fei-tanzania1_S3HiC_R1.fastq.gz,fei-tanzania2_S3HiC_R1.fastq.gz \
    --h2 fei-tanzania1_S3HiC_R2.fastq.gz,fei-tanzania2_S3HiC_R2.fastq.gz \
    SRR29949524_subreads.fastq.gz
```

接下来，合并来自两个单倍型（`hap1` 和 `hap2`）的组装结果和组装图，以便进行后续的挂载：

```bash
# 合并 fasta 文件
cat hifi.asm.hic.hap1.p_ctg.fasta hifi.asm.hic.hap2.p_ctg.fasta > haps.p_ctg.fasta 

# 合并 GFA 文件（不包含序列信息）
cat hifi.asm.hic.hap1.p_ctg.noseq.gfa hifi.asm.hic.hap2.p_ctg.noseq.gfa > haps.p_ctg.noseq.gfa
```

---

## 2. 使用 C-Phasing 进行基因组挂载

运行 `cphasing pipeline` 对合并后的单倍体组装结果进行挂载。由于甘薯是六倍体，染色体基数为 15，因此最终预期的染色体群组（Groups）数量为 90（$15 \times 6$）。

我们在此设置参数 `-n 15:6` 激活单倍型聚类，并启用 `--collapsed-rescue` 模块来挽救和解析压缩（collapsed）区域：

```bash
cphasing pipeline \
    -f haps.p_ctg.fasta \
    --hic1 fei-tanzania1_S3HiC_R1.fastq.gz \
    --hic1 fei-tanzania2_S3HiC_R1.fastq.gz \
    --hic2 fei-tanzania1_S3HiC_R2.fastq.gz \
    --hic2 fei-tanzania2_S3HiC_R2.fastq.gz \
    -t 100 \
    -n 15:6 \
    -e 0 --split-length 0 \
    --collapsed-rescue \
    --gfa haps.p_ctg.noseq.gfa \
    -hcr \
    -p GATC \
    -o cphasing_output
```