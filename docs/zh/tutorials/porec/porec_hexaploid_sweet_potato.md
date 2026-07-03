# 六倍体甘薯组装与挂载教程（Pore-C 模式）

本教程演示了利用 PacBio HiFi 和 Pore-C 测序数据对六倍体甘薯 (*Ipomoea batatas*, 2n = 6x = 90) 进行基因组组装与挂载（Scaffolding）的工作流程。

---

## 1. 使用 Hifiasm 进行单倍型分辨的 Contig 组装

### 步骤 1.1：生成初始 Unitigs
首先，在默认模式下使用 PacBio HiFi 测序数据运行 `hifiasm` 以生成初始的 unitigs，并将 GFA 图结构转换为 FASTA 格式：

```bash
~/software/hifiasm-0.25.0/hifiasm -t 100 \
    -o hifi.asm \
    hifi.fastq.gz

gfatools gfa2fa hifi.asm.p_utg.gfa > hifi.asm.p_utg.fasta
```

### 步骤 1.2：将 Pore-C 测序数据比对回 Unitigs
将您的 Pore-C 数据比对回前面生成的 unitigs 上：

```bash
cphasing mapper hifi.asm.p_utg.fasta porec_reads.fastq.gz -t 100
```

### 步骤 1.3：将 Pore-C 转换为伪 Hi-C 测序数据
将多互作的 Pore-C 比对结果转换为虚拟的双端“伪 Hi-C”（pseudo-Hi-C）数据，使其与 Hifiasm 的单倍型分相（phasing）算法兼容：

```bash
cphasing-rs porec2 porec_reads.porec.gz -l 0 -o porec2hic
```

### 步骤 1.4：使用 Hifiasm 进行分相组装（Hi-C 模式）
使用 PacBio HiFi 测序数据和新生成的双端伪 Hi-C 数据运行 Hi-C 模式下的 `hifiasm`，以获得高连续性的分相 contigs：

!!! note 
    在与步骤 1.1 相同的目录下运行此命令将重用缓存文件，从而允许 `hifiasm` 直接跳过耗时的 HiFi 纠错（correction）和重叠（overlapping）阶段。

```bash
~/software/hifiasm-0.25.0/hifiasm -t 100 \
    -o hifi.asm \
    hifi.fastq.gz \
    --h1 porec2hic_R1.fa.gz \
    --h2 porec2hic_R2.fa.gz 
```

### 步骤 1.5：合并单倍型组装结果
合并来自两个单倍型（`hap1` 和 `hap2`）的组装序列和组装图，为下游的挂载做准备：

```bash
# 合并 FASTA 文件
cat hifi.asm.hic.hap1.p_ctg.fasta hifi.asm.hic.hap2.p_ctg.fasta > haps.p_ctg.fasta 

# 合并 GFA 文件（不包含序列信息）
cat hifi.asm.hic.hap1.p_ctg.noseq.gfa hifi.asm.hic.hap2.p_ctg.noseq.gfa > haps.p_ctg.noseq.gfa
```

---

## 2. 使用 C-Phasing 进行基因组挂载

使用 `cphasing pipeline` 对合并后的单倍体组装结果进行挂载。

由于甘薯是六倍体，染色体基数为 15（$2n = 6x = 90$），因此最终预期的染色体群组数量为 90。我们将染色体单倍型聚类参数设置为 `-n 15:6`（15 个同源群，每个同源群含 6 个单倍型），相应调整 Pore-C 的解析阈值，并启用 `--collapsed-rescue` 模块来挽救和解析压缩（collapsed）区域：

```bash
cphasing pipeline \
    -f haps.p_ctg.fasta \
    -pcd porec_reads.fastq.gz \
    -t 100 \
    -n 15:6 \
    -e 0 --split-length 2m \
    --collapsed-rescue \
    --gfa haps.p_ctg.noseq.gfa \
    -hcr \
    -p GATC \
    -o cphasing_output
```