# Pore-C 教程：基于 hifiasm 单倍型组装结果的单倍型感知挂载（以四倍体为例）

本教程介绍如何使用 **Hifiasm（结合 Pore-C 转换的 Pseudo-Hi-C 数据）** 构建四倍体基因组的**预分离单倍型（pre-separated haplotypes）**，然后使用 **C-Phasing** 的**单倍型感知**模式（`--mode hapaware`）对其进行染色体水平的挂载（Scaffolding）。

该工作流的目标是在保持单倍型一致性的同时，产生染色体水平的挂载结果。

---

## 1. 单倍型解析的 Contig 组装（Pore-C 转 Pseudo-Hi-C）

在挂载之前，我们通过将 Pore-C 数据转换为虚拟的、双端（paired-end）的“pseudo-Hi-C”数据来引导 Hifiasm 的分型算法，从而生成单倍型解析的 contigs。

### 步骤 1.1：生成初始 Unitigs
首先，在默认模式下使用 PacBio HiFi 测序数据运行 `hifiasm` 以生成初级 unitigs，并将 GFA 图结构转换为 FASTA 格式：

```bash
~/software/hifiasm-0.25.0/hifiasm -t 100 \
    -o hifi.asm \
    hifi.fastq.gz

gfatools gfa2fa hifi.asm.p_utg.gfa > hifi.asm.p_utg.fasta
```

### 步骤 1.2：将 Pore-C 数据比对回 Unitigs
将 Pore-C 数据比对回生成的 unitigs 以捕获多维相互作用（multi-way contact）信息：

```bash
cphasing mapper hifi.asm.p_utg.fasta porec_reads.fastq.gz -t 100 -o porec_reads.porec.gz
```

### 步骤 1.3：将 Pore-C 转换为 Pseudo-Hi-C 数据
将多维相互作用的 Pore-C 比对结果转换为虚拟的、标准双端 "pseudo-Hi-C" 数据（即 `porec2hic_R1.fa.gz` 和 `porec2hic_R2.fa.gz`）：

```bash
cphasing-rs porec2reads porec_reads.porec.gz -l 0 -o porec2hic
```

### 步骤 1.4：使用 Hifiasm 进行组装解析（Hi-C 模式）
使用 PacBio HiFi 数据以及新生成的 pseudo-Hi-C 双端数据运行 `hifiasm` 的 Hi-C 模式，以获得分型后的 contig 数据集：

!!! note 提示
    在与步骤 1.1 相同的目录下执行此命令将复用旧的缓存文件，从而允许 `hifiasm` 跳过最耗时的 HiFi 读长纠错（correction）和重叠（overlapping）阶段。

```bash
~/software/hifiasm-0.25.0/hifiasm -t 100 \
    -o hifi.asm \
    hifi.fastq.gz \
    --h1 porec2hic_R1.fa.gz \
    --h2 porec2hic_R2.fa.gz 
```

这将会产生不同单倍型的组装结果（例如，`hifi.asm.hic.hap1.p_ctg.fasta` 到 `hifi.asm.hic.hap4.p_ctg.fasta`）。

---

## 2. 准备合并的 FASTA 文件

C-Phasing 需要输入单个组装 FASTA 文件作为待挂载的 contig 集合。  
对于单倍型感知挂载，需要将刚才输出的各个单倍型组装结果合并为一个单独的 FASTA 文件：

```bash
cat hifi.asm.hic.hap1.p_ctg.fasta \
    hifi.asm.hic.hap2.p_ctg.fasta \
    hifi.asm.hic.hap3.p_ctg.fasta \
    hifi.asm.hic.hap4.p_ctg.fasta > haps.concat.fa
```

### 重要提示：确保 Contig 名称唯一
请确保不同单倍型之间的 Contig ID 是唯一的。如果 hifiasm 默认已经为每个单倍型的 contig 添加了前缀（默认即是如此），则直接合并即可。否则，请重命名 contigs 以避免命名冲突。

---

## 3. 提供初始单倍型聚类文件（`-fc`）

`--mode hapaware` 专为**已预先分离**的单倍型设计。  
提供初始聚类文件有助于 C-Phasing 在挂载时维持并保护这些单倍型的一致性。

```shell
cphasing utils fa2cluster haps.concat.fa -o haps.clusters.txt 
```

生成的 `haps.clusters.txt` 应包含两列内容：

1. `group_id` (单倍型组 ID)
2. `contig_id` (Contig ID)

示例（示意）：

```text
h1   h1tg000001l
h1   h1tg000002l
h2   h2tg000001l
h2   h2tg000002l
h3   h3tg000001l
h4   h4tg000001l
```

> 提示：在单倍型感知模式（hapaware）下，`-fc haps.clusters.txt` 是引入单倍型先验信息的关键输入。

---

## 4. 运行 C-Phasing（Pore-C + hapaware）

结合原始的 Pore-C 数据并开启 `--mode hapaware` 运行主流程。这里我们挂载一个基本染色体数为 12 的四倍体（$2n = 4x = 48$，因此 `-n 4:12`）：

```bash
cphasing pipeline \
  -f haps.concat.fa \
  -pcd porec_reads.fastq.gz \
  --mode hapaware \
  -fc haps.clusters.txt \
  -o cphasing_hapaware_out \
  -t 32 -hcr -p AAGCTT -n 4:12
```

!!! note 提示

    - 如果需要，可以为第二步的分群（partition）设置比对质量阈值：
        ```bash
        cphasing pipeline ... -q2 0 
        ```