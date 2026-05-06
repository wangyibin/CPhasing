# Hi-C 教程：基于 hifiasm hap 组装的 haplotype-aware scaffolding（四倍体示例）

本教程展示如何使用 **C-Phasing** 的 **haplotype-aware** 模式（`--mode hapaware`），对 **hifiasm 已经分好的 hap 组装**进行 **Hi-C scaffolding**。  
目标是在构建染色体级 scaffold 的同时，尽量保持单倍型（haplotype）一致性，减少 hap 之间的错误拼接/混合。

---

## 0. 准备工作

### 输入数据
- **hifiasm hap 组装**（FASTA），例如：
  - `hap1.fa`
  - `hap2.fa`
  - `hap3.fa`
  - `hap4.fa`
- **Hi-C 双端 reads**（FASTQ），例如：
  - `hic_R1.fastq.gz`
  - `hic_R2.fastq.gz`

> 说明  
> - 本教程假设你已经有 hifiasm 输出的 **hap 分离结果**（每个 hap 一个 FASTA）。  
> - 若你拿到的是 hifiasm 的 `*.hap*.p_ctg.fasta` 等文件，选择你希望用于 scaffolding 的那几份即可。

---

## 1. 合并 hap FASTA（推荐）

C-Phasing 的 pipeline 以一个 FASTA 作为 contig 集合输入。  
因此建议把多个 hap FASTA 合并成一个：

```bash
cat hap1.fa hap2.fa hap3.fa hap4.fa > haps.concat.fa
```

### 重要：contig 名字必须唯一
确保不同 hap 中 contig ID 不会重名。  
如果 hifiasm 已经在 contig 名里带了 hap 前缀（常见情况），通常没问题；否则需要重命名（例如加前缀 `h1_ / h2_ / h3_ / h4_`），避免后续分组/统计混乱。

---

## 2. 提供初始 hap 分组文件（`-fc`）

`--mode hapaware` 适用于 **hap 已经分离**的场景。  
你需要用 `-fc/--first-cluster` 提供一个“先验分组”（每个 hap 的 contig 属于哪个 group），以指导 hap-aware 的聚类/拼接。

你可以直接用工具从 FASTA 头信息生成：

```bash
cphasing utils fa2cluster haps.concat.fa -o haps.clusters.txt
```

`haps.clusters.txt` 是一个 **两列**文件：

1) `group_id`（组名/组编号；可用字符串）  
2) `contig_id`

示例（仅示意）：

```text
h1   h1tg000001l
h1   h1tg000002l
h2   h2tg000001l
h2   h2tg000002l
h3   h3tg000001l
h4   h4tg000001l
```

> 提示  
> - hap-aware 模式下，`-fc` 是关键输入：它告诉程序“哪些 contig 本来就是同一个 hap”。

---

## 3. 运行 C-Phasing（Hi-C + hapaware）

使用 Hi-C reads 并开启 `--mode hapaware`：

```bash
cphasing pipeline \
  -f haps.concat.fa \
  -hic1 hic_R1.fastq.gz \
  -hic2 hic_R2.fastq.gz \
  --mode hapaware \
  -fc haps.clusters.txt \
  -o cphasing_hapaware_out \
  -t 32 \
  -hcr -p AAGCTT
```

参数说明（简要）：
- `--mode hapaware`：haplotype-aware scaffolding（不依赖 alleletable；也不会在 hyperpartition 阶段要求额外提供 fasta 来生成 alleles）。
- `-fc haps.clusters.txt`：输入 hap 先验分组。
- `-hcr`：仅使用高置信区域（常用于降低错误连接，尤其在重复较高/多倍体时）。
- `-p AAGCTT`：限制性内切酶识别位点（按你的文库酶切类型填写；不确定可不填/用默认）。

---



## 命令

```bash
cat hap1.fa hap2.fa hap3.fa hap4.fa > haps.concat.fa
cphasing utils fa2cluster haps.concat.fa -o haps.clusters.txt

cphasing pipeline \
  -f haps.concat.fa \
  -hic1 hic_R1.fastq.gz -hic2 hic_R2.fastq.gz \
  --mode hapaware \
  -fc haps.clusters.txt \
  -o cphasing_hapaware_out \
  -t 32 \
  -hcr -p AAGCTT
```