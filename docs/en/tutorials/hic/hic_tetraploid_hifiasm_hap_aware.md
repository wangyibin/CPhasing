# Hi-C tutorial: haplotype-aware scaffolding from hifiasm hap assemblies (tetraploid example)

This tutorial shows how to scaffold **pre-separated haplotypes** produced by **hifiasm** using **C-Phasing** in **haplotype-aware** mode (`--mode hapaware`) with **Hi-C** data.

The goal is to keep haplotypes consistent while producing chromosome-scale scaffolds.

---

## 0. Prerequisites

### Input data
- **hifiasm haplotype assemblies** (FASTA), for example:
  - `hap1.fa`
  - `hap2.fa`
  - `hap3.fa`
  - `hap4.fa`
- **Hi-C read pairs** (FASTQ), for example:
  - `hic_R1.fastq.gz`
  - `hic_R2.fastq.gz`

> Notes
> - This tutorial assumes you already have haplotype-separated assemblies (one FASTA per haplotype).
> - If you have hifiasm `*.hap*.p_ctg.fasta` outputs, convert/select the hap FASTA files you want to scaffold.

---

## 1. Prepare a combined FASTA (recommended)

C-Phasing expects a single assembly FASTA as the contig set to scaffold.  
For hap-aware scaffolding, concatenate the hap assemblies into one FASTA.

```bash
cat hap1.fa hap2.fa hap3.fa hap4.fa > haps.concat.fa
```

### Important: unique contig names
Ensure contig IDs are unique across haplotypes. If hifiasm already prefixes contig names per hap, you are fine. Otherwise, rename contigs to avoid collisions.

A simple pattern is to add a prefix per haplotype, e.g. `h1`, `h2`, `h3`, `h4`.

---

## 2. Provide an initial haplotype clustering file (`-fc`)

`--mode hapaware` is designed for cases where haplotypes are **already separated**.  
Providing an initial cluster file helps C-Phasing keep haplotypes consistent.

```shell
cphasing utils fa2cluster  haps.concat.fa -o haps.clusters.txt 
```

Create a `haps.clusters.txt` with **two columns**:

1. `group_id`
2. `contig_id`

Example (illustrative):

```text
h1   h1tg000001l
h1   h1tg000002l
h2   h2tg000001l
h2   h2tg000002l
h3   h3tg000001l
h4   h4tg000001l
```


> Tip: In hap-aware mode, `-fc first.clusters.txt` is the key input that encodes haplotype priors.

---

## 3. Run C-Phasing (Hi-C + hapaware)

Run the main pipeline with Hi-C read pairs and `--mode hapaware`.

```bash
cphasing pipeline \
  -f haps.concat.fa \
  -hic1 hic_R1.fastq.gz \
  -hic2 hic_R2.fastq.gz \
  --mode hapaware \
  -fc haps.clusters.txt \
  -o cphasing_hapaware_out \
  -t 32 -hcr -p AAGCTT
```

!!! node 

    -  set mapping quality for 2nd partition:
        ```bash
        cphasing pipeline ... -q2 0 
        ```