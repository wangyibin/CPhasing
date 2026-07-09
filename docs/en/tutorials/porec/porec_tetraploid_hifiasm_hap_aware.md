# Pore-C tutorial: haplotype-aware scaffolding from hifiasm hap assemblies (tetraploid example)

This tutorial shows how to construct **pre-separated haplotypes** of a tetraploid genome using **Hifiasm (with Pore-C converted Pseudo-Hi-C data)**, and then scaffold them using **C-Phasing** in **haplotype-aware** mode (`--mode hapaware`).

The goal is to keep haplotypes consistent while producing chromosome-scale scaffolds.

---

## 1. Haplotype-Resolved Contig Assembly (Pore-C to Pseudo-Hi-C)

Before scaffolding, we generate haplotype-resolved contigs by converting Pore-C reads into virtual paired-end "pseudo-Hi-C" reads to guide Hifiasm's phasing.

### Step 1.1: Generate Initial Unitigs
First, run `hifiasm` in default mode using PacBio HiFi reads to generate the primary unitigs, and convert the GFA graph structure to FASTA format:

```bash
~/software/hifiasm-0.25.0/hifiasm -t 100 \
    -o hifi.asm \
    hifi.fastq.gz

gfatools gfa2fa hifi.asm.p_utg.gfa > hifi.asm.p_utg.fasta
```

### Step 1.2: Map Pore-C Reads to Unitigs
Map your Pore-C reads back to the generated unitigs to capture multi-way contact information:

```bash
cphasing mapper hifi.asm.p_utg.fasta porec_reads.fastq.gz -t 100 -o porec_reads.porec.gz
```

### Step 1.3: Convert Pore-C to Pseudo-Hi-C Reads
Convert the multi-contact Pore-C alignments into virtual paired-end "pseudo-Hi-C" reads (`porec2hic_R1.fa.gz` and `porec2hic_R2.fa.gz`):

```bash
cphasing-rs porec2 porec_reads.porec.gz -l 0 -o porec2hic
```

### Step 1.4: Resolve Assembly with Hifiasm (Hi-C Mode)
Run `hifiasm` in Hi-C mode using the PacBio HiFi reads along with the newly generated pseudo-Hi-C paired-end reads to obtain phased contig sets:

!!! note 
    Executing this command in the same directory as Step 1.1 will reuse the cached files, allowing `hifiasm` to bypass the time-consuming HiFi read correction and overlapping phases.

```bash
~/software/hifiasm-0.25.0/hifiasm -t 100 \
    -o hifi.asm \
    hifi.fastq.gz \
    --h1 porec2hic_R1.fa.gz \
    --h2 porec2hic_R2.fa.gz 
```

This will produce phased assemblies for haplotypes (e.g., `hifi.asm.hic.hap1.p_ctg.fasta` to `hifi.asm.hic.hap4.p_ctg.fasta`).

---

## 2. Prepare a combined FASTA

C-Phasing expects a single assembly FASTA as the contig set to scaffold.  
For hap-aware scaffolding, concatenate your output haplotype assemblies into a single FASTA file.

```bash
cat hifi.asm.hic.hap1.p_ctg.fasta \
    hifi.asm.hic.hap2.p_ctg.fasta \
    hifi.asm.hic.hap3.p_ctg.fasta \
    hifi.asm.hic.hap4.p_ctg.fasta > haps.concat.fa
```

### Important: unique contig names
Ensure contig IDs are unique across haplotypes. If hifiasm already prefixes contig names per hap (which it does by default), you are fine. Otherwise, rename contigs to avoid collisions.

---

## 3. Provide an initial haplotype clustering file (`-fc`)

`--mode hapaware` is designed for cases where haplotypes are **already separated**.  
Providing an initial cluster file helps C-Phasing keep haplotypes consistent.

```shell
cphasing utils fa2cluster haps.concat.fa -o haps.clusters.txt 
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

> Tip: In hap-aware mode, `-fc haps.clusters.txt` is the key input that encodes haplotype priors.

---

## 4. Run C-Phasing (Pore-C + hapaware)

Run the main pipeline with original Pore-C reads and `--mode hapaware`. Here we scaffold a tetraploid with a basic chromosome number of 12 ($2n = 4x = 48$, so `-n 4:12`):

```bash
cphasing pipeline \
  -f haps.concat.fa \
  -pcd porec_reads.fastq.gz \
  --mode hapaware \
  -fc haps.clusters.txt \
  -o cphasing_hapaware_out \
  -t 32 -hcr -p AAGCTT -n 4:12
```

!!! note 

    - Set mapping quality for the second partition step if needed:
        ```bash
        cphasing pipeline ... -q2 0 
        ```