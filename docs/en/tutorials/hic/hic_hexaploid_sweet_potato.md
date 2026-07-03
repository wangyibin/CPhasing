
# Hexaploid Sweet Potato Assembly & Scaffolding
This tutorial demonstrates the assembly and scaffolding workflow for the hexaploid sweet potato (*Ipomoea batatas*, 2n = 6x = 90).

## 1. Contig Assembly using Hifiasm (Hi-C Mode)

First, run `hifiasm` in Hi-C mode using PacBio HiFi reads and paired-end Hi-C data to generate haplotype-resolved unitigs:

```bash
~/software/hifiasm-0.25.0/hifiasm -t 190 \
    -o hifi.asm \
    --h1 fei-tanzania1_S3HiC_R1.fastq.gz,fei-tanzania2_S3HiC_R1.fastq.gz \
    --h2 fei-tanzania1_S3HiC_R2.fastq.gz,fei-tanzania2_S3HiC_R2.fastq.gz \
    SRR29949524_subreads.fastq.gz
```
Next, merge the assemblies and assembly graphs from both haplotypes (`hap1` and `hap2`) for downstream scaffolding:

```bash
# Merge fasta files
cat hifi.asm.hic.hap1.p_ctg.fasta hifi.asm.hic.hap2.p_ctg.fasta > haps.p_ctg.fasta 

# Merge GFA files (without sequence)
cat hifi.asm.hic.hap1.p_ctg.noseq.gfa hifi.asm.hic.hap2.p_ctg.noseq.gfa > haps.p_ctg.noseq.gfa
```

## 2. Scaffolding with C-Phasing

Run the `cphasing pipeline` to scaffold the merged haploid assemblies. Since sweet potato is hexaploid with 15 chromosomes per monoploid set, the expected number of final groups is 90 ($15 \times 6$). We configure `-n 15:6` and enable the `--collapsed-rescue` module to resolve collapsed regions:

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
