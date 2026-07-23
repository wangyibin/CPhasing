# Hexaploid Sweet Potato Assembly & Scaffolding (Pore-C Mode)

This tutorial demonstrates the assembly and scaffolding workflow for the hexaploid sweet potato (*Ipomoea batatas*, 2n = 6x = 90) utilizing PacBio HiFi and Pore-C reads.

---

## 1. Haplotype-Resolved Contig Assembly using Hifiasm

### Step 1.1: Generate Initial Unitigs
First, run `hifiasm` in default mode using PacBio HiFi reads to generate the primary unitigs, and convert the GFA graph structure to FASTA format:

```bash
~/software/hifiasm-0.25.0/hifiasm -t 100 \
    -o hifi.asm \
    hifi.fastq.gz

gfatools gfa2fa hifi.asm.p_utg.gfa > hifi.asm.p_utg.fasta
```

### Step 1.2: Map Pore-C Reads to Unitigs
Map your Pore-C reads back to the generated unitigs:

```bash
cphasing mapper hifi.asm.p_utg.fasta porec_reads.fastq.gz -t 100
```

### Step 1.3: Convert Pore-C to Pseudo-Hi-C Reads
Convert the multi-contact Pore-C alignments into virtual paired-end "pseudo-Hi-C" reads to make them compatible with Hifiasm's phasing algorithm:

```bash
cphasing-rs porec2reads porec_reads.porec.gz -l 0 -o porec2hic
```

### Step 1.4: Resolve Assembly with Hifiasm (Hi-C Mode)
Run `hifiasm` in Hi-C mode using the PacBio HiFi reads along with the newly generated pseudo-Hi-C paired-end reads to obtain highly contiguous, phased contigs:

!!! note 
    Executing this command in the same directory as Step 1.1 will reuse the cached files, allowing `hifiasm` to bypass the time-consuming HiFi read correction and overlapping phases.

```bash
~/software/hifiasm-0.25.0/hifiasm -t 100 \
    -o hifi.asm \
    hifi.fastq.gz \
    --h1 porec2hic_R1.fa.gz \
    --h2 porec2hic_R2.fa.gz 
```


### Step 1.5: Merge the Haplotype Assemblies
Merge the assemblies and assembly graphs from both haplotypes (`hap1` and `hap2`) to prepare them for downstream scaffolding:

```bash
# Merge FASTA files
cat hifi.asm.hic.hap1.p_ctg.fasta hifi.asm.hic.hap2.p_ctg.fasta > haps.p_ctg.fasta 

# Merge GFA files (without sequence)
cat hifi.asm.hic.hap1.p_ctg.noseq.gfa hifi.asm.hic.hap2.p_ctg.noseq.gfa > haps.p_ctg.noseq.gfa
```

---

## 2. Scaffolding with C-Phasing

Scaffold the merged haplotype assembly using the `cphasing pipeline`. 

Since sweet potato is a hexaploid organism with a basic chromosome number of 15 ($2n = 6x = 90$), the expected number of final chromosome-level groups is 90. We set the phasing parameters to `-n 15:6` (15 homology groups, each with 6 haplotypes), adjust the Pore-C resolution thresholds, and enable the `--collapsed-rescue` module to resolve collapsed regions:

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