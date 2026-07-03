# Prepartition

The `prepartition` command is designed to pre-cluster contigs by mapping them to a closely related reference genome. 

When a high-quality reference genome of a closely related species (or a monoploid reference) is available, this step provides a reliable initial grouping of contigs based on synteny. This pre-partitioned result can then be passed to the main CPhasing pipeline to guide and accelerate the clustering process.

## Basic Usage

Run `prepartition` by providing the reference genome and your draft contigs:

```bash
cphasing prepartition monoploid.fasta contigs.fasta -t 40 -o prepartition.clusters.txt
```




## Integration with CPhasing Pipeline
Once the prepartition.clusters.txt file is generated, you can use it in the main cphasing pipeline by specifying the `-fc (first cluster)` parameter. This allows the pipeline to skip the initial clustering and directly use the reference-guided groups for subsequent phasing.
```bash
cphasing pipeline -f contigs.fasta \
                  -hic1 hic_R1.fastq.gz \
                  -hic2 hic_R2.fastq.gz \
                  -t 40 \
                  -n 12:4 \
                  -hcr \
                  -p GATC \
                  -fc prepartition.clusters.txt
```
