# Prepartition

The `prepartition` step implements the first-round clustering by mapping contigs to a monoploid genome. This process groups the contigs into different homologous groups prior to Hi-C/Pore-C phasing.

## 1. Run Prepartition

Cluster the assembled contigs based on a reference monoploid genome.

```shell
cphasing prepartition \
  monoploid.fasta \
  contigs.fasta \
  -t 60 \
  -o prepartition.clusters.txt
```


## 2. Pass to the Pipeline

Once you obtained the prepartition result, pass it to the main `cphasing pipeline` using the `-fc` (first cluster) parameter.

```shell
cphasing pipeline \
  -f contigs.fasta \
  -pcd porec_reads.fastq.gz \
  -fc prepartition.clusters.txt \
  -t 60 
```

**Key Parameter:**
* `-fc / --first-cluster`: The cluster file generated in the prepartition step.