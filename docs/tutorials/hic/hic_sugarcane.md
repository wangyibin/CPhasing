# Assemble hybrid sugarcane by Hi-C


Aneuploidy (2n=114), estimate genome size ~ 10 Gb 

## Hi-C data processing
Due to the huge porec data, we mapped each cell to fasta by submitting this command to cluster, respectively. And merged them into one file by `cphasing pairs-merge`.
!!!note 
    Please submit one first; wait until the index is successfully created and the remaining jobs are enabled to submit to the cluster.

```shell
cphasing hic mapper -f sh_hifi.bp.p_utg.fasta -hic1 hic-1_R1.fastq.gz -hic1 hic-1_R1.fastq.gz -t 40 -k 27 -w 14 
```

```shell
cphasing hic mapper -f sh_hifi.bp.p_utg.fasta -hic1 hic-2_R1.fastq.gz -hic1 hic-2_R1.fastq.gz -t 40 -k 27 -w 14 
cphasing hic mapper -f sh_hifi.bp.p_utg.fasta -hic1 hic-3_R1.fastq.gz -hic1 hic-3_R1.fastq.gz -t 40 -k 27 -w 14 
cphasing hic mapper -f sh_hifi.bp.p_utg.fasta -hic1 hic-4_R1.fastq.gz -hic1 hic-4_R1.fastq.gz -t 40 -k 27 -w 14 
```

```shell
cphasing pairs-merge hic-*.pairs.pqs -o hic.merge.pairs.pqs 
```


## Assembling by `cphasing pipeline`
Modern hybrid sugarcane is an aneuploid, which contains an unequal number of chromosomes in each homologous group, so we set `-n 0:0` to automatically output cluster numbers.
```shell
cphasing pipeline -f sh_hifi.bp.p_utg.fasta -pct hic.mrege.pairs.pqs -t 40 -n 0:0 -hcr -p AAGCTT 
```

