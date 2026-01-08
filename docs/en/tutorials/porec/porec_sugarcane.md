# Assemble hybrid sugarcane by Pore-C


Aneuploidy (2n=114), estimate genome size ~ 10 Gb 

## Pore-C data processing
Due to the huge porec data, we mapped each cell to fasta by submitting this command to cluster, respectively. And merged them into one file by `cphasing porec-merge`
```shell
cphasing mapper sh_hifi.bp.p_utg.fasta porec-1.fastq.gz -t 40 
cphasing mapper sh_hifi.bp.p_utg.fasta porec-2.fastq.gz -t 40 
cphasing mapper sh_hifi.bp.p_utg.fasta porec-3.fastq.gz -t 40 
cphasing mapper sh_hifi.bp.p_utg.fasta porec-4.fastq.gz -t 40 
cphasing mapper sh_hifi.bp.p_utg.fasta porec-5.fastq.gz -t 40 
```

```shell
cphasing porec-merge porec-*.porec.gz -o porec.merge.porec.gz 
```


## Assembling by `cphasing pipeline`
- Modern hybrid sugarcane is an aneuploid, which contains an unequal number of chromosomes in each homologous group, so we set `-n 10:0` to automatically output cluster numbers.
```shell
cphasing pipeline -f sh_hifi.bp.p_utg.fasta -pct porec.mrege.porec.gz -t 40 -n 10:0 -hcr -p AAGCTT 
```
