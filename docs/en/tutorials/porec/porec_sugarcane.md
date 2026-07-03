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
>Modern hybrid sugarcane is an aneuploid genome, containing unequal numbers of chromosomes across homologous groups. Therefore, we set `-n 10:0` to automatically determine the number of clusters.

> To rescue collapsed regions, enable CollapseRescue by providing the assembly graph with `--collapsed-rescue --gfa sh_hifi.bp.p_utg.no_seq.gfa`, allowing collapsed unitigs to be identified and reassigned.

> Optionally, the `--refine` parameter can be enabled to further reduce inter-homologous misassignments and improve phasing accuracy.
> Optionally, the `-x very-sensitive` parameter can be enabled for fragmented assemblies (HiFi-only) or low-depth sequencing data.

```shell 
cphasing pipeline -f sh_hifi.bp.p_utg.fasta -pct porec.mrege.porec.gz -t 40 -n 10:0 -hcr -p AAGCTT --collapsed-rescue --gfa sh_hifi.bp.p_utg.no_seq.gfa
```
