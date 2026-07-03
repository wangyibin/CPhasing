# Assemble hybrid sugarcane by Hi-C


Aneuploidy (2n=114), estimate genome size ~ 10 Gb 

## Hi-C data processing
Due to the huge porec data, we mapped each cell to fasta by submitting this command to cluster, respectively. And merged them into one file by `cphasing pairs-merge`.
!!!note 
    Please submit one first; wait until the index is successfully created and the remaining jobs are enabled to submit to the cluster.

```shell
cphasing hic mapper -f sh_hifi.bp.p_utg.fasta -hic1 hic-1_R1.fastq.gz -hic2 hic-1_R2.fastq.gz -t 40 -k 27 -w 14 
```

```shell
cphasing hic mapper -f sh_hifi.bp.p_utg.fasta -hic1 hic-2_R1.fastq.gz -hic2 hic-2_R2.fastq.gz -t 40 -k 27 -w 14 
cphasing hic mapper -f sh_hifi.bp.p_utg.fasta -hic1 hic-3_R1.fastq.gz -hic2 hic-3_R2.fastq.gz -t 40 -k 27 -w 14 
cphasing hic mapper -f sh_hifi.bp.p_utg.fasta -hic1 hic-4_R1.fastq.gz -hic2 hic-4_R2.fastq.gz -t 40 -k 27 -w 14 
```

```shell
cphasing pairs-merge hic-*.pairs.pqs -o hic.merge.pairs.pqs 
```


## Assembling by `cphasing pipeline`
>Modern hybrid sugarcane is an aneuploid genome, containing unequal numbers of chromosomes across homologous groups. Therefore, we set `-n 10:0` to automatically determine the number of clusters.

> To rescue collapsed regions, enable CollapseRescue by providing the assembly graph with `--collapsed-rescue --gfa sh_hifi.bp.p_utg.no_seq.gfa`, allowing collapsed unitigs to be identified and reassigned.

> Optionally, the `--refine` parameter can be enabled to further reduce inter-homologous misassignments and improve phasing accuracy.

```shell
cphasing pipeline -f sh_hifi.bp.p_utg.fasta -pct hic.mrege.pairs.pqs -t 40 -n 10:0 -hcr -p AAGCTT --collapsed-rescue --gfa sh_hifi.bp.p_utg.no_seq.gfa
```
