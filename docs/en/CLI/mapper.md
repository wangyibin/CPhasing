
`cphasing mapper` is designed for process **Pore-C** or **HiFi-C** data, its output two files:   
   (1) porec table (`.porec.gz`) which contain high-order contacts  
   (2) 4DN pairs (`.pairs.pqs`) which only retain VPCs.
!!!note
    The output of HiFi-C, still named `.porec.gz`, contains the same results whether Pore-C or HiFi-C data.
## Examples

### Process one cell Pore-C data
```shell
cphasing mapper draft.contigs.fasta sample.porec.fastq.gz -t 40 
```

### Process multiple cells Pore-C data
#### Process together 
```shell
cphasing mapper draft.contigs.fasta sample1.porec.fastq.gz sample2.porec.fastq.gz -t 40 
```
or 
```shell
cphasing mapper draft.contigs.fasta sample*.porec.fastq.gz -t 40 
```
!!! note
    It will output results using the `sample1.porec` as the prefix


#### Submit each cell to the cluster 

- use multiple scripts to run mapper for each sample
```shell title="run_sample1.sh"
cphasing mapper draft.contigs.fasta sample1.porec.fastq.gz -t 40
```
```shell title="run_sample2.sh"
cphasing mapper draft.contigs.fasta sample2.porec.fastq.gz -t 40
```
```shell title="run_sample3.sh"
cphasing mapper draft.contigs.fasta sample3.porec.fastq.gz -t 40
```

- merge results 
```shell
cphasing porec-merge sample1.porec.porec.gz sample2.porec.porec.gz sample3.porec.porec.gz -o sample.merge.porec.gz
cphasing pairs-merge sample1.porec.pairs.pqs sample2.porec.pairs.pqs sample3.porec.pairs.pqs -o sample.merge.pairs.pqs
```


### Process one cell HiFi-C data
```shell
cphasing mapper draft.contigs.fasta sample.porec.fastq.gz -t 40 --mm2-params "-x map-hifi"
```
### Process multiple cells HiFi-C data
All the steps of this same with processing in the Pore-C data, with the parameter of `--mm2-params "-x map-hifi"`.

