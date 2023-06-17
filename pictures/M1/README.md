```bash
minimap2 -x map-ont -c M1_hifi_ul_50k.p_utg.fasta M1-2.fastq.gz | pigz -p 10 -c > M1-2.paf.gz 
...
cphasing hyperpartition M1-2.edges M1_hifi_ul_50k.p_utg.contigsizes norm.sorted.default.clusters.txt -pt norm.sorted.prune.table -k 8:4 -t 5 -inc
```
`cphasing`, v0.0.22
