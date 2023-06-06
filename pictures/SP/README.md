```bash
minimap2 -x map-ont -k 27 -w 14 -c SweetPo_nhap6_30x.p_utg.fasta GS-3.fastq.gz | pigz -p 10 -c > GS-3.paf.gz 
...
cphasing hyperpartition GS-3.edges SweetPo_nhap6_30x.p_utg.contigsizes min_weight-1.5.norm.sorted.clusters.txt --prune test.norm.sorted.contig.table -k 15:6 -r2 1.5 -t 5 -inc
```