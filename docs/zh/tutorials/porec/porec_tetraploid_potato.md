

## Provide a known homologous chromosome group 

### Prepartition by monoploid or related specices 

```shell
cphasing prepartition DM.chrom.fasta hifi.p_utg.fasta -t 60 -o prepartition.clusters.txt 
```

### Run cphasing pipeline 
```shell
cphasing pipeline -f hifi.p_utg.fasta -paf porec.reads.paf.gz -hcr -p GATC --collapsed-rescue --first-cluster prepartition.clusters.txt -t 30 
```


