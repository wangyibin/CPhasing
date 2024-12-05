
`cphasing mapper` is designed for process Pore-C data, it output two file: (1) porec table (`.porec.gz`) which contain high-order contacts and (2) 4DN pairs (`pairs.gz`) which only retain VPCs.
```shell title="cphasing mapper -h"
 Usage: cphasing mapper [OPTIONS] REFERENCE FASTQ...   
 
 Mapper for pore-c reads.                                                  
 REFERENCE: Path of reference                                              
                                                                           
 FASTQ: Path of pore-c reads, multiple file enabled, the prefix of output  
 default only use sample 1.                                                
                                                                           
╭─ Arguments ─────────────────────────────────────────────────────────────╮
│ *  REFERENCE    (PATH) [required]                                       │
│ *  FASTQ        (PATH) [required]                                       │
╰─────────────────────────────────────────────────────────────────────────╯
╭─ Options ───────────────────────────────────────────────────────────────╮
│ --enzyme        -e        Restrict site pattern, use comma to separate  │
│                           multiple patterns.                            │
│                           (STR)                                         │
│                 -k        kmer size for mapping.                        │
│                           (INT)                                         │
│                           [default: 15]                                 │
│                 -w        minimizer window size for mapping.            │
│                           (INT)                                         │
│                           [default: 10]                                 │
│ --mm2-params              additional parameters for minimap2            │
│                           (STR)                                         │
│                           [default: -x map-ont]                         │
│ --mapq          -q        Minimum quality of mapping [0, 60].           │
│                           (INT)                                         │
│                           [default: 0; 0<=x<=60]                        │
│ --min-identity  -p        Minimum percentage identity of alignments [0, │
│                           1.0].                                         │
│                           (FLOAT)                                       │
│                           [default: 0.8; 0.0<=x<=1.0]                   │
│ --min-length    -l        Minimum length of fragments.                  │
│                           (INT)                                         │
│                           [default: 150]                                │
│ --max-edge      -me       Maximum length of fragment located in the     │
│                           edge of contigs.                              │
│                           (INT)                                         │
│                           [default: 2000]                               │
│ --force         -f        Force run all the command, ignore existing    │
│                           results. The index file also will be removed. │
│ --outprefix     -o        output prefix, if none use the prefix of      │
│                           fastq                                         │
│                           (TEXT)                                        │
│ --threads       -t        Number of threads.                            │
│                           (INT)                                         │
│                           [default: 4]                                  │
│ --help          -h,-help  Show this message and exit.                   │
╰─────────────────────────────────────────────────────────────────────────╯
```

## Examples

### process one cell porec data
```shell
cphasing mapper draft.contigs.fasta sample.porec.fastq.gz -t 40 
```

### process multiple cell porec data
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
cphasing mapper draft.contigs.fasta sample1.poref.fastq.gz -t 40
```
```shell title="run_sample2.sh"
cphasing mapper draft.contigs.fasta sample2.poref.fastq.gz -t 40
```
```shell title="run_sample3.sh"
cphasing mapper draft.contigs.fasta sample3.poref.fastq.gz -t 40
```

- merge results 
```shell
cphasing-rs porec-merge sample1.porec.porec.gz sample2.porec.porec.gz sample3.porec.porec.gz -o sample.merge.porec.gz
cphasing-rs pairs-merge sample1.porec.pairs.gz sample2.porec.pairs.gz sample3.porec.pairs.gz -o sample.merge.pairs.gz
```
