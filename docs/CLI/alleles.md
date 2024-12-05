
Alleles aims to identify the allelic contig pairs by pairwise comparison.
```shell title="cphasing alleles -h"
 Usage: cphasing alleles [OPTIONS]                                              
                                                                                
 Build allele table by kmer similarity.                                         
                                                                                
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --fasta    -f        polyploid contig-level fasta.                        │
│                         (FILE)                                               │
│                         [required]                                           │
│    --output   -o        path of output allele table [default:                │
│                         fasta_prefix.allele.table]                           │
│                         (PATH)                                               │
│               -k        kmer size for similarity calculation.                │
│                         (INT)                                                │
│                         [default: 19]                                        │
│               -w        minimizer window size for similarity calculation.    │
│                         (INT)                                                │
│                         [default: 19]                                        │
│               -c        max occurance of allelic pairs                       │
│                         (INT)                                                │
│                         [default: 60]                                        │
│               -n        minimum chain end                                    │
│                         (INT)                                                │
│                         [default: 5]                                         │
│               -m        minimum k-mer similarity for similarity calculation. │
│                         (FLOAT)                                              │
│                         [default: 0.5; 0<=x<=1]                              │
│               -d        minimum different threshold between contig pairs.    │
│                         (FLOAT RANGE)                                        │
│                         [default: 0.1; 0<=x<=1]                              │
│    --threads  -t        Number of threads.                                   │
│                         (INT)                                                │
│                         [default: 4]                                         │
│    --help     -h,-help  Show this message and exit.                          │
╰──────────────────────────────────────────────────────────────────────────────╯
```



## Examples
### Total length of genome < 10 Gb
```shell
cphasing alleles -f draft.asm.fasta
```
### Total length of genome >= 10 Gb
```shell
cphasing alleles -f draft.asm.fasta -k 27 -w 14
```

!!! warning
    Now, `alleles` do not support contig with a length larger than 135 Mb. 
