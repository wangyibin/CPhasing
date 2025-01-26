# HiC Mapper
## one cell
```shell
cphasing hic mapper -f draft.asm.fasta -1 hic_R1.fastq.gz -2 hic_R2.fastq.gz -t 40
```
!!!note
    If the total length of your input genome is larger than 8 Gb, the `-k 27 -w 14` should be specified, to avoid the error of `chromap`. 


## mutiple cells

#### Submit each cell to the cluster 
- use multiple scripts to run mapper for each sample  
```shell  title="run_sample1.sh"
cphasing hic mapper -f draft.asm.fasta -1 hic-1_R1.fastq.gz -2 hic-1_R2.fastq.gz -t 40
```

    ```shell  title="run_sample2.sh"
    cphasing hic mapper -f draft.asm.fasta -1 hic-2_R1.fastq.gz -2 hic-2_R2.fastq.gz -t 40
    ```

    ```shell  title="run_sample3.sh"
    cphasing hic mapper -f draft.asm.fasta -1 hic-3_R1.fastq.gz -2 hic-3_R2.fastq.gz -t 40
    ```
!!!note
    Please submit one first; wait until the index is successfully created and the remaining jobs are enabled to submit to the cluster.

- merge results
```shell
cphasing pairs-merge hic-*.pairs.gz -o hic.merge.pairs.gz
```


## Parameters
```shell
                                                                                                  
 Usage: cphasing hic mapper [OPTIONS]                                                             
                                                                                                  
 Mapper for reads mapping.                                                                        
                                                                                                  
╭─ Options ──────────────────────────────────────────────────────────────────────────────────────╮
│ *  --fasta,--reference  -f,-r     Path of reference fasta file.                                │
│                                   (FILE)                                                       │
│                                   [required]                                                   │
│ *  --read1              -1        Path of read 1.                                              │
│                                   (FILE)                                                       │
│                                   [required]                                                   │
│ *  --read2              -2        Path of read 2.                                              │
│                                   (FILE)                                                       │
│                                   [required]                                                   │
│                         -k        kmer size for mapping.                                       │
│                                   (INT)                                                        │
│                                   [default: 17]                                                │
│                         -w        minimizer window size for mapping.                           │
│                                   (INT)                                                        │
│                                   [default: 7]                                                 │
│    --mapq               -q        Minimum quality of mapping [0, 60].                          │
│                                   (INT)                                                        │
│                                   [default: 0; 0<=x<=60]                                       │
│    --aligner            -a        Aligner executable. _chromap is the modifed version in       │
│                                   C-Phasing, if you want to use the offical version you can    │
│                                   set aligner to chromap                                       │
│                                   (_chromap|chromap)                                           │
│                                   [default: _chromap]                                          │
│    --threads            -t        Number of threads.                                           │
│                                   (INT)                                                        │
│                                   [default: 4]                                                 │
│    --help               -h,-help  Show this message and exit.                                  │
╰────────────────────────────────────────────────────────────────────────────────────────────────╯
```