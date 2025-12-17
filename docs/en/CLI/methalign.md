# Methalign 

To improve mapping accuracy by utilizing information from allele-specific 5mC sites, we developed a module called **`Methalign`**. This module comprises a pipeline for correcting alignments by the fifth base (5mCG). This module is designed for ultra-complex polyploids, which contain many high-similarity homologous regions and are hard to partition (unstable).

!!!info 
    If you want to process Pore-C/HiFi-C data without 5mC information, please use the [:octicons-arrow-right-24:Mapper](mapper.md)

!!! note
    The input bam should contain the `MM/ML` tags

## Activate the environment of methalign
```shell
source activate_cphasing methalign
```
!!! info
    The network should be accessible if the methalign environment is first activated.

## Calculate the 5mCG sites of contig assembly

=== "HiFi reads"
    
    Align the HiFi reads by pbmm2

    ```shell
    pbmm2 index --preset CCS contigs.fasta index.mmi 
    pbmm2 align --preset CCS index.mmi HiFi_reads.bam | samtools view - -b -o HiFi.align.bam
    samtools sort HiFi.align.bam -o HiFi.align.sorted.bam
    samtools index HiFi.align.sorted.bam
    ```

    Calculate the 5mC sites by pb-cpg-tools
    ```shell
    aligned_bam_to_cpg_scores --bam HiFi.align.sorted.bam --ref contigs.fasta --model pileup_calling_model.v1.tflite --modsites-mode reference -q 2 
    awk '{print $1,$2,$3,$9}' OFS='\t' HiFi_align_sorted_bam_to_cpg_scores.combined.bed > output.methyl.bg
    ```

=== "ONT reads" 
    Align the ONT reads by dorado

    ```shell
    dorado aligner contigs.fasta ont_reads.bam -t 20 | samtools view -q 2 -@ 10 -b | samtools sort -@ 10 > ont.align.sorted.bam 
    samtools index ont.align.sorted.bam 
    ```

    Calculate the 5mC sites by modkit

    ```shell
    modkit pileup ont.align.sorted.bam output.methyl.bed --log-filepath pileup.log --cpg --ref contigs.fasta -t 60 --combine-strands --ignore h
    awk '{print $1,$2,$3,$11}' OFS='\t' output.methyl.bed > output.methyl.bg
    ```

## Align the candidate reads to contig assembly
```shell
dorado aligner contigs.fasta porec_reads.bam --secondary=yes > porec.align.bam
```
!!! note
    Replace `--secondary=yes` to `--mm2-opts "--secondary=yes"` when using dorado >= 0.8.0

## Refine alignments

- Refine the alignments by methylation information
```shell
methalign -t 20 contigs.fasta output.methyl.bg porec.reads*.bam -o porec.align.refined.paf.gz
```
!!! note
    This step will output `methalign.refined.paf.gz` and `methalign.refined.porec.gz`

- After refine
Input `porec.align.refined.porec.gz` 
```shell
cphasing pipeline -f contigs.fasta -pct methalign.refined.porec.gz -t 10 -n 8:4 -hcr -p AAGCTT 
```

!!! tips
    If you input a large aligned bam you can use `cphasing-rs split-bam` to split bam into several parts and run `methalign`, respectively. 
    - Split bam to speed up the refine step
    ```shell
    cphasing-rs split-bam porec.align.bam -o split_outdir/output.split
    cd split_outdir
    ```