Collapsed contigs are commonly observed in polyploid hybrids due to the presence of highly similar homologous regions introduced through hybridization. These regions pose significant challenges for direct *de novo* assembly using current computational approaches (e.g. hifiasm). To address this limitation, we developed a strategy comprising two steps:     
    **(1)** Collapsed contigs detection: Identification of candidate contigs (copy number ≥2) through integrated analysis of sequencing depth profiles from HiFi, ONT, and Pore-C data.    
    **(2)** Collapsed contigs rescuing: Duplicating and putting collapsed contigs into correctly groups.
!!! warning    
    This methodology’s particular efficacy in resolving localized collapsed regions can not resolve the collapsed regions that nearly whole chromosomes.

## Collapsed contigs detection
=== "HiFi data"
    - Custom mapping
    ```shell
    minimap2 -cx map-hifi -I 16g --secondary=no -t 40 draft.asm.fasta hifi.fastq.gz | pigz -p 10 -c > hifi.align.paf.gz

    cphasing-rs paf2depth hifi.align.paf.gz -w 5000 -s 1000 -o hifi.align.depth
    cphasing collapse from-depth hifi.align.depth
    ```

    - Directly use the hitig results    
    `output.collapsed.contigs.list`

    
=== "ONT data"
    - Custom mapping    
    ```shell
    minimap2 -cx map-ont -I 16g -t 40 --secondary=no draft.asm.fasta ont.fastq.gz | pigz -p 10 -c > ont.align.paf.gz

    cphasing-rs paf2depth hifi.align.paf.gz -w 5000 -s 1000 -o ont.align.depth
    cphasing collapse from-depth ont.align.depth
    ```
    - Directly use the hitig results    
    `output.collapsed.contigs.list`

=== "Pore-C data"
    The `porec.align.paf.gz` generated from `cphasing mapper`.
    ```shell
    cphasing-rs paf2depth porec.align.paf.gz -w 5000 -s 1000 -o porec.align.depth
    cphasing collapse from-depth porec.align.depth
    ```

=== "GFA"
    The GFA file generated from `hifiasm` record the read number, which can be used to identify the collapsed contigs.
    ```shell
    cphasing collapse from-gfa hifi.p_utg.noseq.gfa.gz 
    ```

## Collapsed contigs rescuing
=== "Pore-C/CiFi"

    ```shell
        cphasing collapse rescue sample.porec.gz draft.asm.contigsizes 3.hyperpartition/output.clusters.txt contigs.collapsed.contig.list -n 4 -at 3.hyperpartition/draft.asm.allele.table
    ```

=== "Hi-C"

    ```shell
        cphasing collapse rescue sample.hic.pairs.pqs draft.asm.contigsizes 3.hyperpartition/output.clusters.txt contigs.collapsed.contig.list -n 4 -at 3.hyperpartition/draft.asm.allele.table
    ```

!!! note
    Currently, we output the rescue result into a `collapsed.rescue.clusters.txt` file, which the user executes the next step of `4.scaffolding` manually. And the input file can directly use the previous.


## After scaffolding

- Generate new contig-level fasta and agp, which renamed the duplicated contigs (e.g. utg000001l -> utg000001_d2)
The name of collapsed contigs in agp is need to rename for subsequence processes.
```bash
cphasing agp2fasta groups.rescued.agp draft.asm.fasta --contigs > draft.dup.fasta
cphasing collapse agp-dup groups.rescued.agp > groups.dup.agp
```

- Generate a new `pairs.gz` or `pairs.pqs` file 
```shell
cphasing collapse pairs-dup sample.pairs.pqs collapsed.rescue.contigs.list -o sample.dup.pairs.pqs 
```

- Curation by juicebox
```shell
cphasing utils agp2assembly groups.dup.agp > groups.dup.assembly
cphasing pairs2mnd sample.dup.pairs.pqs -o sample.dup.mnd.txt
```

- Rename
```shell
cphasing rename -r ref.fa -f draft.dup.fasta -a groups.dup.agp -t 40 
```

- Plot
```shell
cphasing pairs2cool sample.dup.pairs.pqs sample.dup.pairs.pqs/_contigsizes sample.10k.cool 
cphasing plot -a groups.dup.agp -m sample.10k.cool -o sample.500k.wg.png 
```