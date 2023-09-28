# **ontig**

## Examples
- `split-reads`
```
cphasing ontig split-reads -i 20220608-UNL303-P6-PAK13867.pass.fastq.gz -t 4
```
- `correct-alignments`
> mapping split reads to draft assembly and correct the alignments 
```bash
cphasing ontig correct-alignments -f draft.contig.fasta -i split_fastq/20220608-UNL303-P6-PAK13867.pass.part_001_5.0k.fastq,split_fastq/20220608-UNL303-P6-PAK13867.pass.part_002_5.0k.fastq,split_fastq/20220608-UNL303-P6-PAK13867.pass.part_003_5.0k.fastq,split_fastq/20220608-UNL303-P6-PAK13867.pass.part_004_5.0k.fastq,split_fastq/20220608-UNL303-P6-PAK13867.pass.part_005_5.0k.fastq,split_fastq/20220608-UNL303-P6-PAK13867.pass.part_006_5.0k.fastq,split_fastq/20220608-UNL303-P6-PAK13867.pass.part_007_5.0k.fastq,split_fastq/20220608-UNL303-P6-PAK13867.pass.part_008_5.0k.fastq,split_fastq/20220608-UNL303-P6-PAK13867.pass.part_009_5.0k.fastq,split_fastq/20220608-UNL303-P6-PAK13867.pass.part_010_5.0k.fastq -t 20
```
- `find-chimeric`
> find split-alignments and idendify the chimeric contigs
```bash
cphasing ontig find-chimeric -p output.corrected.paf -l outputLIS.gtf -f draft.contig.fasta
```
- `hcr`
> Identifing high confidence region.
```bash
cphasing ontig hcr -l outputLIS.gtf -sa output.mergedSplitAlign.txt -d output.depth
```

