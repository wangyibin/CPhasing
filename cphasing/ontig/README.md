# **ontig**
**ONTig**: Correct con**tig**s or identify high confidince regions by **ONT** data

|         |                                                                    |
| ------- | ------------------------------------------------------------------ |
| Main author | Jiaxin Yu ([Yujiaxin419](http://github.com/Yujiaxin419))       |
| Co-authors  | Yibin Wang ([wangyibin](http://github.com/wangyibin))          |
|             | Xingtan Zhang ([tangerzhang](https://github.com/tangerzhang/)) |

## Examples
- `split-reads`
> split fastq into several parts to and split each read by window.
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
cphasing ontig find-chimeric -p output.corrected.paf -l outputmapq.LIS.gtf -f draft.contig.fasta
```
- `hcr`
> Identifing high confidence region.
```bash
## result is output.hcr_all.bed
cphasing ontig hcr -l outputmapq.LIS.gtf -sa output.mergedSplitAlign.txt -d output.depth -b output.breakPos.txt -c draft.contig.contigsizes
```

## How to use hcr results in phasing steps
1. use corrected fasta as draft.asm.fasta
```bash
cphasing mapper corrected.fasta sample.fastq.gz
```
2. use hcr.bed to filter the mapping results 
```bash
cphasing alignments porec-intersection sample.porec.gz output.hcr_all.bed sample_hcr.porec.gz
cphasing-rs porec2pairs sample_hcr.porec.gz corrected.contigsizes -o sample_hcr.pairs.gz 
```
And then, use `sample_hcr.porec.gz` and `sample_hcr.pairs.gz` to run subsequence steps.

