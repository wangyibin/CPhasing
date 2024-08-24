# **hitig**
**hitig**: Get **hi**gh quality con**tig**s or regions (**Hi**gh confidence regions) by ONT/HiFi sequencing data

|         |                                                                          |
| ------- | -------------------------------------------------------------------------|
| Main author | Jiaxin Yu ([Yujiaxin419](http://github.com/Yujiaxin419))             |
| Co-authors  | Yibin Wang ([wangyibin](http://github.com/wangyibin))                |
|             | Xingtan Zhang ([tangerzhang](https://github.com/tangerzhang/))       |

## Examples
### One command 

```bash
hitig pipeline -f draft.contig.fasta -i sample.fastq.gz -n ploidy-level
```

### step by step
- `correct-alignments`
> mapping split reads to draft assembly and correct the alignments 
```bash
hitig correct-alignments -f draft.contig.fasta -i sample.fastq.gz -t 20 -n ploidy-level
```
- `find-chimeric`
> find split-alignments and idendify the chimeric contigs
```bash
hitig find-chimeric -p output.corrected.paf -l output.mapq.LIS.gtf -f draft.contig.fasta
```
- `hcr`
> Identifing high confidence region.
```bash
## result is output.hcr_all.bed
hitig hcr \
    -l output.mapq.LIS.gtf \
    -sa output.mergedSplitAlign.txt \
    -d output.depth -b output.breakPos.txt \
    -f draft.contig.fasta \ 
    -p output.paf
```

## How to use hitig results in phasing steps
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

