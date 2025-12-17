
Alleles aims to identify the allelic contig pairs by pairwise comparison.

## Examples
### Total length of genome < 10 Gb
```shell
cphasing alleles -f draft.asm.fasta
```
### Total length of genome >= 10 Gb
```shell
cphasing alleles -f draft.asm.fasta -k 27 -w 14
```
