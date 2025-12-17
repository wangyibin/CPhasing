# Rename
Rename and reorient the chromosome, according to a monoploid or closely relative species's genome.


## Rename phased assembly from `phasing` mode
The `phasing` mode of `C-Phasing` has been clustering the chromsomes into the same homologous and oriented them to parallel.
```shell
cphasing rename -r mono.fa -f draft.asm.fasta -a groups.agp -t 20 
```
![phasing](https://pic.superbed.cc/item/67960ddcfa9f77b4dc145ac9.png)


## Rename all of the chromosomes
Rename assembly from `basal` or `basal_withprune` mode, with the `--unphased` parameter.
```shell
cphasing rename -r mono.fa -f draft.asm.fasta -a groups.agp -t 20 --unphased 
```

![uphased](https://pic.superbed.cc/item/67960d93fa9f77b4dc1458f2.png)

## Use the uppercase letter as the suffix of the chromosome name
> `Chr01A Chr01B Chr01C ...`
```shell
cphasing rename -r mono.fa -f draft.asm.fasta -a groups.agp -t 20 -s upperletter
```

![upperletter](https://pic.superbed.cc/item/67960d32fa9f77b4dc1455ab.png)
