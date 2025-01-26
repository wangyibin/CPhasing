# Rename
以单倍型或者近缘物种的参考基因组作为参考，将当前组装的多倍体进行重命名和染色体方向调整.

## 调整`phasing` mode的组装结果
默认情况下`phasing` mode已经将同源染色体聚成一组，并且不同同源染色体将的方向已经调整至平行，为了节省时间开销，我们仅将第一个单倍型即`Chr??.g1`比对至参考基因组上，剩下的单倍型命名和方向跟随`g1`进行调整。如果不符合以上说明，则需加`--unphased`参数对所有染色体进行比对。
```shell
cphasing rename -r mono.fa -f draft.asm.fasta -a groups.agp -t 20 
```
![phasing](https://pic.superbed.cc/item/67960ddcfa9f77b4dc145ac9.png)


## 对所有染色体进行重命名
`basal`或者`basal_withprune` mode，需要对所有染色体进行比对，然后重命名和重新调整方向，需要加`--unphased`参数。
```shell
cphasing rename -r mono.fa -f draft.asm.fasta -a groups.agp -t 20 --unphased 
```

![uphased](https://pic.superbed.cc/item/67960d93fa9f77b4dc1458f2.png)

## 用字母作为不同单倍型的后缀
> `Chr01A Chr01B Chr01C ...`
```shell
cphasing rename -r mono.fa -f draft.asm.fasta -a groups.agp -t 20 -s upperletter
```

![upperletter](https://pic.superbed.cc/item/67960d32fa9f77b4dc1455ab.png)

## Parameters
```shell
 Usage: cphasing rename [OPTIONS]                                                                      
                                                                                                       
 Rename and orientation the groups according to a refernce.                                            
 To speed up this function, we only align the first haplotype to reference, because the orientation    
 among different haplotypes has been paralleled in scaffolding. If you want to orient all of the       
 haplotypes you can specify the --unphased parameters.                                                 
                                                                                                       
╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --ref           -r         Genome of monoploid or closely related species.                       │
│                               (REFERENCE FASTA)                                                     │
│                               [required]                                                            │
│ *  --fasta         -f         Path to draft assembly                                                │
│                               (DRAFT FASTA)                                                         │
│                               [required]                                                            │
│ *  --agp           -a         agp file.                                                             │
│                               (AGP)                                                                 │
│                               [required]                                                            │
│    --unphased      -unphased  Haplotypes are not parallel aligned. Rename all of the chromosomes.   │
│    --suffix-style  -s         The suffix style of renamed chromosome among different haplotypes.    │
│                                 number: g1, g2, g3, ... (default)                                   │
│                                 upperletter: A, B, C, ...                                           │
│                                 lowerletter: a, b, c, ...                                           │
│                               (STR)                                                                 │
│                               [default: number]                                                     │
│    --output        -o         Output path of renamed agp file                                       │
│                               (OUTPUT_AGP)                                                          │
│    --threads       -t         Number of threads.                                                    │
│                               (INT)                                                                 │
│                               [default: 8]                                                          │
│    --help          -h,-help   Show this message and exit.                                           │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────╯
```