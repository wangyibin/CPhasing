## Simulation


### Chimeric and collapsed simulation
**Note: These scripts depend on chromosome name and contig name**
> Default support for `1A.ctg1, 1A.ctg2 ... 1B.ctg1 ... 2B.ctg1`.  
> Need input contig-level assembly
- chimeric contigs  
```bash
simulate_chimeric.py ploidy-12.n500k.fasta --chimeric_ratio 0.2 > simulate_chimeric.0.2.fasta
```
- collapsed contigs  
```bash
simulate_collapse.py ploidy-2.2.n500k.fasta --ratio 0.2 > collapsed_0.2.fasta
```