### The results of the first round partition are unsatisfactory.
In our two-round partition algorithm, the first round partition depends on the h-trans errors between homologous chromosomes; if you input a contig assembly with low level switch errors or input a high accuracy pore-c data, the h-trans will be not enough to cluster all contigs to correct homologous groups, resulting in unsatisfactory results. You can set the `-q1 0` for `hyperpartition` to increase the rate of h-trans errors. However, this parameter may raise error of `out of memory` when you input huge pore-c data in porec table or hic contacts in pairs file. 

### How to set the `-n` parameter when assembling an aneuploid genome. 
The aneuploid genome, such as modern cultivated sugarcane, contains unequal homologous chromosomes. The `-n` parameter can be set to zero (`-n 10:0`) to automatically partition contigs into different chromosomes within a homologous group.     
However, we also allow the user to input a file with two columns: the first column is the index(1-base) of the first round partition, and the second column is the chromosome number of each homologous. And then specified the `-n 10:second.number.tsv` in `cphasing pipeline` or `cphasing hyperpartition`.
- `second.number.tsv`
```text
1    13
2    12
3    12
4    11
5    10
6    12
7    12
8    10
9    12
10    12
```