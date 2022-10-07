# C-Phasing
phasing and scaffolding polyploid genomes based on Pore-C data

## Dependencies
- [chromap](https://github.com/haowenz/chromap)
- [gmap]
- hisat2
- samtools
- bedtools

## Installation
```
git clone https://github.com/wangyibin/CPhasing.git
export PATH=/path/to/CPhasing/bin:$PATH
export PYTHONPATH=/path/to/CPhasing:$PYTHONPATH
```

## Citiing

If you use the hic data to scaffold genomes please cite ALLHiC as while.  
Zhang, X. ,  Zhang, S. ,  Zhao, Q. ,  Ming, R. , &  Tang, H. . (2019). Assembly of allele-aware, chromosomal-scale autopolyploid genomes based on hi-c data. Nature Plants, 5(5). doi: [10.1038/s41477-019-0487-8](https://doi.org/10.1038/s41477-019-0487-8)

```bibtex
@article{2019Assembly,
        title={Assembly of allele-aware, chromosomal-scale autopolyploid genomes based on Hi-C data},
        author={ Zhang, Xingtan  and  Zhang, Shengcheng  and  Zhao, Qian  and  Ming, Ray  and  Tang, Haibao },
        journal={Nature Plants},
        volume={5},
        number={5},
        year={2019},
}
```