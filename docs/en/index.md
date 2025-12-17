<img src="https://raw.githubusercontent.com/wangyibin/CPhasing/main/docs/logo/logo.png" alt="C-Phasing logo" width="140px" align="left" />
<h1 align="center"><b>C</b>-Phasing</h1>
<p align="center"> <b>Phasing</b> and scaffolding polyploid genomes based on Pore-<b>C</b>, HiFi-<b>C</b>(<b>C</b>iFi), Ultra-long, or Hi-<b>C</b> data</p>.

## Introduction
One of the major problems with Hi-C scaffolding of polyploid genomes is a large proportion of ambiguous short-read mapping, leading to a high-level of switched or chimeric assemblies. Now, the long-read-based chromosome conformation capture technology, e.g., **Pore-C**, **HiFi-C**, provides an effective way to overcome this problem. Here, we developed a new pipeline, namely `C-Phasing`, which is specifically tailored for polyploid phasing by leveraging the advantage of Pore-C or HiFi-C data.  It also works on **Hi-C** (or other 3C technologies, e.g. Omni-C and Micro-C) data and **diploid** genome assembly.  
  
The advantages of `C-Phasing`:   
    - High speed.   
    - High anchor rate of genome.  
    - High accuracy of polyploid phasing.   

![Summary_of_CPhasing](https://raw.githubusercontent.com/wangyibin/CPhasing/main/pictures/Summary_of_CPhasing.png)



<div class="grid cards" markdown>
-  __Basic Usage__

    ---
    [:lucide-hard-drive-download: Installation](installation.md){ .md-button }

    [:material-clock-fast: Getting started](basic_usage.md){ .md-button .md-button--primary }

- __Tutorials__
    
    ---

    Tutorials to help you know how to run `cphasing`.  
    [:lucide-book-open-text: Tutorials](tutorials/porec/porec_decaploid.md){ .md-button }

    &emsp; [:octicons-arrow-right-24: Pore-C](tutorials/porec/porec_decaploid.md) 

    &emsp; [:octicons-arrow-right-24: Hi-C](tutorials/hic/hic_decaploid.md)

</div>


## Citation
[:octicons-arrow-right-24:Citation](citation.md)