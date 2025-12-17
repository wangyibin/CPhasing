<img src="https://raw.githubusercontent.com/wangyibin/CPhasing/main/docs/logo/logo.png" alt="C-Phasing logo" width="140px" align="left" />
<h1 align="center"><b>C</b>-Phasing</h1>
<p align="center"> <b>Phasing</b> and scaffolding polyploid genomes based on Pore-<b>C</b>, HiFi-<b>C</b>, Ultra-long, or Hi-<b>C</b> data</p>.

## 介绍
使用Hi-C数据实现多倍体基因组分型挂载，存在主要问题之一是大量不明确的短读长比对，这容易导致高水平的交换或嵌合组装错误。现在，基于长读长的染色体构象捕获技术，如**Pore-C**、**HiFi-C**，为克服这一问题提供了有效途径。在这里，我们开发了一个新的流程，即“C-Phasing”，它是专门为多倍体分型组装量身定制的，旨在充分利用Pore-C或HiFi-C数据的优势。此外，它也可用于**Hi-C**或**Micro-C**等基于二代的染色质构象捕获技术产生的数据和二倍体基因组组装。

`C-Phasing`的优势:   
    - 快速.   
    - 更高的挂载率.  
    - 更高的分型组装质量.
![Summary_of_CPhasing](https://raw.githubusercontent.com/wangyibin/CPhasing/main/pictures/Summary_of_CPhasing.png)




<div class="grid cards" markdown>
-  __基本用法__

    ---
    [:lucide-hard-drive-download: Installation](installation.md){ .md-button }

    [:material-clock-fast: Getting started](basic_usage.md){ .md-button .md-button--primary }


- __教程__
    
    ---

    一些案例教你如何运行 `cphasing`.  
    [:lucide-book-open-text: Tutorials](tutorials/porec/porec_decaploid.md){ .md-button }

</div>


[:octicons-arrow-right-24:Citation](citation.md)