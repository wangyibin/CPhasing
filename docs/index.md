<img src="logo/logo.png" alt="C-Phasing logo" width="140px" align="left" />
<h1 align="center"><b>C</b>-Phasing</h1>
<p align="center"> <b>Phasing</b> and scaffolding polyploid genomes based on Pore-<b>C</b>, Ultra-long, or Hi-<b>C</b> data</p>.

## Introduction
One of the major problems with Hi-C scaffolding of polyploid genomes is a large proportion of ambiguous short-read mapping, leading to a high-level of switched or chimeric assemblies. Now, the long-read-based chromosome conformation capture technology, e.g., **Pore-C**, provides an effective way to overcome this problem. Here, we developed a new pipeline, namely `C-Phasing`, which is specifically tailored for polyploid phasing by leveraging the advantage of Pore-C data. It also works on **Hi-C** data and diploid genome assembly.  
  
The advantages of `C-Phasing`:   
    - High speed.   
    - High anchor rate of genome.  
    - High accuracy of polyploid phasing.   

## Installation

### Via activate_cphasing (Recommended)
```shell
## Download C-Phasing and install all dependencies
git clone https://github.com/wangyibin/CPhasing.git

## activate environment (For the first configuration, run it when the network is accessible.)
./CPhasing/bin/activate_cphasing
```
!!! tip
    If you do not have direct access to the anaconda repository, you can load the image in the following way.