<img src="logo/logo.png" alt="C-Phasing logo" width="140px" align="left" />
<h1 align="center"><b>C</b>-Phasing</h1>
<p align="center"> <b>Phasing</b> and scaffolding polyploid genomes based on Pore-<b>C</b>, HiFi-<b>C</b>(<b>C</b>iFi), Ultra-long, or Hi-<b>C</b> data</p>.

## Introduction
One of the major problems with Hi-C scaffolding of polyploid genomes is a large proportion of ambiguous short-read mapping, leading to a high-level of switched or chimeric assemblies. Now, the long-read-based chromosome conformation capture technology, e.g., **Pore-C**, **HiFi-C**, provides an effective way to overcome this problem. Here, we developed a new pipeline, namely `C-Phasing`, which is specifically tailored for polyploid phasing by leveraging the advantage of Pore-C or HiFi-C data.  It also works on **Hi-C** (or other 3C technologies, e.g. Omni-C and Micro-C) data and **diploid** genome assembly.  
  
The advantages of `C-Phasing`:   
    - High speed.   
    - High anchor rate of genome.  
    - High accuracy of polyploid phasing.   

![Summary_of_CPhasing](pictures/Summary_of_CPhasing.png)

## Installation

=== "recommended" 
    Download the CPhasing and activate environments by the `activate_cphasing`
    ```shell
    git clone https://github.com/wangyibin/CPhasing.git 
    
    ### activate environment
    source ./CPhasing/bin/activate_cphasing

    ### deactivate
    source deactivate 
    ```
    !!! note
        For the first configuration, run it when the network is accessible.
    
    !!! note
        If you do not have direct access to the anaconda repository, you can set the mirror for pixi.
    !!! note
        For the platform of **`linux-aarch64`**, please download from github release.
    
=== "conda"
    Download the CPhasing and install environment by conda
    ```shell
    git clone https://github.com/wangyibin/CPhasing.git

    cd CPhasing
    conda env create -f environment.yml
    conda activate cphasing

    ```
    Add the following to the `~/.bash_profile`
    ```bash title="~/.bash_profile"
    export PATH=/path/to/CPhasing/bin:$PATH
    export PYTHONPATH=/path/to/CPhasing:$PYTHONPATH
    ```
    !!! note
        The hic pipeline require GLIBCXX_3.4.29, or you can add the following to the start of cphasing execute script, e.g.: `run.sh`
        ```bash title="run.sh"
        export LD_LIBRARY_PATH=/path/to/anaconda3/envs/cphasing/lib:$LD_LIBRARY_PATH
        ```

=== "custom"
    #### Install C-Phasing
    ```bash
    ## Download C-Phasing and install python dependencies
    git clone https://github.com/wangyibin/CPhasing.git
    cd CPhasing
    pip install .

    ```
    Add following to the `~/.bash_profile` or `~/.bashrc`
    ```bash title="~/.bash_profile" 
    export PATH=/path/to/CPhasing/bin:$PATH
    ```
    #### Dependencies
    1. For core function  
        - [bedtools](https://bedtools.readthedocs.io/en/latest/)  
        - [seqkit](https://bioinf.shenwei.me/seqkit/)  
        - [pigz](https://github.com/madler/pigz)  
    2. For Pore-C pipeline   
        - [minimap2](https://github.com/lh3/minimap2)(>= v2.24)  
    


<div class="grid cards" markdown>
-  __Basic Usage__

    ---

    Use `pipeline` to run `cphasing`.  

    [:octicons-arrow-right-24:Basic Usage](basic_usage)

- __Tutorials__
    
    ---

    Tutorials to help you know how to run `cphasing`.  
    [:octicons-arrow-right-24:Tutorials](tutorials/porec/porec_decaploid.md)  
    &emsp; [:octicons-arrow-right-24: Pore-C](tutorials/porec/porec_decaploid.md)  
    &emsp; [:octicons-arrow-right-24: Hi-C](tutorials/hic/hic_decaploid.md)

</div>


## Citation
[:octicons-arrow-right-24:Citation](citation.md)