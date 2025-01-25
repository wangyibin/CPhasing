<img src="logo/logo.png" alt="C-Phasing logo" width="140px" align="left" />
<h1 align="center"><b>C</b>-Phasing</h1>
<p align="center"> <b>Phasing</b> and scaffolding polyploid genomes based on Pore-<b>C</b>, HiFi-<b>C</b>, Ultra-long, or Hi-<b>C</b> data</p>.

## 介绍
使用Hi-C数据实现多倍体基因组分型挂载，存在主要问题之一是大量不明确的短读长比对，这容易导致高水平的交换或嵌合组装错误。现在，基于长读长的染色体构象捕获技术，如**Pore-C**、**HiFi-C**，为克服这一问题提供了有效途径。在这里，我们开发了一个新的流程，即“C-Phasing”，它是专门为多倍体分型组装量身定制的，旨在充分利用Pore-C或HiFi-C数据的优势。此外，它也可用于**Hi-C**数据和二倍体基因组组装。

`C-Phasing`的优势:   
    - 快速.   
    - 更高的挂载率.  
    - 更高的分型组装质量.

## 安装

=== "recommended" 
    从GitHub官网下载`C-Phasing`软件，并加载`activate_cphasing`激活环境。
    ```shell
    git clone https://github.com/wangyibin/CPhasing.git
    
    source ./CPhasing/bin/activate_cphasing
    ```
    !!! note
        第一次配置需要在有网络的情况下运行`./CPhasing/bin/activate_cphasing`.
    
    !!! note
        如果你的服务器不能访问anaconda的仓库，可以配置镜像源，如配置浙江大学的镜像源：  
        1. 打开配置文件
        ```bash
        ~/.pixi/bin/pixi config edit
        ```
        2. 写入以下内容
        ```yaml
        [mirrors]
        "https://conda.anaconda.org/conda-forge" = [
            "https://mirrors.zju.edu.cn/anaconda/cloud/conda-forge"
        ]

        "https://conda.anaconda.org/bioconda" = [
            "https://mirrors.zju.edu.cn/anaconda/cloud/bioconda"
        ]
        ```
        3. 按`ctrl + X`，然后再按`Y`, 保存退出
        
        

=== "conda"
    从GitHub官网下载`C-Phasing`软件，使用anaconda或者miniconda配置环境和激活环境。
    ```shell
    git clone https://github.com/wangyibin/CPhasing.git

    cd CPhasing
    conda env create -f environment.yml
    conda activate cphasing

    ```
    接着配置以下命令道环境变量 `.bash_profile`或者`.bashrc`
    ```bash title="~/.bash_profile"
    export PATH=/path/to/CPhasing/bin:$PATH
    export PYTHONPATH=/path/to/CPhasing:$PYTHONPATH
    ```
    !!! note
        在比对Hi-C数据的时候，有用户遇到找不到`GLIBCXX_3.4.29`的报错，这种情况可以在运行`cphasing`命令前，配置以下命令：
        ```bash title="run.sh"
        export LD_LIBRARY_PATH=/path/to/anaconda3/envs/cphasing/lib:$LD_LIBRARY_PATH
        ```

=== "custom"
    #### 从GitHub下载`CPhasing`, 然后使用`pip`安装。
    ```bash
    ## Download C-Phasing and install python dependencies
    git clone https://github.com/wangyibin/CPhasing.git
    cd CPhasing
    pip install .

    ```
    接着配置以下命令道环境变量 `~/.bash_profile` 或者`~/.bashrc`
    ```bash title="~/.bash_profile" 
    export PATH=/path/to/CPhasing/bin:$PATH
    ```
    #### 第三方依赖包
    1. 核心功能依赖：  
        - [bedtools](https://bedtools.readthedocs.io/en/latest/)  
        - [seqkit](https://bioinf.shenwei.me/seqkit/)  
        - [pigz](https://github.com/madler/pigz)  
    2. Pore-C流程依赖：   
        - [minimap2](https://github.com/lh3/minimap2)(>= v2.24)  
    


<div class="grid cards" markdown>
-  __基本用法__

    ---

    使用 `pipeline`命令开始运行`cphasing`.  

    [:octicons-arrow-right-24:基本用法](basic_usage)

- __教程__
    
    ---

    一些案例教你如何运行 `cphasing`.  
    [:octicons-arrow-right-24:Tutorials](tutorials/porec/porec_decaploid.md)

</div>