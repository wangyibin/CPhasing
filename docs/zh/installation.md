---
icon: lucide/hard-drive-download
tags: 
    - Setup
---


=== "recommended" 
    从GitHub官网下载`C-Phasing`软件，并加载`activate_cphasing`激活环境。
    ```shell
    git clone https://github.com/wangyibin/CPhasing.git
    
    source ./CPhasing/bin/activate_cphasing
    ```
    !!! note
        第一次配置需要在有网络的情况下运行`./CPhasing/bin/activate_cphasing`.
    
    !!! note
        大陆用户由于网络问题可能没办法安装`pixi`这个软件，可以直接自行从github release快速下载通道下载：
        例如：
        ```bash
        mkdir -p ~/.pixi/bin
        cd ~/.pixi/bin
        wget https://gh.zwy.one/https://github.com/prefix-dev/pixi/releases/download/v0.60.0/pixi-x86_64-unknown-linux-musl.tar.gz
        tar xzvf pixi-x86_64-unknown-linux-musl.tar.gz

        ```

        之后再次运行`activate_cphasing`安装`CPhasing`依赖的包

        
        

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
        3. ++ctrl+x+y++ : 按`ctrl + X`，然后再按`Y`, 保存退出
        
    !!! note
        如果你想运行在arm(**aarch64**)平台上，请从github release页面下载，git仓库里默认存放x86-64版本
        

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
    
