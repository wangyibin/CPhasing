---
icon: lucide/hard-drive-download
tags: 
    - Setup
---


=== "recommended" 
    从GitHub release下载`C-Phasing`软件，并加载`activate_cphasing`激活环境。
    ```shell
    LATEST_URL=$(curl -s https://api.github.com/repos/wangyibin/CPhasing/releases/latest | grep "browser_download_url.*linux-64.tar.gz" | cut -d '"' -f 4)
    wget $LATEST_URL
    tar xzvf CPhasing*.tar.gz
    
    source ./CPhasing*/bin/activate_cphasing
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
        如果你想运行在arm(**aarch64**)平台上，请从[github release](https://github.com/wangyibin/CPhasing/releases)页面下载，git仓库里默认存放x86-64版本
        

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
    

## 系统要求与软件依赖 (System Requirements & Dependencies)

<details>
<summary><b>1. 操作系统与硬件环境要求</b> (点击展开)</summary>

- **支持的操作系统**: 在 **Linux** 系统上（如 Ubuntu >= 18.04, CentOS >= 7, Debian >= 10）进行了完整测试和支持。Windows 用户推荐通过 WSL2 运行。
- **系统架构**: 支持 `linux-64` 和 `linux-aarch64` 架构。
- **基础系统库**: 
  - Linux 内核版本 `>= 3.10.0`
  - GNU C Library (glibc) `libc >= 2.17`
</details>

<details>
<summary><b>2. 核心生物信息学软件</b> (点击展开)</summary>

- **序列比对与三维基因组构建**:
  - `minimap2` `>= 2.28, < 3`
  - `chromap` `>= 0.3.2, < 0.4`
  - `bwa-mem2` `>= 2.3, < 3`
  - `minibwa` `>= 0.1, < 0.4`
  - `minigraph` `>= 0.21, < 0.22`
  - `wfmash` `0.17.0.*`
- **序列处理与文件操作**:
  - `samtools` `>= 1.20, < 1.21`
  - `bedtools` `>= 2.31.1, < 3`
  - `seqkit` `>= 2.9.0, < 3`
  - `samblaster` `>= 0.1.26, < 0.2`
  - `pigz` `>= 2.8, < 3`
  - `crabz` `>= 0.10.0, < 0.11`
</details>

<details>
<summary><b>3. Python 运行时与库</b> (点击展开)</summary>

- **Python 运行时**: `3.12.0.*`
- **生物信息分析相关库**:
  - `biopython` `>= 1.84, < 2`
  - `pysam` `>= 0.22.1, < 0.23`
  - `cooler` `>= 0.10.2, < 0.11`
  - `hicmatrix` `>= 17.2, < 18`
  - `pyranges` `>= 0.1.2, < 0.2`
  - `ncls` `>= 0.0.68, < 0.0.69`
  - `needletail` `>= 0.7.1, < 0.8`
- **数据处理与计算性能**:
  - `pandas` `>= 2.2.3, < 3`
  - `numpy` `>= 1.26.4, < 2`
  - `polars` `>= 1.17.1, < 1.18.0`
  - `pyarrow` `>= 18.1.0, < 19`
  - `dask` `>= 2024.11.2, < 2025`
  - `joblib` `>= 1.4.2, < 2`
  - `pandarallel` `>= 1.6.5, < 2`
  - `scikit-learn` `>= 1.5.2, < 2`
  - `sparse_dot_mkl` `>= 0.9.10, < 0.10` *(仅限 linux-64 平台)*
- **图论与网络分析**:
  - `networkx` `>= 3.4.2, < 3.5`
  - `python-igraph` `>= 0.11.8, < 0.12`
  - `cdlib` `>= 0.4.0, < 0.5`
  - `graph-tool` `>= 2.97, < 3`
- **数据可视化库**:
  - `matplotlib` `>= 3.9.3, < 4`
  - `seaborn` `>= 0.13.2, < 0.14`
  - `plotly` `>= 6.2.0, < 7`
  - `plotnine` `>= 0.15.3, < 0.16`
  - `patchworklib` `>= 0.6.3, < 0.7`
  - `colormaps` `>= 0.4.2, < 0.5`
- **命令行界面及公共组件**:
  - `click` `>= 8.1.8, < 8.2`
  - `rich-click` `1.9.7.*`
  - `click-didyoumean` `>= 0.3.1, < 0.4`
</details>

<details>
<summary><b>4. 特定分析流程与可选依赖</b> (点击展开)</summary>

#### 甲基化比对环境 (`methalign`)
- `ont-modkit` `>= 0.4.3, < 0.5`
- `bammap2` `>= 0.1.7, < 0.2`
- `pbmm2` `>= 1.16.99, < 2` *(仅限 linux-64 平台)*
- `pb-cpg-tools` `>= 3.0.0, < 4` *(仅限 linux-64 平台)*

#### 基因组评估与共线性环境 (`eval` & `eval2`)
- **Python**: `3.8.*` *(eval专属环境)*
- `syri` `>= 1.6.3, < 2`
- `plotsr` `>= 1.1.1, < 2`
- `nucflag` `== 1.0.0a2`
</details>