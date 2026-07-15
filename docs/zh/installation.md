---
icon: lucide/hard-drive-download
tags: 
    - Setup
---

# 安装 (Installation)

## 系统要求 (System Requirements)

CPhasing 使用现代包管理器极大简化了依赖项管理。所有必备的软件工具和库文件都可以被自动解析与安装。

### 1. 操作系统与系统架构
- **支持的操作系统**: 原生 **Linux**（例如 Ubuntu >= 18.04, CentOS >= 7, Debian >= 10）。Windows 用户可以通过 WSL2 运行。macOS 用户可以通过 Docker 运行。
- **系统架构**: 支持 `linux-64` 和 `linux-aarch64`。
- **基础系统库**: 
  - Linux 内核版本 `>= 3.10.0`
  - GNU C Library (glibc) `libc >= 2.17`

### 2. 环境管理器（必需）
要安装 CPhasing，你需要准备好以下环境管理器之一：
- **Pixi（推荐）**: 极速、跨平台的包管理器。Pixi 会自动下载并隔离运行环境中的所有外部生物信息学工具和 Python 包。
- **Conda / Miniconda**: 支持作为传统管理环境的备选方案。

!!! note "自动依赖解析"
    你**不需要**手动安装具体的第三方工具（如 `minimap2`、`samtools`、`bedtools`）或 Python 库。在通过 **Pixi** 或 **Conda** 进行安装时，所有这一切都预配置好并会自动被解析。

---

## 安装指南 (Installation Instructions)

=== "推荐方式 (Pixi)" 
    下载预编译的 CPhasing 组件并使用 `activate_cphasing` 脚本激活环境：
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
        

=== "Conda"
    从 GitHub 仓库下载 C-Phasing 代码，并使用 Anaconda 或 Miniconda 导入环境依赖：
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

=== "自定义安装 (Custom)"
    如果你希望手动管理环境和路径，可以通过 pip 直接从源码构建：
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
    

## 软件具体依赖包明细列表 (Detailed Software & Library Dependencies)

<details>
<summary><b>1. 核心生物信息学工具</b> (点击展开)</summary>

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
<summary><b>2. Python 运行时与核心库</b> (点击展开)</summary>

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
- **数据可视化**:
  - `matplotlib` `>= 3.9.3, < 4`
  - `seaborn` `>= 0.13.2, < 0.14`
  - `plotly` `>= 6.2.0, < 7`
  - `plotnine` `>= 0.15.3, < 0.16`
  - `patchworklib` `>= 0.6.3, < 0.7`
  - `colormaps` `>= 0.4.2, < 0.5`
</details>

<details>
<summary><b>3. 特定分析流程与可选依赖</b> (点击展开)</summary>

#### 甲基化比对分析流程 (`methalign`)
- `ont-modkit` `>= 0.4.3, < 0.5`
- `bammap2` `>= 0.1.7, < 0.2`
- `pbmm2` `>= 1.16.99, < 2` *(仅限 linux-64 平台)*
- `pb-cpg-tools` `>= 3.0.0, < 4` *(仅限 linux-64 平台)*

#### 基因组评估与共线性比对流程 (`eval` & `eval2`)
- **Python**: `3.8.*` *(该评估子工具中指定的环境)*
- `syri` `>= 1.6.3, < 2`
- `plotsr` `>= 1.1.1, < 2`
- `nucflag` `== 1.0.0a2`
</details>

---

## 验证安装与运行示例 (Verifying Installation)

为了验证 CPhasing 及其所需的所有外部三方软件是否正确配置，我们提供了一套预包装的轻量化示例数据集和一键式自动测试脚本。详情亦可参阅完整的 [examples/README.md](../../examples/README.md) 说明文档。

### 第一步：下载测试数据集 (Toy Dataset)
导航到 CPhasing 的 `examples`（示例代码）目录，进行下载和解压：
```bash
# 如果是用推荐的一键式 tar.gz 方式安装的：
cd CPhasing*/examples

# 如果您是通过 Git 方式拉取的代码：
cd CPhasing/examples
```

下载并提取小体积（Toy）训练集：
```bash
# 下载轻量玩具测试数据集
wget https://github.com/wangyibin/CPhasing/releases/download/v0.3.2/example_data.tar.gz

# 提取并解压（这会自动生成 data/ 文件夹）
tar -xzvf example_data.tar.gz
```

### 第二步：激活 CPhasing 环境
在使用测试组件之前，确保你的应用隔离环境已被正确激活：
* **如果是 Pixi 一键激活（推荐手法）:**
  ```bash
  source ../bin/activate_cphasing
  ```
* **如果是 Conda 环境:**
  ```bash
  conda activate cphasing
  ```

### 第三步：一键运行示例脚本
在 `examples/` 目录下执行自动验证脚本：
```bash
bash run_example.sh
```

### 预期结果与运行耗时
脚本会将 **Pore-C pipeline** 和 **Hi-C pipeline** 流程全部运行并验证完毕：
* **Pore-C 输出文件夹:** `examples/cphasing_output_porec/`
* **Hi-C 输出文件夹:** `examples/cphasing_output_hic/`
* **预期总运行时间:** 在主流 8 核处理器工作站上耗时约 **5 分钟**。

在无报错的情况下干净利落地运行完毕，代表您的 `CPhasing` 状态已被完全正确部署，即可放心用于您大型的真实基因组数据集上。