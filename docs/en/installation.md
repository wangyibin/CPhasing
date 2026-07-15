---
icon: lucide/hard-drive-download
tags: 
    - Setup
---



# Installation

## System Requirements

CPhasing simplifies dependency management using modern package managers. All required software tools and libraries can be resolved automatically.

### 1. Operating Systems & Architecture
- **Operating Systems**: Native **Linux** (e.g., Ubuntu >= 18.04, CentOS >= 7, Debian >= 10). Windows users can run it via WSL2. MacOS users can run it via Docker.  
- **Architectures**: Supported on `linux-64` and `linux-aarch64`.  
- **System Libraries**:  
  - Linux kernel `>= 3.10.0`  
  - GNU C Library (glibc) `libc >= 2.17`  

### 2. Environment Managers (Required)
To install CPhasing, you should have one of the following environment managers ready:
- **Pixi (Recommended)**: A fast, multi-platform package manager. Pixi automatically downloads and isolates all external bioinformatics tools and Python packages inside the work environment.  
- **Conda / Miniconda**: Supported as a traditional alternative for managing environments.  

!!! note "Automatic Dependency Resolution"
    You do **not** need to manually install individual third-party tools (such as `minimap2`, `samtools`, `bedtools`) or PyPI libraries. All of them are pre-configured and will be resolved automatically when you install via **Pixi** or **Conda**.

---


=== "recommended" 
    Download the CPhasing and activate environments by the `activate_cphasing`
    ```shell
    LATEST_URL=$(curl -s https://api.github.com/repos/wangyibin/CPhasing/releases/latest | grep "browser_download_url.*linux-64.tar.gz" | cut -d '"' -f 4)
    wget $LATEST_URL
    
    tar xzvf CPhasing*.tar.gz
    
    ### activate environment
    source ./CPhasing*/bin/activate_cphasing

    ### deactivate
    source deactivate 
    ```
    !!! note
        For the first configuration, run it when the network is accessible.
    
    !!! note
        If you can not download the `pixi`, you can download it from the github release page:
        ```bash
        mkdir -p ~/.pixi/bin
        cd ~/.pixi/bin
        wget https://github.com/prefix-dev/pixi/releases/download/v0.60.0/pixi-x86_64-unknown-linux-musl.tar.gz
        tar xzvf pixi-x86_64-unknown-linux-musl.tar.gz
        ```

        And rerun the `activate_cphasing` to install the dependencies of `CPhasing`

    !!! note
        If you do not have direct access to the anaconda repository, you can set the mirror for pixi.
    !!! note
        For the platform of **`linux-aarch64`**, please download from [github release](https://github.com/wangyibin/CPhasing/releases).
    
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
    

## Detailed Software & Library Dependencies

<details>
<summary><b>1. Core Bioinformatics Tools</b> (Click to expand)</summary>

```plaintext
# Alignment & Mapping
minimap2 >= 2.28, < 3
chromap >= 0.3.2, < 0.4
bwa-mem2 >= 2.3, < 3
minibwa >= 0.1, < 0.4
minigraph >= 0.21, < 0.22
wfmash 0.17.0.*

# Sequence Processing & Manipulation
samtools >= 1.20, < 1.21
bedtools >= 2.31.1, < 3
seqkit >= 2.9.0, < 3
samblaster >= 0.1.26, < 0.2
pigz >= 2.8, < 3
crabz >= 0.10.0, < 0.11
```

</details>

<details>
<summary><b>2. Python Runtime & Core Libraries</b> (Click to expand)</summary>

```plaintext
# Python Runtime
python 3.12.0.*

# Bioinformatics Libraries
biopython >= 1.84, < 2
pysam >= 0.22.1, < 0.23
cooler >= 0.10.2, < 0.11
hicmatrix >= 17.2, < 18
pyranges >= 0.1.2, < 0.2
ncls >= 0.0.68, < 0.0.69
needletail >= 0.7.1, < 0.8

# Data Manipulation & Performance
pandas >= 2.2.3, < 3
numpy >= 1.26.4, < 2
polars >= 1.17.1, < 1.18.0
pyarrow >= 18.1.0, < 19
dask >= 2024.11.2, < 2025
joblib >= 1.4.2, < 2
pandarallel >= 1.6.5, < 2
scikit-learn >= 1.5.2, < 2
sparse_dot_mkl >= 0.9.10, < 0.10  (linux-64 only)

# Graph & Networks
networkx >= 3.4.2, < 3.5
python-igraph >= 0.11.8, < 0.12
cdlib >= 0.4.0, < 0.5
graph-tool >= 2.97, < 3

# Data Visualization
matplotlib >= 3.9.3, < 4
seaborn >= 0.13.2, < 0.14
plotly >= 6.2.0, < 7
plotnine >= 0.15.3, < 0.16
patchworklib >= 0.6.3, < 0.7
colormaps >= 0.4.2, < 0.5
```

</details>

<details>
<summary><b>3. Optional & Environment-Specific Dependencies</b> (Click to expand)</summary>

```plaintext
# Methylation Alignment Flow (methalign)
ont-modkit >= 0.4.3, < 0.5
bammap2 >= 0.1.7, < 0.2
pbmm2 >= 1.16.99, < 2  (linux-64 only)
pb-cpg-tools >= 3.0.0, < 4  (linux-64 only)

# Genome Evaluation & Assembly Comparison (eval & eval2)
python 3.8.*  (specific runtimes inside evaluation workflows)
syri >= 1.6.3, < 2
plotsr >= 1.1.1, < 2
nucflag == 1.0.0a2
```

</details>


# Verifying Installation (Example Run)

To verify that CPhasing and all its required external dependencies are correctly configured, we provide a pre-packaged lightweight example dataset and an automated verification script. For more info, see the full [examples/README.md](../../examples/README.md) guide.

### Step 1: Download the Example Dataset
Navigate to the `examples` directory of your CPhasing installation:
```bash
# If installed via the recommended precompiled tar.gz:
cd CPhasing*/examples

# Or if you cloned the repository via Git:
cd CPhasing/examples
```


Download the test dataset and extract it:
```bash
# Download the lightweight example package
wget https://github.com/wangyibin/CPhasing/releases/download/v0.3.2/example_data.tar.gz

# Extract the archive (this will create the data/ folder)
tar -xzvf example_data.tar.gz
```

### Step 2: Activate Your CPhasing Environment
Ensure your installation environment is active:
* **Via Pixi (Recommended):**
  ```bash
  source ../bin/activate_cphasing
  ```
* **Via Conda:**
  ```bash
  conda activate cphasing
  ```

### Step 3: Run the Example Script
Execute the test automation script in the `examples/` directory:
```bash
bash run_example.sh
```

### Expected Results & Runtime
The script will run both the **Pore-C pipeline** and the **Hi-C pipeline** validation.
* **Pore-C Output Directory:** `examples/cphasing_output_porec/`
* **Hi-C Output Directory:** `examples/cphasing_output_hic/`
* **Expected Runtime:** ~5 minutes on an 8-core CPU workstation.

Successful runs that exit without errors indicate your installation and runtime dependencies are ready for real-world datasets.