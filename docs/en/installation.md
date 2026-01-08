---
icon: lucide/hard-drive-download
tags: 
    - Setup
---

=== "recommended" 
    Download the CPhasing and activate environments by the `activate_cphasing`
    ```shell
    LATEST_URL=$(curl -s https://api.github.com/repos/wangyibin/CPhasing/releases/latest | grep "browser_download_url.*linux-x86.tar.gz" | cut -d '"' -f 4)
    wget $LATEST_URL
    
    tar xzvf CPhasing*.tar.gz
    
    source ./CPhasing*/bin/activate_cphasing
    
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
        For the platform of **`linux-aarch64`**, please download from [github release](https://github.com/wangyibin/CPhasing/releases/tag/v0.2.8).
    
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
    
