# CPhasing Examples

This directory contains a minimal Example dataset and an automated test script to verify that `CPhasing` is correctly installed and all its dependencies are working as expected.

## Directory Structure
```text
examples/
├── run_example.sh          # Automated testing script
├── README.md            # This instruction file
└── data/                # Lightweight Example dataset (minimal contigs & reads)
    ├── chroms.fasta.gz
    ├── contigs.fasta.gz
    ├── porec_reads.fasta.gz
    ├── hic_R1.fasta.gz
    └── hic_R2.fasta.gz
```

---

## How to Run the Example

### Step 1: Activate CPhasing Environment
Make sure you have activated the environment before running the test script.

* **If you installed via the recommended method (Pixi):**
  ```bash
  source /path/to/CPhasing/bin/activate_cphasing
  ```

* **If you installed via Conda:**
  ```bash
  conda activate cphasing
  ```

### Step 2: Download the Example Dataset
Download the lightweight test dataset from the GitHub release `v0.3.2` and extract it into the `examples/` directory:
```bash
# Navigate to the examples directory
cd examples

# Download the Example dataset archive
wget https://github.com/wangyibin/CPhasing/releases/download/v0.3.2/example_data.tar.gz

# Extract the dataset (which creates the data/ folder)
tar -xzvf example_data.tar.gz
```


### Step 3: Run the Automated Test Script
Execute the test script in the `examples/` folder:
```bash
bash run_example.sh
```
---

## What Happens During the Test?

The `run_example.sh` script automatically executes two representative C-Phasing workflows.

1. **Pore-C Workflow**

   Runs the long-read chromatin conformation pipeline (`cphasing pipeline`) using the provided example Pore-C reads.

   **Output:** Results are written to `examples/cphasing_output_porec/`, including `4.scaffolding/groups.agp`, `4.scaffolding/groups.asm.fasta`, `5.plot/groups.q1.100k.wg.png`, and log files.

2. **Hi-C Workflow**

   Runs the short-read chromatin conformation pipeline (`cphasing pipeline`) using the provided paired-end Hi-C example reads.

   **Output:** Results are written to `examples/cphasing_output_hic/`, including `4.scaffolding/groups.agp`, `4.scaffolding/groups.asm.fasta`, `5.plot/groups.q1.100k.wg.png`, and log files.

### Expected Runtime
Typically completes within 5 minutes on an 8-core workstation using the example dataset.
Successful completion without errors indicates that C-Phasing and its required dependencies have been installed correctly.