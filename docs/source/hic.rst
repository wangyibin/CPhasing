.. toctree::
    :maxdepth: 2
    :hidden:


Scaffolding use Hi-C data
==========================


Step by step 
-------------

1. mapping Hi-C reads to draft assembly

.. code-block:: bash

    cphasing hic mapper -r draft.asm.fasta -1 Lib_R1.fastq.gz -2 Lib_R2.fastq.gz -t 10 

2. build Hi-C contact matrix

.. code-block:: bash 

    cphasing-rs chromsizes draft.asm.fasta > draft.asm.contigsizes
    cphasing prepare pairs2cool Lib.pairs.gz draft.asm.contigsizes Lib.10000.cool 

3. extract Hi-C contact matrix 

.. code-block:: bash 

    cphasing extract Lib.pairs.gz draft.asm.contigsizes Lib.edges --pairs

4. partition

.. code-block:: bash 

    cphasing hyperpartition Lib.edges draft.asm.contigsizes out.clusters.txt 


    