.. C-Phasing documentation master file, created by
   sphinx-quickstart on Tue May 16 14:44:08 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Overview

   self


Getting started
****************


**C-Phasing**: a new pipeline for scaffolding and phasing polyploid genome by Pore-C or Hi-C data.

The major problem of scaffolding polyploid genome by Hi-C is that the lower unique mapping. 
Now, the long reads based chromosome conformation capture technology, called Pore-C, provide a fine way to solve this problem. 
So, we developed a new pipeline, called **C-Phasing**, specificially tailored to the polyploid phasing and scaffolding by Pore-C data.
Also it can be used to scaffold by Hi-C data, but will be slowly. 

Installation
============
Requirements
------------
- Python 3.7+
- Scientific Python packages

Install using conda
-------------------

Compile and install `C-Phasing` and its Python dependencies

.. code-block:: bash 
   
   git clone https://github.com/wangyibin/CPhasing.git
   cd CPhasing
   conda env create -f environment.yml

   export PATH=/path/to/CPhasing/bin:$PATH
   export PYTHONPATH=/path/to/CPhasing:$PYTHONPATH
   
   source activate cphasing


.. toctree::
   :maxdepth: 1
   :caption: Tutorials
   
   pore-c
   hic
   manual_adjust 


.. toctree::
   :maxdepth: 1
   :caption: Reference

   releases