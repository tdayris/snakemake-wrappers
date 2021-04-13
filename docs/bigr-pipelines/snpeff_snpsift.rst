.. _`SnpEff_SnpSift`:

SNPEFF_SNPSIFT
==============

Annotate VCF files with SnpEff and SNpSift

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your directory containing VCF files

  cd /path/to/VCF/dir

  # Copy/paste the following line for **HG38**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/snpeff_snpsift/run.sh

  # Copy/paste the following line for **HG19**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/snpeff_snpsift/run.sh hg19


Input/Output
------------


**Input:**

 
  
* VCF files
  
 


**Output:**

 
  
* Annotated VCF file
  
 






Used wrappers
-------------

The following individual wrappers are used in this pipeline:


* :ref:`master/bio/snpeff/annotate`

* :ref:`master/bio/snpsift/varType`

* :ref:`master/bio/snpsift/genesets`

* :ref:`master/bio/snpsift/annotate`

* :ref:`master/bio/snpsift/gwascat`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.





Authors
-------


* Thibault Dayris
