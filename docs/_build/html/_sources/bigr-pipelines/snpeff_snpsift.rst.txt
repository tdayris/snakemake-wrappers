.. _`SnpEff_SnpSift`:

SNPEFF_SNPSIFT
==============

Annotate VCF with SnpEff and SNpSift

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your directory containing VCF files

  cd /path/to/CEL/dir

  # Copy/paster the following line

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/snpeff_snpsift/run.sh


Input/Output
------------


**Input:**

 
  
* VCF files
  
 


**Output:**

 
  
* Annotated VCF file
  
 





Used meta-wrappers
------------------

The following individual meta-wrappers are used in this pipeline:


* :ref:`master/bio/snpeff/annotate`

* :ref:`master/bio/snpsift/varType`

* :ref:`master/bio/snpsift/genesets`

* :ref:`master/bio/snpsift/annotate`

* :ref:`master/bio/snpsift/gwascat`


Please refer to each meta-wrapper in above list for additional configuration parameters and information about the executed code.






Authors
-------


* Thibault Dayris
