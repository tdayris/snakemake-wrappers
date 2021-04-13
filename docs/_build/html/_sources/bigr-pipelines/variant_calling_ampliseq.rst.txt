.. _`Variant_Calling_Ampliseq`:

VARIANT_CALLING_AMPLISEQ
========================

Perform Variant calling on Ampliseq

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Copy/paste the following line for **HG19**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/snpeff_snpsift/run.sh hg19


Input/Output
------------


**Input:**

 
  
* VCF files
  
 


**Output:**

 
  
* Annotated VCF file
  
 





Used meta-wrappers
------------------

The following individual meta-wrappers are used in this pipeline:


* :ref:`master/bio/bwa_fixmate`

* :ref:`master/bio/gatk_bqsr`

* :ref:`master/bio/varscan2_calling`


Please refer to each meta-wrapper in above list for additional configuration parameters and information about the executed code.




Used wrappers
-------------

The following individual wrappers are used in this pipeline:


* :ref:`master/bio/bigr/copy`

* :ref:`master/bio/fastp`

* :ref:`master/bio/snpeff/annotate`

* :ref:`master/bio/snpsift/varType`

* :ref:`master/bio/snpsift/genesets`

* :ref:`master/bio/snpsift/annotate`

* :ref:`master/bio/snpsift/gwascat`

* :ref:`master/bio/compress/pbgzip`

* :ref:`master/bio/tabix`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.




Notes
-----

The only difference with a classic WES pipeline is the absence of duplicates removal.

Prerequisites:

* A TSV formatted design file with the following columns:

..list-table:: Desgin file format
  :widths: 33 33 33
  :header-rows: 1

  * - Sample_id
    - Upstream_fastq
    - Downstream_fastq
  * - Name of the Sample1
    - Path to upstream fastq file
    - Path to downstream fastq file




Authors
-------


* Thibault Dayris
