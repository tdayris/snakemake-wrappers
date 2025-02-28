.. _`WES_Germline`:

WES_GERMLINE
============

Perform Variant calling on WES (germline)

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your working directory

  cd /path/to/my/working/directory

  # Build a design file (see below)

  # Copy/paste the following line for **HG38**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/WES_germline/run.sh

  # Copy/paste the following line for **HG19**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/WES_germline/run.sh hg19


Input/Output
------------


**Input:**

 
  
* Fastq files
  
 


**Output:**

 
  
* Annotated VCF file
  
 





Used meta-wrappers
------------------

The following individual meta-wrappers are used in this pipeline:


* :ref:`meta/bio/bwa_fixmate`

* :ref:`meta/bio/gatk_bqsr`

* :ref:`meta/bio/varscan2_germline`

* :ref:`meta/bio/mutect2_germline`

* :ref:`meta/bio/meta_caller_germline`

* :ref:`bigr_pipelines/snpeff_snpsift`

* :ref:`bigr_pipelines/fastqc_multiqc`


Please refer to each meta-wrapper in above list for additional configuration parameters and information about the executed code.




Used wrappers
-------------

The following individual wrappers are used in this pipeline:


* :ref:`bio/picard/markduplicates`

* :ref:`bio/fastp`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.




Notes
-----

The only difference with a classic WES pipeline is the absence of duplicates removal.

Prerequisites:

* A TSV formatted design file, *named 'design.tsv'* with the following columns:

.. list-table:: Desgin file format
  :widths: 33 33 33
  :header-rows: 1

  * - Sample_id
    - Upstream_fastq
    - Downstream_fastq
  * - Name of the Sample1
    - Path to upstream fastq file
    - Path to downstream fastq file
  * - Name of the Sample2
    - Path to upstream fastq file
    - Path to downstream fastq file
  * - ...
    - ...
    - ...

If this design file is missing, Fastq pairs of files are assumed to be available within the working directory.




Authors
-------


* Thibault Dayris

* M boyba Diop

* Marc Deloger
