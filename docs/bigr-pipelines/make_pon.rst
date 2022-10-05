.. _`make_pon`:

MAKE_PON
========

Build Panel of Normal over a list of Fastq files

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  cd /path/to/project

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/make_pon/run.sh


Input/Output
------------


**Input:**

 
  
* Fastq files
  
 


**Output:**

 
  
* PoN
  
 









Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    include: "rules/000.commons.smk"
    include: "rules/001.index_reference.smk"
    include: "rules/002.copy_trim.smk"
    include: "rules/003.bwa.smk"
    include: "rules/004.correct_mappings.smk"
    include: "rules/005.mutect2.smk"
    include: "rules/006.PoN.smk"
    include: "rules/007.qc.smk"


    onstart:
        shell("rm --force --verbose ERROR DONE && touch ON_GOING")


    onerror:
        shell("rm --force --verbose ON_GOING DONE && touch ERROR")


    onsuccess:
        shell("rm --force --verbose ERROR ON_GOING && touch DONE")


    rule target:
        input:
            "data_output/PoN.vcf.gz",
            "data_output/PoN.gdb",




Authors
-------


* Thibault Dayris
