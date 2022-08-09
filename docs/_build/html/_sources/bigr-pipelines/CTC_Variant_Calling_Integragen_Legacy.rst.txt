.. _`CTC_Variant_Calling_Integragen_Legacy`:

CTC_VARIANT_CALLING_INTEGRAGEN_LEGACY
=====================================

Perform analyses alike mat&met available in Integragen reports

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  Do not run


Input/Output
------------


**Input:**

 
  
* BAM files
  
 


**Output:**

 
  
* Called variants
  
 









Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    include: "rules/001.picard.samtofq.smk"
    include: "rules/002.cutadapt.smk"
    include: "rules/003.bwa.mem.smk"
    include: "rules/004.sambamba.smk"
    include: "rules/005.gatk.germline.smk"
    include: "rules/006.grep.smk"
    include: "rules/007.gatk.baseline.smk"
    include: "rules/008.grep.baseline.smk"
    include: "rules/009.gatk.somatic.smk"
    include: "rules/010.bam_readcount.smk"
    include: "rules/011.bcr2vep.smk"

    import pandas

    configfile: config.yaml
    design = pandas.read_csv(config["design"], sep="\t", header=0)

    rule target:
        input:
            expand(
                "gatk/mutect2/{sample}.vcf.gz",
                sample=design["Sample_id"]
            ),





Authors
-------


* Thibault Dayris
