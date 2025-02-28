.. _`bs_ampliseq`:

BS_AMPLISEQ
===========

Analyse bisulfite ampliseq with Bismark

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  cd /path/to/project

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/bs_ampliseq/run.sh


Input/Output
------------


**Input:**

 
  
* Fastq files
  
 


**Output:**

 
  
* Bismark mapping + Methylation
  
 









Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    include: "rules/000.common.smk"
    include: "rules/001.bigr_copy.smk"
    include: "rules/002.trimming.smk"
    include: "rules/003.bismark.mapping.smk"
    include: "rules/004.bismark.extract.smk"
    include: "rules/005.bismark.report.smk"
    include: "rules/006.multiqc.smk"


    rule target:
        input:
            multiqc="data_output/MultiQC/Bismark.html",
            mbias_r1=expand(
                "data_output/Bismark/{sample}.M-bias_R1.png", sample=sample_list
            ),
            mbias_r2=expand(
                "data_output/Bismark/{sample}.M-bias_R2.png", sample=sample_list
            ),
            mbias_report=expand("bismark/meth/{sample}.M-bias.txt", sample=sample_list),
            splitting_report=expand(
                "bismark/meth/{sample}_splitting_report.txt", sample=sample_list
            ),
            methylome_CpG_cov=expand(
                "data_output/Bismark/{sample}.cpg.cov.gz", sample=sample_list
            ),
            methylome_CpG_mlevel_bedGraph=expand(
                "data_output/Bismark/{sample}.cpg.bedGraph.gz",
                sample=sample_list,
            ),
            read_base_meth_state_cpg=expand(
                "data_output/Bismark/CpG_context_{sample}.txt.gz", sample=sample_list
            ),
            read_base_meth_state_chg=expand(
                "data_output/Bismark/CHG_context_{sample}.txt.gz", sample=sample_list
            ),
            read_base_meth_state_chh=expand(
                "data_output/Bismark/CHH_context_{sample}.txt.gz", sample=sample_list
            ),




Authors
-------


* Thibault Dayris
