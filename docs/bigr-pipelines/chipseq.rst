.. _`Chipseq`:

CHIPSEQ
=======

Perform various analyses on Chipseq

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your directory containing VCF files

  cd /path/to/project

  # Copy/paste the following line for **HG38**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/chipseq/run.sh


Input/Output
------------


**Input:**

 
  
* Fastq files
  
 


**Output:**

 
  
* Peak called
  
 









Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    include: "rules/000.common.smk"
    include: "rules/001.bigr_copy.smk"
    include: "rules/002.trimming.smk"
    include: "rules/003.bowtie2.smk"
    include: "rules/004.sambamba.smk"
    include: "rules/005.samtools.smk"
    include: "rules/006.bedtools_bam2bed.smk"
    include: "rules/008.genomecov.smk"
    include: "rules/010.macs2.smk"
    include: "rules/011.deeptools.smk"
    include: "rules/012.mapping_qc.smk"
    include: "rules/013.multiqc.smk"
    include: "rules/009.seacr.smk"


    rule target:
        input:
            expand("samtools/view/{sample}.bam", sample=sample_list),
            expand("samtools/view/{sample}.bam.bai", sample=sample_list),
            expand(
                "macs2/callpeak/{peaktype}/{sample}_peaks.{peaktype}.bed",
                peaktype=peak_types,
                sample=sample_list,
            ),
            expand("deeptools/bamcoverage/{sample}.bw", sample=sample_list),
            "multiqc/Chipseq.html",
            expand("seacr/{sample}.{mode}.bed", sample=sample_list, mode=mode_list)




Authors
-------


* Thibault Dayris
