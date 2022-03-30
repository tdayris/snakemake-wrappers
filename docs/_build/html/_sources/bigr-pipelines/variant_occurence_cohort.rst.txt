.. _`variant_occurence_cohort (Under development)`:

VARIANT_OCCURENCE_COHORT (UNDER DEVELOPMENT)
============================================

Compute variant occurences over a given cohort

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your working directory

  cd /path/to/my/working/directory

  # Build a design file (see below)

  # Copy/paste the following line

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_occurence_cohort/run.sh


Input/Output
------------


**Input:**

 
  
* (Annotated) VCF files (PBGZIP + TBI)
  
 


**Output:**

 
  
* Annotated VCF files (raw vcf)
  
 






Used wrappers
-------------

The following individual wrappers are used in this pipeline:


* :ref:`bio/snpsift/rminfo`

* :ref:`bio/variantoccurence/chromosomes`

* :ref:`bio/variantoccurence/annotate`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.




Notes
-----

Prerequisites:

* A TSV formatted design file, *named 'design.tsv'* with the following columns:

.. list-table:: Desgin file format
  :widths: 33 33 33
  :header-rows: 1

  * - Sample_id
    - Upstream_file
  * - Name of the Sample1
    - Path to upstream vcf file
  * - Name of the Sample2
    - Path to upstream vcf file
  * - ...
    - ...





Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    import datetime
    import logging
    import os
    import pandas
    import sys
    from pathlib import Path

    worflow_source_dir = Path(next(iter(workflow.get_sources()))).absolute().parent
    common = str(worflow_source_dir / "../common/python")
    sys.path.append(common)

    from file_manager import *
    from files_linker import *
    from write_yaml import *
    from messages import *

    from snakemake.utils import min_version
    min_version("6.0")

    logging.basicConfig(
        filename="snakemake.variant_calling_somatic.log",
        filemode="w",
        level=logging.DEBUG
    )

    container: "docker://continuumio/miniconda3:4.4.10"
    localrules: bigr_copy

    design = get_design(os.getcwd(), search_vcf)


    rule target:
       input:
          vcf = expand()



    """
    Cleaning annotation field in order to avoid double
    VarOcc field in final VCF.
    """
    rule clean_var_occ_from_vcf:
       input:
          call="data_input/calls/{sample}.vcf.gz"
       output:
          call=temp("snpsift/rminfo/{sample}.vcf")
       threads: 2
       resources:
          mem_mb=lambda wildcards, attempt: attempt * 1024 * 5,
          time_min=lambda wildcards, attempt: attempt * 35,
          tmpdir="tmp"
       log:
          "snpsfit/rminfo/{sample}.log"
       params:
          extra = "VarOcc"
       wrapper:
          "bio/snpsift/rminfo"


    """
    We use the input data to compute occurence since its annotation does not impact results
    """
    rule variant_occurence_per_chr:
       input:
          calls=expand(
             "data_input/calls/{sample}.vcf",
             sample=design["Sample_id"]
          )
       output:
          txt=temp("bigr/occurence/{chr}.txt")
       threads: 7
       resources:
          mem_mb=lambda wildcards, attempt: attempt * 1024,
          time_min=lambda wildcards, attempt: attempt * 45,
          tmpdir="tmp"
       log:
          "logs/bigr/variant_occurence/{chr}.log"
       wrapper:
          "bio/variantoccurence/chromosomes"


    rule variant_occurence_annotate:
        input:
            calls = ["snpsift/rminfo/{sample}.vcf"],
            occurence = "bigr/occurences/all_chroms.txt"
        output:
            calls = [temp("bigr/occurence_annotated/{sample}.vcf")]
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/variant_occurence/uncompress/{sample}.log"
        wrapper:
            "bio/variantoccurence/annotate"


    rule concatenate_per_chr_information:
        input:
            expand(
                "bigr/occurence/{chr}.txt",
                chr=config["params"]["chr"]
            )
        output:
            temp("bigr/occurences/all_chroms.txt")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/variant_occurence/all.log"
        shell:
            "for i in {input}; do sed '1d' ${{i}}; done > {output} 2> {log}"




Authors
-------


* Thibault Dayris

* Marc Deloger
