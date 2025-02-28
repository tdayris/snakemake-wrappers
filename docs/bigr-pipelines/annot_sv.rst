.. _`annot_sv`:

ANNOT_SV
========

Annotate Facets VCF with AnnotSV

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your working directory

  cd /path/to/my/working/directory

  # Build a design file (see below)

  # Copy/paste the following line for **HG38**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/annotsv/run.sh


Input/Output
------------


**Input:**

 
  
* VCF files
  
 


**Output:**

 
  
* Annotated TSV files with Census data
  
 







Notes
-----

Prerequisites:

* A TSV formatted design file, *named 'design.tsv'* with the following columns:

.. list-table:: Desgin file format
  :widths: 33 33 33
  :header-rows: 1

  * - Sample_id
    - VCF_file
  * - Name of the Sample1
    - Path to Facets VCF file
  * - ...
    - ...





Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    import logging
    import os
    import pandas
    import sys
    from pathlib import Path

    logging.basicConfig(
        filename="snakemake.snpeff_snpsift.log",
        filemode="w",
        level=logging.DEBUG
    )

    worflow_source_dir = Path(snakemake.workflow.srcdir("."))
    common = str(worflow_source_dir / ".." / "common" / "python")
    sys.path.append(common)

    from file_manager import *
    from files_linker import *
    from write_yaml import *
    from messages import *
    from snakemake.utils import min_version
    min_version("6.0")

    default_config = read_yaml(worflow_source_dir / "config.hg38.yaml")
    configfile: get_config(default_config)
    design = get_design(os.getcwd(), search_vcf_files)
    design["Sample_id"] = design["Sample_id"].str.replace("-", "_")

    container: "docker://continuumio/miniconda3:4.4.10"
    localrules: bigr_copy

    ruleorder: bigr_copy_idx > index_vcf

    samples_list = design["Sample_id"]

    wildcard_constraints:
        sample = r"|".join(samples_list)

    vcf_links = link_vcf(
        design.Sample_id,
        design.Upstream_file
    )


    ##################
    ### Flag rules ###
    ##################

    onsuccess:
        shell("touch DONE && rm --force --verbose ON_GOING ERROR")

    onerror:
        shell("touch ERROR && rm --force --verbose ON_GOING DONE")

    onstart:
        shell("touch ON_GOING && rm --force --verbose ERROR DONE")

    rule target:
        input:
            expand("annot_sv/{sample}.census.tsv", sample=samples_list)

    #################################################
    ### Gather files from iRODS or mounting point ###
    #################################################

    rule bigr_copy_vcf:
        output:
            "data_input/calls/{sample}.vcf.gz"
        message:
            "Getting {wildcards.sample} VCF file"
        threads: 1
        resources:
          mem_mb=lambda wildcards, attempt: min(attempt * 1024, 2048),
          time_min=lambda wildcards, attempt: attempt * 45,
        params:
            input=lambda wildcards, output: vcf_links[output[0]]
        group:
            "get_data"
        log:
            "logs/bigr_copy/{sample}.log"
        wrapper:
            "bio/BiGR/copy"


    rule bigr_copy_idx:
        output:
            "data_input/calls/{sample}.vcf.gz.tbi"
        message:
            "Getting {wildcards.sample} TBI file"
        threads: 1
        resources:
          mem_mb=lambda wildcards, attempt: min(attempt * 1024, 2048),
          time_min=lambda wildcards, attempt: attempt * 45,
        params:
            input=lambda wildcards, output: vcf_links[output[0][:-4]] + ".tbi"
        group:
            "get_data"
        log:
            "logs/bigr_copy/{sample}.log"
        wrapper:
            "bio/BiGR/copy"


    rule index_vcf:
        input:
            "data_input/calls/{sample}.vcf.gz"
        output:
            "data_input/calls/{sample}.vcf.gz.tbi"
        message: "Indexing facet vcf for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=512,
            time_min=10,
            tmpdir="tmp"
        group:
            "get_data"
        log:
            "logs/tabix/{sample}.log"
        params:
            "-p vcf"
        wrapper:
            "bio/tabix"


    ###################
    ### Run AnnotSV ###
    ###################

    rule annot_sv:
        input:
            vcf="data_input/calls/{sample}.vcf.gz",
            vcf_idx="data_input/calls/{sample}.vcf.gz.tbi"
        output:
            vcf=temp("annot_sv/{sample}.tsv")
        message: "Annotating Facet VCF for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
            time_min=lambda wildcards, attempt: attempt * 30,
            tmpdir="tmp"
        log:
            "logs/annot_sv/raw/{sample}.log"
        params:
            # install_readme=workflow.source_path("source/AnnotSV/README.md"),
            install_dir="/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/annot_sv/source/AnnotSV/",
            extra="-annotationMode both -snvIndelPASS 1 -tx ENSEMBL -SVinputInfo 0",
            genome=f"-genomeBuild {config['ref'].get('genome', 'GRCh38')}",

        conda:
            "envs/annot_sv.yaml"
        shell:
            """
    export ANNOTSV="{params.install_dir}" &&
    "{params.install_dir}/bin/AnnotSV" {params.genome} -SVinputFile {input.vcf} -outputFile {output.vcf} > {log} 2>&1
            """


    rule add_census:
        input:
            annot_sv="annot_sv/{sample}.tsv",
            census=config["ref"]["cancer_census"]
        output:
            census="annot_sv/{sample}.census.tsv"
        message:
            "Adding census data to {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 3,
            time_min=lambda wildcards, attempt: attempt * 5,
            tmpdir="tmp"
        log:
            "logs/annot_sv/census/{sample}.log"
        script:
            "scripts/add_census.py"




Authors
-------


* Thibault Dayris
