.. _`retrotransposons`:

RETROTRANSPOSONS
================

Quantify and search differentially expressed retrotransposons from RNASeq

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 


Input/Output
------------


**Input:**

 
  
* Pairs of fastq files
  
 


**Output:**

 
  
* DESeq2 results
  
 







Notes
-----

This is a variation of the initial RNASeq differential expression pipeline





Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-

    from snakemake.utils import min_version
    from pathlib import Path
    from yaml import dump
    min_version("6.0")

    import sys

    worflow_source_dir = Path(next(iter(workflow.get_sources()))).absolute().parent
    common = str(worflow_source_dir / "../common/python")
    sys.path.append(common)

    from dataframes import *
    from file_manager import *
    from files_linker import *
    from graphics import *
    from write_yaml import *
    from messages import message


    #################
    ### Preambule ###
    #################

    logging.basicConfig(
        filename="snakemake.retrotransposons.log",
        filemode="w",
        level=logging.DEBUG
    )

    default_config = read_yaml(worflow_source_dir / "config.hg38.yaml")
    configfile: get_config(default_config)
    design = get_design(os.getcwd(), search_fastq_pairs)

    fastq_links = link_fq(
        design.Sample_id,
        design.Upstream_file,
        design.Downstream_file
    )

    sample_list = design.Sample_id.to_list()


    rule target:
        input:
            "homer/repeats.tsv",
            "multiqc/MultiQC.html"


    ########################
    ### Quality controls ###
    ########################


    rule multiqc:
        input:
            summary=expand(
                "picard/alignment_summary/{sample}.summary.txt",
                sample=sample_list
            ),
            html=expand("fastp/html/pe/{sample}.fastp.html",
                sample=sample_list
            ),
            json=expand("fastp/json/pe/{sample}.fastp.json",
                sample=sample_list
            )
        output:
            "multiqc/MultiQC.html"
        message:
            "Aggregating quality reports from Fastp and Salmon"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35,
            tmpdir="tmp"
        log:
            "logs/multiqc.log"
        wrapper:
            "bio/multiqc"


    #############
    ### Homer ###
    #############

    rule analyze_repeats:
        input:
            fasta = config["ref"]["fasta"],
            tag_directories = expand(
                "homer/tag_directories/{sample}",
                sample=sample_list
            )
        output:
            "homer/repeats.tsv"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        params:
            extra = config["params"].get(
                "analyzerepeats_extra", "-format sam"
            )
        log:
            "logs/homer/analyzerepeats.log"
        wrapper:
            "bio/homer/analyzeRepeats"


    rule make_tag_directory:
        input:
            bam="star/{sample}/Aligned.out.sam"
        output:
            directory("homer/tag_directories/{sample}")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        params:
            extra=config["params"]["homer_tag_dirs"]
        log:
            "logs/homer/makeTagDirectory/{sample}.log"
        wrapper:
            "bio/homer/makeTagDirectory"


    ###################################
    ### STAR mapping and correcting ###
    ###################################

    star_mapping_config = {
        "fasta": config["ref"]["fasta"],
        "gtf": config["ref"]["gtf"],
        "samtools_view_extra": config["params"]["samtools_view"],
        "star_mapping_extra": config["params"]["star_mapping"],
        "star_index_extra": config["params"]["star_index_extra"],
        "sjdbOverhang": config["params"]["sjdbOverhang"]
    }

    module star_mapping:
        snakefile: "../../meta/bio/star_mapping/test/Snakefile"
        config: star_mapping_config


    use rule * from star_mapping



    ###################
    ## Bam indexing ###
    ###################

    index_dataset_config = {}

    module index_datasets:
        snakefile: "../../meta/bio/index_datasets/test/Snakefile"
        config: index_dataset_config

    use rule samtools_index_bam from index_datasets



    ############################
    ### FASTP FASTQ CLEANING ###
    ############################

    rule fastp_clean:
        input:
            sample=expand(
                "reads/{sample}.{stream}.fq.gz",
                stream=["1", "2"],
                allow_missing=True
            ),
        output:
            trimmed=temp(expand(
                "fastp/trimmed/pe/{sample}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True
            )),
            html="fastp/html/pe/{sample}.fastp.html",
            json="fastp/json/pe/{sample}.fastp.json"
        message: "Cleaning {wildcards.sample} with Fastp"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 4096, 15360),
            time_min=lambda wildcard, attempt: attempt * 45,
            tmpdir="tmp"
        params:
            adapters=config["params"].get("fastp_adapters", None),
            extra=config["params"].get("fastp_extra", "")
        log:
            "logs/fastp/{sample}.log"
        wrapper:
            "bio/fastp"


    #################################################
    ### Gather files from iRODS or mounting point ###
    #################################################

    rule bigr_copy:
        output:
            "reads/{sample}.{stream}.fq.gz"
        message:
            "Gathering {wildcards.sample} fastq file ({wildcards.stream})"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
            time_min=lambda wildcard, attempt: attempt * 45
        params:
            input=lambda wildcards, output: fastq_links[output[0]]
        log:
            "logs/bigr_copy/{sample}.{stream}.log"
        wrapper:
            "bio/BiGR/copy"




Authors
-------


* Thibault Dayris
