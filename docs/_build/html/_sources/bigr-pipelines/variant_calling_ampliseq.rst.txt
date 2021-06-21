.. _`Variant_Calling_Ampliseq (under development)`:

VARIANT_CALLING_AMPLISEQ (UNDER DEVELOPMENT)
============================================

Perform Variant calling on Ampliseq

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your working directory

  cd /path/to/my/working/directory

  # Build a design file (see below)

  # Copy/paste the following line for **HG19**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_ampliseq/run.sh

  # Copy/paste the following line for **HG38**

  {'\# WARNING': 'HG38 remains not tested (possible wrong capture kit)'}

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_ampliseq/run.sh hg38


Input/Output
------------


**Input:**

 
  
* VCF files
  
 


**Output:**

 
  
* Annotated VCF file
  
 





Used meta-wrappers
------------------

The following individual meta-wrappers are used in this pipeline:


* :ref:`meta/bio/bwa_fixmate`

* :ref:`meta/bio/gatk_bqsr`

* :ref:`meta/bio/varscan2_calling`

* :ref:`bigr_pipelines/snpeff_snpsift`


Please refer to each meta-wrapper in above list for additional configuration parameters and information about the executed code.




Used wrappers
-------------

The following individual wrappers are used in this pipeline:


* :ref:`bio/bigr/copy`

* :ref:`bio/fastp`

* :ref:`bio/compress/pbgzip`

* :ref:`bio/tabix`

* :ref:`bio/multiqc`

* :ref:`bio/picard/collectalignmentsummarymetrics`

* :ref:`bio/fastq_screen`


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





Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    import logging
    import os
    import pandas
    import sys

    sys.path.append("/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/python/")

    from file_manager import *
    from files_linker import link_fq
    from write_yaml import read_yaml
    from pathlib import Path
    from messages import CustomFormatter
    from snakemake.utils import min_version
    min_version("6.0")

    logging.basicConfig(
        filename="snakemake.variant_calling_ampliseq.log",
        filemode="w",
        level=logging.DEBUG
    )

    container: "docker://continuumio/miniconda3:4.4.10"
    localrules: bigr_copy
    ruleorder: bwa_mem > bwa_fixmate_meta_bwa_mem


    default_config = read_yaml("/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_ampliseq/config.hg38.yaml")
    configfile: get_config(default_config)
    design = get_design(os.getcwd(), search_fastq_pairs)


    wildcard_constraints:
        sample = r"|".join(design["Sample_id"]),
        stream = r"1|2|R1|R2"


    fastq_links = link_fq(
        design.Sample_id,
        design.Upstream_file,
        design.Downstream_file
    )

    rule all:
        input:
            calls=expand(
                "snpsift/dbnsfp/{sample}.vcf.gz{index}",
                sample=design["Sample_id"].tolist(),
                index=["", ".tbi"]
            ),
            html="multiqc/variant_calling_ampliseq.html"
        message:
            "Finishing the Ampliseq variant calling"


    #################
    ### Gather QC ###
    #################

    rule multiqc:
        input:
            html=expand(
                "fastp/html/pe/{sample}.fastp.html",
                sample=design["Sample_id"]
            ),
            json=expand(
                "fastp/json/pe/{sample}.fastp.json",
                sample=design["Sample_id"]
            ),
            picard=expand(
                "picard/alignment_summary/{sample}.summary.txt",
                sample=design["Sample_id"]
            ),
            fastq_screen=expand(
                "fastq_screen/{sample}.{stream}.fastq_screen.{ext}",
                sample=design["Sample_id"],
                stream=["1", "2"],
                ext=["txt", "png"]
            )
        output:
            report(
                "multiqc/variant_calling_ampliseq.html",
                caption="../common/reports/multiqc.rst",
                category="Quality Controls"
            )
        message:
            "Aggregating quality reports from SnpEff"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35
        log:
            "logs/multiqc.log"
        wrapper:
            "/bio/multiqc"


    rule alignment_summary:
        input:
            bam="samtools/sort/{sample}.bam",
            bam_index="samtools/sort/{sample}.bam.bai",
            ref=config['ref']['fasta'],
            ref_idx=get_fai(config['ref']['fasta']),
            ref_dict=get_dict(config['ref']['fasta']),
        output:
            temp("picard/alignment_summary/{sample}.summary.txt")
        message:
            "Collecting alignment metrics on GATK recalibrated {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1020,
            time_min=lambda wildcards, attempt: attempt * 45
        log:
            "logs/picard/alignment_summary/{sample}.log"
        params:
            "VALIDATION_STRINGENCY=LENIENT "
            "METRIC_ACCUMULATION_LEVEL=null "
            "METRIC_ACCUMULATION_LEVEL=SAMPLE"
        wrapper:
            "/bio/picard/collectalignmentsummarymetrics"


    rule fastq_screen:
        input:
            "reads/{sample}.{stream}.fq.gz"
        output:
            txt=temp("fastq_screen/{sample}.{stream}.fastq_screen.txt"),
            png=temp("fastq_screen/{sample}.{stream}.fastq_screen.png")
        message:
            "Assessing quality of {wildcards.sample}, {wildcards.stream}"
        threads: config.get("threads", 20)
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 4096, 8192),
            time_min=lambda wildcard, attempt: attempt * 50
        params:
            fastq_screen_config=config["fastq_screen"],
            subset=100000,
            aligner='bowtie2'
        log:
            "logs/fastqc/{sample}.{stream}.log"
        wrapper:
            "/bio/fastq_screen"


    #################################
    ### FINAL VCF FILE INDEXATION ###
    #################################

    module compress_index_vcf_meta:
        snakefile: "../../meta/bio/compress_index_vcf/test/Snakefile"
        config: config

    use rule * from compress_index_vcf_meta as compress_index_vcf_*

    ######################
    ### VCF annotation ###
    ######################


    module snpeff_meta:
        snakefile: "../../meta/bio/snpeff_annotate/test/Snakefile"
        config: config

    use rule snpeff from snpeff_meta with:
        input:
            calls="meta_caller/calls/{sample}.vcf.gz",
            calls_index=get_tbi("meta_caller/calls/{sample}.vcf.gz"),
            db=config["ref"]["snpeff"]


    module snpsift:
        snakefile: "../../meta/bio/snpsift/test/Snakefile"
        config: config

    use rule * from snpsift as snpsift_*



    #####################################
    ### Merge variant calling results ###
    #####################################

    module metacaller_germline_meta:
        snakefile: "../../meta/bio/meta_caller_germline/test/Snakefile"
        config: {"genome": config["ref"]["fasta"], "bed": config["ref"]["capture_kit_bed"]}


    use rule * from metacaller_germline_meta as *


    ############################################################################
    ### Correcting Mutect2 :                                                 ###
    ### AS_FilterStatus: Number=1 and not Number=A which violates VCF format ###
    ### AD becomes ADM: AD is reserved for Allele Depth, Mutect2 stores      ###
    ###                 multiple information under "AD" field.               ###
    ############################################################################

    rule correct_mutect2_vcf:
        input:
            "mutect2/filter_reheaded/{sample}.vcf.gz"
        output:
            temp("mutect2/corrected/{sample}.vcf")
        message:
            "Renaming reserved AD field and fixing AS_FilterStrand format error"
            " on {wildcards.sample}"
        threads: 3
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 256,
            time_min=lambda wildcards, attempt: attempt * 20
        log:
            "logs/mutect2/correct_fields/{sample}.log"
        params:
            rename_ad="'s/=AD;/=ADM;/g'",
            rename_ad_format="'s/:AD:/:ADM:/g'",
            fix_as_filterstatus="'s/ID=AS_FilterStatus,Number=A/ID=AS_FilterStatus,Number=1/g'"
        shell:
            "(gunzip -c {input} | "
            "sed {params.rename_ad} | "
            "sed {params.rename_ad_format} | "
            "sed {params.fix_as_filterstatus}) "
            "> {output} 2> {log}"

    ###############################
    ### Variant calling Mutect2 ###
    ###############################


    module gatk_mutect2_germline_meta:
        snakefile: "../../meta/bio/mutect2_germline/test/Snakefile"
        config: {"genome": config["ref"]["fasta"], "known": config["ref"]["af_only"], "bed": config["ref"]["capture_kit_bed"], "dbsnp": config["ref"]["dbsnp"]}

    use rule * from gatk_mutect2_germline_meta as gatk_mutect2_germline_*

    ################################
    ### Variant Calling Varscan2 ###
    ################################

    module varscan2_germline_meta:
        snakefile: "../../meta/bio/varscan2_germline/test/Snakefile"
        config: {"genome": config["ref"]["fasta"], "bed": config["ref"]["capture_kit_bed"]}

    use rule * from varscan2_germline_meta as *

    use rule samtools_mpilup from varscan2_germline_meta with:
        input:
            bam="gatk/recal_bam/{sample}.bam",
            reference_genome=config['ref']['fasta'],
            reference_genome_idx=get_fai(config['ref']['fasta']),


    ##############################
    ### GATK BAM RECALIBRATION ###
    ##############################

    module gatk_bqsr_meta:
        snakefile: "../../meta/bio/gatk_bqsr/test/Snakefile"
        config: {"threads": config["threads"], "genome": config["ref"]["fasta"], "dbsnp": config["ref"]["dbsnp"]}


    use rule gatk_apply_baserecalibrator from gatk_bqsr_meta with:
        input:
            bam="samtools/sort/{sample}.bam",
            bam_index="samtools/sort/{sample}.bam.bai",
            ref=config['ref']['fasta'],
            ref_idx=get_fai(config['ref']['fasta']),
            ref_dict=get_dict(config['ref']['fasta']),
            recal_table="gatk/recal_data_table/{sample}.grp"


    use rule gatk_compute_baserecalibration_table from gatk_bqsr_meta with:
        input:
            bam="samtools/sort/{sample}.bam",
            bam_index="samtools/sort/{sample}.bam.bai",
            ref=config['ref']['fasta'],
            ref_idx=get_fai(config['ref']['fasta']),
            ref_dict=get_dict(config['ref']['fasta']),
            known=config['ref']['dbsnp'],
            known_idx=get_tbi(config['ref']['dbsnp'])


    ###################
    ### BWA MAPPING ###
    ###################

    module bwa_fixmate_meta:
        snakefile: "../../meta/bio/bwa_fixmate/test/Snakefile"
        config: {"threads": config["threads"], "genome": config["ref"]["fasta"]}

    use rule * from bwa_fixmate_meta as bwa_fixmate_meta_*

    use rule bwa_mem from bwa_fixmate_meta with:
        input:
            reads=expand(
                "fastp/trimmed/pe/{sample}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True
            ),
            index=multiext(
                "bwa_mem2/index/genome", ".0123", ".amb", ".ann", ".pac"
            )


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
            trimmed=expand(
                "fastp/trimmed/pe/{sample}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True
            ),
            html="fastp/html/pe/{sample}.fastp.html",
            json=temp("fastp/json/pe/{sample}.fastp.json")
        message: "Cleaning {wildcards.sample} with Fastp"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 4096, 15360),
            time_min=lambda wildcard, attempt: attempt * 45
        params:
            adapters=config.get("fastp_adapters", None),
            extra=config.get("fastp_extra", "")
        log:
            "logs/fastp/{sample}.log"
        wrapper:
            "/bio/fastp"


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
            "/bio/BiGR/copy"




Authors
-------


* Thibault Dayris

* M boyba Diop

* Marc Deloger
