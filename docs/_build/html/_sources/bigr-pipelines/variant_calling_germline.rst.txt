.. _`Variant_Calling_Germline (under development)`:

VARIANT_CALLING_GERMLINE (UNDER DEVELOPMENT)
============================================

Perform Variant calling on Germline

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your working directory

  cd /path/to/my/working/directory

  # Build a design file (see below)

  # Copy/paste the following line for **HG19**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_germline/run.sh

  # Copy/paste the following line for **HG38**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_germline/run.sh hg38


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
    from pathlib import Path

    worflow_source_dir = Path(next(iter(workflow.get_sources()))).absolute().parent
    common = str(worflow_source_dir / "../common/python")
    sys.path.append(common)

    from file_manager import *
    from files_linker import *
    from write_yaml import *
    from messages import *
    from snakemake.utils import min_version
    from pathlib import Path
    min_version("6.0")

    logging.basicConfig(
        filename="snakemake.variant_calling_germline.log",
        filemode="w",
        level=logging.DEBUG
    )

    container: "docker://continuumio/miniconda3:4.4.10"
    localrules: bigr_copy


    default_config = read_yaml(worflow_source_dir / "config.hg38.yaml")
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
            html="multiqc/variant_calling_germline.html"
        message:
            "Finishing the WES Germline Variant Calling pipeline"


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
            ),
            picards_metrics=expand(
                "picard/markduplicates/metrics/{sample}.picard.metrics.txt",
                sample=design["Sample_id"]
            ),
            snpeff_stats=expand(
                "snpeff/report/{sample}.html",
                sample=design["Sample_id"]
            ),
            snpeff_csvstats=expand(
                "snpeff/csvstats/{sample}.csv",
                sample=design["Sample_id"]
            )
        output:
            report(
                "multiqc/variant_calling_germline.html",
                caption="../common/reports/multiqc.rst",
                category="Quality Controls"
            )
        message:
            "Aggregating quality reports from SnpEff"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35,
            tmpdir=lambda wildcards, input: f"tmp/multiqc.tmp"
        log:
            "logs/multiqc.log"
        wrapper:
            "bio/multiqc"


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
            time_min=lambda wildcards, attempt: attempt * 45,
            tmpdir=lambda wildcards, input: f"tmp/{wildcards.sample}.tmp"
        log:
            "logs/picard/alignment_summary/{sample}.log"
        params:
            "VALIDATION_STRINGENCY=LENIENT "
            "METRIC_ACCUMULATION_LEVEL=null "
            "METRIC_ACCUMULATION_LEVEL=SAMPLE"
        wrapper:
            "bio/picard/collectalignmentsummarymetrics"


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
            time_min=lambda wildcard, attempt: attempt * 50,
            tmpdir=lambda wildcards, input: "tmp/{}.{}.tmp".format(
              wildcards.sample, wildcards.stream
            )
        params:
            fastq_screen_config=config["fastq_screen"],
            subset=100000,
            aligner='bowtie2'
        log:
            "logs/fastqc/{sample}.{stream}.log"
        wrapper:
            "bio/fastq_screen"


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

    use rule * from snpsift


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
            time_min=lambda wildcards, attempt: attempt * 20,
            tmpdir=lambda wildcards, input: f"tmp/{wildcards.sample}.tmp"
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

    gatk_mutect2_germline_meta_config = {
        "genome": config["ref"]["fasta"],
        "known": config["ref"]["af_only"],
        "bed": config["ref"]["capture_kit_bed"],
        "dbsnp": config["ref"]["dbsnp"]
    }

    module gatk_mutect2_germline_meta:
        snakefile: "../../meta/bio/mutect2_germline/test/Snakefile"
        config: gatk_mutect2_germline_meta_config

    use rule * from gatk_mutect2_germline_meta


    use rule muterc2_filter from gatk_mutect2_germline_meta with:
        input:
            vcf="mutect2/call/{sample}.vcf.gz",
            ref=config["ref"]["fasta"],
            fasta_index=get_fai(config["ref"]["fasta"]),
            fasta_dict=get_dict(config["ref"]["fasta"]),
            contamination="summary/{sample}_calculate_contamination.table",
            bam="picard/markduplicates/mapping/{sample}.bam",
            bam_index=get_bai("picard/markduplicates/mapping/{sample}.bam"),
            f1r2="gatk/artifacts_prior/{sample}.artifacts_prior.tar.gz"


    use rule get_pileup_summaries from gatk_mutect2_germline_meta with:
        input:
            bam="picard/markduplicates/mapping/{sample}.bam",
            bam_index=get_bai("picard/markduplicates/mapping/{sample}.bam"),
            intervals=config["ref"]["capture_kit_bed"],
            variants=config["ref"]["af_only"],
            variants_index=get_tbi(config["ref"]["af_only"])


    use rule mutect2_germline from gatk_mutect2_germline_meta with:
        input:
            fasta=config["ref"]["fasta"],
            fasta_index=get_fai(config["ref"]["fasta"]),
            fasta_dict=get_dict(config["ref"]["fasta"]),
            map="picard/markduplicates/mapping/{sample}.bam",
            map_index=get_bai("picard/markduplicates/mapping/{sample}.bam"),
            germline=config["ref"]["af_only"],
            germline_tbi=get_tbi(config["ref"]["af_only"]),
            intervals=config["ref"]["capture_kit_bed"]

    ################################
    ### Variant Calling Varscan2 ###
    ################################

    varscan2_config = {
        "genome": config["ref"]["fasta"],
        "bed": config["ref"]["capture_kit_bed"]
    }

    module varscan2_meta:
        snakefile: "../../meta/bio/varscan2_germline/test/Snakefile"
        config: varscan2_config

    use rule * from varscan2_meta as varscan2_meta_*


    ###################
    ### BWA MAPPING ###
    ###################

    module bwa_meta:
        snakefile: "../../meta/bio/bwa_fixmate/test/Snakefile"
        config: {"threads": config["threads"], "genome": config["ref"]["fasta"]}

    use rule * from bwa_meta as *

    use rule bwa_mem from bwa_meta with:
        input:
            reads=expand(
                "fastp/trimmed/pe/{sample}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True
            ),
            index=multiext(
                "bwa_mem2/index/genome", ".0123", ".amb", ".ann", ".pac"
            )


    #####################
    ### Deduplicating ###
    #####################

    rule picard_markduplicates:
        input:
            bam="samtools/sort/{sample}.bam"
        output:
            bam=temp("picard/markduplicates/mapping/{sample}.bam"),
            metrics=temp("picard/markduplicates/metrics/{sample}.picard.metrics.txt")
        message:
            "Removing duplicates on {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 5120, 10240),
            time_min=lambda wildcards, attempt: attempt * 45,
            tmpdir=lambda wildcards, input: f"tmp/{wildcards.sample}.tmp"
        log:
            "logs/picard/markduplicates/{sample}.markdup.log"
        params:
            "--ASSUME_SORT_ORDER coordinate --REMOVE_DUPLICATES true"
        wrapper:
            "bio/picard/markduplicates"


    use rule samtools_index from bwa_meta as picard_index with:
        input:
            "picard/markduplicates/mapping/{sample}.bam"
        output:
            temp(get_bai("picard/markduplicates/mapping/{sample}.bam"))
        message:
            "Indexing Picard deduplicated bam for {wildcards.sample}"
        log:
            "logs/picard/markduplicates/{sample}.index.log"


    ##############################
    ### GATK BAM RECALIBRATION ###
    ##############################

    gatk_bqsr_config = {
        "threads": config["threads"],
        "genome": config["ref"]["fasta"],
        "dbsnp": config["ref"]["dbsnp"]
    }

    module gatk_bqsr_meta:
        snakefile: "../../meta/bio/gatk_bqsr/test/Snakefile"
        config: gatk_bqsr_config


    use rule samtools_index from bwa_meta as gatk_index with:
        input:
            "gatk/recal_bam/{sample}.bam"
        output:
            temp(get_bai("gatk/recal_bam/{sample}.bam"))
        message:
            "Indexing GATK recalibrated bam for {wildcards.sample}"
        log:
            "logs/gatk/apply_baserecalibrator/{sample}.index.log"


    use rule gatk_apply_baserecalibrator from gatk_bqsr_meta with:
        input:
            bam="picard/markduplicates/mapping/{sample}.bam",
            bam_index=get_bai("picard/markduplicates/mapping/{sample}.bam"),
            ref=config['ref']['fasta'],
            ref_idx=get_fai(config['ref']['fasta']),
            ref_dict=get_dict(config['ref']['fasta']),
            recal_table="gatk/recal_data_table/{sample}.grp"


    use rule gatk_compute_baserecalibration_table from gatk_bqsr_meta with:
        input:
            bam="picard/markduplicates/mapping/{sample}.bam",
            bam_index=get_bai("picard/markduplicates/mapping/{sample}.bam"),
            ref=config['ref']['fasta'],
            ref_idx=get_fai(config['ref']['fasta']),
            ref_dict=get_dict(config['ref']['fasta']),
            known=config['ref']['dbsnp'],
            known_idx=get_tbi(config['ref']['dbsnp'])


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
            time_min=lambda wildcard, attempt: attempt * 45,
            tmpdir=lambda wildcards, input: f"tmp/{wildcards.sample}.tmp"
        params:
            adapters=config.get("fastp_adapters", None),
            extra=config.get("fastp_extra", "")
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
          mem_mb=lambda wildcards, attempt: min(attempt * 1024, 2048),
          time_min=lambda wildcards, attempt: attempt * 45,
          tmpdir=lambda wildcards, input: "tmp/{}.{}.tmp".format(
            wildcards.sample, wildcards.stream
          )
        params:
            input=lambda wildcards, output: fastq_links[output[0]]
        log:
            "logs/bigr_copy/{sample}.{stream}.log"
        wrapper:
            "bio/BiGR/copy"




Authors
-------


* Thibault Dayris

* M boyba Diop

* Marc Deloger
