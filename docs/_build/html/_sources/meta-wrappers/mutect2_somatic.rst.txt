.. _`mutect2_somatic`:

MUTECT2_SOMATIC
===============

Call variants on BAM files with Mutect2, considering a pair of bam file.


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    default_config_mutect2_germline = {
        # Your genome sequence
        "genome": "reference/genome.fasta",
        # Path to a VCF containing AF fields
        "known": "reference/dbsnp.vcf.gz",
        # Path to a BED containing the kit's catpured regions
        "bed": "reference/regions.bed",
        # Path to dbSNP vcf, its tbi should be aside.
        "dbsnp": "reference/dbsnp.vcf.gz"
    }

    try:
        if config == dict():
            config = default_config_mutect2_germline
    except NameError:
        config = default_config_mutect2_germline


    def get_fai(genome_path: str) -> str:
        return genome_path + ".fai"

    def get_dict(genome_path: str) -> str:
        return ".".join(genome_path.split(".")[:-1]) + ".dict"

    def get_tbi(vcf_path: str) -> str:
        return vcf_path + ".tbi"

    def get_bai(bam_path: str) -> str:
        return bam_path + ".bai"


    ###########################################
    ### Estimate cross-sample contamination ###
    ###########################################


    """
    Estimate possible contaminations
    """
    rule calculate_tumor_contamination:
        input:
            summary="gatk/getpileupsummaries/tumor/{sample}_tumor_getpileupsummaries.table",
            normal=="gatk/getpileupsummaries/normal/{sample}_normal_getpileupsummaries.table",
        output:
            table=temp("summary/{sample}_calculate_contamination.table"),
            segmentation=temp("summary/{sample}_segments.table")
        group:
            "Contamination_Estimate"
        message:
            "Summarizing read support for known variant sites to further "
            "estimate contamination on {wildcards.sample}"
            " (on tumor only)"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 5120, 15360),
            time_min=lambda wildcards, attempt: attempt * 35
        log:
            "logs/gatk/CalculateContamination/{sample}.log"
        wrapper:
            "/bio/gatk/calculatecontamination"


    """
    Summarize the read support over known variants
    """
    rule get_pileup_summaries:
        input:
            bam="samtools/sort/{status}/{sample}_{status}.bam",
            bam_index=get_bai("samtools/sort/{status}/{sample}_{status}.bam"),
            intervals=config["bed"],
            variants=config["known"],
            variants_index=get_tbi(config["known"])
        output:
            table=temp("gatk/getpileupsummaries/{status}/{sample}_{status}_getpileupsummaries.table")
        group:
            "Contamination_Estimate"
        message:
            "Summarizing read support for known variant sites to further "
            "estimate contamination on {wildcards.sample}"
            " ({wildcards.status})"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 5120, 15360),
            time_min=lambda wildcards, attempt: attempt * 35
        log:
            "logs/gatk/GetPileupSummaries/{sample}.{status}.log"
        wrapper:
            "/bio/gatk/getpileupsummaries"


    ######################
    ### Actual Calling ###
    ######################
    """
    This rule calls somatic variants with GATK Mutect2
    """
    rule mutect2_somatic:
        input:
            fasta=config["genome"],
            fasta_index=get_fai(config["genome"]),
            fasta_dict=get_dict(config["genome"]),
            map="samtools/sort/normal/{sample}_normal.bam",
            map_index=get_bai("samtools/sort/normal/{sample}_normal.bam"),
            tumor="samtools/sort/tumor/{sample}_tumor.bam",
            tumor_index=get_bai("samtools/sort/tumor/{sample}_tumor.bam"),
            germline=config["known"],
            germline_tbi=get_tbi(config["known"]),
            intervals=config["bed"]
        output:
            vcf=temp("mutect2/call/{sample}.vcf.gz"),
            f1r2=temp("mutect2/f1r2/{sample}.tar.gz")
        message:
            "Calling variants on {wildcards.sample} with GATK Mutect2"
        threads: 4
        resources:
            time_min=lambda wildcards, attempt: attempt * 45,
            mem_mb=lambda wildcards, attempt: min(attempt * 8192, 20480)
        params:
            extra=(
                "--max-reads-per-alignment-start 0 "
                "--tumor-sample Mutect2_{sample}_tumor "
                "--normal Mutect2_{sample}_normal "
                "--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter "
            )
        log:
            "logs/gatk/mutect2/call/{sample}.log"
        wrapper:
            "/bio/gatk/mutect"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/gatk/filtermutectcalls`

* :ref:`bio/gatk/learnreadorientationmodel`

* :ref:`bio/gatk/calculatecontamination`

* :ref:`bio/gatk/getpileupsummaries`

* :ref:`bio/gatk/mutect`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

Bam are expected to be mate-fixed (see bwa_fixmate meta-wrapper).




Authors
-------


* Thibault Dayris

