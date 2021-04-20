.. _`mutect2_germline`:

MUTECT2_GERMLINE
================

Call variants on BAM files with Mutect2, considering a single bam file.


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    default_config_mutect2_germline = {
        # Your genome sequence
        "genome": "reference/genome.fasta",
        # Path to a VCF containing AF fields
        "known": "reference/germline.vcf"
        # Path to a BED containing the kit's catpured regions
        "bed": "reference/regions.bed",
        # Path to dbSNP vcf, its tbi should be aside.
        "dbsnp": "reference/dbsnp.vcf.gz"
    }

    def get_fai(genome_path: str) -> str:
        return genome_path + ".fai"

    def get_dict(genome_path: str) -> str:
        return ".".join(".".split(genome_path)[:-1]) + ".dict"

    def get_tbi(vcf_path: str) -> str:
        return vcf_path + "tbi"

    def get_bai(bam_path: str) -> str:
        return bam_path + "bai"


    ## Required modules :
    ## This module includes fasta indexes and fasta dictionnaries
    # module index_datasets:
    #     snakefile: "../../index_datasets/test/Snakefile"
    #     config: config
    #
    # use rule * from index_datasets as index_datasets_*
    #
    ## This module includes BWA mapping and Samtools corrections
    # module bwa_fixmate:
    #     snakefile: "../../bwa_fixmate/test/Snakefile"
    #     config: config
    #
    # use rule * from bwa_fixmate as bwa_fixmate_*

    #####################
    ### Apply filters ###
    #####################

    """
    Filter over estimated contaminations
    """
    rule muterc2_filter:
        input:
            vcf="gatk/mutect2/{sample}.vcf.gz",
            ref=config["genome"],
            fasta_index=get_fai(config["genome"]),
            fasta_dict=get_dict(config["genome"]),
            contamination="summary/{sample}_calculate_contamination.table",
            bam="samtools/sort/{sample}.bam",
            bam_index=get_bai("samtools/sort/{sample}.bam"),
            f1r2="gatk/artifacts_prior/{sample}.artifacts_prior.tar.gz"
        output:
            vcf="gatk/mutect2/filter/contamination/{sample}.vcf.gz"
        message:
            "Filtering Mutect2 calls for {wildcards.sample}"
        threads: 1
        resources:
            time_min=lambda wildcards, attempt: attempt * 45,
            mem_mb=lambda wildcards, attempt: min(attempt * 8192, 20480)
        params:
            extra=""
        wrapper:
            "/bio/gatk/filtermutectcalls"


    ################################
    ### Estimate sequencing bias ###
    ################################
    """
    Build orientation model from f1r2 counts made in Mutect2
    """
    rule learn_read_orientation_model:
        input:
            f1r2="gatk/mutect2/f1r2/{sample}.tar.gz"
        output:
            temp("gatk/artifacts_prior/{sample}.artifacts_prior.tar.gz")
        message:
            "Build model over orientation bias on {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 8192, 15360),
            time_min=lambda wildcards, attempt: attempt * 45
        log:
            "logs/gatk/learnreadorientationmodel/{sample}.log"
        wrapper:
            "/bio/gatk/learnreadorientationmodel"


    ###########################################
    ### Estimate cross-sample contamination ###
    ###########################################


    """
    Estimate possible contaminations
    """
    rule calculate_contamination:
        input:
            summary="{sample}_tumor_pileup_summary.table""
        output:
            table=temp("summary/{sample}_calculate_contamination.table")
        group:
            "Contamination_Estimate"
        message:
            "Summarizing read support for known variant sites to further "
            "estimate contamination on {wildcards.sample}"
        threads: 1
        resource:
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
            bam="samtools/sort/{sample}.bam",
            bam_index=get_bai("samtools/sort/{sample}.bam"),
            interval=config["bed"],
            variants=config["known"],
            variants_index=get_tbi(config["known"])
        output:
            table=temp("gatk/getpileupsummaries/{sample}_getpileupsummaries.table")
        group:
            "Contamination_Estimate"
        message:
            "Summarizing read support for known variant sites to further "
            "estimate contamination on {wildcards.sample}"
        threads: 1
        resource:
            mem_mb=lambda wildcards, attempt: min(attempt * 5120, 15360),
            time_min=lambda wildcards, attempt: attempt * 35
        log:
            "logs/gatk/GetPileupSummaries/{sample}.log"
        wrapper:
            "/bio/gatk/getpileupsummaries"


    ######################
    ### Actual Calling ###
    ######################

    """
    This rule calls germline variants with GATK Mutect2
    """
    rule mutect2_germline:
        input:
            fasta=config["genome"],
            fasta_index=get_fai(config["genome"]),
            fasta_dict=get_dict(config["genome"]),
            map="samtools/sort/{sample}.bam",
            map_index=get_bai("samtools/sort/{sample}.bam"),
            germline=config["known"],
            germline_tbi=get_tbi(config["known"]),
            intervals=config["bed"]
        output:
            vcf=temp("gatk/mutect2/call/{sample}.vcf.gz"),
            f1r2=temp("gatk/mutect2/f1r2/{sample}.tar.gz")
        message:
            "Calling variants on {wildcards.sample} with GATK Mutect2"
        threads: 4
        resources:
            time_min=lambda wildcards, attempt: attempt * 45,
            mem_mb=lambda wildcards, attempt: min(attempt * 8192, 20480)
        params:
            extra=(
                "--max-reads-per-alignment-start 0 "
                "--tumor-sample {sample} "
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

