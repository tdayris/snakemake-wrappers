.. _`varscan2_somatic`:

VARSCAN2_SOMATIC
================

Call variants on BAM files with Varscan2, considering a normal background and a tumor state.


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    default_config_varscan2_somatic = {
        # Your genome sequence
        "genome": "reference/genome.fasta",
        # Path to a BED containing the kit's catpured regions
        "bed": "reference/regions.bed"
    }

    def get_fai(genome_path: str) -> str:
        return genome_path + ".fai"

    def get_bai(bam_path: str) -> str:
        return bam_path + "bai"

    def get_mpileup_input(config: dict[str, str]) -> dict[str, str]:
        """
        If user provides a bed file in the config, then it should be used as a
        region file through samtools mpileup.
        """
        if config["bed"] is None:
            return {
                bam=[
                    "gatk/recal_bam/{sample}_tumor.bam",
                    "gatk/recal_bam/{sample}_normal.bam"
                ]
                bam_index=get_bai("samtools/sort/{sample}.bam",),
                reference_genome=config["genome"],
                reference_genome_idx=get_fai(config["genome"])
            }
        return {

        }

    ## Required modules :
    ## This module includes fasta indexes and fasta dictionnaries
    # module index_datasets:
    #     snakefile: "../../index_datasets/test/Snakefile"
    #     config: config
    #
    # use rule * from index_datasets as index_datasets_*

    ## This module includes GATK recalibration
    # module gatk_bqsr:
    #     snakefile: "../../gatk_bqsr/test/Snakefile"
    #     config: config
    #
    # use rule * from gatk_bqsr as gatk_bqsr_*

    ## This module includes BWA mapping and Samtools corrections
    # module bwa_fixmate:
    #     snakefile: "../../bwa_fixmate/test/Snakefile"
    #     config: config
    #
    # use rule * from bwa_fixmate as bwa_fixmate_*


    """
    This rule concats snp and indel callings from Varscan2 in order to produce a
    full VCF file with both kind of variations.
    """
    rule bcftools_concat:
        input:
            calls=expand(
                "vascan2/somatic/{sample}.{content}.vcf",
                content=["snp", "indel"],
                allow_missing=True
            )
        output:
            temp("vascan2/concat/{sample}.vcf.gz")
        threads: 2
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1025, 4096),
            time_min=lambda wildcards, attempt: attempt * 45
        params:
            "--output-type z --remove-duplicates --allow-overlaps"
        log:
            "logs/varscan/pileup2indel/concat/{sample}.log"
        wrapper:
            "/bio/bcftools/concat"


    """
    This rule performs germline calling with Varscan2
    """
    rule varscan2_somatic:
        input:
            mpileup="samtools/mupleup/{sample}.mpileup.gz"
        output:
            snp=temp("vascan2/somatic/{sample}.snp.vcf.gz")
            indel=temp("varscan2/somatic/{sample}.indel.vcf.gz")
        message:
            "Calling variants on {wildcards.sample} with Varscan2 mpileup2cns"
        threads: 2
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 8192, 20480),
            time_min=lambda wildcards, attempt: attempt * 45
        params:
            extra="--p-value 0.05 --variants"
        log:
            "logs/varscan2/mpileup2cns/{sample}.call.log"
        wrapper:
            "/bio/varscan/somatic"


    """
    This rule runs samtools mpileup to list each single difference between mapped
    reads and reference genome
    """
    rule samtools_mpilup:
        input:
            bam=[
                "picard/markduplicates/mapping/{sample}_tumor.bam",
                "picard/markduplicates/mapping/{sample}_normal.bam"
            ],
            bam_index=[
                get_bai("picard/markduplicates/mapping/{sample}_tumor.bam"),
                get_bai("picard/markduplicates/mapping/{sample}_normal.bam")
            ],
            reference_genome=config["genome"],
            reference_genome_idx=get_fai(config["genome"]),
            bed=config["bed"]
        output:
            temp("samtools/mpileup/{sample}.mpileup.gz")
        message:
            "Building mpilup on {wildcards.sample} with samtools (tumor/normal)"
        threads: 2
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 4096, 20480),
            time_min=lambda wildcards, attempt: attempt * 120
        log:
            "logs/samtools/mpileup/{sample}.log"
        params:
            extra="--count-orphans --no-BAQ"
        wrapper:
            "/bio/samtools/mpileup"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/bcftools/concat`

* :ref:`bio/varscan/somatic`

* :ref:`bio/samtools/mpileup`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

Bam are expected to be mate-fixed (see bwa_fixmate meta-wrapper), and recalibrated (see gatk_bqsr meta-wrapper).




Authors
-------


* Thibault Dayris

