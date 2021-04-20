.. _`varscan2_germline`:

VARSCAN2_GERMLINE
=================

Call variants on BAM files with Varscan2, considering a single bam file.


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    default_config_varscan2_germline = {
        # Your genome sequence
        "genome": "reference/genome.fasta",
        # Path to a BED containing the kit's catpured regions
        "bed": "reference/regions.bed"
    }

    def get_fai(genome_path: str) -> str:
        return genome_path + ".fai"

    def get_bai(bam_path: str) -> str:
        return bam_path + "bai"

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
                "varscan2/calling/{kind}/{sample}.vcf",
                kind=["snp", "indel"],
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
    This rule calls small indel variants with Varscan2. Results are set to unzipped
    VCF because BCFTools merging will be faster that way. Thus, unzipped VCF are
    temporary files.
    """
    rule varscan2_indel_calling:
        input:
            pileup="samtools/mupleup/{sample}.mpileup.gz",
            sample_list="varscan2/mpileup2cns/{sample}.sample.list"
        output:
            temp("varscan2/calling/indel/{sample}.vcf")
        message: "Calling Indels on {wildcards.sample} with Varscan2"
        threads: 2
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 4096, 15360),
            time_min=lambda wildcards, attempt: attempt * 45
        params:
            extra="--p-value 0.05 --variants"
        log:
            "logs/varscan/pileup2indel/call/{sample}.log"
        wrapper:
            "/bio/varscan/mpileup2indel"


    """
    This rule calls snp variants with Varscan2. Results are set to unzipped
    VCF because BCFTools merging will be faster that way. Thus, unzipped VCF are
    temporary files.
    """
    rule varscan2_snp_calling:
        input:
            pileup="samtools/mupleup/{sample}.mpileup.gz",
            sample_list="varscan2/mpileup2cns/{sample}.sample.list"
        output:
            temp("varscan2/calling/snp/{sample}.vcf")
        message: "Calling SNP on {wildcards.sample} with Varscan2"
        threads: 2
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 4096, 15360),
            time_min=lambda wildcards, attempt: attempt * 45
        params:
            extra="--p-value 0.05 --variants"
        log:
            "logs/varscan/pileup2snp/call/{sample}.log"
        wrapper:
            "/bio/varscan/mpileup2snp"


    """
    This rule provides a list of sample names (here, only one) for Varscan2
    """
    rule varscan2_sample_list:
        output:
            temp("varscan2/mpileup2cns/{sample}.sample.list")
        message:
            "Building sample list for Varscan2 mpileup2cns"
        threads: 1
        resources:
            mem_mb=128,
            time_min=2
        params:
            '"{sample}"'
        log:
            "logs/varscan2/mpileup2cns/{sample}.list.log"
        shell:
            "echo {params} > {output} 2> {log}"


    """
    This rule runs samtools mpileup to list each single difference between mapped
    reads and reference genome
    """
    rule samtools_mpilup:
        input:
            bam="gatk/recal_bam/{sample}.bam",
            bam_index=get_bai("gatk/recal_bam/{sample}.bam"),
            reference_genome=config["genome"],
            reference_genome_idx=get_fai(config["genome"]),
            bed=config["bed"]
        output:
            temp("samtools/mpileup/{sample}.mpileup.gz")
        message: "Building mpilup on {wildcards.sample} with samtools"
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


* :ref:`bio/samtools/mpileup`

* :ref:`bio/varscan/mpileup2snp`

* :ref:`bio/varscan/mpileup2indel`

* :ref:`bio/bcftools/concat`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

Bam are expected to be mate-fixed (see bwa_fixmate meta-wrapper), and recalibrated (see gatk_bqsr meta-wrapper).




Authors
-------


* Thibault Dayris

