.. _`varscan2_calling`:

VARSCAN2_CALLING
================

Call variants on BAM files with Varscan2


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    default_config_varscan2_calling = {
        "genome": "reference/genome.fasta",
    }

    def get_fasta_index_from_genome_path(genome_path: str) -> str:
        return genome_path + ".fai"


    try:
        if config == dict():
            config = default_config_varscan2_calling
    except NameError:
        config = default_config_varscan2_calling


    """
    This rule indexes the gzipped VCF as almost every downstream tool require
    an indexed VCF file.
    """
    rule tabix_index:
        input:
            "bcftools/{sample}.vcf.gz"
        output:
            "bcftools/{sample}.vcf.gz.tbi"
        message: "Indexing {wildcards.sample} calls with Tabix"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1024, 10240),
            time_min=lambda wildcards, attempt: attempt * 60
        params:
            "-p vcf"
        log:
            "logs/tabix/index/{sample}.log"
        wrapper:
            "0.72.0-558-g97f44f15f/bio/tabix"


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
            "bcftools/{sample}.vcf.gz"
        threads: 2
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1025, 4096),
            time_min=lambda wildcards, attempt: attempt * 45
        params:
            "--output-type z --threads 2"
        log:
            "logs/bcftools/concat/{sample}.log"
        wrapper:
            "0.72.0-558-g97f44f15f/bio/bcftools/concat"


    """
    This rule calls small indel variants with Varscan2. Results are set to unzipped
    VCF because BCFTools merging will be faster that way. Thus, unzipped VCF are
    temporary files.
    """
    rule varscan2_indel_calling:
        input:
            "samtools/mpileup/{sample}.mpileup.gz"
        output:
            temp("varscan2/calling/indel/{sample}.vcf")
        message: "Calling Indels on {wildcards.sample} with Varscan2"
        threads: 2
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 4096, 1536),
            time_min=lambda wildcards, attempt: attempt * 45
        params:
            extra="--output-vcf 1"
        log:
            "logs/varscan/pileup2indel/call/{sample}.log"
        wrapper:
            "0.72.0-558-g97f44f15f/bio/varscan/mpileup2indel"


    """
    This rule calls snp variants with Varscan2. Results are set to unzipped
    VCF because BCFTools merging will be faster that way. Thus, unzipped VCF are
    temporary files.
    """
    rule varscan2_snp_calling:
        input:
            "samtools/mpileup/{sample}.mpileup.gz"
        output:
            temp("varscan2/calling/snp/{sample}.vcf")
        message: "Calling SNP on {wildcards.sample} with Varscan2"
        threads: 2
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 4096, 1536),
            time_min=lambda wildcards, attempt: attempt * 45
        params:
            extra="--output-vcf 1"
        log:
            "logs/varscan/pileup2snp/call/{sample}.log"
        wrapper:
            "0.72.0-558-g97f44f15f/bio/varscan/mpileup2snp"


    """
    This rule runs samtools mpileup to list each single difference between mapped
    reads and reference genome
    """
    rule samtools_mpilup:
        input:
            bam="mapped/{sample}.bam",
            reference_genome=config["genome"]
            reference_genome_idx=get_fasta_index_from_genome_path(config["genome"]),
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
            extra=""
        wrapper:
            "0.72.0-558-g97f44f15f/bio/samtools/mpileup"


    """
    This rule indexes the input genome sequence with Samtools. It is not
    explicitely requested by Samtools, but it will crash if the genome sequence
    is not indexed.

    This rule is cached since it should be used only once per reference sequence
    """
    rule samtools_faidx:
        input:
            config["genome"]
        output:
            get_fasta_index_from_genome_path(config["genome"])
        message: "Indexing reference fasta with Samtools"
        cache: True
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1024, 4098),
            time_min=lambda wildcards, attempt: attempt * 45
        params:
            ""
        log:
            "logs/samtools/faidx/{genome}.log"
        wrapper:
            "0.72.0-558-g97f44f15f/bio/samtools/faidx"

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

* :ref:`bio/samtools/faidx`

* :ref:`bio/samtools/mpileup`

* :ref:`bio/tabix`

* :ref:`bio/varscan/mpileup2indel`

* :ref:`bio/varscan/mpileup2snp`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

Bam are expected to be mate-fixed, and recalibrated.




Authors
-------


* Thibault Dayris

